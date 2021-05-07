"""A class defining an airfoil."""

import math
import json
import copy
import os
import operator
import warnings

import matplotlib.pyplot as plt
import scipy.interpolate as interp
import subprocess as sp
import numpy as np

from .poly_fits import multivariablePolynomialFit, multivariablePolynomialFunction, autoPolyFit, multivariableRMS, multivariableR2
from .exceptions import DatabaseBoundsError, CamberSolverNotConvergedError, PolyFitBoundsError


class Airfoil:
    """A class defining an airfoil. If the airfoil geometry is defined using outline
    points, then when this class is initialized, a solver will be automatically
    run to determine the camber line and thickness distribution of the airfoil.
    The parameters "camber_relaxation", "le_loc", and "camber_termination_tol" then
    have bearing on this solver.

    When using the member methods get_CL, get_CD, get_Cm, get_CLa, get_CLM, get_CLRe,
    and get_aL0, the default parameters are dependent upon the type of airfoil.

    For all airfoil types except 'functional', 'alpha', 'trailing_flap_deflection',
    and 'trailing_flap_fraction' default to 0.0.

    For 'functional' airfoils, the defaults are specified by the user.

    For 'linear' airfoils, 'Mach' and 'Rey' have no effect on computations.

    For 'database' airfoils, 'Mach' and 'Rey' default to the average value in the database,
    if the database is dependent on that variable. Otherwise, they default to the value
    specified as constant when the database was generated.

    For 'poly_fit' airfoils, 'Mach' and 'Rey' default to the values given in the fit file.
    If export_polynomial_fits() is used, these values are the same as those for the 
    'database' type.

    Parameters
    ----------
    name : str
        Name of the airfoil.

    airfoil_input : dict or str
        Dictionary or path to JSON object describing the airfoil.

    verbose : bool, optional
        Whether to display information on the progress of parameterizing the geometry
        for an airfoil defined by a set of outline points. Defaults to False.

    camber_relaxation : float, optional
        A value between 0.0 and 1.0 that defines how much of the update at each
        iteration of the camber line solver should be accepted. Helpful for some
        poorly behaved cases. Defaults to 1.0 (full update).
    
    le_loc : list, optional
        Gives the location of the leading edge relative to the given points for an airfoil
        defined by a set of outline points. If this is given, the camber line will be forced
        to intersect this point. THIS POINT SHOULD LIE ON THE AIRFOIL OUTLINE. If not given,
        the camber line solver will try to iteratively find the leading edge (the point where
        the camber line intersects the front of the profile).

    camber_termination_tol : float, optional
        The tolerance below which the maximum approximate error in the camber line estimate
        must fall in order for the camber line solver to terminate. Defaults to 1e-10.

    max_iterations : int, optional
        Maximum number of iterations for the camber line solver. Defaults to 100.
    """

    def __init__(self, name, airfoil_input, **kwargs):
        
        self.name = name
        self._load_params(airfoil_input)
        if "camber_solver_kwargs" in list(self._input_dict.keys()):
            kwarg_dict = self._input_dict["camber_solver_kwargs"]
            self._verbose = kwarg_dict.get("verbose", False)
            self._camber_relaxation = kwarg_dict.get("camber_relaxation", 1.0)
            self._le_loc = kwarg_dict.get("le_loc", None)
            self._camber_termination_tol = kwarg_dict.get("camber_termination_tol", 1e-10)
            self._max_iterations = kwarg_dict.get("max_iterations", 100)
        else:
            self._verbose = kwargs.get("verbose", False)
            self._camber_relaxation = kwargs.get("camber_relaxation", 1.0)
            self._le_loc = kwargs.get("le_loc", None)
            self._camber_termination_tol = kwargs.get("camber_termination_tol", 1e-10)
            self._max_iterations = kwargs.get("max_iterations", 100)

        # Load flaps
        self._load_flaps()

        # Store undeformed outlines
        self._initialize_geometry()

        # Specify database DOF parameters
        self._allowable_dofs = ["alpha", "Rey", "Mach", "trailing_flap_deflection", "trailing_flap_fraction"]
        self._dof_defaults = {
            "alpha" : 0.0,
            "Rey" : 1e6,
            "Mach" : 0.0,
            "trailing_flap_deflection" : 0.0,
            "trailing_flap_fraction" : 0.0
        }

        # Load input file
        input_file = self._input_dict.get("input_file", None)
        if input_file is not None:
            if self._type == "database":
                self.import_database(filename=input_file)
            elif self._type == "poly_fit":
                self.import_polynomial_fits(filename=input_file)


        self._raise_poly_bounds_error = True


    def set_err_state(self, **kwargs):
        """Sets the error state for the airfoil. Each may be specified as 'raise' or 'ignore'.

        Parameters
        ----------
        poly_fit_bounds : str, optional
            How to handle PolyFitBoundsError. Defaults to 'raise'.
        """

        # Polynomial fit bounds
        if kwargs.get("poly_fit_bounds", "raise") == "raise":
            self._raise_poly_bounds_error = True
        else:
            self._raise_poly_bounds_error = False


    def set_verbosity(self, verbosity):
        """Sets the verbosity of the airfoil."""
        self._verbose = verbosity


    def set_type(self, database_type):
        """Determines how the aerodynamic coefficients will be calculated.

        Parameters
        ----------
        database_type: str
            "linear", "functional", "database", or "poly_fit". Airfoil
            will automatically check if it has the necessary to perform
            the given type of computation and throw a warning if it does not.

        """

        # Check for proper type spec
        if database_type not in ["linear", "database", "poly_fit", "functional"]:
            raise IOError("{0} is not a valid type specification.".format(database_type))

        # Check for linear data
        if database_type == "linear":
            if not hasattr(self, "_CLa"):
                raise RuntimeWarning("Airfoil {0} does not have linear coefficients specified. Reverting to type '{1}' for computations.".format(self.name, self._type))
            else:
                self._type = database_type

        # Check for database
        if database_type == "database":
            if not hasattr(self, "_data"):
                raise RuntimeWarning("Airfoil {0} does not have a database of coefficients. Reverting to type '{1}' for computations.".format(self.name, self._type))
            else:
                self._type = database_type

            # Set up data normalization
            self._data_norms = np.zeros((1,self._num_dofs))
            for i in range(self._num_dofs):
                self._data_norms[0,i] = np.max(np.abs(self._data[:,i]))

            # Make sure we don't divide by zero
            self._data_norms[np.where(self._data_norms==0.0)] = 1.0

            # Normalize independent vars
            self._normed_ind_vars = self._data[:,:self._num_dofs]/self._data_norms

            # Determine default Mach and Reynolds number
            if "Rey" in list(self._dof_db_cols.keys()):
                i = self._dof_db_cols["Rey"]
                Re_min = np.min(self._data[:,i])
                Re_max = np.max(self._data[:,i])
                self._dof_defaults["Rey"] = 0.5*(Re_max+Re_min)

            if "Mach" in list(self._dof_db_cols.keys()):
                i = self._dof_db_cols["Mach"]
                M_min = np.min(self._data[:,i])
                M_max = np.max(self._data[:,i])
                self._dof_defaults["Mach"] = 0.5*(M_max+M_min)

        # Check for polynomial fits
        if database_type == "poly_fit":
            if not hasattr(self, "_CL_poly_coefs"):
                raise RuntimeWarning("Airfoil {0} does not have a set of polynomial fits. Reverting to type '{1}' for computations.".format(self.name, self._type))
            else:
                self._type = database_type

        # Check for functional definition
        if database_type == "functional":
            if not hasattr(self, "_CL"):
                raise RuntimeWarning("Airfoil {0} does not have functional definitions of coefficients. Reverting to type '{1}' for computations.".format(self.name, self._type))
            else:
                self._type = database_type


    def _load_params(self, airfoil_input):

        # Load input dict
        if isinstance(airfoil_input, str):
            # Load JSON object
            with open(airfoil_input, 'r') as json_handle:
                self._input_dict = json.load(json_handle)
        elif isinstance(airfoil_input, dict):
            self._input_dict = airfoil_input
        else:
            raise IOError("{0} is not an allowed airfoil definition. Must be path or dictionary.".format(airfoil_input))

        # Store type
        self._type = self._input_dict.get("type", "linear")

        # If linear, store coefficients
        if self._type == "linear":
            self._aL0 = self._input_dict.get("aL0", 0.0)
            self._CLa = self._input_dict.get("CLa", 2*np.pi)
            self._Cma = self._input_dict.get("Cma", 0.0)
            self._CD0 = self._input_dict.get("CD0", 0.0)
            self._CD1 = self._input_dict.get("CD1", 0.0)
            self._CD2 = self._input_dict.get("CD2", 0.0)
            self._CL_max = self._input_dict.get("CL_max", np.inf)
            if abs(self._CL_max) < 1e-10:
                warnings.warn("You have specified a maximum lift coefficient of 0. Are you sure you want to do this?...")
            self._CmL0 = self._input_dict.get("CmL0", 0.0)

        # For functional, store functions
        elif self._type == "functional":
            self._CL = self._input_dict["CL"]
            self._CD = self._input_dict["CD"]
            self._Cm = self._input_dict["Cm"]


    def _load_flaps(self):
        # Loads flaps based on the input dict
        self._trailing_flap_type = self._input_dict.get("trailing_flap_type", None)
        self._trailing_flap_hinge_height = self._input_dict.get("trailing_flap_hinge_height", 0.0)


    def _initialize_geometry(self):
        # Initialize geometry based on whether points or a NACA designation were given

        # Check that there's only one geometry definition
        geom_dict = self._input_dict.get("geometry", {})
        outline_points = geom_dict.get("outline_points", None)
        NACA = geom_dict.get("NACA", None)
        if outline_points is not None and NACA is not None:
            raise IOError("Outline points and a NACA designation may not be both specified for airfoil {0}.".format(self.name))

        # Check for user-given points
        if outline_points is not None:
            self.geom_specification = "points"

            # Array given
            if isinstance(outline_points, np.ndarray):
                self._raw_outline = outline_points
            
            # File
            else:

                # Import
                with open(outline_points, 'r') as input_handle:
                    self._raw_outline = np.genfromtxt(input_handle)

                # Check for comma-delimited
                if np.isnan(self._raw_outline).any():
                    with open(outline_points, 'r') as input_handle:
                        self._raw_outline = np.genfromtxt(input_handle, delimiter=',')

            # Get number of points
            self._N_orig = self._raw_outline.shape[0]

            # Rearrange the coordinates if necessary to be top first then bottom
            top_ind = np.argmax(self._raw_outline[:,1])
            bot_ind = np.argmin(self._raw_outline[:,1])
            if bot_ind < top_ind: # Bottom came first
                self._raw_outline = self._raw_outline[::-1]

            # Calculate camber line and thickness
            self._calc_geometry_from_points()

        # NACA definition
        elif NACA is not None:
            self.geom_specification = "NACA"
            self._naca_closed_te = geom_dict.get("NACA_closed_te", False)
            self._naca_des = NACA

            self._calc_geometry_from_NACA()

        # No geometry given
        else:
            self.geom_specification = "none"
            self._max_camber = geom_dict.get("max_camber", 0.0)
            self._max_thickness = geom_dict.get("max_thickness", 0.0)
            return

    
    def _calc_geometry_from_NACA(self):
        # Creates thickness, camber, camber derivative, and outline splines base on the NACA equations

        # 4-digit series
        if len(self._naca_des) == 4:

            # Stores the camber and thickness getters based on the NACA designation of the airfoil
            self._m = float(self._naca_des[0])/100
            self._p = float(self._naca_des[1])/10
            self._t = float(self._naca_des[2:])/100
            self._max_camber = self._m
            self._max_thickness = self._t

            # Camber line
            def camber(x):
                if self._p != 0.0:
                    return np.where(x<self._p, self._m/(self._p*self._p)*(2*self._p*x-x*x), self._m/((1-self._p)*(1-self._p))*(1-2*self._p+2*self._p*x-x*x))
                else:
                    return np.zeros_like(x)

            self._camber_line = camber

            # Camber line derivative
            def camber_deriv(x):
                if abs(self._m)<1e-10 or abs(self._p)<1e-10: # Symmetric
                    return np.zeros_like(x)
                else:
                    return np.where(x<self._p, 2*self._m/(self._p*self._p)*(self._p-x), 2*self._m/((1-self._p)*(1-self._p))*(self._p-x))

            self._camber_deriv = camber_deriv

            # Thickness
            def thickness(x):
                if self._naca_closed_te:
                    return 5.0*self._t*(0.2969*np.sqrt(x)-0.1260*x-0.3516*x*x+0.2843*x*x*x-0.1036*x*x*x*x)
                else:
                    return 5.0*self._t*(0.2969*np.sqrt(x)-0.1260*x-0.3516*x*x+0.2843*x*x*x-0.1015*x*x*x*x)

            self._thickness = thickness

        else:
            raise IOError("Only 4-digit NACA specifications may be given. {0} is not 4-digit.".format(self._naca_des))

        # Cosine distribution of chord locations
        theta = np.linspace(-np.pi, np.pi, 200)
        x = 0.5*(1-np.cos(theta))

        y_c = self._camber_line(x)
        dy_c_dx = self._camber_deriv(x)
        t = self._thickness(x)

        # Outline points
        X = x-t*np.sin(np.arctan(dy_c_dx))*np.sign(theta)
        Y = y_c+t*np.cos(np.arctan(dy_c_dx))*np.sign(theta)

        outline_points = np.concatenate([X[:,np.newaxis], Y[:,np.newaxis]], axis=1)

        # Create splines defining the outline
        self._x_outline, self._y_outline = self._create_splines_of_s(outline_points)
        _,_,self._s_le = self._get_intersection_point(0.0, 0.0, dy_c_dx[0], "leading_edge")


    def _create_splines_of_s(self, outline_points):
        # Creates x and y coordinate splines which are a function of the normalized distance along the outline

        # Calculate distances
        x_diff = np.diff(outline_points[:,0])
        y_diff = np.diff(outline_points[:,1])
        ds = np.sqrt(x_diff*x_diff+y_diff*y_diff)
        ds = np.insert(ds, 0, 0.0)

        # Create s distribution
        s = np.cumsum(ds)
        s_normed = s/s[-1]

        # Create splines
        x_outline = interp.UnivariateSpline(s_normed, outline_points[:,0], k=5, s=1e-10)
        y_outline = interp.UnivariateSpline(s_normed, outline_points[:,1], k=5, s=1e-10)

        return x_outline, y_outline


    def _calc_geometry_from_points(self):
        # Calculates the camber line and thickness distribution from the outline points.
        # Also scales the airfoil, places the leading edge at the origin, and rotates it so it is at zero angle of attack.
        if self._verbose:
            print("Calculating camber line...")
            print("{0:<20}{1:<20}".format("Iteration", "Max Approx Error"))
            print("".join(["-"]*40))

        # Camber line estimate parameters
        num_camber_points = self._N_orig//5 # As recommended by Cody Cummings to avoid discontinuities in dyc/dx
        camber_deriv_edge_order = 2

        # Create splines defining the outline
        self._x_outline, self._y_outline = self._create_splines_of_s(self._raw_outline)

        # Get distribution of x locations
        le_ind = np.argmin(np.abs(self._raw_outline[:,0]))
        te = (self._raw_outline[0]+self._raw_outline[-1])*0.5
        if self._le_loc is None:
            x_tip = np.min(self._raw_outline[:,0])
            x_c = np.linspace(x_tip+0.0001, te[0], num_camber_points)
        else:
            x_c = np.linspace(self._le_loc[0], te[0], num_camber_points)

        # Get initial estimate for the camber line
        y_t = np.interp(x_c, self._raw_outline[le_ind::-1,0], self._raw_outline[le_ind::-1,1])
        y_b = np.interp(x_c, self._raw_outline[le_ind:,0], self._raw_outline[le_ind:,1])
        y_c = 0.5*(y_t+y_b)
        if self._le_loc is not None:
            y_c[0] = self._le_loc[1]

        # A useful distribution
        x_space = np.linspace(0.0, 1.0, 1000)

        # Show
        if self._verbose:
            plt.figure()
            plt.plot(self._raw_outline[:,0], self._raw_outline[:,1], 'b-', label='Outline Data')
            plt.plot(self._x_outline(x_space), self._y_outline(x_space), 'k--', label='Outline Spline')
            plt.plot(x_c, y_c, 'r--', label='Initial Camber Line Estimate')
            plt.plot(x_c, np.gradient(y_c, x_c, edge_order=camber_deriv_edge_order), 'g--', label='Camber Line Derivative')
            plt.legend()
            plt.gca().set_aspect('equal', adjustable='box')
            plt.show()

        # Initialize
        x_t = np.zeros(num_camber_points)
        x_b = np.zeros(num_camber_points)
        x_c_new = np.zeros(num_camber_points)
        y_c_new = np.zeros(num_camber_points)

        # Check for symmetric airfoil
        if np.allclose(y_c, 0.0):
            camber_error = 0.0
        else:
            camber_error = 1.0

        # Iterate until convergence
        iteration = 0
        while camber_error > self._camber_termination_tol:
            iteration += 1

            # Determine camber line slope
            dyc_dx = np.gradient(y_c, x_c, edge_order=camber_deriv_edge_order)
            
            # Determine slope of lines perpendicular to the camber
            with np.errstate(divide="ignore"):
                b = -1.0/dyc_dx

            # Loop through points on the camber line to find where their normal intersects the outline
            for i in range(num_camber_points):
                if self._le_loc is not None and i == 0:
                    continue

                # Get point information
                xc = x_c[i]
                yc = y_c[i]
                bi = b[i]

                # Estimate the intersection points
                x_t[i], y_t[i], _ = self._get_intersection_point(xc, yc, bi, "top")
                x_b[i], y_b[i], _ = self._get_intersection_point(xc, yc, bi, "bottom")

            # Calculate new camber line points
            x_c_new = 0.5*(x_t+x_b)
            y_c_new = 0.5*(y_t+y_b)
            if self._le_loc is not None:
                x_c_new[0] = self._le_loc[0]
                y_c_new[0] = self._le_loc[1]

            # Plot new and old estimate
            if False:
                plt.figure()
                plt.plot(self._raw_outline[:,0], self._raw_outline[:,1], 'b-', label='Outline Data')
                plt.plot(x_c, y_c, 'r--', label='Old Camber Line Estimate')
                plt.plot(x_c_new, y_c_new, 'g--', label='New Camber Line Estimate')
                plt.plot(x_t, y_t, 'rx')
                plt.plot(x_b, y_b, 'ro')
                plt.legend()
                plt.gca().set_aspect('equal', adjustable='box')
                plt.show()

            # Approximate error
            x_diff = x_c_new-x_c
            y_diff = y_c_new-y_c
            camber_error = np.max(np.sqrt(x_diff*x_diff+y_diff*y_diff))
            if self._verbose:
                print("{0:<20}{1:<20}".format(iteration, camber_error))

            # Sort, just in case things got messed up
            sorted_ind = np.argsort(x_c_new)
            x_c_new = x_c_new[sorted_ind]
            y_c_new = y_c_new[sorted_ind]

            # Update for next iteration
            x_c = self._camber_relaxation*x_c_new+(1.0-self._camber_relaxation)*x_c
            y_c = self._camber_relaxation*y_c_new+(1.0-self._camber_relaxation)*y_c

            # Check iteration limit
            if iteration == self._max_iterations:
                raise CamberSolverNotConvergedError(self.name, camber_error)

        if self._verbose: print("Camber line solver converged.")

        # Calculate where the camber line intersects the outline to find the leading edge
        if self._le_loc is None:
            dyc_dx = np.gradient(y_c, x_c, edge_order=2)
            b = dyc_dx[0]
            x_le, y_le, self._s_le = self._get_intersection_point(x_c[0], y_c[0], b, "leading_edge")
            le = np.array([x_le, y_le])
            x_c = np.insert(x_c, 0, le[0])
            y_c = np.insert(y_c, 0, le[1])
        else:
            le = self._le_loc

        if self._verbose: print("Leading edge: {0}".format(le))

        if self._verbose:
            # Plot
            plt.figure()
            plt.plot(self._raw_outline[:,0], self._raw_outline[:,1], 'b-', label='Outline Data')
            plt.plot(x_c, y_c, 'g--', label='Final Camber Line Estimate')
            plt.legend()
            plt.gca().set_aspect('equal', adjustable='box')
            plt.show()

        # Get trailing edge
        te = (self._raw_outline[0]+self._raw_outline[-1])*0.5
        if self._verbose: print("Trailing edge: {0}".format(te))

        # Renormalize using new leading and trailing edge
        self._normalize_points(self._raw_outline, le, te)
        self._x_outline, self._y_outline = self._create_splines_of_s(self._raw_outline)

        # Normalize camber line points
        camber_points = np.concatenate([x_c[:,np.newaxis], y_c[:,np.newaxis]], axis=1)
        self._normalize_points(camber_points, le, te)

        # Store camber and max camber
        self._camber_line = interp.UnivariateSpline(camber_points[:,0], camber_points[:,1], k=5, s=1e-10)
        self._max_camber = np.max(camber_points[:,1])

        # Store camber line derivative
        y_c = self._camber_line(x_space)
        dyc_dx= np.gradient(y_c, x_space)
        self._camber_deriv = interp.UnivariateSpline(x_space, dyc_dx, k=5, s=1e-10)
        with np.errstate(divide="ignore"):
            b = -1.0/dyc_dx

        # Find points on the surface to determine the thickness
        if self._verbose: print("Calculating thickness distribution...", end="")
        x_t_t = np.zeros_like(x_space)
        x_b_t = np.zeros_like(x_space)
        y_t_t = np.zeros_like(x_space)
        y_b_t = np.zeros_like(x_space)
        for i, xc in enumerate(x_space):
            if i == 0: continue # This point intersects the outline by definition and so will have zero thickness
            yc = y_c[i]
            bi = b[i]

            x_t_t[i], y_t_t[i],_ = self._get_intersection_point(xc, yc, bi, "top")
            x_b_t[i], y_b_t[i],_ = self._get_intersection_point(xc, yc, bi, "bottom")

        # Store thickness distribution
        t = 0.5*np.sqrt((x_t_t-x_b_t)*(x_t_t-x_b_t)+(y_t_t-y_b_t)*(y_t_t-y_b_t))
        self._thickness = interp.UnivariateSpline(x_space, t, k=5, s=1e-10)
        self._max_thickness = np.max(t)
        if self._verbose: print("Done")

        # Calculate estimated top and bottom points
        y_c_pred = self._camber_line(x_space)
        dyc_dx = np.gradient(y_c_pred, x_space)
        t_pred = self._thickness(x_space)
        theta = np.arctan(dyc_dx)
        S_theta = np.sin(theta)
        C_theta = np.cos(theta)

        x_b = x_space+t_pred*S_theta
        y_b = y_c_pred-t_pred*C_theta
        x_t = x_space-t_pred*S_theta
        y_t = y_c_pred+t_pred*C_theta
        
        # Plot
        if self._verbose:
            plt.plot(self._raw_outline[:,0], self._raw_outline[:,1], 'b-', label='Original Data')
            plt.plot(x_t, y_t, 'r--', label='Top Fit')
            plt.plot(x_b, y_b, 'r--', label='Bottom Fit')
            plt.legend()
            plt.gca().set_aspect('equal', adjustable='box')
            plt.show()


    def _normalize_points(self, points, le_coords, te_coords):
        # Takes the given points, translates them to the origin, rotates them to be at zero angle of attack, and scales so the chord length is unity

        # Translate
        points[:,0] -= le_coords[0]
        points[:,1] -= le_coords[1]

        # Rotate
        x_diff = te_coords[0]-le_coords[0]
        y_diff = te_coords[1]-le_coords[1]
        theta = -np.arctan2(y_diff, x_diff)
        R = np.array([[np.cos(theta), -np.sin(theta)],
                      [np.sin(theta),  np.cos(theta)]])
        points = np.matmul(R, points.T)

        #Scale
        c_unscaled = np.sqrt(x_diff*x_diff+y_diff*y_diff)
        points = points/c_unscaled


    def _get_intersection_point(self, xc, yc, b, surface, plot=False):
        # Calculates the point on the surface where the line extending from (xc, yc) with slope of b intersects
        # Uses the secant method to converge in s

        # Start s in the middle of the respective surface
        if surface == "top":
            s0 = 0.24
            s1 = 0.26
        elif surface == "bottom":
            s0 = 0.74
            s1 = 0.76
        elif surface == "leading_edge":
            s0 = 0.49
            s1 = 0.51

        # Initial points on outline
        x0 = self._x_outline(s0)
        y0 = self._y_outline(s0)
        x1 = self._x_outline(s1)
        y1 = self._y_outline(s1)

        # Initial distances
        d0 = self._distance(x0, y0, xc, yc, b)
        d1 = self._distance(x1, y1, xc, yc, b)

        # Secant method
        while abs(d1) > 1e-10:
            
            # Get new estimate in s
            if d1 > 0.2 or (s1 > 0.35 and s1 < 0.65): # Apply some relaxation when we're far away or near the leading edge (to keep from shooting to the other surface)
                s2 = s1-0.2*d1*(s0-s1)/(d0-d1)
            else:
                s2 = s1-d1*(s0-s1)/(d0-d1)

            # Make sure we're in bounds
            if s2 > 1.1:
                s2 = 1-0.01*s2
            if s2 < -0.1:
                s2 = -0.01*s2

            # Get new point
            x2 = self._x_outline(s2)
            y2 = self._y_outline(s2)

            # Get new distance
            d2 = self._distance(x2, y2, xc, yc, b)

            # Plot
            if plot:
                plt.figure()
                s_space = np.linspace(0.0, 1.0, 10000)
                plt.plot(self._x_outline(s_space), self._y_outline(s_space), 'b')
                plt.plot(xc, yc, 'or')
                plt.plot(x0, y0, 'bx')
                plt.plot(x1, y1, 'gx')
                plt.plot(x2, y2, 'rx')
                plt.gca().set_aspect('equal', adjustable='box')
                plt.show()

            # Update for next iteration
            s0 = s1
            d0 = d1
            x0 = x1
            y0 = y1

            s1 = s2
            d1 = d2
            x1 = x2
            y1 = y2

        return x1, y1, s1


    def _distance(self, x0, y0, x, y, b):
        if not np.isnan(b) and not np.isinf(b):
            return (-b*x0+y0+b*x-y)/np.sqrt(b*b+1)
        else:
            return abs(x0-x)


    def check_against_NACA(self, naca_des):
        """Checks the error in the camber and thickness against that predicted by the NACA equations. This is recommended as a check for the user
        if unusual geometries are being imported. Checks against the open trailing edge formulation of the NACA equations.

        Parameters
        ----------
        naca_des : str
            NACA designation of the airfoil to compare against as a string. May only be 4-digit series.
        """

        print("Checking estimated thickness and camber against NACA equations for NACA {0}...".format(naca_des))

        # 4-digit series
        if len(naca_des) == 4:

            # Generate true coordinates
            x_space = np.linspace(0.0, 1.0, 10000)
            m = float(naca_des[0])/100
            p = float(naca_des[1])/10
            if p != 0.0:
                y_c_true = np.where(x_space<p, m/p**2*(2*p*x_space-x_space**2), m/(1-p)**2*(1-2*p+2*p*x_space-x_space**2))
            else:
                y_c_true = np.zeros_like(x_space)

            # True thickness
            t = float(naca_des[2:])/100
            t_true = 5*t*(0.2969*np.sqrt(x_space)-0.1260*x_space-0.3516*x_space**2+0.2843*x_space**3-0.1015*x_space**4)

        else:
            raise IOError("Can only compare to NACA 4-digit series. Cannot be {0}.".format(naca_des))

        # Estimate camber
        y_c_est = self._camber_line(x_space)

        # Plot camber
        plt.figure()
        plt.plot(x_space, y_c_true, 'g--', label='True Camber')
        plt.plot(x_space, y_c_est, 'r--', label='Estimated Camber')
        plt.legend()
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

        # Estimate thickness
        t_est = self._thickness(x_space)

        # Plot thickness
        plt.figure()
        plt.plot(x_space, t_true, 'g--', label='True Thickness')
        plt.plot(x_space, t_est, 'r--', label='Estimated Thickness')
        plt.legend()
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

        with np.errstate(divide='ignore', invalid='ignore'):
            # Camber error
            camber_error = np.abs((y_c_est-y_c_true)/y_c_true)
            max_camber_error = np.max(np.where(np.isfinite(camber_error), camber_error, 0.0))
            print("Max camber error: {0:.5f}%".format(max_camber_error*100))

            # Thickness error
            thickness_error = np.abs((t_est-t_true)/t_true)
            max_thickness_error = np.max(np.where(np.isfinite(thickness_error), thickness_error, 0.0))
            print("Max thickness error: {0:.5f}%".format(max_thickness_error*100))


    def _get_flap_influence(self, c_f, d_f):
        # Returns the influence of the flap on the effective angle of attack
        # From Phillips *Mechanics of Flight* Ch. 1, Sec. 7
        theta_f = np.arccos(2.0*c_f-1.0)
        eps_flap_ideal = 1.0-(theta_f-np.sin(theta_f))/np.pi
        hinge_eff = 3.9598*np.arctan((c_f+0.006527)*89.2574+4.898015)-5.18786
        flap_eff = np.where(np.abs(d_f)>0.19198621771937624, -0.4995016675499485*np.abs(d_f)+1.09589743589744, 1.0)
        return hinge_eff*flap_eff*eps_flap_ideal*d_f


    def _get_database_data(self, data_index, **kwargs):
        # Returns an interpolated data point from the database.
        # data_index: 0 = CL, 1 = CD, 2 = Cm

        # Determine size of query
        max_size = 1
        for dof in self._dof_db_order:
            param = kwargs.get(dof, self._dof_defaults[dof])
            if isinstance(param, np.ndarray):
                max_size = max(max_size, param.shape[0])

        # Get params
        param_vals = np.zeros((max_size,self._num_dofs))
        for i, dof in enumerate(self._dof_db_order):
            param_vals[:,i] = kwargs.get(dof, self._dof_defaults[dof])

        # Interpolate
        return_val = interp.griddata(self._normed_ind_vars, self._data[:,self._num_dofs+data_index].flatten(), param_vals/self._data_norms, method='linear').flatten()
        #return_val = interp.griddata(self._data[:,:self._num_dofs], self._data[:,self._num_dofs+data_index].flatten(), param_vals, method='linear').flatten()

        # Check for going out of bounds
        if np.isnan(return_val).any():
            raise DatabaseBoundsError(self.name, np.argwhere(np.isnan(return_val)).flatten(), kwargs)

        # Return
        if max_size == 1:
            return return_val.item()
        else:
            return return_val


    def get_CL(self, **kwargs):
        """Returns the coefficient of lift. Note: all parameters can be given as numpy arrays, in which case a numpy array of the coefficient will be returned.
        To do this, all parameter arrays must have only one dimension and must have the same length.

        Parameters
        ----------
        alpha : float, optional
            Angle of attack in radians. Defaults to 0.0.

        Rey : float, optional
            Reynolds number.

        Mach : float, optional
            Mach number.

        trailing_flap_deflection : float, optional
            Trailing flap deflection in radians. Defaults to 0.

        trailing_flap_fraction : float, optional
            Trailing flap fraction of the chord length. Defaults to 0.

        Returns
        -------
        float or ndarray
            Lift coefficient
        """
        
        # Linearized model
        if self._type == "linear":

            # Get params
            alpha = kwargs.get("alpha", 0.0)
            d_f = kwargs.get("trailing_flap_deflection", 0.0)
            c_f = kwargs.get("trailing_flap_fraction", 0.0)

            # Calculate lift coefficient
            if np.array(d_f == 0.0).all() or np.array(c_f == 0.0).all():
                CL = self._CLa*(alpha-self._aL0)
            else:
                delta_a_d_f = self._get_flap_influence(c_f, d_f)
                CL = self._CLa*(alpha-self._aL0+delta_a_d_f)
            
            # Saturate
            with np.errstate(invalid='ignore'):
                if isinstance(CL, np.ndarray):
                    CL = np.where((CL > self._CL_max) | (CL < -self._CL_max), np.sign(CL)*self._CL_max, CL)
                elif CL > self._CL_max or CL < -self._CL_max:
                    CL = np.sign(CL)*self._CL_max

        # Functional model
        elif self._type == "functional":
            CL = self._CL(**kwargs)

        # Generated/imported database
        elif self._type == "database":
            CL = self._get_database_data(0, **kwargs)

        # Fits
        elif self._type == "poly_fit":
            CL = self._get_polynomial_data(0, **kwargs)

        return CL


    def get_CD(self, **kwargs):
        """Returns the coefficient of drag. note: all parameters can be given as numpy arrays, in which case a numpy array of the coefficient will be returned.
        to do this, all parameter arrays must have only one dimension and must have the same length.

        Parameters
        ----------
        alpha : float, optional
            Angle of attack in radians. Defaults to 0.0.

        Rey : float, optional
            Reynolds number.

        Mach : float, optional
            Mach number.

        trailing_flap_deflection : float, optional
            Trailing flap deflection in radians. Defaults to 0.

        trailing_flap_fraction : float, optional
            Trailing flap fraction of the chord length. Defaults to 0.

        Returns
        -------
        float or ndarray
            Drag coefficient
        """

        # Linear type
        if self._type == "linear":
            d_f = kwargs.pop("trailing_flap_deflection", 0.0)
            CL = self.get_CL(**kwargs)
            CD_flap = 0.002*np.abs(np.degrees(d_f)) # A *very* rough estimate for flaps
            CD = self._CD0+self._CD1*CL+self._CD2*CL**2+CD_flap

        # Functional model
        elif self._type == "functional":
            CD = self._CD(**kwargs)

        # Generated/imported database
        elif self._type == "database":
            CD = self._get_database_data(1, **kwargs)

        # Fits
        elif self._type == "poly_fit":
            CD = self._get_polynomial_data(1, **kwargs)

        return CD


    def get_Cm(self, **kwargs):
        """Returns the moment coefficient. note: all parameters can be given as numpy arrays, in which case a numpy array of the coefficient will be returned.
        to do this, all parameter arrays must have only one dimension and must have the same length.

        Parameters
        ----------
        alpha : float, optional
            Angle of attack in radians. Defaults to 0.0.

        Rey : float, optional
            Reynolds number.

        Mach : float, optional
            Mach number.

        trailing_flap_deflection : float, optional
            Trailing flap deflection in radians. Defaults to 0.

        trailing_flap_moment_deriv : float or ndarray, optional
            Change in section moment with respect to trailing flap deflection. Defaults to 0.

        trailing_flap_fraction : float, optional
            Trailing flap fraction of the chord length. Defaults to 0.

        Returns
        -------
        float or ndarray
            Moment coefficient
        """

        # Linear type
        if self._type == "linear":

            # Get parameters
            alpha = kwargs.get("alpha", 0.0)
            d_f = kwargs.get("trailing_flap_deflection", 0.0)
            c_f = kwargs.get("trailing_flap_fraction", 0.0)

            # No control deflection
            if np.array(d_f == 0.0).all() or np.array(c_f == 0.0).all():
                Cm = self._CmL0+self._Cma*(alpha-self._aL0)
            else:
                theta_f = np.arccos(2.0*c_f-1.0)
                Cm_df = 0.25*(np.sin(2.0*theta_f)-2.0*np.sin(theta_f))
                Cm = self._CmL0+self._Cma*(alpha-self._aL0)+Cm_df*d_f

        # Functional model
        elif self._type == "functional":
            Cm = self._Cm(**kwargs)

        # Generated/imported database
        elif self._type == "database":
            Cm = self._get_database_data(2, **kwargs)

        # Fits
        elif self._type == "poly_fit":
            Cm = self._get_polynomial_data(2, **kwargs)

        return Cm


    def get_aL0(self, **kwargs):
        """Returns the zero-lift angle of attack, taking flap deflection into account.

        Parameters
        ----------
        Rey : float, optional
            Reynolds number.

        Mach : float, optional
            Mach number.

        trailing_flap_deflection : float, optional
            Trailing flap deflection in radians. Defaults to 0.

        trailing_flap_fraction : float, optional
            Trailing flap fraction of the chord length. Defaults to 0.

        Returns
        -------
        float
            Zero-lift angle of attack
        """

        # Linear airfoil model
        if self._type == "linear":

            # Get params
            d_f = kwargs.get("trailing_flap_deflection", 0.0)
            c_f = kwargs.get("trailing_flap_fraction", 0.0)

            # Calculate lift coefficient
            if np.array(d_f == 0.0).all() or np.array(c_f == 0.0).all():
                aL0 = self._aL0
            else:
                delta_a_d_f = self._get_flap_influence(c_f, d_f)
                aL0 = self._aL0-delta_a_d_f
            
            return aL0

        # Database/poly fit/functional
        # Use secant method in alpha to find a_L0
        elif self._type == "database" or self._type == "poly_fit" or self._type == "functional":

            # Remove alpha from kwargs
            kwargs.pop('alpha', None)

            # Initialize secant method
            a0 = 0.0
            CL0 = self.get_CL(alpha=a0, **kwargs)
            a0 = np.zeros_like(CL0)
            a1 = np.zeros_like(CL0)+0.01
            CL1 = self.get_CL(alpha=a1, **kwargs)
            a2 = np.zeros_like(CL1)

            # If we're outside the domain of the database, aL0 should be nan
            if a2.size == 1:
                if np.isnan(CL1):
                    a2 = np.nan
            else:
                a2[np.where(np.isnan(CL1))] = np.nan
            
            # Iterate
            np.seterr(invalid='ignore')
            not_converged = np.where(np.array(np.abs(CL1)>1e-10))[0]
            while not_converged.size>0:
                
                # Update estimate
                if a2.size == 1:
                    a2 = (a1-CL1*(a0-a1)/(CL0-CL1))
                else:
                    a2[not_converged] = (a1-CL1*(a0-a1)/(CL0-CL1))[not_converged]
                CL2 = self.get_CL(alpha=a2, **kwargs)

                # Update for next iteration
                a0 = np.copy(a1)
                CL0 = np.copy(CL1)
                a1 = np.copy(a2)
                CL1 = np.copy(CL2)

                # Check convergence
                not_converged = np.where(np.array(np.abs(CL1)>1e-10))[0]

            np.seterr()
            return a2


    def get_CLM(self, **kwargs):
        """Returns the lift slope with respect to Mach number using a forward-difference approximation.
        Simply returns 0 for a type 'linear' airfoil.

        Parameters
        ----------
        alpha : float, optional
            Angle of attack in radians. Defaults to 0.0.

        Rey : float, optional
            Reynolds number.

        Mach : float, optional
            Mach number.

        trailing_flap_deflection : float, optional
            Trailing flap deflection in radians. Defaults to 0.

        trailing_flap_fraction : float, optional
            Trailing flap fraction of the chord length. Defaults to 0.

        dx : float
            Step size for finite-difference equation. Defaults to 0.05.

        Returns
        -------
        float or ndarray
            Lift slope with respect to Mach number
        """

        # Linear model
        if self._type == "linear":
            return 0.0

        # Database
        elif self._type == "database" or self._type == "poly_fit" or self._type == "functional":

            # Check the database is dependent on Mach
            if (self._type == "database" or self._type == "poly_fit") and "Mach" not in self._dof_db_order:
                return 0.0

            # Get center Mach value
            dx = kwargs.get("dx", 0.05)
            Mach = kwargs.pop("Mach", 0.0)
            
            # Calculate forward and center points (since we'll often be at M=0 and negative M doesn't work)
            CL1 = self.get_CL(Mach=Mach+dx, **kwargs)
            CL0 = self.get_CL(Mach=Mach, **kwargs)

            return (CL1-CL0)/dx


    def get_CLRe(self, **kwargs):
        """Returns the lift slope with respect to Reynolds number using a central-difference approximation.
        Simply returns 0 for a type 'linear' airfoil.

        Parameters
        ----------
        alpha : float, optional
            Angle of attack in radians. Defaults to 0.0.

        Rey : float, optional
            Reynolds number.

        Mach : float, optional
            Mach number.

        trailing_flap_deflection : float, optional
            Trailing flap deflection in radians. Defaults to 0.

        trailing_flap_fraction : float, optional
            Trailing flap fraction of the chord length. Defaults to 0.

        dx : float
            Step size for finite-difference equation. Defaults to 1000.

        Returns
        -------
        float or ndarray
            Lift slope with respect to Reynolds number
        """

        # Linear model
        if self._type == "linear":
            return 0.0

        # Database
        elif self._type == "database" or self._type == "poly_fit" or self._type == "functional":

            # Check the database is dependent on Re
            if (self._type == "database" or self._type == "poly_fit") and "Rey" not in self._dof_db_order:
                return 0.0

            # Get center Re value
            dx = kwargs.get("dx", 1000)
            Rey = kwargs.pop("Rey", 1000000)
            
            # Calculate forward and backward points
            CL1 = self.get_CL(Rey=Rey+dx, **kwargs)
            CL0 = self.get_CL(Rey=Rey-dx, **kwargs)

            return (CL1-CL0)/(2*dx)


    def get_CLa(self, **kwargs):
        """Returns the lift slope using a central-difference approximation.

        Parameters
        ----------
        alpha : float, optional
            Angle of attack in radians. Defaults to 0.0.

        Rey : float, optional
            Reynolds number.

        Mach : float, optional
            Mach number.

        trailing_flap_deflection : float, optional
            Trailing flap deflection in radians. Defaults to 0.

        trailing_flap_fraction : float, optional
            Trailing flap fraction of the chord length. Defaults to 0.

        dx : float, optional
            Step size for finite-difference equation. Defaults to 0.001 radians.

        Returns
        -------
        float or ndarray
            Lift slope
        """

        # Linear model
        if self._type == "linear":
            return self._CLa

        # Database
        elif self._type == "database" or self._type == "poly_fit" or self._type == "functional":

            # Check the database is dependent on alpha
            if (self._type == "database" or self._type == "poly_fit") and "alpha" not in self._dof_db_order:
                return 0.0

            # Get center alpha value
            dx = kwargs.get("dx", 0.001)
            alpha = kwargs.pop("alpha", 0.0)
            
            # Calculate forward and backward points
            CL1 = self.get_CL(alpha=alpha+dx, **kwargs)
            CL0 = self.get_CL(alpha=alpha-dx, **kwargs)

            return (CL1-CL0)/(2*dx)


    def get_thickness(self, x):
        """Returns the thickness normal to the camber line at the specified x location(s).

        Parameters
        ----------
        x : float or ndarray
            x location(s) at which to get the thickness.

        Returns
        -------
        float or ndarray
            Thickness as a percentage of the chord.
        """
        if isinstance(x, float):
            return self._thickness(x).item()
        else:
            return self._thickness(x)


    def get_max_thickness(self):
        """Returns the maximum thickness of the airfoil, divided by the chord length.
        """
        return self._max_thickness


    def get_camber(self, x):
        """Returns the y coordinate(s) of the camber line at the specified x location(s).

        Parameters
        ----------
        x : float or ndarray
            x location(s) at which to get the thickness.

        Returns
        -------
        float or ndarray
            Camber line y coordinate(s) as a percentage of the chord.
        """
        if isinstance(x, float):
            return self._camber_line(x).item()
        else:
            return self._camber_line(x)


    def get_max_camber(self):
        """Returns the maximum camber of the airfoil, divided by the chord length.
        """
        return self._max_camber


    def get_outline_points(self, **kwargs):
        """Returns an array of outline points showing the geometry of the airfoil.

        Parameters
        ----------
        N : int, optional
            The number of outline points to return. This function will not always return exactly this many
            but never more. Defaults to 200.

        cluster : bool, optional
            Whether to cluster points about the leading and trailing edges. Defaults to True.

        trailing_flap_deflection : float, optional
            Trailing flap deflection in radians (positive down). Defaults to zero.

        trailing_flap_fraction : float, optional
            Trailing flap fraction of the chord. Defaults to zero.

        export : str, optional
            If specified, the outline points will be saved to a file. Defaults to no file.

        top_first : bool, optional
            The order of the coordinates when exported. Defaults to going from the trailing edge along the top and
            the around to the bottom.

        close_te : bool, optional
            Whether the top and bottom trailing edge points should be forced to be equal. Defaults to False

        plot : bool, optional
            Whether a plot of the outline points should be displayed once computed. Defaults to False.

        original_points : bool, optional
            Whether this should simply return the original inputted points. If set to True, all other kwargs relating
            to the geometry will be ignored. Defaults to False.

        Returns
        -------
        ndarray
            Outline points in airfoil coordinates.
        """

        # Check the geometry has been defined
        if self.geom_specification != "none":

            # Get kwargs
            N = kwargs.get("N", 200)
            cluster = kwargs.get("cluster", True)
            trailing_flap_deflection = kwargs.get("trailing_flap_deflection", 0.0)
            trailing_flap_fraction = kwargs.get("trailing_flap_fraction", 0.0)
            close_te = kwargs.get("close_te", False)

            if not kwargs.get("original_points", False):
            
                # Zach's method of determining NACA airfoils with deflected parabolic flaps (actually works really well for all of them...)

                # Get flap parameters
                d_f = trailing_flap_deflection
                x_f = 1.0-trailing_flap_fraction
                y_f = self._input_dict.get("trailing_flap_hinge_height", self._camber_line(x_f))
            
                # Calculate distributions
                theta_c = np.linspace(np.pi, -np.pi, N) # Loop over entire airfoil from TE->top->LE->bottom->TE
                if cluster:
                    x_c = 0.5*(1.0-np.cos(theta_c))
                else:
                    x_c = abs(np.linspace(-1.0, 1.0, N))
            
                # Calculate undeformed camber line and thickness
                y_c = self._camber_line(x_c)
                t = self._thickness(x_c)
                dyc_dx = self._camber_deriv(x_c)
            
                if d_f != 0.0 and x_f < 1.0:
                    # Determine which camber points belong to the flap
                    flap_ind = np.where(x_c>x_f)

                    # Linear flap
                    if self._trailing_flap_type == "linear":

                        # Calculate deflected camber line Eqs. (8-10) in "Geometry and Aerodynamic Performance of Parabolic..." by Hunsaker, et al. 2018
                        r = np.sqrt((y_c-y_f)*(y_c-y_f)+(x_c-x_f)*(x_c-x_f))
                        psi = np.arctan((y_c-y_f)/(x_c-x_f))
                        x_c[flap_ind] = x_f+(r*np.cos(d_f-psi))[flap_ind]
                        y_c[flap_ind] = y_f-(r*np.sin(d_f-psi))[flap_ind]
                        dyc_dx[flap_ind] = np.gradient(y_c[flap_ind], x_c[flap_ind], edge_order=2)

                    else:
                        # Parabolic flap from "Geometry and Aerodynamic Performance of Parabolic..." by Hunsaker, et al. 2018

                        # Calculate the neutral line parameters
                        l_n = np.sqrt(y_f*y_f+(1-x_f)*(1-x_f))
                        phi_n = -np.arctan2(y_f, 1-x_f)

                        # Calculate the location of the deflected trailing edge
                        tan_df = np.tan(d_f)
                        R = np.sqrt(4*tan_df*tan_df+1)+np.arcsinh(2*tan_df)/(2*tan_df)
                        E_te = 2.0*l_n/R

                        # Find E_p using secant method

                        # Constants
                        R_tan_df = R*tan_df
                        R_tan_df2 = R_tan_df*R_tan_df
                        E_0 = ((x_c-x_f)/(1-x_f)*l_n)[flap_ind]

                        # Initial guesses
                        E_p0 = E_te*E_0/l_n
                        R0 = E_p0/2*np.sqrt(E_p0**2/l_n**2*R_tan_df2+1)+l_n/(2*R_tan_df)*np.arcsinh(E_p0/l_n*R_tan_df)-E_0
                        E_p1 = E_te*E_0/l_n+0.001
                        R1 = E_p1/2*np.sqrt(E_p1**2/l_n**2*R_tan_df2+1)+l_n/(2*R_tan_df)*np.arcsinh(E_p1/l_n*R_tan_df)-E_0

                        # Suppress warnings because an error will often occur within the np.where that has no effect on computation
                        with np.errstate(invalid='ignore'):
                        
                            # Iterate
                            while (abs(R1)>1e-10).any():

                                # Update value
                                E_p2 = np.where(np.abs(R0-R1) != 0.0, E_p1-R1*(E_p0-E_p1)/(R0-R1), E_p1)

                                # Get residual
                                R2 = E_p2/2*np.sqrt(E_p2**2/l_n**2*R_tan_df2+1)+l_n/(2*R_tan_df)*np.arcsinh(E_p2/l_n*R_tan_df)-E_0

                                # Update for next iteration
                                E_p0 = E_p1
                                R0 = R1
                                E_p1 = E_p2
                                R1 = R2

                        # Store final result
                        E_p = E_p1
                        n_p = -E_p*E_p/E_te*tan_df

                        # Calculate deflected neutral line
                        x_p = x_f+E_p*np.cos(phi_n)-n_p*np.sin(phi_n)
                        y_p = y_f+E_p*np.sin(phi_n)+n_p*np.cos(phi_n)
                        y_nl = y_f*(1-(x_c-x_f)/(1.0-x_f))
                        dy_c = (y_c-y_nl)[flap_ind]

                        # Calculate deflected camber line
                        C = np.arctan(2*E_p/E_te*tan_df)
                        x_c[flap_ind] = x_p+dy_c*np.sin(C)
                        y_c[flap_ind] = y_p+dy_c*np.cos(C)

                        dyc_dx[flap_ind] = (dyc_dx[flap_ind] - 2*E_p*tan_df/E_te) / (1 + 2*E_p*tan_df/E_te*dyc_dx[flap_ind])

                # Outline points
                X = x_c-t*np.sin(np.arctan(dyc_dx))*np.sign(theta_c)
                Y = y_c+t*np.cos(np.arctan(dyc_dx))*np.sign(theta_c)

                # Trim overlapping points from linear flap deflection
                if self._trailing_flap_type == "linear" and d_f != 0.0 and x_f < 1.0:

                    # Get indices of top and bottom points
                    top_ind = theta_c > 0.0
                    bot_ind = theta_c <= 0.0

                    # Split points into top and bottom
                    X_t = X[top_ind]
                    Y_t = Y[top_ind]
                    X_b = X[bot_ind]
                    Y_b = Y[bot_ind]

                    # Find the point on the surface where the hinge breaks
                    # These equations determine the slope of the line bisecting the angle between the camber line before the flap
                    # and the camber line of the flap
                    psi_prime = np.pi+np.arctan(self._camber_deriv(x_f)) # Angle to the camber line before the flap
                    d_f_prime = -d_f+np.arctan(self._camber_deriv(x_f)) # Angle to the camber line after the flap
                    phi = 0.5*(psi_prime+d_f_prime)
                    b = np.tan(phi)
                    x_t_break, y_t_break, _ = self._get_intersection_point(x_f, y_f, b, "top")
                    x_b_break, y_b_break, _ = self._get_intersection_point(x_f, y_f, b, "bottom")

                    # Trim overlapping points off of both surfaces
                    X_b, Y_b, num_trimmed_from_bot = self._trim_surface(X_b, Y_b, "forward")
                    X_t, Y_t, num_trimmed_from_top = self._trim_surface(X_t, Y_t, "backward")

                    # Check if we need to fill anything in so the number of outline points remains the same
                    if num_trimmed_from_bot>0:
                        r_t = self._get_cart_dist(x_t_break, y_t_break, x_f, y_f)
                        X_t, Y_t = self._fill_surface(X_t, Y_t, x_f, y_f, r_t, x_t_break, num_trimmed_from_bot, "backward")

                    if num_trimmed_from_top>0:
                        r_b = self._get_cart_dist(x_b_break, y_b_break, x_f, y_f)
                        X_b, Y_b = self._fill_surface(X_b, Y_b, x_f, y_f, r_b, x_b_break, num_trimmed_from_top, "forward")

                    # Concatenate top and bottom points
                    X = np.concatenate([X_t, X_b])
                    Y = np.concatenate([Y_t, Y_b])

                # Make sure the trailing edge is sealed
                if close_te:
                    x_te = 0.5*(X[0]+X[-1])
                    y_te = 0.5*(Y[0]+Y[-1])
                    X[0] = x_te
                    Y[0] = y_te
                    X[-1] = x_te
                    Y[-1] = y_te

            # Just use the original points
            else:
                X = self._raw_outline[:,0]
                Y = self._raw_outline[:,1]

            # Plot result
            if kwargs.get("plot", False):
                plt.figure()
                plt.plot(X, Y)
                #plt.plot(x_t_break, y_t_break, 'rx')
                #plt.plot(x_b_break, y_b_break, 'rx')
                #plt.plot(x_f, y_f, 'rx')
                plt.plot(x_c, y_c, 'r--')
                plt.xlabel('x')
                plt.ylabel('y')
                plt.title(self.name)
                plt.gca().set_aspect('equal', adjustable='box')
                plt.show()

            # Concatenate x and y
            outline_points = np.concatenate([X[:,np.newaxis], Y[:,np.newaxis]], axis=1)
            if not kwargs.get("top_first", True):
                outline_points = outline_points[::-1,:]
                    
            # Save to file
            export = kwargs.get("export", None)
            if export is not None:
                np.savetxt(export, outline_points, fmt='%10.8f')

            return outline_points

        else:
            raise RuntimeError("The geometry has not been defined for airfoil {0}. Outline points cannot be generated.".format(self.name))


    def _trim_surface(self, X, Y, direction):
        # Trims any points within a region of doubling back in x

        # Loop through points
        trim_indices = []
        rev = False
        fwd = False

        # Determine direction
        indices = range(X.shape[0])
        if direction == "backward":
            indices = indices[::-1]
        
        # Loop through points
        trim_indices = []
        i_prev = 0
        for i in indices:
            
            # Skip first third of the section so we don't end up trimming off the nose
            if X[i] < 0.3:
                continue

            # Check if we've started going backwards
            if not rev and X[i] < X[i_prev]:
                rev = i_prev

            # Check if we've started going forward again
            if rev and not fwd and X[i] > X[i_prev]:
                fwd = i_prev
                break

            i_prev = i
        
        # Trim and insert
        if rev and fwd:

            # Determine where to trim
            if direction == "forward":
                trim_indices = list(range(rev, fwd+1))
            else:
                trim_indices = list(range(fwd, rev+1))

            # Trim
            X = np.delete(X, trim_indices)
            Y = np.delete(Y, trim_indices)

            return X, Y, len(trim_indices)
        
        else:
            return X, Y, 0



    def _fill_surface(self, X, Y, x_h, y_h, r, x_break, num_points, direction):
        # Fills in num_points along the arc defined by x_h, y_h, and r

        # Determine direction
        indices = range(X.shape[0])
        if direction == "backward":
            indices = indices[::-1]

        # Find the two points we need to fill in between
        for i in indices:
            
            # Skip first third of the section
            if X[i] < 0.3:
                continue

            if X[i] > x_break:
                fill_start = i-1
                fill_stop = i
                break

        # No filling needed, apparently...
        else:
            return X, Y

        # Get angles from the hinge point to the start and stop points
        theta0 = math.atan2(Y[fill_start]-y_h, X[fill_start]-x_h)
        theta1 = math.atan2(Y[fill_stop]-y_h, X[fill_stop]-x_h)

        # Find fill in points
        theta_fill = np.linspace(theta0, theta1, num_points+2)[1:-1]
        x_fill = x_h+r*np.cos(theta_fill)
        y_fill = y_h+r*np.sin(theta_fill)

        return np.insert(X, fill_stop, x_fill), np.insert(Y, fill_stop, y_fill)


    def _get_closest_point_on_surface(self, x, y, surface):
        # Finds the point on the surface of the airfoil which is closest to the given point using Golden section search
        # Returns the coordinates and the radius from the given point

        # Decide where to start in s
        if surface == "top":
            s0 = 0.0
            s3 = 0.5
        else:
            s0 = 0.5
            s3 = 1.0

        # Get Golden ratio
        R = 2.0/(1.0+math.sqrt(5.0))

        # Loop until interval converges
        while abs(s0-s3)>1e-10:
            
            # Get interior points
            diff = s3-s0
            s1 = s3-diff*R
            s2 = s0+diff*R

            # Calculate points
            x1 = self._x_outline(s1)
            y1 = self._y_outline(s1)
            x2 = self._x_outline(s2)
            y2 = self._y_outline(s2)
    
            # Calculate distances
            d1 = self._get_cart_dist(x, y, x1, y1)
            d2 = self._get_cart_dist(x, y, x2, y2)

            # Check
            if d1 < d2:
                s3 = s2
            else:
                s0 = s1

        return x1, y1, d1


    def _get_cart_dist(self, x0, y0, x1, y1):
        x_diff = x0-x1
        y_diff = y0-y1
        return math.sqrt(x_diff*x_diff+y_diff*y_diff)


    def generate_database(self, **kwargs):
        """Makes calls to Xfoil to calculate CL, CD, and Cm as a function of each given degree of freedom.

        Parameters
        ----------
        degrees_of_freedom : dict
            A dict specifying which degrees of freedom the database should perturb. Allowable degrees of 
            freedom are "alpha", "Rey", "Mach", "trailing_flap_deflection", and "trailing_flap_fraction".

            Each key should be one of these degrees of freedom. To specify a range for the degree of freedom,
            a dictionary with the following keys should be given:

                "range" : list
                    The lower and upper limits for this DOF.

                "steps" : int
                    The number of points in the range to interrogate.

                "index" : int
                    Index of the column for this degree of freedom in the database.

                "log_step" : bool, optional
                    Whether the steps in this dof should be spaced linearly (False) or logarithmically
                    (True). Defaults to False.

            If instead of perturbing a variable, you want the database to be evaluated at a constant value of
            that variable, a float should be given instead of a dictionary. That variable will then not be considered
            a degree of freedom of the database and will not appear in the database.
            
            An example is shown below:

                dofs = {
                    "alpha" : {
                        "range" : [math.radians(-15.0), math.radians(15.0)],
                        "steps" : 21,
                        "index" : 1
                    },
                    "Rey" : 2900000.0,
                    "trailing_flap_deflection" : {
                        "range" : [math.radians(-20.0), math.radians(20.0)],
                        "steps" : 1,
                        "index" : 0
                    },
                    "trailing_flap_fraction" : 0.25
                }

            The above input will run the airfoil through angles of attack from -15 to 15 degrees at a Reynolds number
            of 2900000.0 and through flap deflections from -20 to 20 degrees with a chord fraction of 25%.

            If a float, the degree of freedom is assumed to be constant at that value.
            
            If not specified, each degree of freedom defaults to the following:

                "alpha" : 0.0

                "Rey" : 1000000.0

                "Mach" : 0.0

                "trailing_flap_deflection" : 0.0

                "trailing_flap_fraction" : 0.0

            Please note that all angular degreees of freedom are in radians, rather than degrees.

            Currently, the number of steps in Reynolds number multiplied by the number of steps in Mach number
            may not exceed 12. This is due to internal limitations in Xfoil. However, due to the weak dependence
            of airfoil properties on Reynolds number, we do not expect this to be a great hinderance.

        N : int, optional
            Number of panel nodes for Xfoil to use. Defaults to 200.

        max_iter : int, optional
            Maximum iterations for Xfoil. Defaults to 100.

        x_trip : float or list, optional
            x location, non-dimensionalized by the chord length, of the boundary layer trip position. This is 
            specified for the top and bottom of the airfoil. If a float, the value is the same for the top
            and the bottom. If a list, the first list element is the top trip location and the second list element
            is the bottom trip location. Defaults to 1.0 for both.

        N_crit : float or list, optional
            Critical amplification exponent for the boundary layer in Xfoil. If a float, the value is the same for
            the top and the bottom. If a list, the first list element is for the top and the second list element is
            for the bottom. Defaults to 9.0 for both.

        update_type : bool, optional
            Whether to update the airfoil to use the newly computed database for calculations. Defaults to True.

        show_xfoil_output : bool, optional
            Display whatever Xfoil prints out. Defaults to False.

        show_xfoil_plots : bool, optional
            Display Xfoil plots. Defaults to True.

        resize_xfoil_window : float, optional
            resizes the xfoil window to screen size fraction. Xfoil defaults to 0.8 window/screen size.
            This variable defaults to None. Has no effect if show_xfoil_plots is False.

        CD_type : str, optional
            Which drag coefficient to read in. May be 'total', 'friction', or 'pressure'.
            Defaults to 'total'.

        verbose : bool, optional
            Defaults to True
        """

        # Set up lists of independent vars
        xfoil_args = {}
        self._dof_db_cols = {}
        num_total_runs = 1
        for dof, params in kwargs.pop("degrees_of_freedom", {}).items():
            if dof not in self._allowable_dofs:
                raise IOError("{0} is not an allowable DOF.".format(dof))
            vals, column_index = self._setup_ind_var(params)
            xfoil_args[dof] = vals
            num_total_runs *= len(vals)
            if column_index is not None:
                self._dof_db_cols[dof] = column_index

        # Get coefficients
        CL, CD, Cm = self.run_xfoil(**xfoil_args, **kwargs)

        # Determine the rows and cols in the database; each independent var and coefficient is a column to be iterpolated using scipy.interpolate.griddata
        self._num_dofs = len(list(self._dof_db_cols.keys()))
        num_cols = 3+self._num_dofs
        num_rows = CL.size-np.count_nonzero(np.isnan(CL))
        dof_sorted = sorted(self._dof_db_cols.items(), key=operator.itemgetter(1))
        self._dof_db_order = [x[0] for x in dof_sorted]

        # Arrange into 2D database
        self._data = np.zeros((num_rows, num_cols))
        database_row = 0
        coef_shape = CL.shape

        for i in range(coef_shape[0]):
            for j in range(coef_shape[1]):
                for k in range(coef_shape[2]):
                    for l in range(coef_shape[3]):
                        for m in range(coef_shape[4]):

                            # Check for nan
                            if np.isnan(CL[i,j,k,l,m]):
                                continue

                            # Append independent vars to database
                            for n, dof in enumerate(self._dof_db_order):
                                if dof == "alpha":
                                    ind = i
                                elif dof == "Rey":
                                    ind = j
                                elif dof == "Mach":
                                    ind = k
                                elif dof == "trailing_flap_deflection":
                                    ind = l
                                else:
                                    ind = m
                                self._data[database_row,n] = xfoil_args[dof][ind]
                        
                            # Append coefficients
                            self._data[database_row,self._num_dofs] = CL[i,j,k,l,m]
                            self._data[database_row,self._num_dofs+1] = CD[i,j,k,l,m]
                            self._data[database_row,self._num_dofs+2] = Cm[i,j,k,l,m]

                            database_row += 1

        # Sort by columns so the first column is perfectly in order
        dtype = ",".join(['i8' for i in range(num_cols)])
        for i in range(self._num_dofs,-1,-1):
            self._data = self._data[self._data[:,i].argsort(axis=0, kind='stable')] # 'stable' option is necessary to maintain ordering of columns not being actively sorted

        # Let the user know how much of the design space we actually got results for
        if kwargs.get("verbose", True):
            percent_success = round(self._data.shape[0]/num_total_runs*100, 2)
            print("\nDatabase generation complete.")
            print("Convergent results obtained from Xfoil for {0}% of the requested points.".format(percent_success))

        # Update type
        if kwargs.get("update_type", True):
            self.set_type("database")


    def _setup_ind_var(self, input_dict):
        # Sets up a range of independent variables

        # Constant value
        if isinstance(input_dict, float):
            lower = input_dict
            upper = input_dict
            N = 1
            index = None
            log_step = False

        # Range
        else:
            limits = input_dict.get("range")
            lower = limits[0]
            upper = limits[1]
            N = input_dict.get("steps")
            index = input_dict.get("index", None)
            log_step = input_dict.get("log_step", False)

            # Check that the range will actually have multiple points
            if N == 1:
                raise IOError("A range with only one step may not be specified for a degree of freedom. Please give a single float instead.")
        
        if log_step:
            return list(np.logspace(np.log10(lower), np.log10(upper), N)), index
        else:
            return list(np.linspace(lower, upper, N)), index


    def export_database(self, **kwargs):
        """Exports the database generated by generate_database().

        Parameters
        ----------
        filename : str
            File to export the database to.
        """

        filename = kwargs.get("filename")

        # Check the database is there
        if hasattr(self, "_data"):

            # Create header
            header = []

            # Add degrees of freedom
            for dof in sorted(self._dof_db_cols.items(), key=operator.itemgetter(1)):
                header.append("{:<25s}".format(dof[0]))

            # Add coefficients
            header.append("{:<25s}".format('CL'))
            header.append("{:<25s}".format('CD'))
            header.append("{:<25s}".format('Cm'))
            header = " ".join(header)

            # Export
            with open(filename, 'w') as db_file:
                np.savetxt(db_file, self._data, '%25.10E', header=header)
        else:
            raise RuntimeError("No database has been generated for airfoil {0}. Please create a database before exporting.".format(self.name))


    def import_database(self, **kwargs):
        """Imports the specified database. Please note that if you have generated your own database not
        using AirfoilDatabase, angle of attack should be stored in radians, rather than degrees.

        Parameters
        ----------
        filename : str
            File to import the database from

        update_type : bool, optional
            Whether to update the airfoil to use the newly imported database for calculations. Defaults to True.
        """

        filename = kwargs.get("filename")

        # Load data from file
        with open(filename, 'r') as db_file:
            self._data = np.loadtxt(db_file)

        # Determine the column indices
        with open(filename, 'r') as db_file:
            header = db_file.readline().strip('#')
        self._dof_db_cols = {}
        self._num_dofs = 0
        for i, col_name in enumerate(header.split()):

            # Stop once we get to coefficient columns
            if col_name == "CL":
                break

            # Add
            self._dof_db_cols[col_name] = i
            self._num_dofs += 1

        # Figure out the order of the columns in the database
        dof_sorted = sorted(self._dof_db_cols.items(), key=operator.itemgetter(1))
        self._dof_db_order = [x[0] for x in dof_sorted]

        # Update type
        if kwargs.get("update_type", True):
            self.set_type("database")


    def run_xfoil(self, **kwargs):
        """Calls Xfoil and extracts the aerodynamic coefficients at the given state.

        Parameters
        ----------
        alpha : float or list of float
            Angle(s) of attack to calculate the coefficients at in radians. Defaults to 0.0.

        Rey : float or list of float
            Reynolds number(s) to calculate the coefficients at. Defaults to 1000000.

        Mach : float or list of float
            Mach number(s) to calculate the coefficients at. Defaults to 0.0.

        trailing_flap_deflection : float or list of float
            Flap deflection(s) to calculate the coefficients at in radians. Defaults to 0.0.

        trailing_flap_fraction : float or list of float
            Flap fraction(s) to calculate the coefficients at in radians. Defaults to 0.0.

        N : int, optional
            Number of panels for Xfoil to use. Defaults to 200.

        max_iter : int, optional
            Maximum iterations for Xfoil. Defaults to 100.

        x_trip : float or list, optional
            x location, non-dimensionalized by the chord length, of the boundary layer trip position. This is 
            specified for the top and bottom of the airfoil. If a float, the value is the same for the top
            and the bottom. If a list, the first list element is the top trip location and the second list element
            is the bottom trip location. Defaults to 1.0 for both.

        xycm : list, optional
            x-y coordinates, non-dimensionalized by the chord length, of the reference point for determining the
            moment coefficient. Defaults to the quater-chord.

        N_crit : float or list, optional
            Critical amplification exponent for the boundary layer in Xfoil. If a float, the value is the same for
            the top and the bottom. If a list, the first list element is for the top and the second list element is
            for the bottom. Defaults to 9.0 for both.

        show_xfoil_output : bool, optional
            Display whatever Xfoil outputs from the command line interface. Defaults to False.

        show_xfoil_plots : bool, optional
            Display Xfoil plots. Defaults to True.
        
        resize_xfoil_window : float, optional
            resizes the xfoil window to screen size fraction. Xfoil defaults to 0.8 window/screen size.
            This variable defaults to None. Has no effect if show_xfoil_plots is False.

        CD_type : str, optional
            Which drag coefficient to read in. May be 'total', 'friction', or 'pressure'.
            Defaults to 'total'.

        verbose : bool, optional

        Returns
        -------
        CL : ndarray
            Coefficient of lift. First dimension will match the length of alpha, second will match Rey, etc.

        CD : ndarray
            Coefficient of drag. Dimensions same as CL.

        Cm : ndarray
            Moment coefficient. Dimensions same as CL.
        """
        N = kwargs.get("N", 200)
        max_iter = kwargs.get("max_iter", 100)
        verbose = kwargs.get("verbose", True)
        show_xfoil_output = kwargs.get("show_xfoil_output", False)
        x_trip = kwargs.get("x_trip", [1.0, 1.0])
        xycm = kwargs.get("xycm", [0.25, 0.0])
        if isinstance(x_trip, float):
            x_trip = [x_trip, x_trip]
        N_crit = kwargs.get("N_crit", [9.0, 9.0])
        if isinstance(N_crit, float):
            N_crit = [N_crit, N_crit]

        # Get states
        # Angle of attack
        alphas = kwargs.get("alpha", [0.0])
        if isinstance(alphas, float):
            alphas = [alphas]
        first_dim = len(alphas)
    
        # Reynolds number
        Reys = kwargs.get("Rey", [1000000.0])
        if isinstance(Reys, float):
            Reys = [Reys]
        second_dim = len(Reys)

        # Mach number
        Machs = kwargs.get("Mach", [0.0])
        if isinstance(Machs, float):
            Machs = [Machs]
        third_dim = len(Machs)

        # Flap deflections
        delta_fts = kwargs.get("trailing_flap_deflection", [0.0])
        if isinstance(delta_fts, float):
            delta_fts = [delta_fts]
        fourth_dim = len(delta_fts)

        # Flap fractions
        c_fts = kwargs.get("trailing_flap_fraction", [0.0])
        if isinstance(c_fts, float):
            c_fts = [c_fts]
        fifth_dim = len(c_fts)

        # xfoil window resizing
        resize_xfoil_window = kwargs.get('resize_xfoil_window', None)

        # Initialize coefficient arrays
        CL = np.empty((first_dim, second_dim, third_dim, fourth_dim, fifth_dim))
        CD = np.empty((first_dim, second_dim, third_dim, fourth_dim, fifth_dim))
        Cm = np.empty((first_dim, second_dim, third_dim, fourth_dim, fifth_dim))
        CL[:] = np.nan
        CD[:] = np.nan
        Cm[:] = np.nan

        # Clean up from previous iterations
        dir_list = os.listdir()
        for item in dir_list:
            if os.path.isfile(item) and ".pacc" in item or ".geom" in item:
                os.remove(item)

        # Loop through flap deflections and fractions
        if verbose:
            print("Running Xfoil...")
            print("{0:>25}{1:>25}{2:>25}".format("Percent Complete", "Flap Deflection [deg]", "Flap Fraction"))
            print(''.join(['-']*75))
        num_xfoil_runs = fourth_dim*fifth_dim
        for l, delta_ft in enumerate(delta_fts):
            for m, c_ft in enumerate(c_fts):
                
                # Clear pacc file list
                pacc_files = []
            
                # Export geometry
                outline_points = "a_{0:1.6f}_{1:1.6f}.geom".format(delta_ft, c_ft)
                #outline_points = os.path.abspath("xfoil_geom_{0:1.6f}.geom".format(delta_ft))
                self.get_outline_points(N=N, trailing_flap_deflection=delta_ft, trailing_flap_fraction=c_ft, export=outline_points, close_te=False)

                # Display update
                if verbose:
                    percent_complete = round((l*fifth_dim+m)/(num_xfoil_runs)*100)
                    print("{0:>24}%{1:>25}{2:>25}".format(percent_complete, math.degrees(delta_ft), c_ft))

                # Loop through Reynolds number
                for Re in Reys:

                    # Initialize xfoil execution
                    with sp.Popen(['xfoil'], stdin=sp.PIPE, stdout=sp.PIPE) as xfoil_process:

                        commands = []

                        # Turn off plots
                        if not kwargs.get("show_xfoil_plots", True):
                            commands += ['PLOP',
                                         'G',
                                         '']
                        elif resize_xfoil_window != None:
                            commands += ['PLOP',
                                         'W {}'.format(resize_xfoil_window),
                                         '']

                        # Read in geometry
                        commands += ['LOAD {0}'.format(outline_points),
                                     '{0}'.format(self.name)]

                        # Set panelling ratio and let Xfoil makes its own panels
                        commands += ['PPAR',
                                     'N',
                                     '{0}'.format(N),
                                     'T',
                                     '1',
                                     '',
                                     '']

                        # Set moment reference point
                        commands += ['XYCM',
                                     str(xycm[0]),
                                     str(xycm[1])]

                        # Set viscous mode
                        commands += ['OPER',
                                     'VISC',
                                     '']
                        pacc_index = 0

                        # Set boundary layer parameters
                        commands += ['VPAR',
                                     'Xtr',
                                     str(x_trip[0]),
                                     str(x_trip[1]),
                                     'NT',
                                     str(N_crit[0]),
                                     'NB',
                                     str(N_crit[1]),
                                     '',
                                     '']

                        # Loop through Mach numbers
                        for M in Machs:

                            # Polar accumulation file
                            file_id = str(np.random.randint(0, 10000))
                            pacc_file = "xfoil_results_{0}.pacc".format(file_id)
                            pacc_files.append(pacc_file)

                            # Set Mach, Reynolds number, iteration limit, and polar accumulation
                            commands += ['OPER',
                                        'RE',
                                        str(Re),
                                        'MACH',
                                        str(M),
                                        'ITER {0}'.format(max_iter),
                                        'PACC',
                                        pacc_file,
                                        '']

                            # Sweep from 0 aoa up
                            zero_ind = np.argmin(np.abs(alphas))
                            if zero_ind != len(alphas)-1:
                                for a in alphas[zero_ind:]:
                                    commands.append('ALFA {0:1.6f}'.format(math.degrees(a)))
                            
                            # Reset solver
                            commands.append('INIT')

                            # Sweep from 0 aoa down
                            if zero_ind != 0:
                                for a in alphas[zero_ind-1::-1]:
                                    commands.append('ALFA {0:1.6f}'.format(math.degrees(a)))
                            
                            # Reset solver
                            commands.append('INIT')

                            # End polar accumulation
                            commands += ['PACC {0}'.format(pacc_index),
                                         '']
                            pacc_index += 1

                        # Finish commands
                        commands += ['',
                                     'QUIT']

                        # Run Xfoil
                        xfoil_input = '\r'.join(commands).encode('utf-8')
                        response = xfoil_process.communicate(xfoil_input)

                        # Show output
                        if show_xfoil_output:
                            print(response[0].decode('utf-8'))
                            if response[1] is not None:
                                print(response[1].decode('utf-8'))

                # Clean up geometry
                os.remove(outline_points)

                # Read in files and store arrays
                for filename in pacc_files:

                    # Read in file
                    try:
                        alpha_i, CL_i, CD_i, Cm_i, Re_i, M_i = self.read_pacc_file(filename, CD_type=kwargs.get('CD_type', 'total'))
                    except FileNotFoundError:
                        warnings.warn("Couldn't find results file {0}. Usually an indication of Xfoil crashing.".format(filename))
                        continue

                    # Determine the Reynolds and Mach indices
                    j = min(range(len(Reys)), key=lambda i: abs(Reys[i]-Re_i))

                    k = min(range(len(Machs)), key=lambda i: abs(Machs[i]-M_i))

                    # Loop through alphas
                    i_true = 0
                    for i_iter, alpha in enumerate(alpha_i):

                        # Line up with our original independent alpha, as Xfoil does not output a non-converged result
                        i_true = min(range(len(alphas)), key=lambda i: abs(alphas[i]-alpha))

                        CL[i_true,j,k,l,m] = CL_i[i_iter]
                        CD[i_true,j,k,l,m] = CD_i[i_iter]
                        Cm[i_true,j,k,l,m] = Cm_i[i_iter]

                    # Interpolate missing values
                    for i, alpha in enumerate(alpha_i):
                        if np.isnan(CL[i,j,k,l,m]): # Result did not converge

                            # Mid-value
                            if i != 0 and i != len(alpha_i)-1:
                                weight = (alpha_i[i+1]-alpha)/(alpha_i[i+1]-alpha_i[i-1])
                                CL[i,j,k,l,m] = CL[i-1,j,k,l,m]*(1-weight)+CL[i+1,j,k,l,m]*weight
                                CD[i,j,k,l,m] = CD[i-1,j,k,l,m]*(1-weight)+CD[i+1,j,k,l,m]*weight
                                Cm[i,j,k,l,m] = Cm[i-1,j,k,l,m]*(1-weight)+Cm[i+1,j,k,l,m]*weight

                    # Clean up polar files
                    os.remove(filename)

        return CL, CD, Cm


    def read_pacc_file(self, filename, CD_type='total'):
        """Reads in and formats an Xfoil polar accumulation file.

        Parameters
        ----------
        filename : str
            File to read in.

        CD_type : str, optional
            Which drag coefficient to read in. May be 'total', 'friction', or 'pressure'.
            Defaults to 'total'.

        Returns
        -------
        alpha : list
            List of angles of attack from the accumulation polar.

        CL : list
            Coefficient of lift at each alpha.

        CD : list
            Coefficient of drag at each alpha.
        
        Cm : list
            Moment coefficient at each alpha.

        Re : float
            Reynolds number for the polar.

        M : float
            Mach number for the polar.
        """

        # Open file
        with open(filename, 'r') as file_handle:

            # Read in line by line
            lines = []
            line = file_handle.readline()
            while line:
                lines.append(line)
                line = file_handle.readline()

        # Find Mach and Reynolds number
        mr_line = lines[8].split()
        M = float(mr_line[2])
        Re = float(''.join(mr_line[5:8]))

        # Collect alpha and coefficients
        alpha = []
        CL = []
        CD = []
        Cm = []
        for line in lines[12:]:
            split_line = line.split()
            alpha.append(math.radians(float(split_line[0])))
            CL.append(float(split_line[1]))
            if CD_type=='total':
                CD.append(float(split_line[2]))
            elif CD_type=='pressure':
                CD.append(float(split_line[3]))
            elif CD_type=='friction':
                CD.append(float(split_line[2])-float(split_line[3]))
            Cm.append(float(split_line[4]))

        # Sort in alpha
        sorted_indices = np.argsort(alpha)
        alpha = np.array(alpha)[sorted_indices]
        CL = np.array(CL)[sorted_indices]
        CD = np.array(CD)[sorted_indices]
        Cm = np.array(Cm)[sorted_indices]

        return alpha, CL, CD, Cm, Re, M


    def _create_filled_database(self):
        # Fills in missing values in the database to make it a consistent d-dimensional array.

        # Initialize storage for independent vars
        self._dof_filled = []
        for i in range(self._num_dofs):
            self._dof_filled.append([])

        # Gather independent vars
        for i in range(self._data.shape[0]):
            for j in range(self._num_dofs):

                # Check if the value is already there
                if not self._data[i,j] in self._dof_filled[j]:
                    self._dof_filled[j].append(self._data[i,j])

        # Sort independent vars and initialize filled array
        shape = []
        for dof in self._dof_filled:
            dof.sort()
            shape.append(len(dof))
        filled_CL = np.zeros(tuple(shape))
        filled_CD = np.zeros(tuple(shape))
        filled_Cm = np.zeros(tuple(shape))
        N = filled_CL.size

        # Create grid of independent vars
        raw_grid = np.meshgrid(*self._dof_filled)

        # Parse independent vars for passing to getters
        params = {}
        for i, dof in enumerate(self._dof_db_order):
            params[dof] = raw_grid[i].flatten()

        # Fill in values
        filled_CL_view = filled_CL.reshape(N)
        filled_CL_view[:] = self.get_CL(**params)
        filled_CD_view = filled_CD.reshape(N)
        filled_CD_view[:] = self.get_CD(**params)
        filled_Cm_view = filled_Cm.reshape(N)
        filled_Cm_view[:] = self.get_Cm(**params)

        # Huh. Turns out I don't actually need this, so I'm not going to bother developing it further. But I'll keep it here in case it becomes useful


    def generate_polynomial_fit(self, **kwargs):
        """Generates a set of multivariable polynomials using least-squares regression to approximate the database.
        Note: This airfoil must have a database already for fits to be created.

        Parameters
        ----------
        CL_degrees : dict or str, optional
            If dict, order of fit polynomial for the coefficient of lift for each degree of freedom,
            formatted as

            {
                "<DOF1_NAME>" : <ORDER>,
                "<DOF2_NAME>" : <ORDER>,
                ...
            }

            Orders must be integers. Defaults to 1 for any not specified.

            Can also be specified as "auto". In this case, the fit degrees will be determined automatically
            using the method described in Morelli, "Global Nonlinear Aerodynamic Modeling Using
            Multivariate Orthogonal Functions," Journal of Aircraft, 1995. The algorithm tried to minimize
            the RMS error between the data and the prediction while also minimizing the degree of the fit. 
            This will automatically determine the fit order for each degree of freedom.

        CD_degrees : dict, optional
            Same as CL_degrees.

        Cm_degrees : dict, optional
            Same as CL_degrees.

        CL_kwargs : dict, optional
            keyword arguments sent to the CL polynomial fit function
            
            When CL_degrees is specified as "auto" then CL_kwargs can be
            
                "max_order" : int, optional
                    gives the max order of polynomial for any one of the independent varialbes
                    to try. Defaults to 6.
                "tol" : float, optional
                    Gives the cut-off value for any polynomial coefficient to not be included
                    in the final results. If a coefficient has an absolute value below tol,
                    it won't be included. Defaults to 1e-12.
                "sigma" : float, optional
                    value used to determine the trade off between how good of a fit to perform
                    and how many terms to keep. Defaults to None, which causes the function to
                    calculate sigma automatically using the mean squared of the difference of
                    the independent variable values with respect to the mean independent
                    variable value of the dataset
                "sigma_multiplier" : float, optional
                    term multiplied onto sigma to change it's value. Allows using a multiple
                    of the automatically determined sigma value. Defaults to 1.
            
            Otherwise CL_kwargs could be
            
                "interaction" : boolean, optional
                    value with default set to True. This variable determines whether or not
                    interaction terms are included in the fit function. If set to True,
                    interaction terms up the max order for each independent variable are
                    included, i.e. if Nvec = [3,2] then the highest interaction term included
                    is x_1^3*x_2^2. Specific interaction terms can be omitted using the
                    constraints input
                "sym" : list, optional
                    Defaults to an empty list. If used, the length should be V and each element
                    should contain a boolean, True or False. The ith element determines if the
                    ith independent variable is symmetric either even or odd, which is
                    determined by the order given in Nvec. This will also remove the
                    cooresponding interaction terms if they are enabled.
                "sym_same" : list, optional
                    Defaults as an empty list. If used, the entries in the list should be
                    tuples with two integers. The integers represent the independent variables
                    that the "same" symmetry condition will be applied. The "same" symmetry
                    ensures all interaction terms with exponents of the two independent
                    variables that are either odd-odd or even-even to be forced to zero
                "sym_diff" : = list, optional
                    Defaults as an empty list. Similar to "sym_same" except it enforces the
                    "diff" symmetry condition which ensures all interaction terms with exponents
                    of the two independent variables that are either odd-even or even-odd to be
                    forced to zero
                "zeroConstraints" : list, optional
                    Defaults as an empty list. Entries in the list contain integer tuples of
                    length V. The integer values represent the powers of the independent variables
                    whose coefficient will be forced to 0 before the best fit calculations are
                    performed, allowing the user to omit specific interaction terms or regular
                    polynomial terms
                "constraints" : list, optional
                    Defaults to an empty list. Entries in the list contain tuples of length 2.
                    The first entry is a list of integers that represent the powers of the
                    independent variables whose coefficient will then be forced to be equal to the
                    second entry in the tuple, which should be a float.
                "percent" : boolean, optional
                    Default set to False. When set to True the least squares is performed on the
                    percent error squared. This option should not be used if y contains any zero
                    or near zero values, as this might cause a divide by zero error.
                "weighting" : function, optional
                    Defaults to None. If given, weighting should be a function that takes as
                    arguments x, y, and p where x and y are the independent and dependent
                    variables defined above and p is the index representing a certain data point.
                    weighting should return a 'weighting factor' that determines how important
                    that datapoint is. Returning a '1' weights the datapoint normally.
        
        CD_kwargs : dict, optional
            Same as CL_kwargs

        Cm_kwargs : dict, optional
            Same as CL_kwargs

        update_type : bool, optional
            Whether to update the airfoil to use the newly computed polynomial fits for calculations. Defaults to True.

        verbose : bool, optional
        """

        # Check for database
        if not hasattr(self, "_data"):
            raise RuntimeError("No database found! Please generate or import a database before trying to create polynomial fits.")

        # Determine what the maximum fit order is for autoPolyFit
        CL_degrees = kwargs.pop("CL_degrees", {})
        CD_degrees = kwargs.pop("CD_degrees", {})
        Cm_degrees = kwargs.pop("Cm_degrees", {})
        CL_kwargs  = kwargs.pop("CL_kwargs", {})
        CD_kwargs  = kwargs.pop("CD_kwargs", {})
        Cm_kwargs  = kwargs.pop("Cm_kwargs", {})
        verbose = kwargs.pop('verbose', True)
        if verbose: print('Generating Polynomial Fits for airfoil {}'.format(self.name))
        CL_kwargs['verbose'] = verbose
        CD_kwargs['verbose'] = verbose
        Cm_kwargs['verbose'] = verbose

        # CL
        if verbose: print('Performing CL curve fit')
        if CL_degrees=="auto":
            self._CL_poly_coefs, self._CL_degrees, self._CLfit_R2 = autoPolyFit(self._data[:,:self._num_dofs], self._data[:, self._num_dofs], **CL_kwargs)

        elif isinstance(CL_degrees, dict):

            # Sort fit degrees
            self._CL_degrees = []
            for dof in self._dof_db_order:
                self._CL_degrees.append(CL_degrees.get(dof, 1))

            # Generate
            self._CL_poly_coefs, self._CLfit_R2 = multivariablePolynomialFit(self._CL_degrees, self._data[:,:self._num_dofs], self._data[:, self._num_dofs], **CL_kwargs)

        else:
            raise IOError("Fit degree specification must be 'auto' or type(dict). Got {0} type {1}.".format(CL_degrees, type(CL_degrees)))
        
        # CD
        if verbose: print('Performing CD curve fit')
        if CD_degrees=="auto":
            self._CD_poly_coefs, self._CD_degrees, self._CDfit_R2 = autoPolyFit(self._data[:,:self._num_dofs], self._data[:,self._num_dofs+1], **CD_kwargs)

        elif isinstance(CL_degrees, dict):

            # Sort fit degrees
            self._CD_degrees = []
            for dof in self._dof_db_order:
                self._CD_degrees.append(CD_degrees.get(dof, 1))

            # Generate
            self._CD_poly_coefs, self._CDfit_R2 = multivariablePolynomialFit(self._CD_degrees, self._data[:,:self._num_dofs], self._data[:,self._num_dofs+1], **CD_kwargs)

        else:
            raise IOError("Fit degree specification must be 'auto' or type(dict). Got {0} type {1}.".format(CL_degrees, type(CL_degrees)))
        
        # Cm
        if verbose: print('Performing Cm curve fit')
        if Cm_degrees=="auto":
            self._Cm_poly_coefs, self._Cm_degrees, self._Cmfit_R2 = autoPolyFit(self._data[:,:self._num_dofs], self._data[:,self._num_dofs+2], **Cm_kwargs)

        elif isinstance(Cm_degrees, dict):

            # Sort fit degrees
            self._Cm_degrees = []
            for dof in self._dof_db_order:
                self._Cm_degrees.append(Cm_degrees.get(dof, 1))

            # Generate polynomial fit
            self._Cm_poly_coefs, self._Cmfit_R2 = multivariablePolynomialFit(self._Cm_degrees, self._data[:,:self._num_dofs], self._data[:, self._num_dofs+2], **Cm_kwargs)

        # Store limits
        self._dof_limits = []
        for i in range(self._num_dofs):
            self._dof_limits.append([np.min(self._data[:,i]), np.max(self._data[:,i])])

        # Update type
        if kwargs.get("update_type", True):
            self.set_type("poly_fit")
        
        self._CLfit_RMS, self._CLfit_RMSN = multivariableRMS(self._data[:,:self._num_dofs], self._data[:,self._num_dofs], self._CL_poly_coefs, self._CL_degrees, verbose=verbose)
        self._CDfit_RMS, self._CDfit_RMSN = multivariableRMS(self._data[:,:self._num_dofs], self._data[:,self._num_dofs+1], self._CD_poly_coefs, self._CD_degrees, verbose=verbose)
        self._Cmfit_RMS, self._Cmfit_RMSN = multivariableRMS(self._data[:,:self._num_dofs], self._data[:,self._num_dofs+2], self._Cm_poly_coefs, self._Cm_degrees, verbose=verbose)
        
        if verbose:
            print('\nCL fits\n'+'='*20)
            print('R^2 : {}'.format(self._CLfit_R2))
            print('RMS : {}'.format(self._CLfit_RMS))
            print('RMSN: {}\n'.format(self._CLfit_RMSN))
            print('CD fits\n'+'='*20)
            print('R^2 : {}'.format(self._CDfit_R2))
            print('RMS : {}'.format(self._CDfit_RMS))
            print('RMSN: {}\n'.format(self._CDfit_RMSN))
            print('Cm fits\n'+'='*20)
            print('R^2 : {}'.format(self._Cmfit_R2))
            print('RMS : {}'.format(self._Cmfit_RMS))
            print('RMSN: {}\n'.format(self._Cmfit_RMSN))


    def export_polynomial_fits(self, **kwargs):
        """Save the polynomial fit to a JSON object.

        Parameters
        ----------
        filename : str
            JSON object to write polynomial fit data to.

        write_limits : bool, optional
            Whether to limit the polynomial fits based on the original range of data.
            Defaults to True.

        """

        # Get filename
        filename = kwargs.get("filename")

        # Create export dictionary
        export = {}
        export["tag"] = "Polynomial fit to database for {0} airfoil.".format(self.name)
        export["degrees_of_freedom"] = self._dof_db_order
        if kwargs.get("write_limits"):
            export["limits"] = self._dof_limits
        export["defaults"] = self._dof_defaults
        export["fit_degrees"] = {}
        export["fit_degrees"]["CL"] = self._CL_degrees
        export["fit_degrees"]["CD"] = self._CD_degrees
        export["fit_degrees"]["Cm"] = self._Cm_degrees
        export['fit_error'] = {}
        export['fit_error']['CL'] = {}
        export['fit_error']['CL']['R^2'] = self._CLfit_R2
        export['fit_error']['CL']['RMS'] = self._CLfit_RMS
        export['fit_error']['CL']['RMSN'] = self._CLfit_RMSN
        export['fit_error']['CD'] = {}
        export['fit_error']['CD']['R^2'] = self._CDfit_R2
        export['fit_error']['CD']['RMS'] = self._CDfit_RMS
        export['fit_error']['CD']['RMSN'] = self._CDfit_RMSN
        export['fit_error']['Cm'] = {}
        export['fit_error']['Cm']['R^2'] = self._Cmfit_R2
        export['fit_error']['Cm']['RMS'] = self._Cmfit_RMS
        export['fit_error']['Cm']['RMSN'] = self._Cmfit_RMSN
        export["fit_coefs"] = {}
        export["fit_coefs"]["CL"] = list(self._CL_poly_coefs)
        export["fit_coefs"]["CD"] = list(self._CD_poly_coefs)
        export["fit_coefs"]["Cm"] = list(self._Cm_poly_coefs)

        # Export data
        with open(filename, 'w') as export_file_handle:
            json.dump(export, export_file_handle, indent=4)


    def import_polynomial_fits(self, **kwargs):
        """Read in polynomial fit data from a JSON object.

        Parameters
        ----------
        filename : str
            JSON object to read polynomial fit data from.

        update_type : bool, optional
            Whether to update the airfoil to use the newly imported polynomial fits for calculations. Defaults to True.

        """

        # Get filename
        filename = kwargs.get("filename")

        # Read in data
        with open(filename, 'r') as import_file_handle:
            input_dict = json.load(import_file_handle)

        # Parse input dict
        self._dof_db_order = input_dict["degrees_of_freedom"]
        self._num_dofs = len(self._dof_db_order)
        self._dof_limits = input_dict.get("limits", [[-np.inf, np.inf]]*self._num_dofs)
        self._CL_degrees = input_dict["fit_degrees"]["CL"]
        self._CD_degrees = input_dict["fit_degrees"]["CD"]
        self._Cm_degrees = input_dict["fit_degrees"]["Cm"]
        self._CL_poly_coefs = np.array(input_dict["fit_coefs"]["CL"])
        self._CD_poly_coefs = np.array(input_dict["fit_coefs"]["CD"])
        self._Cm_poly_coefs = np.array(input_dict["fit_coefs"]["Cm"])
        if isinstance(input_dict.get('fit_error', None), dict):
            if isinstance(input_dict['fit_error'].get('CL', None), dict):
                self._CLfit_R2 = input_dict['fit_error']['CL'].get('R^2', None)
                self._CLfit_RMS = input_dict['fit_error']['CL'].get('RMS', None)
                self._CLfit_RMSN = input_dict['fit_error']['CL'].get('RMSN', None)
            if isinstance(input_dict['fit_error'].get('CD', None), dict):
                self._CDfit_R2 = input_dict['fit_error']['CD'].get('R^2', None)
                self._CDfit_RMS = input_dict['fit_error']['CD'].get('RMS', None)
                self._CDfit_RMSN = input_dict['fit_error']['CD'].get('RMSN', None)
            if isinstance(input_dict['fit_error'].get('Cm', None), dict):
                self._Cmfit_R2 = input_dict['fit_error']['Cm'].get('R^2', None)
                self._Cmfit_RMS = input_dict['fit_error']['Cm'].get('RMS', None)
                self._Cmfit_RMSN = input_dict['fit_error']['Cm'].get('RMSN', None)

        # Update type
        if kwargs.get("update_type", True):
            self.set_type("poly_fit")


    def _get_polynomial_data(self, coef_index, **kwargs):
        # Determines the value of the given coefficient from the polynomial fit
        # coef_index | coef
        #   0        |  CL
        #   1        |  CD
        #   2        |  Cm

        # Stack up independent vars
        ind_vars = []
        N = 1
        for dof in self._dof_db_order:

            # Get independent variable values
            ind_var = kwargs.get(dof, self._dof_defaults[dof])
            
            # Check for array of size 1
            if not np.isscalar(ind_var) and len(ind_var) == 1:
                ind_vars.append(np.asscalar(ind_var))
            else:
                ind_vars.append(ind_var)

            # Update max size
            if not np.isscalar(ind_var):
                N = max(N, len(ind_var))

        # Fill any values
        if N > 1:
            for i, ind_var in enumerate(ind_vars):
                if np.isscalar(ind_var):
                    ind_vars[i] = np.full(N, ind_var)

        # Setup input matrix
        if N == 1:
            x = np.array(ind_vars)[:,np.newaxis]
        else:
            x = np.array(ind_vars)

        # Get data
        if coef_index == 0:
            coef = multivariablePolynomialFunction(self._CL_poly_coefs, self._CL_degrees, x)
        elif coef_index == 1:
            coef = multivariablePolynomialFunction(self._CD_poly_coefs, self._CD_degrees, x)
        elif coef_index == 2:
            coef = multivariablePolynomialFunction(self._Cm_poly_coefs, self._Cm_degrees, x)

        # Check limits
        if self._raise_poly_bounds_error:
            for i in range(N):
                for j in range(self._num_dofs):
                    if x[j,i] > self._dof_limits[j][1] or x[j,i] < self._dof_limits[j][0]:
                        coef[i] = np.nan

            # Check for going out of bounds
            if np.isnan(coef).any():
                raise PolyFitBoundsError(self.name, np.argwhere(np.isnan(coef)).flatten(), kwargs)

        return coef


    def generate_linear_model(self, **kwargs):
        """Creates a linearized model of the airfoil coefficients in alpha. Pulls from the
        database or poly_fit information to generate this model. Cannot be used on a type
        "linear" airfoil.
        
        linear_limits : list, optional
            Limits in alpha for which the behavior of the airfoil can be considered linear. If not
            given, the user will be prompted to graphically select this region. Given in degrees.

        update_type : bool, optional
            Whether to change the type of the airfoil to "linear" once the model is determined.
            Defaults to True.

        plot_model : bool, optional
            Whether to display a polar of the linear region determined. Defaults to False.

        Rey : float, optional
            Reynolds number at which to evaluate the model.

        Mach : float, optional
            Mach number at which to evaluate the model.
        """

        # Check for correct type
        if self._type == "linear":
            raise RuntimeError("generate_linear_model() cannot be called on a 'linear' type airfoil.")

        # Determine range of alpha
        linear_limits = kwargs.get("linear_limits", None)
        if linear_limits is None:
            alpha = np.linspace(-20.0, 20.0, 30)
        else:
            alpha = np.linspace(linear_limits[0], linear_limits[1], 30)

        # Generate dataset
        CL = np.zeros(30)
        Cm = np.zeros(30)
        CD = np.zeros(30)
        for i, a in enumerate(alpha):
            try:
                CL[i] = self.get_CL(alpha=np.radians(a), **kwargs)
                Cm[i] = self.get_Cm(alpha=np.radians(a), **kwargs)
                CD[i] = self.get_CD(alpha=np.radians(a), **kwargs)
            except DatabaseBoundsError:
                continue

        # Plot dataset for user to select linear region
        if linear_limits is None:

            # Define picker
            self._selected_ind = []
            def on_pick(event):

                # Get index
                i = int(event.ind[0])
                self._selected_ind.append(i)

                # Add x
                fig.axes[0].plot(alpha[i], CL[i], 'rx')

                # Check if we have enough
                if len(self._selected_ind) >= 2:
                    plt.close(fig)

            # Display
            plt.ion()
            fig, ax = plt.subplots()
            ax.plot(alpha, CL, 'bo', picker=3)
            plt.xlabel("Alpha [deg]")
            plt.ylabel("Lift Coefficient")
            plt.title("To generate a linear model, select the lower and upper limits of the linear region.")
            fig.canvas.mpl_connect('pick_event', on_pick)
            plt.show(block=True)
            plt.ioff()

            # Trim arrays
            self._selected_ind = sorted(self._selected_ind)
            if len(self._selected_ind) != 2:
                raise RuntimeError("No points were selected to determine the linear region.")
            alpha = np.radians(alpha[self._selected_ind[0]:self._selected_ind[1]+1])
            CL = CL[self._selected_ind[0]:self._selected_ind[1]+1]
            Cm = Cm[self._selected_ind[0]:self._selected_ind[1]+1]
            CD = CD[self._selected_ind[0]:self._selected_ind[1]+1]

        else:
            alpha = np.radians(alpha)

        # Get CL model
        coef_array = np.polyfit(alpha, CL, 1)
        self._aL0 = coef_array[1]
        self._CLa = coef_array[0]
        self._CL_max = np.max(CL)

        # Get Cm model
        coef_array = np.polyfit(alpha, Cm, 1)
        self._Cma = coef_array[0]
        self._CmL0 = self._Cma*self._aL0+coef_array[1]

        # Get CD model
        coef_array = np.polyfit(CL, CD, 2)
        self._CD0 = coef_array[2]
        self._CD1 = coef_array[1]
        self._CD2 = coef_array[0]

        # Plot model within linear region
        if kwargs.get("plot_model", False):

            # CL
            plt.close('all')
            fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3)
            fig.suptitle("Linear Model for {0}".format(self.name))
            ax0.set_title("CL")
            ax0.plot(alpha, CL, "gx", label="Data")
            ax0.plot(alpha, alpha*self._CLa+self._aL0, "g-", label="Model")
            ax0.legend()
            ax0.set_xlabel("Angle of Attack [rad]")
            ax0.set_ylabel("Lift Coefficient")

            # Cm
            ax1.set_title("Cm")
            ax1.plot(alpha, Cm, 'bx', label="Data")
            ax1.plot(alpha, (alpha-self._aL0)*self._Cma+self._CmL0, "b-", label="Model")
            ax1.legend()
            ax1.set_xlabel("Angle of Attack [rad]")
            ax1.set_ylabel("Moment Coefficient")

            # CD
            ax2.set_title("CD")
            ax2.plot(CL, CD, 'rx', label="Data")
            ax2.plot(CL, self._CD0+self._CD1*CL+self._CD2*CL*CL, 'r-', label="Model")
            ax2.legend()
            ax2.set_xlabel("Lift Coefficient")
            ax2.set_ylabel("Drag Coefficient")
            plt.show()

        # Update type
        if kwargs.get("update_type", True):
            self._type = "linear"


    def export_linear_model(self, **kwargs):
        """Exports the linear coefficients used to predict the behavior of the airfoil.

        Parameters
        ----------
        filename : str
            JSON file to export the model data to.
        """

        # Check there is a model to export
        if not hasattr(self, "_aL0"):
            raise RuntimeError("A linear model for {0} could not be exported because it has none.".format(self.name))

        # Parse coefficients
        model_dict = {
            "aL0" : self._aL0,
            "CLa" : self._CLa,
            "CmL0" : self._CmL0,
            "Cma" : self._Cma,
            "CD0" : self._CD0,
            "CD1" : self._CD1,
            "CD2" : self._CD2,
            "CL_max" : self._CL_max
        }

        # Export
        filename = kwargs.get("filename")
        with open(filename, 'w') as export_file_handle:
            json.dump(model_dict, export_file_handle, indent=4)
