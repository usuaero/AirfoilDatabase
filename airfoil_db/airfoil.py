"""A class defining an airfoil."""

import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import scipy.signal as sig
import subprocess as sp
import json
import copy
import os
import operator
from .poly_fits import multivariablePolynomialFit, multivariablePolynomialFunction, autoPolyFit
import io
import sys
import warnings


class Airfoil:
    """A class defining an airfoil.

    Parameters
    ----------
    name : str
        Name of the airfoil.

    airfoil_input : dict or str
        Dictionary or path to JSON object describing the airfoil.

    verbose : bool
    """

    def __init__(self, name, airfoil_input, verbose=False):
        
        self.name = name
        self._load_params(airfoil_input)
        self._verbose = verbose

        # Load flaps
        self._load_flaps()

        # Store undeformed outlines
        self._initialize_geometry()

        # Specify database DOF parameters
        self._allowable_dofs = ["alpha", "Rey", "Mach", "trailing_flap_deflection", "trailing_flap_fraction"]
        self._dof_defaults = {
            "alpha" : 0.0,
            "Rey" : 100000.0,
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


    def set_verbosity(self, verbosity):
        """Sets the verbosity of the airfoil."""
        self._verbose = verbosity


    def set_type(self, database_type):
        """Determines how the aerodynamic coefficients will be calculated.

        Parameters
        ----------
        database_type: str
            "linear", "database", or "poly_fit". Airfoil will automatically
            check if it has the necessary to perform the given type of 
            computation and throw a warning if it does not.

        """

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

        # Check for polynomial fits
        if database_type == "poly_fit":
            if not hasattr(self, "_CL_poly_coefs"):
                raise RuntimeWarning("Airfoil {0} does not have a set of polynomial fits. Reverting to type '{1}' for computations.".format(self.name, self._type))
            else:
                self._type = database_type


    def _load_params(self, airfoil_input):

        # Load input dict
        if isinstance(airfoil_input, str):
            # Load JSON object
            with open(airfoil_input) as json_handle:
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
            self._CmL0 = self._input_dict.get("CmL0", 0.0)
            self._Cma = self._input_dict.get("Cma", 0.0)
            self._CD0 = self._input_dict.get("CD0", 0.0)
            self._CD1 = self._input_dict.get("CD1", 0.0)
            self._CD2 = self._input_dict.get("CD2", 0.0)
            self._CL_max = self._input_dict.get("CL_max", np.inf)

            self._CLM = self._input_dict.get("CLM", 0.0)
            self._CLRe = self._input_dict.get("CLRe", 0.0)


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
            self._naca_des = NACA

            self._calc_geometry_from_NACA()

        # No geometry given
        else:
            self.geom_specification = "none"
            return

    
    def _calc_geometry_from_NACA(self):
        # Creates thickness, camber, camber derivative, and outline splines base on the NACA equations

        # 4-digit series
        if len(self._naca_des) == 4:

            # Stores the camber and thickness getters based on the NACA designation of the airfoil
            self._m = float(self._naca_des[0])/100
            self._p = float(self._naca_des[1])/10
            self._t = float(self._naca_des[2:])/100

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
                    return np.where(x<self._p, 2*self._m/(self._p*self._p)*(self._p-x), 2*self._m/(1-self._p*self._p)*(self._p-x))

            self._camber_deriv = camber_deriv

            # Thickness
            def thickness(x):
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
        if self._verbose: print("Calculating camber line...")

        # Find the pseudo leading and trailing edge positions
        te = (self._raw_outline[0]+self._raw_outline[-1])*0.5
        le_ind = np.argmin(self._raw_outline[:,0])
        le = self._raw_outline[le_ind]

        # Translate to origin, rotate, and scale
        self._normalize_points(self._raw_outline, le, te)
        le_ind = np.argmin(np.abs(self._raw_outline[:,0]))
        le = [0.0, 0.0]
        te = (self._raw_outline[0]+self._raw_outline[-1])*0.5

        # Create splines defining the outline
        self._x_outline, self._y_outline = self._create_splines_of_s(self._raw_outline)

        # Camber line estimate parameters
        num_camber_points = self._N_orig//5 # As recommended by Cody Cummings to avoid discontinuities in dyc/dx
        camber_deriv_edge_order = 2
        self._cosine_cluster = False
        le_offset = 0.0001

        # Use cosine clustering to space points along the camber line
        if self._cosine_cluster:
            theta = np.linspace(-np.pi, np.pi, num_camber_points)
            x_c = 0.5*(1-np.cos(theta))*(te[0]-le[0])
        else:
            x_c = np.linspace(le[0]+le_offset, te[0], num_camber_points)


        # Get initial estimate for the camber line and thickness distributions
        y_t = np.interp(x_c, self._raw_outline[le_ind::-1,0], self._raw_outline[le_ind::-1, 1])
        y_b = np.interp(x_c, self._raw_outline[le_ind:,0], self._raw_outline[le_ind:, 1])
        y_c = 0.5*(y_t+y_b)

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
        x_c_new= np.zeros(num_camber_points)
        y_c_new= np.zeros(num_camber_points)

        # Iterate through camber line estimates
        camber_error = 1

        # Check for symmetric airfoil
        if np.allclose(y_c, 0.0):
            camber_error = 0.0

        # Iterate until convergence
        while camber_error > 1e-10:

            # Determine camber line slope
            dyc_dx = np.gradient(y_c, x_c, edge_order=camber_deriv_edge_order)
            
            # Determine slope of lines perpendicular to the camber
            b = -1.0/dyc_dx

            # Loop through points on the camber line to find where their normal intersects the outline
            for i in range(num_camber_points):
                xc = x_c[i]
                yc = y_c[i]
                bi = b[i]

                x_t[i], y_t[i], _ = self._get_intersection_point(xc, yc, bi, "top")
                x_b[i], y_b[i], _ = self._get_intersection_point(xc, yc, bi, "bottom")

            # Calculate new camber line points
            x_c_new = 0.5*(x_t+x_b)
            y_c_new = 0.5*(y_t+y_b)

            # Plot new and old estimate
            if False:
                plt.figure()
                plt.plot(self._raw_outline[:,0], self._raw_outline[:,1], 'b-', label='Outline Data')
                plt.plot(x_c, y_c, 'r--', label='Old Camber Line Estimate')
                plt.plot(x_c_new, y_c_new, 'g--', label='New Camber Line Estimate')
                plt.legend()
                plt.gca().set_aspect('equal', adjustable='box')
                plt.show()

            # Approximate error
            x_diff = x_c_new-x_c
            y_diff = y_c_new-y_c
            camber_error = np.max(np.sqrt(x_diff*x_diff+y_diff*y_diff))
            if self._verbose: print("Camber error: {0}".format(camber_error))

            # Sort, just in case things got messed up
            sorted_ind = np.argsort(x_c_new)
            x_c_new = x_c_new[sorted_ind]
            y_c_new = y_c_new[sorted_ind]

            # Update for next iteration
            x_c = x_c_new
            y_c = y_c_new

        # Calculate where the camber line intersects the outline to find the leading edge
        dyc_dx = np.gradient(y_c, x_c, edge_order=2)
        b = dyc_dx[0]
        x_le, y_le, self._s_le = self._get_intersection_point(x_c[0], y_c[0], b, "leading_edge")
        le = np.array([x_le, y_le])
        x_c = np.insert(x_c, 0, le[0])
        y_c = np.insert(y_c, 0, le[1])
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

        # Store camber
        self._camber_line = interp.UnivariateSpline(camber_points[:,0], camber_points[:,1], k=5, s=1e-10)

        # Calculate thickness
        y_c = self._camber_line(x_space)
        dyc_dx= np.gradient(y_c, x_space)
        b = -1.0/dyc_dx

        # Find points on the surface to determine the thickness
        x_t_t = np.zeros(x_space.shape)
        x_b_t = np.zeros(x_space.shape)
        y_t_t = np.zeros(x_space.shape)
        y_b_t = np.zeros(x_space.shape)
        for i, xc in enumerate(x_space):
            if i == 0: continue # This point intersects the outline by definition and so will have zero thickness
            yc = y_c[i]
            bi = b[i]

            x_t_t[i], y_t_t[i],_ = self._get_intersection_point(xc, yc, bi, "top")
            x_b_t[i], y_b_t[i],_ = self._get_intersection_point(xc, yc, bi, "bottom")

        t = 0.5*np.sqrt((x_t_t-x_b_t)*(x_t_t-x_b_t)+(y_t_t-y_b_t)*(y_t_t-y_b_t))
        self._thickness = interp.UnivariateSpline(x_space, t, k=5, s=1e-10)

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
            s0 = 0.25
            s1 = 0.26
        elif surface == "bottom":
            s0 = 0.75
            s1 = 0.76
        elif surface == "leading_edge":
            s0 = 0.50
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
            if d1 > 0.2: # Apply some relaxation when we're far away
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
        return (-b*x0+y0+b*x-y)/np.sqrt(b*b+1)


    def check_against_NACA(self, naca_des):
        """Checks the error in the camber and thickness against that predicted by the NACA equations. This is recommended as a check for the user
        if unusual geometries are being imported.

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


    def get_CL(self, **kwargs):
        """Returns the coefficient of lift. Note: all parameters can be given as numpy arrays, in which case a numpy array of the coefficient will be returned.
        To do this, all parameter arrays must have only one dimension and must have the same length.

        Parameters
        ----------
        alpha : float, optional
            Angle of attack in radians. Defaults to 0.0.

        Rey : float, optional
            Reynolds number. Defaults to 100000.

        Mach : float, optional
            Mach number. Defaults to 0.

        trailing_flap_deflection : float, optional
            Trailing flap deflection in radians. Defaults to 0.

        trailing_flap_fraction : float, optional
            Trailing flap fraction of the chord length. Defaults to 0.

        trailing_flap_efficiency : float, optional
            Trailing flap efficiency. Defaults to 1.0.

        Returns
        -------
        float or ndarray
            Lift coefficient
        """
        
        # Linearized model
        if self._type == "linear":

            # Get params
            alpha = kwargs.get("alpha", 0.0)
            trailing_flap = kwargs.get("trailing_flap_deflection", 0.0)
            trailing_flap_efficiency = kwargs.get("trailing_flap_efficiency", 1.0)

            # Calculate lift coefficient
            CL = self._CLa*(alpha-self._aL0+trailing_flap*trailing_flap_efficiency)
            
            # Saturate
            if isinstance(CL, np.ndarray):
                CL = np.where((CL > self._CL_max) | (CL < -self._CL_max), np.sign(CL)*self._CL_max, CL)
            elif CL > self._CL_max or CL < -self._CL_max:
                CL = np.sign(CL)*self._CL_max

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
            Reynolds number. Defaults to 100000.

        Mach : float, optional
            Mach number. Defaults to 0.

        trailing_flap_deflection : float, optional
            Trailing flap deflection in radians. Defaults to 0.

        trailing_flap_fraction : float, optional
            Trailing flap fraction of the chord length. Defaults to 0.

        trailing_flap_efficiency : float, optional
            Trailing flap efficiency. Defaults to 1.0.

        Returns
        -------
        float or ndarray
            Drag coefficient
        """
        if self._type == "linear":
            df = kwargs.pop("trailing_flap", 0.0)
            kwargs.pop("trailing_flap_efficiency", None)
            CL = self.get_CL(**kwargs)
            CD_flap = 0.002*np.abs(df)*180/np.pi # A rough estimate for flaps
            CD = self._CD0+self._CD1*CL+self._CD2*CL**2+CD_flap

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
            Reynolds number. Defaults to 100000.

        Mach : float, optional
            Mach number. Defaults to 0.

        trailing_flap_deflection : float, optional
            Trailing flap deflection in radians. Defaults to 0.

        trailing_flap_fraction : float, optional
            Trailing flap fraction of the chord length. Defaults to 0.

        trailing_flap_efficiency : float, optional
            Trailing flap efficiency. Defaults to 1.0.

        Returns
        -------
        float or ndarray
            Moment coefficient
        """
        if self._type == "linear":
            alpha = kwargs.get("alpha", 0.0)
            df = kwargs.get("trailing_flap", 0.0)
            flap_eff = kwargs.get("trailing_flap_efficiency", 1.0)
            Cm =  self._Cma*alpha+self._CmL0+df*flap_eff

        # Generated/imported database
        elif self._type == "database":
            Cm = self._get_database_data(2, **kwargs)

        # Fits
        elif self._type == "poly_fit":
            Cm = self._get_polynomial_data(2, **kwargs)

        return Cm


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
        if max_size == 1:
            return interp.griddata(self._data[:,:self._num_dofs], self._data[:,self._num_dofs+data_index].flatten(), param_vals, method='linear').item()
        else:
            return interp.griddata(self._data[:,:self._num_dofs], self._data[:,self._num_dofs+data_index].flatten(), param_vals, method='linear')


    def get_aL0(self, **kwargs):
        """Returns the zero-lift angle of attack

        Parameters
        ----------
        Rey : float, optional
            Reynolds number. Defaults to 100000.

        Mach : float, optional
            Mach number. Defaults to 0.

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
            return self._aL0

        # Database
        elif self._type == "database" or self._type == "poly_fit":

            # Use secant method in alpha to find a_L0
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

        Parameters
        ----------
        alpha : float, optional
            Angle of attack in radians. Defaults to 0.0.

        Rey : float, optional
            Reynolds number. Defaults to 100000.

        Mach : float, optional
            Mach number. Defaults to 0.

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
            return self._CLM

        # Database
        elif self._type == "database" or self._type == "poly_fit":

            # Check the database is dependent on Mach
            if "Mach" not in self._dof_db_order:
                return 0.0

            # Get center Mach value
            dx = kwargs.get("dx", 0.05)
            Mach = kwargs.pop("Mach", 0.0)
            
            # Calculate forward and center points (since we'll often be at M=0 and negative M doesn't work)
            CL1 = self.get_CL(Mach=Mach+dx, **kwargs)
            CL0 = self.get_CL(Mach=Mach, **kwargs)

            return (CL1-CL0)/dx


    def get_CLRe(self, **kwargs):
        """Returns the lift slope with respect to Reynolds number using a backward-difference approximation.

        Parameters
        ----------
        alpha : float, optional
            Angle of attack in radians. Defaults to 0.0.

        Rey : float, optional
            Reynolds number. Defaults to 100000.

        Mach : float, optional
            Mach number. Defaults to 0.

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
            return self._CLRe

        # Database
        elif self._type == "database" or self._type == "poly_fit":

            # Check the database is dependent on Re
            if "Rey" not in self._dof_db_order:
                return 0.0

            # Get center Re value
            dx = kwargs.get("dx", 1000)
            Rey = kwargs.pop("Rey", 100000)
            
            # Calculate forward and backward points
            CL1 = self.get_CL(Rey=Rey+dx, **kwargs)
            CL0 = self.get_CL(Rey=Rey-dx, **kwargs)

            return (CL1-CL0)/(2*dx)


    def get_CLa(self, **kwargs):
        """Returns the lift slope using a backward-difference approximation.

        Parameters
        ----------
        alpha : float, optional
            Angle of attack in radians. Defaults to 0.0.

        Rey : float, optional
            Reynolds number. Defaults to 100000.

        Mach : float, optional
            Mach number. Defaults to 0.

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
        elif self._type == "database" or self._type == "poly_fit":

            # Check the database is dependent on alpha
            if "alpha" not in self._dof_db_order:
                return 0.0

            # Get center alpha value
            dx = kwargs.get("dx", 0.001)
            alpha = kwargs.pop("alpha", 0.0)
            
            # Calculate forward and backward points
            CL1 = self.get_CL(alpha=alpha+dx, **kwargs)
            CL0 = self.get_CL(alpha=alpha-dx, **kwargs)

            return (CL1-CL0)/(2*dx)


    def get_outline_points(self, N=200, cluster=True, trailing_flap_deflection=0.0, trailing_flap_fraction=0.0, export=None, top_first=True, close_te=True):
        """Returns an array of outline points showing the geometry of the airfoil.

        Parameters
        ----------
        N : int
            The number of outline points to return. This function will not always return exactly this many
            but never more. Defaults to 200.

        cluster : bool
            Whether to cluster points about the leading and trailing edges. Defaults to True.

        trailing_flap_deflection : float, optional
            Trailing flap deflection in radians (positive down). Defaults to zero.

        trailing_flap_fraction : float, optional
            Trailing flap fraction of the chord. Defaults to zero.

        export : str
            If specified, the outline points will be saved to a file. Defaults to no file.

        top_first : bool
            The order of the coordinates when exported. Defaults to going from the trailing edge along the top and
            the around to the bottom.

        clos_te : bool
            Whether the top and bottom trailing edge points should be forced to be equal.

        Returns
        -------
        ndarray
            Outline points in airfoil coordinates.
        """

        # Check the geometry has been defined
        if self.geom_specification != "none":

            # Case with no deflection or flap at all
            if trailing_flap_deflection == 0.0 or trailing_flap_fraction == 0.0:

                # Determine spacing of points
                if cluster:
                    # Divide points between top and bottom
                    N_t = int(N*self._s_le)
                    N_b = N-N_t
                    
                    # Create distributions using cosine clustering
                    theta_t = np.linspace(0.0, np.pi, N_t)
                    s_t = 0.5*(1-np.cos(theta_t))*self._s_le
                    theta_b = np.linspace(0.0, np.pi, N_b)
                    s_b = 0.5*(1-np.cos(theta_b))*(1-self._s_le-0.0001)+self._s_le+0.0001
                    # I tack on the miniscule offset to keep it from generating 2 points at exactly the
                    # leading edge. Xfoil seems to not like that...
                    s = np.concatenate([s_t, s_b])
                else:
                    s = np.linspace(0.0, 1.0, N)

                # Get outline
                X = self._x_outline(s)
                Y = self._y_outline(s)

            # Trailing flap deflection
            else:

                # Get flap parameters
                df = trailing_flap_deflection
                x_f = 1.0-trailing_flap_fraction
                y_f = self._trailing_flap_hinge_height

                # Get undeformed camber points
                if cluster:
                    theta_c = np.linspace(0.0, np.pi, N//2)
                    x_c = 0.5*(1-np.cos(theta_c))
                else:
                    x_c = np.linspace(0.0, 1.0, N//2)
                
                y_c = self._camber_line(x_c)

                # Get thickness before we deform the camber
                t = self._thickness(x_c)

                # Determine which camber points belong to the flap
                flap_ind = np.where(x_c>x_f)

                # Linear flap
                if self._trailing_flap_type == "linear":
                    
                    # Calculate deflected camber line Eqs. (8-10) in "Geometry and Aerodynamic Performance of Parabolic..." by Hunsaker, et al. 2018
                    r = np.sqrt((y_c-y_f)*(y_c-y_f)+(x_c-x_f)*(x_c-x_f))
                    psi = np.arctan((y_c-y_f)/(x_c-x_f))
                    x_c[flap_ind] = x_f+(r*np.cos(df-psi))[flap_ind]
                    y_c[flap_ind] = y_f-(r*np.sin(df-psi))[flap_ind]

                # Parabolic flap from "Geometry and Aerodynamic Performance of Parabolic..." by Hunsaker, et al. 2018
                elif self._trailing_flap_type == "parabolic":

                    # Calculate the neutral line parameters
                    l_n = np.sqrt(y_f*y_f+(1-x_f)*(1-x_f))
                    phi_n = -np.arctan2(y_f, 1-x_f)

                    # Calculate the location of the deflected trailing edge
                    tan_df = np.tan(df)
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

                    # Iterate
                    while (abs(R1)>1e-10).any():

                        # Update value
                        # Suppress warnings because an error will often occur within the np.where that has no effect on computation
                        np.seterr(invalid='ignore')
                        E_p2 = np.where(np.abs(R0-R1) != 0.0, E_p1-R1*(E_p0-E_p1)/(R0-R1), E_p1)
                        np.seterr(invalid='warn')

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

                else:
                    raise IOError("{0} is not an allowable flap type.".format(self._trailing_flap_type))

                # Calculate camber line gradient
                dyc_dx = np.gradient(y_c, x_c, edge_order=2)

                # Outline points
                X_t = x_c-t*np.sin(np.arctan(dyc_dx))
                Y_t = y_c+t*np.cos(np.arctan(dyc_dx))
                X_b = x_c+t*np.sin(np.arctan(dyc_dx))
                Y_b = y_c-t*np.cos(np.arctan(dyc_dx))

                # Find the point on the surface where the hinge breaks
                _, y_h_t, r_t = self._get_closest_point_on_surface(x_f, y_f, "top")
                _, y_h_b, r_b = self._get_closest_point_on_surface(x_f, y_f, "bottom")

                # Trim overlapping points off of both surfaces
                X_b, Y_b, num_trimmed_from_bot = self._trim_surface(X_b, Y_b, "forward")
                X_t, Y_t, num_trimmed_from_top = self._trim_surface(X_t, Y_t, "forward")

                # Check if we need to fill anything in
                fill_top = (df > 0 and y_f < y_h_t) or (df < 0 and y_f > y_h_t)
                if fill_top:
                    X_t, Y_t = self._fill_surface(X_t, Y_t, x_f, y_f, r_t, num_trimmed_from_top)

                fill_bot = (df > 0 and y_f < y_h_b) or (df < 0 and y_f > y_h_b)
                if fill_bot:
                    X_b, Y_b = self._fill_surface(X_b, Y_b, x_f, y_f, r_b, num_trimmed_from_bot)

                # Make sure the trailing edge is sealed
                if close_te and self._get_cart_dist(X_t[-1], Y_t[-1], X_b[-1], Y_b[-1]) > 1e-3:
                    X_t = np.concatenate([X_t, X_b[-1][np.newaxis]])
                    Y_t = np.concatenate([Y_t, Y_b[-1][np.newaxis]])

                # Concatenate top and bottom points
                X = np.concatenate([X_t[::-1], X_b])
                Y = np.concatenate([Y_t[::-1], Y_b])

            # Plot result
            if False:
                plt.figure()
                plt.plot(X, Y)
                #plt.plot(x_f, y_f, 'rx')
                plt.gca().set_aspect('equal', adjustable='box')
                plt.show()

            # Concatenate x and y
            outline_points = np.concatenate([X[:,np.newaxis], Y[:,np.newaxis]], axis=1)
            if not top_first:
                outline_points = outline_points[::-1,:]
                    
            # Save to file
            if export is not None:
                np.savetxt(export, outline_points, fmt='%10.8f')

            return outline_points

        else:
            raise RuntimeError("The geometry has not been defined for airfoil {0}.".format(self.name))


    def _trim_surface(self, X, Y, direction):
        # Trims any points that reverse direction in x

        # Loop through points
        trim_indices = []
        last_x_before_break = None

        # Determine direction
        if direction == "forward":
            indices = range(X.shape[0])
        
        for i in indices:
            
            # Skip first point
            if i == 0:
                continue

            # Check if we've passed the break point and store break point
            if last_x_before_break is None and direction == "forward" and X[i] < X[i-1]:
                last_x_before_break = X[i-1]

            # Check if we've now gone backwards
            if last_x_before_break is not None and X[i] < last_x_before_break:
                trim_indices.append(i)

        return np.delete(X, trim_indices), np.delete(Y, trim_indices), len(trim_indices)


    def _fill_surface(self, X, Y, x_h, y_h, r, num_points):
        # Fills in num_points along the arc defined by x_h, y_h, and r

        # Find the two points we need to fill in between
        for i in range(X.shape[0]):
            if X[i] > x_h:
                fill_start = i-1
                fill_stop = i
                break

        # No filling needed, apparently...
        else:
            return X, Y

        # Get angles from the hinge point to the start and stop points
        theta0 = math.atan2(Y[fill_start]-y_h, X[fill_start]-x_h)
        theta1 = math.atan2(Y[fill_stop]-y_h, X[fill_stop]-x_h)

        # Make sure we're actually going to fill in points
        if num_points < 1:
            num_points = 1

        # Find fill in points
        theta_fill = np.linspace(theta0, theta1, num_points+2)[1:-2]
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
            freedom are:

                "alpha"
                "Rey"
                "Mach"
                "trailing_flap_deflection"
                "trailing_flap_fraction"

            Each key should be one of these degrees of freedom. The value is a dictionary or a float
            describing how that DOF should be perturbed. If a dictionary, the following keys must be
            specified:

                "range" : list
                    The lower and upper limits for this DOF.

                "steps" : int
                    The number of points in the range to interrogate.

                "index" : int
                    Index of the column for this degree of freedom in the database.

            If a float, the degree of freedom is assumed to be constant at that value.
            
            If not specified, the above degrees of freedom default to the following:

                "alpha" : 0.0
                "Rey" : 100000.0
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
            print("Convergent results obtained from Xfoil for {0}% of the design space.".format(percent_success))

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

        # Range
        else:
            limits = input_dict.get("range")
            lower = limits[0]
            upper = limits[1]
            N = input_dict.get("steps")
            index = input_dict.get("index", None)
        
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
        """Imports the specified database.

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
            Reynolds number(s) to calculate the coefficients at. Defaults to 100000.

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

        N_crit : float or list, optional
            Critical amplification exponent for the boundary layer in Xfoil. If a float, the value is the same for
            the top and the bottom. If a list, the first list element is for the top and the second list element is
            for the bottom. Defaults to 9.0 for both.

        show_xfoil_output : bool, optional
            Display whatever Xfoil outputs from the command line interface. Defaults to False.

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
        if isinstance(x_trip, float):
            x_trip = [x_trip, x_trip]
        N_crit = kwargs.get("N_crit", [9.0, 9.0])
        if isinstance(N_crit, float):
            N_crit = [N_crit, N_crit]

        # Get states
        # Angle of attack
        alphas = kwargs.get("alpha", [self._dof_defaults["alpha"]])
        if isinstance(alphas, float):
            alphas = [alphas]
        first_dim = len(alphas)
    
        # Reynolds number
        Reys = kwargs.get("Rey", [self._dof_defaults["Rey"]])
        if isinstance(Reys, float):
            Reys = [Reys]
        second_dim = len(Reys)

        # Mach number
        Machs = kwargs.get("Mach", [self._dof_defaults["Mach"]])
        if isinstance(Machs, float):
            Machs = [Machs]
        third_dim = len(Machs)

        # Flap deflections
        delta_fts = kwargs.get("trailing_flap_deflection", [self._dof_defaults["trailing_flap_deflection"]])
        if isinstance(delta_fts, float):
            delta_fts = [delta_fts]
        fourth_dim = len(delta_fts)

        # Flap fractions
        c_fts = kwargs.get("trailing_flap_fraction", [self._dof_defaults["trailing_flap_fraction"]])
        if isinstance(c_fts, float):
            c_fts = [c_fts]
        fifth_dim = len(c_fts)

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
                sp.call(['rm', item])


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

                            # Loop through alphas
                            zero_ind = np.argmin(np.abs(alphas))
                            if zero_ind != len(alphas)-1:
                                for a in alphas[zero_ind:]:
                                    commands.append('ALFA {0:1.6f}'.format(math.degrees(a)))
                            if zero_ind != 0:
                                for a in alphas[zero_ind-1::-1]:
                                    commands.append('ALFA {0:1.6f}'.format(math.degrees(a)))

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
                sp.call(['rm', outline_points])

                # Read in files and store arrays
                for filename in pacc_files:

                    # Read in file
                    try:
                        alpha_i, CL_i, CD_i, Cm_i, Re_i, M_i = self.read_pacc_file(filename)
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
                    sp.call(['rm', filename])

        return CL, CD, Cm


    def read_pacc_file(self, filename):
        """Reads in and formats an Xfoil polar accumulation file.

        Parameters
        ----------
        filename : str
            File to read in.

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
            CD.append(float(split_line[2]))
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

        interaction : bool, optional
            Whether to include interaction terms in the polynomial fit (i.e.
            x^2*y). Defaults to False.

        update_type : bool, optional
            Whether to update the airfoil to use the newly computed polynomial fits for calculations. Defaults to True.

        max_order : int, optional
            Gives the max order of polynomial for any one of the independent varialbes to try. Defaults to the largest number
            of independent variable values present in the database, minus 1.

        sigma : float, optional
            Value used to determine the trade off between how good of a fit to perform and how many terms to keep.
            Defaults to None, which causes the function to calculate sigma automatically using the mean squared of the
            difference of the independent variable values with respect to the mean independent variable value of the dataset
            
        sigma_multiplier : float, optional
            Term multiplied onto sigma to change it's value. Allows using a multiple of the automatically determined sigma
            value. Defaults to 1.

        verbose : bool, optional
        """

        # Check for database
        if not hasattr(self, "_data"):
            raise RuntimeError("No database found! Please generate or import a database before trying to create polynomial fits.")

        # Determine what the maximum fit order is for autoPolyFit
        CL_degrees = kwargs.pop("CL_degrees", {})
        CD_degrees = kwargs.pop("CD_degrees", {})
        Cm_degrees = kwargs.pop("Cm_degrees", {})
        max_order = kwargs.pop("max_order", None)
        if "auto" in [CL_degrees, CD_degrees, Cm_degrees] and max_order is None:
            max_order = 1
            for i in range(self._num_dofs):
                dof_points = np.unique(self._data[:,i])
                dof_max_order = len(dof_points)-1
                max_order = max(max_order, dof_max_order)

        # CL
        if CL_degrees=="auto":
            self._CL_poly_coefs, self._CL_degrees, R2_CL = autoPolyFit(self._data[:,:self._num_dofs], self._data[:, self._num_dofs], max_order=max_order, **kwargs)

        elif isinstance(CL_degrees, dict):

            # Sort fit degrees
            self._CL_degrees = []
            for dof in self._dof_db_order:
                self._CL_degrees.append(CL_degrees.get(dof, 1))

            # Generate
            self._CL_poly_coefs, R2_CL = multivariablePolynomialFit(self._CL_degrees, self._data[:,:self._num_dofs], self._data[:, self._num_dofs], **kwargs)

        else:
            raise IOError("Fit degree specification must be 'auto' or type(dict). Got {0} type {1}.".format(CL_degrees, type(CL_degrees)))
        
        # CD
        if CD_degrees=="auto":
            self._CD_poly_coefs, self._CD_degrees, R2_CD = autoPolyFit(self._data[:,:self._num_dofs], self._data[:,self._num_dofs+1], max_order=max_order, **kwargs)

        elif isinstance(CL_degrees, dict):

            # Sort fit degrees
            self._CD_degrees = []
            for dof in self._dof_db_order:
                self._CD_degrees.append(CD_degrees.get(dof, 1))

            # Generate
            self._CD_poly_coefs, R2_CD = multivariablePolynomialFit(self._CD_degrees, self._data[:,:self._num_dofs], self._data[:,self._num_dofs+1], **kwargs)

        else:
            raise IOError("Fit degree specification must be 'auto' or type(dict). Got {0} type {1}.".format(CL_degrees, type(CL_degrees)))
        
        # Cm
        if Cm_degrees=="auto":
            self._Cm_poly_coefs, self._Cm_degrees, R2_Cm = autoPolyFit(self._data[:,:self._num_dofs], self._data[:,self._num_dofs+1], max_order=max_order, **kwargs)

        elif isinstance(Cm_degrees, dict):

            # Sort fit degrees
            self._Cm_degrees = []
            for dof in self._dof_db_order:
                self._Cm_degrees.append(Cm_degrees.get(dof, 1))

            # Generate polynomial fit
            self._Cm_poly_coefs, R2_Cm = multivariablePolynomialFit(self._Cm_degrees, self._data[:,:self._num_dofs], self._data[:, self._num_dofs+2], **kwargs)

        # Store limits
        self._dof_limits = []
        for i in range(self._num_dofs):
            self._dof_limits.append([np.min(self._data[:,i]), np.max(self._data[:,i])])

        # Update type
        if kwargs.get("update_type", True):
            self.set_type("poly_fit")


    def export_polynomial_fits(self, **kwargs):
        """Save the polynomial fit to a JSON object.

        Parameters
        ----------
        filename : str
            JSON object to write polynomial fit data to.

        """

        # Get filename
        filename = kwargs.get("filename")

        # Create export dictionary
        export = {}
        export["tag"] = "Polynomial fit to database for {0} airfoil.".format(self.name)
        export["degrees_of_freedom"] = self._dof_db_order
        export["limits"] = self._dof_limits
        export["fit_degrees"] = {}
        export["fit_degrees"]["CL"] = self._CL_degrees
        export["fit_degrees"]["CD"] = self._CD_degrees
        export["fit_degrees"]["Cm"] = self._Cm_degrees
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
        self._dof_limits = input_dict["limits"]
        self._CL_degrees = input_dict["fit_degrees"]["CL"]
        self._CD_degrees = input_dict["fit_degrees"]["CD"]
        self._Cm_degrees = input_dict["fit_degrees"]["Cm"]
        self._CL_poly_coefs = np.array(input_dict["fit_coefs"]["CL"])
        self._CD_poly_coefs = np.array(input_dict["fit_coefs"]["CD"])
        self._Cm_poly_coefs = np.array(input_dict["fit_coefs"]["Cm"])

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
        for i in range(N):
            for j in range(self._num_dofs):
                if x[j,i] > self._dof_limits[j][1] or x[j,i] < self._dof_limits[j][0]:
                    coef[i] = np.nan

        return coef