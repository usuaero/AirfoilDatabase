"""A class defining an airfoil."""

import numpy as np
import math as m
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import scipy.signal as sig
import subprocess as sp
import json
import copy
import os

class Flap:
    """A class defining a flap.

    Parameters
    ----------
    input_dict : dict
        A dictionary describing the flap.

    loc : str
        Can be "trailing" or "leading".
    """

    def __init__(self, input_dict, loc):
        self.loc = loc
        self.type = input_dict.get("type", "linear")
        if self.loc == "trailing":
            def_x = 1.0
        else:
            def_x = 0.0
        self.x = input_dict.get("x", def_x) # Default behavior is no flap
        self.y = input_dict.get("y", 0.0)


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


    def set_verbosity(self, verbosity):
        """Sets the verbosity of the airfoil."""
        self._verbose = verbosity


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
        self._trailing_flap = Flap(self._input_dict.get("trailing_flap", {}), "trailing")


    def _initialize_geometry(self):
        # Initialize geometry based on whether points or a NACA designation were given

        # Check that there's only one geometry definition
        geom_dict = self._input_dict.get("geometry", {})
        geom_file = geom_dict.get("outline_points", None)
        NACA = geom_dict.get("NACA", None)
        if geom_file is not None and NACA is not None:
            raise IOError("Outline points and a NACA designation may not be both specified for airfoil {0}.".format(self.name))

        # Check for user-given points
        if geom_file is not None:
            self.geom_specification = "points"

            # Import
            with open(geom_file, 'r') as input_handle:
                self._raw_outline = np.genfromtxt(input_handle)

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
                return np.where(x<self._p, self._m/(self._p*self._p)*(2*self._p*x-x*x), self._m/((1-self._p)*(1-self._p))*(1-2*self._p+2*self._p*x-x*x))

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


    def get_CL(self, inputs):
        """Returns the coefficient of lift.

        Parameters
        ----------
        inputs : ndarray
            Parameters which can affect the airfoil coefficients. The first
            three are always alpha, Reynolds number, and Mach number. Fourth 
            is flap efficiency and fifth is flap deflection.

        Returns
        -------
        float
            Lift coefficient
        """
        if self._type == "linear":
            CL = self._CLa*(inputs[0]-self._aL0+inputs[3]*inputs[4])
            if CL > self._CL_max or CL < -self._CL_max:
                CL = np.sign(CL)*self._CL_max
            return CL


    def get_CD(self, inputs):
        """Returns the coefficient of drag

        Parameters
        ----------
        inputs : ndarray
            Parameters which can affect the airfoil coefficients. The first
            three are always alpha, Reynolds number, and Mach number. Fourth 
            is flap efficiency and fifth is flap deflection.

        Returns
        -------
        float
            Drag coefficient
        """
        if self._type == "linear":
            delta_flap = inputs[4]
            inputs_wo_flap = copy.copy(inputs)
            inputs_wo_flap[3:] = 0.0
            CL = self.get_CL(inputs_wo_flap)
            CD_flap = 0.002*np.abs(delta_flap)*180/np.pi # A rough estimate for flaps
            return self._CD0+self._CD1*CL+self._CD2*CL**2+CD_flap


    def get_Cm(self, inputs):
        """Returns the moment coefficient

        Parameters
        ----------
        inputs : ndarray
            Parameters which can affect the airfoil coefficients. The first
            three are always alpha, Reynolds number, and Mach number. Fourth 
            is flap efficiency and fifth is flap deflection.

        Returns
        -------
        float
            Moment coefficient
        """
        if self._type == "linear":
            return self._Cma*inputs[0]+self._CmL0+inputs[3]*inputs[4]


    def get_aL0(self, inputs):
        """Returns the zero-lift angle of attack

        Parameters
        ----------
        inputs : ndarray
            Parameters which can affect the airfoil coefficients. The first
            three are always alpha, Reynolds number, and Mach number. Fourth 
            is flap efficiency and fifth is flap deflection.

        Returns
        -------
        float
            Zero-lift angle of attack
        """
        if self._type == "linear":
            return self._aL0


    def get_CLM(self, inputs):
        """Returns the lift slope with respect to Mach number

        Parameters
        ----------
        inputs : ndarray
            Parameters which can affect the airfoil coefficients. The first
            three are always alpha, Reynolds number, and Mach number. Fourth 
            is flap efficiency and fifth is flap deflection.

        Returns
        -------
        float
            Lift slope with respect to Mach number
        """
        if self._type == "linear":
            return self._CLM


    def get_CLRe(self, inputs):
        """Returns the lift slope with respect to Reynolds number

        Parameters
        ----------
        inputs : ndarray
            Parameters which can affect the airfoil coefficients. The first
            three are always alpha, Reynolds number, and Mach number. Fourth 
            is flap efficiency and fifth is flap deflection.

        Returns
        -------
        float
            Lift slope with respect to Reynolds number
        """
        if self._type == "linear":
            return self._CLRe


    def get_CLa(self, inputs):
        """Returns the lift slope

        Parameters
        ----------
        inputs : ndarray
            Parameters which can affect the airfoil coefficients. The first
            three are always alpha, Reynolds number, and Mach number. Fourth 
            is flap efficiency and fifth is flap deflection.

        Returns
        -------
        float
            Lift slope
        """
        if self._type == "linear":
            return self._CLa


    def get_outline_points(self, N=200, cluster=True, trailing_flap_deflection=0.0, export=None, top_first=True):
        """Returns an array of outline points showing the geometry of the airfoil.

        Parameters
        ----------
        N : int
            The number of outline points to return. Defaults to 200.

        cluster : bool
            Whether to cluster points about the leading and trailing edges. Defaults to True.

        trailing_flap_deflection : float, optional
            Trailing flap deflection in degrees (positive down). Defaults to zero.

        export : str
            If specified, the outline points will be saved to a file. Defaults to no file.

        top_first : bool
            The order of the coordinates when exported. Defaults to going from the trailing edge along the top and
            the around to the bottom.

        Returns
        -------
        ndarray
            Outline points in airfoil coordinates.
        """

        # Check the geometry has been defined
        if self.geom_specification != "none":

            # Case with no deflection
            if trailing_flap_deflection == 0.0:

                # Determine spacing of points
                if cluster:
                    # Divide points between top and bottom
                    N_t = int(N*self._s_le)
                    N_b = N-N_t
                    
                    # Create distributions using cosine clustering
                    theta_t = np.linspace(0.0, np.pi, N_t)
                    s_t = 0.5*(1-np.cos(theta_t))*self._s_le
                    theta_b = np.linspace(0.0, np.pi, N_b)
                    s_b = 0.5*(1-np.cos(theta_b))*(1-self._s_le)+self._s_le
                    s = np.concatenate([s_t, s_b])
                else:
                    s = np.linspace(0.0, 1.0, N)

                # Get outline
                X = self._x_outline(s)
                Y = self._y_outline(s)

            # Trailing flap deflection
            else:
                
                # Linear flap
                if self._trailing_flap.type == "linear":

                    # Get undeformed camber points
                    if cluster:
                        theta_c = np.linspace(0.0, np.pi, N//2)
                        x_c = 0.5*(1-np.cos(theta_c))
                    else:
                        x_c = np.linspace(0.0, 1.0, N//2)

                    y_c = self._camber_line(x_c)

                    # Get thickness before we deform the camber
                    t = self._thickness(x_c)
                    
                    # Calculate deflected camber line Eqs. (8-10) in "Geometry and Aerodynamic Performance of Parabolic..." by Hunsaker, et al. 2018
                    x_f = self._trailing_flap.x
                    y_f = self._trailing_flap.y
                    r = np.sqrt((y_c-y_f)*(y_c-y_f)+(x_c-x_f)*(x_c-x_f))
                    psi = np.arctan((y_c-y_f)/(x_c-x_f))
                    x_c = np.where(x_c<x_f, x_c, x_f+r*np.cos(np.radians(trailing_flap_deflection)-psi))
                    y_c = np.where(x_c<x_f, y_c, y_f-r*np.sin(np.radians(trailing_flap_deflection)-psi))

                    # Calculate outline
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
                    X_b, Y_b = self._trim_surface(X_b, Y_b, "forward")
                    X_t, Y_t = self._trim_surface(X_t, Y_t, "forward")

                    # Check if we need to fill anything in
                    fill_top = (trailing_flap_deflection > 0 and y_f < y_h_t) or (trailing_flap_deflection < 0 and y_f > y_h_t)
                    if fill_top:
                        X_t, Y_t = self._fill_surface(X_t, Y_t, x_f, y_f, r_t)

                    fill_bot = (trailing_flap_deflection > 0 and y_f < y_h_b) or (trailing_flap_deflection < 0 and y_f > y_h_b)
                    if fill_bot:
                        X_b, Y_b = self._fill_surface(X_b, Y_b, x_f, y_f, r_b)

                    # Concatenate top and bottom points
                    X = np.concatenate([X_t[::-1], X_b])
                    Y = np.concatenate([Y_t[::-1], Y_b])

                    # Plot result
                    if True:
                        plt.figure()
                        plt.plot(X, Y)
                        plt.plot(x_f, y_f, 'rx')
                        plt.gca().set_aspect('equal', adjustable='box')
                        plt.show()

                # Parabolic flap
                else:
                    pass

            # Concatenate x and y
            outline_points = np.concatenate([X[:,np.newaxis], Y[:,np.newaxis]], axis=1)
            if not top_first:
                outline_points = outline_points[::-1,:]
                    
            # Save to file
            if export is not None:
                np.savetxt(export, outline_points, fmt='%10.5f')

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

        return np.delete(X, trim_indices), np.delete(Y, trim_indices)


    def _fill_surface(self, X, Y, x_h, y_h, r):
        # Fills in points along the arc defined by x_h, y_h, and r

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
        theta0 = m.atan2(Y[fill_start]-y_h, X[fill_start]-x_h)
        theta1 = m.atan2(Y[fill_stop]-y_h, X[fill_stop]-x_h)

        # Find fill in points
        theta_fill = (theta1+theta0)/2.0
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
        R = 2.0/(1.0+m.sqrt(5.0))

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
        return m.sqrt(x_diff*x_diff+y_diff*y_diff)


    def generate_database(self, degrees_of_freedom):
        """Makes calls to Xfoil to calculate CL, CD, and Cm as a function of each given degree of freedom.

        Parameters
        ----------
        degrees_of_freedom : dict
            A dict specifying which degrees of freedom the database should perturb. Allowable degrees of 
            freedom are:

                "angle_of_attack"
                "reynolds_number"
                "mach_number"
                "trailing_flap_deflection"

            Each key should be one of these degrees of freedom. The value should then be a dictionary 
            describing how that DOF should be perturbed. The following keys must be specified:

                "range" : list
                    The upper and lower limits for this DOF.

                "steps" : int
                    The number of points in the range to interrogate.
        """
        pass


    def run_xfoil(self, **kwargs):
        """Calls Xfoil and extracts the aerodynamic coefficients at the given state.

        Parameters
        ----------
        alpha : float or list of float
            Angle(s) of attack to calculate the coefficients at. Defaults to 0.0.

        Rey : float or list of float
            Reynolds number(s) to calculate the coefficients at. Defaults to 10000.

        Mach : float or list of float
            Mach number(s) to calculate the coefficients at. Defaults to 0.0.

        trailing_flap_deflection : float or list of float
            Flap deflection(s) to calculate the coefficients at. Defaults to 0.0.

        N : int, optional
            Number of panel nodes for Xfoil to use. Defaults to 200.

        max_iter : int, optional
            Maximum iterations for Xfoil. Defaults to 100.
        
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
        max_iter = kwargs.get("max_iter", 50)

        # Get states
        # Angle of attack
        alphas = kwargs.get("alpha", [0.0])
        if isinstance(alphas, float):
            alphas = [alphas]
        first_dim = len(alphas)
    
        # Reynolds number
        Reys = kwargs.get("Rey", [0.0])
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

        # Initialize coefficient arrays
        CL = np.zeros((first_dim, second_dim, third_dim, fourth_dim))
        CD = np.zeros((first_dim, second_dim, third_dim, fourth_dim))
        Cm = np.zeros((first_dim, second_dim, third_dim, fourth_dim))

        # Clean up from previous iterations
        dir_list = os.listdir()
        for item in dir_list:
            if os.path.isfile(item) and ".pacc" in item or ".geom" in item:
                sp.call(['rm', item])


        # Loop through flap deflections
        for l, delta_ft in enumerate(delta_fts):
            pacc_files = []
            
            # Export geometry
            geom_file = os.path.abspath("xfoil_geom_{0}.geom".format(delta_ft))
            self.get_outline_points(N=N, trailing_flap_deflection=delta_ft, export=geom_file)

            # Initialize xfoil execution
            with sp.Popen(['xfoil'], stdin=sp.PIPE, stdout=sp.PIPE) as xfoil_process:

                commands = []

                # Read in geometry
                commands += ['LOAD {0}'.format(geom_file),
                             '{0}'.format(self.name)]

                # Set viscous mode
                commands += ['OPER',
                             'VISC',
                             '',
                             '']
                pacc_index = 0

                # Loop through Mach and Reynolds number
                for Re in Reys:
                    for M in Machs:

                        # Polar accumulation file
                        file_id = str(np.random.randint(0, 10000))
                        pacc_file = "xfoil_results_{0}.pacc".format(file_id)
                        pacc_files.append(pacc_file)

                        # Set up commands
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
                        for a in alphas:
                            commands.append('ALFA {0}'.format(a))

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
                if self._verbose:
                    print(response[0].decode('utf-8'))
                    if response[1] is not None:
                        print(response[1].decode('utf-8'))

            # Clean up geometry
            sp.call(['rm', geom_file])

            # Read in files and store arrays
            for filename in pacc_files:

                # Read in file
                alpha_i, CL_i, CD_i, Cm_i, Re_i, M_i = self.read_pacc_file(filename)

                # Determine the Reynolds and Mach indices
                j = Reys.index(Re_i)
                k = Machs.index(M_i)

                # Loop through alphas
                for i_iter, alpha in enumerate(alpha_i):
                    i_true = alphas.index(alpha)
                    CL[i_true,j,k,l] = CL_i[i_iter]
                    CD[i_true,j,k,l] = CD_i[i_iter]
                    Cm[i_true,j,k,l] = Cm_i[i_iter]

        # Clean up polar files
        dir_list = os.listdir()
        for item in dir_list:
            if os.path.isfile(item) and ".pacc" in item:
                sp.call(['rm', item])

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
                alpha.append(float(split_line[0]))
                CL.append(float(split_line[1]))
                CD.append(float(split_line[2]))
                Cm.append(float(split_line[4]))

        return alpha, CL, CD, Cm, Re, M
