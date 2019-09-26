"""A class defining an airfoil."""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import scipy.signal as sig

class Airfoil:
    """A class defining an airfoil.

    Parameters
    ----------
    name : str
        User specified name for the airfoil.

    geom_file : str, optional
        Path to a geometry file containing the outline points to the airfoil. Must not have a header. Data should
        be organized in two columns, x and y, using standard airfoil coordinates. x Must go from 1.0 to 0.0 and back
        to 1.0. The top coordinates should be listed first, then the bottom. The leading edge is assumed to be the 
        median of these points.

    compare_to_NACA : str, optional
        The NACA designation for an airfoil you wish to compare the geometry to. Must be 4-digit.

    NACA : str, optional
        NACA designation for the airfoil. May not be specified along with a geometry file. Must be 4-digit.

    verbose : bool
    """

    def __init__(self, **kwargs):
        
        self.name = kwargs.get("name")
        geom_file = kwargs.get("geom_file", None)
        compare_to_NACA = kwargs.get("compare_to_NACA", None)
        self._verbose = kwargs.get("verbose", False)
        NACA = kwargs.get("NACA", None)

        # Store undeformed outlines
        self._initialize_geometry(geom_file, NACA)

        # Check how well the camber line algorithm did
        if self.geom_specification == "points" and compare_to_NACA is not None:
            self._check_against_NACA(compare_to_NACA)


    def _initialize_geometry(self, geom_file, NACA):
        # Initialize geometry based on whether points or a NACA designation were given

        # Check that there's only one geometry definition
        if geom_file is not None and NACA is not None:
            raise IOError("Outline points and a NACA designation may not be both specified for airfoil {0}.".format(self.name))

        # Check for user-given points
        if geom_file is not None:
            self.geom_specification = "points"

            with open(geom_file, 'r') as input_handle:
                self._raw_outline = np.genfromtxt(input_handle)

            self._N_orig = self._raw_outline.shape[0]

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
                return np.where(x<self._p, self._m/self._p*self._p*(2*self._p*x-x*x), self._m/((1-self._p)*(1-self._p))*(1-2*self._p+2*self._p*x-x*x))

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

        # Create splines defining the outline
        self._x_outline, self._y_outline = self._create_splines_of_s(self._raw_outline)

        # Camber line estimate parameters
        num_camber_points = 100
        camber_deriv_edge_order = 2
        self._cosine_cluster = False

        # Use cosine clustering to space points along the camber line
        if self._cosine_cluster:
            theta = np.linspace(-np.pi, np.pi, num_camber_points)
            x_c = 0.5*(1-np.cos(theta))*(te[0]-le[0])
        else:
            x_c = np.linspace(le[0], te[0], num_camber_points)


        # Get initial estimate for the camber line and thickness distributions
        y_t = np.interp(x_c, self._raw_outline[le_ind::-1,0], self._raw_outline[le_ind::-1, 1])
        y_b = np.interp(x_c, self._raw_outline[le_ind:,0], self._raw_outline[le_ind:, 1])
        y_c = 0.5*(y_t+y_b)

        x_space = np.linspace(0.0, 1.0, 10000)

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

        iteration = 0
        while camber_error > 1e-10:
            iteration += 1

            # Determine camber line slope
            dyc_dx = np.gradient(y_c, x_c, edge_order=camber_deriv_edge_order)
            
            # Determine slope of lines perpendicular to the camber
            b = -1.0/dyc_dx

            # Loop through points on the camber line to find where their normal intersects the outline
            for i in range(num_camber_points):
                xc = x_c[i]
                yc = y_c[i]
                bi = b[i]

                x_t[i], y_t[i] = self._get_intersection_point(xc, yc, bi, "top")
                x_b[i], y_b[i] = self._get_intersection_point(xc, yc, bi, "bottom")

            # Calculate new camber line points
            x_c_new = 0.5*(x_t+x_b)
            y_c_new = 0.5*(y_t+y_b)

            # Approximate error
            x_diff = x_c_new-x_c
            y_diff = y_c_new-y_c
            camber_error = np.max(np.sqrt(x_diff*x_diff+y_diff*y_diff))
            if self._verbose: print("Camber error: {0}".format(camber_error))

            # Update for next iteration
            x_c = x_c_new
            y_c = y_c_new

        # Calculate where the camber line intersects the outline to find the leading edge
        dyc_dx = np.gradient(y_c, x_c, edge_order=2)
        b = dyc_dx[0]
        x_le, y_le = self._get_intersection_point(x_c[0], y_c[0], b, "leading_edge")
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

            x_t_t[i], y_t_t[i] = self._get_intersection_point(xc, yc, bi, "top")
            x_b_t[i], y_b_t[i] = self._get_intersection_point(xc, yc, bi, "bottom")

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
                s_space = np.linspace(0.0, 1.0, 1000)
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

        return x1, y1


    def _distance(self, x0, y0, x, y, b):
        return (-b*x0+y0+b*x-y)/np.sqrt(b*b+1)


    def _check_against_NACA(self, naca_des):
        # Checks the error in the camber and thickness against that predicted by the NACA equations
        print("Checking estimated thickness and camber against NACA equations for NACA {0}...".format(naca_des))

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


    def get_outline_points(self, N=200, cluster=True, trailing_flap_deflection=0.0, export=None):
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
                # TODO: Implement cosine clustering
                    s = np.linspace(0.0, 1.0, N)
                else:
                    s = np.linspace(0.0, 1.0, N)

                # Get outline
                X = self._x_outline(s)
                Y = self._y_outline(s)
                outline_points=  np.concatenate([X[:,np.newaxis], Y[:,np.newaxis]], axis=1)

            # TODO: Implement flap equations
            else:
                pass

            # Save to file
            if export is not None:
                np.savetxt(export, outline_points, fmt='%10.5f')

        else:
            raise RuntimeError("The geometry has not been defined for airfoil {0}.".format(self.name))


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

            For "trailing_flap_deflection", two more keys may be specified:

                "type" : str
                    May be "parabolic" or "traditional".

                "chord_fraction" : float
                    The fraction of the total chord occupied by the flap.
        """
        pass


    def run_xfoil(self, **kwargs):
        pass