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

    leading_edge_coords : list, optional
        x,y coordinates of the leading edge for an airfoil specified by geom_file. Defaults to [0.0, 0.0].

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

        # Import geometry points
        if geom_file is not None:
            self.type = "geometric"
            self._raw_outline = np.genfromtxt(geom_file, dtype=float)
            self._N0 = self._raw_outline.shape[0]
            self._x_le, self._y_le = kwargs.get("leading_edge_coords", [0.0, 0.0])

            self._calc_camber_line()

            if compare_to_NACA is not None:
                self._check_against_NACA(compare_to_NACA)

        elif NACA is not None:
            self.type = "NACA"
            self._naca_des = NACA

            self._setup_NACA_equations()

        else:
            raise IOError("Airfoil designation not given.")


    def _calc_camber_line(self):
        # Calculates the camber line and thickness distribution from the outline points
        if self._verbose: print("Calculating camber line...")

        # Create splines defining the outline
        le_ind = np.argmin(self._raw_outline[:,0])
        self._top_surface = interp.UnivariateSpline(self._raw_outline[le_ind::-1,0], self._raw_outline[le_ind::-1, 1], k=5, s=1e-10)
        self._bottom_surface = interp.UnivariateSpline(self._raw_outline[le_ind:,0], self._raw_outline[le_ind:, 1], k=5, s=1e-10)

        # Find points where |dy_dx| is greater than 1 and store these in a separate spline as a function of y
        dy_dx = np.gradient(self._raw_outline[:,1], self._raw_outline[:,0])
        self._leading_surface_indices = np.where(np.abs(dy_dx)>1.0)[0][::-1]
        self._leading_surface = interp.UnivariateSpline(self._raw_outline[self._leading_surface_indices,1], self._raw_outline[self._leading_surface_indices,0], k=5, s=1e-10)

        # Find the leading and trailing edge positions
        self._x_te = np.max(self._raw_outline[:,0])

        # Camber line estimate parameters
        num_camber_points = 100
        le_offset = 0.0
        te_offset = 0.001
        camber_deriv_edge_order = 2
        self._cosine_cluster = False

        # Use cosine clustering to space points along the camber line
        if self._cosine_cluster:
            theta = np.linspace(0, np.pi, num_camber_points)
            x_c = 0.5*(1-np.cos(theta))*(self._x_te-te_offset-self._x_le-le_offset)+le_offset
        else:
            x_c = np.linspace(self._x_le+le_offset, self._x_te-te_offset, num_camber_points)


        # Get initial estimate for the camber line and thickness distributions
        y_t = self._top_surface(x_c)
        y_b = self._bottom_surface(x_c)
        y_c = 0.5*(y_t+y_b)

        plt.figure()
        plt.plot(x_c, y_b)
        plt.plot(x_c, y_t)
        plt.plot(x_c, y_c)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

        x_space = np.linspace(0.0, 1.0, 1000)

        # Show
        if self._verbose:
            plt.figure()
            plt.plot(self._raw_outline[:,0], self._raw_outline[:,1], 'b-', label='Outline Data')

            plt.plot(x_space, self._top_surface(x_space), 'k--', label='Top Outline Spline')
            plt.plot(x_space, self._bottom_surface(x_space), 'k--', label='Bottom Outline Spline')
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
            x_c[0] = self._x_le # Force the camber line to always intersect the leading edge
            y_c[0] = self._y_le

            # Sort in x
            sorted_indices = np.argsort(x_c)
            x_c = x_c[sorted_indices]
            y_c = y_c[sorted_indices]

            #plt.figure()
            #plt.plot(x_b, y_b)
            #plt.plot(x_t, y_t)
            #plt.plot(x_c, y_c)
            #plt.gca().set_aspect('equal', adjustable='box')
            #plt.show()

        if self._verbose:
            # Plot
            plt.figure()
            plt.plot(self._raw_outline[:,0], self._raw_outline[:,1], 'b-', label='Outline Data')
            plt.plot(x_c_new, y_c_new, 'g--', label='Final Camber Line Estimate')
            plt.legend()
            plt.gca().set_aspect('equal', adjustable='box')
            plt.show()

        # Store camber
        self._camber_line = interp.UnivariateSpline(x_c_new, y_c_new, k=5, s=1e-10)

        # Calculate thickness
        y_c = self._camber_line(x_space)
        dyc_dx= np.gradient(y_c, x_space)
        b = -1.0/dyc_dx

        # Do the same thing as before but with smaller discretization
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


    def _get_intersection_point(self, xc, yc, b, surface):
        # Calculates the point on the surface where the line extending from (xc, yc) with slope of b.

        # Decide in which direction to go
        if surface == "top":
            surface_range = range(self._N0-1)
        else:
            surface_range = range(self._N0-1, -1, -1)

        # Loop through synthesized panels on the surface to see which intersects the line perpendicular to the camber.
        for j in surface_range:

            # Determine slope of panel
            x0 = self._raw_outline[j,0]
            y0 = self._raw_outline[j,1]
            if surface == "top":
                x1 = self._raw_outline[j+1,0]
                y1 = self._raw_outline[j+1,1]
            else:
                x1 = self._raw_outline[j-1,0]
                y1 = self._raw_outline[j-1,1]
            d = (y1-y0)/(x1-x0)

            # Find point of intersection
            x = (-d*x0+y0+b*xc-yc)/(b-d)
            y = b*(x-xc)+yc

            # Determine if the point of intersection is on the panel
            sorted_x = sorted([x0, x1])
            sorted_y = sorted([y0, y1])

            if sorted_x[0] <= x <= sorted_x[1] and sorted_y[0] <= y <= sorted_y[1]:
                return x, y


    def _check_against_NACA(self, naca_des):
        # Checks the error in the camber and thickness against that predicted by the NACA equations
        print("Checking estimated thickness and camber against NACA equations for NACA {0}...".format(naca_des))

        if len(naca_des) == 4:

            # Generate true coordinates
            x_space = np.linspace(0.0, 1.0, 10000)
            m = float(naca_des[0])/100
            p = float(naca_des[1])/10
            y_c_true = np.where(x_space<p, m/p**2*(2*p*x_space-x_space**2), m/(1-p)**2*(1-2*p+2*p*x_space-x_space**2))

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


    def _setup_NACA_equations(self):
        # Stores the camber and thickness getters based on the NACA designation of the airfoil
        self._m = float(self._naca_des[0])/100
        self._p = float(self._naca_des[1])/10
        self._t = float(self._naca_des[2:])/100

        # Camber line
        def camber(x):
            return np.where(x<p, self._m/self._p*self._p*(2*self._p*x-x*x), self._m/((1-self._p)*(1-self._p))*(1-2*self._p+2*self._p*x-x*x))

        self._camber_line = camber

        # Thickness
        def thickness(x):
            return 5.0*self._t*(0.2969*np.sqrt(x)-0.1260*x-0.3516*x*x+0.2843*x*x*x-0.1015*x*x*x*x)

        self._thickness = thickness


    def export_outline(self, filename, N=200, cosine_cluster=True, trailing_flap_deflection=0.0):
        """Exports the outline points of the airfoil for running through Xfoil.

        Parameters
        ----------
        filename : str
            Name of the file to output the results to.

        N : int, optional
            Number of outline points to export. Must be even. Defaults to 200.

        cosine_cluster : bool, optional
            Whether to use cosine clustering of outline points. Defaults to True

        trailing_flap_deflection : float, optional
            Trailing flap deflection in degrees (positive down). Defaults to zero.
        """
        
        # Check for even number of outline points
        if N%2 != 0:
            raise IOError("Number of outline points to export must be even.")

        # Case with no deflection
        if trailing_flap_deflection == 0.0:
            
            # Determine spacing of x points
            if cosine_cluster:
                theta = np.linspace(0.0, np.pi, N/2)
                x_c = 0.5*(1-np.cos(theta))
            else:
                x_c = np.linspace(0.0, 1.0, N/2)

            # Get camber and thickness
            y_c = self._camber_line(x_c)
            dyc_dx = np.gradient(y_c, x_c)
            t = self._thickness(x_c)
            theta = np.arctan(dyc_dx)
            S_theta = np.sin(theta)
            C_theta = np.cos(theta)

            x_b = x_c+t*S_theta
            y_b = y_c-t*C_theta
            x_t = x_c-t*S_theta
            y_t = y_c+t*C_theta

            outline_points = np.zeros((N, 2))
            mid = int(N/2)
            outline_points[:mid,0] = x_t[::-1]
            outline_points[:mid,1] = y_t[::-1]
            outline_points[mid:,0] = x_b[:]
            outline_points[mid:,1] = y_b[:]

        np.savetxt(filename, outline_points, fmt='%10.5f')