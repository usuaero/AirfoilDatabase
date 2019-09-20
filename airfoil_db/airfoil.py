"""A class defining an airfoil."""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

class Airfoil(object):
    """A class defining an airfoil.

    Parameters
    ----------
    geom_file : str, optional
        Path to a geometry file containing the outline points to the airfoil. Must not have a header. Data should
        be organized in two columns, x and y, using standard airfoil coordinates. x Must go from 1.0 to 0.0 and back
        to 1.0. The top coordinates should be listed first, then the bottom. The leading edge is assumed to be the 
        median of these points.
    """

    def __init__(self, geom_file=None):

        # Import geometry points
        if geom_file is not None:
            self._raw_outline = np.genfromtxt(geom_file, dtype=float)
            self._N0 = self._raw_outline.shape[0]

            self._calc_camber_line()


    def _calc_camber_line(self):
        # Calculates the camber line and thickness distribution from the outline points

        # Create a spline of the outline
        min_ind = np.argmin(self._raw_outline[:,0])
        le_ind = min_ind
        self._x_t = np.copy(self._raw_outline[le_ind::-1,0])
        self._y_t = np.copy(self._raw_outline[le_ind::-1, 1])
        self._x_b = np.copy(self._raw_outline[le_ind:,0])
        self._y_b = np.copy(self._raw_outline[le_ind:, 1])

        self._top_surface = interp.UnivariateSpline(self._x_t, self._y_t, k=5, s=1e-10)
        self._bottom_surface = interp.UnivariateSpline(self._x_b, self._y_b, k=5, s=1e-10)

        # Get first estimate of the camber line and thickness
        self._x_le = self._raw_outline[le_ind,0]
        self._y_le = self._raw_outline[le_ind,1]
        self._x_te = np.max(self._raw_outline[:,0])

        num_camber_points = 100
        le_offset = 0.0
        theta = np.linspace(0, np.pi, num_camber_points)
        x_c = np.linspace(self._x_le+le_offset, self._x_te, num_camber_points)
        y_c = (self._top_surface(x_c)+self._bottom_surface(x_c))/2.0
        t_t = abs(self._top_surface(x_c)-y_c)
        t_b = abs(self._bottom_surface(x_c)-y_c)

        # Show
        plt.figure()
        plt.plot(self._raw_outline[:,0], self._raw_outline[:,1], 'b-', label='Raw Outline')

        x_space = np.linspace(self._x_le, self._x_te, 10000)
        plt.plot(x_space, self._top_surface(x_space), 'k--', label='Top Spline')
        plt.plot(x_space, self._bottom_surface(x_space), 'k--', label='Bottom Spline')

        plt.plot(x_c, y_c, 'r--', label='Camber Line')
        plt.plot(x_c, np.gradient(y_c, x_c, edge_order=2), 'g--', label='Camber Line Derivative')
        plt.legend()
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

        # Iterate
        camber_error = 1
        while camber_error > 1e-9:
            # Determine camber line slope
            dyc_dx = np.gradient(y_c, x_c, edge_order=2)
            
            # Determine slope of lines perpendicular to the camber
            b = -1.0/dyc_dx

            x_t = np.zeros(num_camber_points)
            y_t = np.zeros(num_camber_points)
            x_b = np.zeros(num_camber_points)
            y_b = np.zeros(num_camber_points)

            # Loop through points on the camber line
            for i in range(num_camber_points):
                xc = x_c[i]
                yc = y_c[i]

                # Loop through panels on the surface to see which intersects the line perpendicular to the camber. Start on top.
                for j in range(self._N0-1):

                    # Determine slope of panel
                    x0 = self._raw_outline[j,0]
                    y0 = self._raw_outline[j,1]
                    x1 = self._raw_outline[j+1,0]
                    y1 = self._raw_outline[j+1,1]
                    d = (y1-y0)/(x1-x0)

                    # Find point of intersection
                    x = (-d*x0+y0+b[i]*xc-yc)/(b[i]-d)
                    y = b[i]*(x-xc)+yc

                    # If the intersection is on the panel, then keep that thickness
                    sorted_x = sorted([x0, x1])
                    sorted_y = sorted([y0, y1])
                    if sorted_x[0] <= x <= sorted_x[1] and sorted_y[0] <= y <= sorted_y[1]:
                        t_t[i] = np.sqrt((x-xc)**2+(y-yc)**2)
                        x_t[i] = x
                        y_t[i] = y
                        break

                # Now go through the bottom
                for j in range(self._N0-1, -1, -1):

                    # Determine slope of panel
                    x0 = self._raw_outline[j,0]
                    y0 = self._raw_outline[j,1]
                    x1 = self._raw_outline[j-1,0]
                    y1 = self._raw_outline[j-1,1]
                    d = (y1-y0)/(x1-x0)

                    # Find point of intersection
                    x = (-d*x0+y0+b[i]*xc-yc)/(b[i]-d)
                    y = b[i]*(x-xc)+yc

                    # If the intersection is on the panel, then keep that thickness
                    sorted_x = sorted([x0, x1])
                    sorted_y = sorted([y0, y1])
                    if sorted_x[0] <= x <= sorted_x[1] and sorted_y[0] <= y <= sorted_y[1]:
                        t_b[i] = np.sqrt((x-xc)**2+(y-yc)**2)
                        x_b[i] = x
                        y_b[i] = y
                        break

            # Calculate new camber line points
            x_c_new = 0.5*(x_t+x_b)
            y_c_new = 0.5*(y_t+y_b)

            # Calculate new thickness
            t_t = 0.5*np.sqrt((x_t-x_b)**2+(y_t-y_b)**2)
            t_b = 0.5*np.sqrt((x_t-x_b)**2+(y_t-y_b)**2)

            # Plot
            plt.figure()
            plt.plot(self._raw_outline[:,0], self._raw_outline[:,1], 'b-', label='Outline')
            plt.plot(x_c, y_c, 'r--', label='Old Camber Line')
            plt.plot(x_c_new, y_c_new, 'g--', label='New Camber Line')
            plt.legend()
            plt.gca().set_aspect('equal', adjustable='box')
            plt.show()

            # Update for next iteration
            x_diff = x_c_new-x_c
            y_diff = y_c_new-y_c
            camber_error = np.max(np.sqrt(x_diff*x_diff+y_diff*y_diff))
            print("Camber error: {0}".format(camber_error))

            x_c = x_c_new
            y_c = y_c_new

        # Plot
        plt.figure()
        plt.plot(self._raw_outline[:,0], self._raw_outline[:,1], 'b-', label='Outline')
        plt.plot(x_c, y_c, 'r--', label='Old Camber Line')
        plt.plot(x_c_new, y_c_new, 'g--', label='New Camber Line')
        plt.legend()
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

        # Store camber and thickness
        self._x_c = x_c_new
        self._y_c = y_c_new
        self._camber_line = interp.CubicSpline(self._x_c, self._y_c, extrapolate=True)
        self._t = t_t
        self._thickness = interp.CubicSpline(self._x_c, self._t, extrapolate=True)

        # Show estimate from camber and thickness splines
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
        
        plt.plot(self._raw_outline[:,0], self._raw_outline[:,1], 'b-', label='Original Data')
        plt.plot(self._x_c, self._y_c, 'g--', label='Camber Line')
        plt.plot(x_t, y_t, 'r--', label='Top Fit')
        plt.plot(x_b, y_b, 'r--', label='Bottom Fit')
        plt.legend()
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    
    def _get_closest_outline_point_vertically(self, x, y):
        # Returns the y value on the outline that is closest to (x, y) and has the same x-coordinate
        y_t = self._top_surface(x)
        y_b = self._bottom_surface(x)

        return np.where(x<self._x_le, self._y_le, np.where(np.abs(y_t-y)<np.abs(y_b-y), y_t, y_b))