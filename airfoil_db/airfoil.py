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
        self._top_surface = interp.CubicSpline(self._raw_outline[le_ind::-1,0], self._raw_outline[le_ind::-1, 1], extrapolate=True)
        self._bottom_surface = interp.CubicSpline(self._raw_outline[le_ind:,0], self._raw_outline[le_ind:, 1], extrapolate=True)

        # Get first estimate of the camber line and thickness
        x_c = np.linspace(self._raw_outline[le_ind,0], np.max(self._raw_outline[:,0]), 100)
        y_c = (self._top_surface(x_c)+self._bottom_surface(x_c))/2.0
        t_t = abs(self._top_surface(x_c)-y_c)
        t_b = abs(self._bottom_surface(x_c)-y_c)

        # Show
        plt.figure()
        plt.plot(self._raw_outline[:,0], self._raw_outline[:,1], 'b-', label='Raw Outline')
        x_space = np.linspace(self._raw_outline[le_ind,0], np.max(self._raw_outline[:,0]), 10000)
        plt.plot(x_space, self._top_surface(x_space), 'k--', label='Top Spline')
        plt.plot(x_space, self._bottom_surface(x_space), 'k--', label='Bottom Spline')
        plt.plot(x_c, y_c, 'r--', label='Camber Line')
        plt.legend()
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

        # Iterate
        camber_error = 1
        while camber_error > 1e-9:
            print(camber_error)

            # Determine camber line slope
            dyc_dx = np.gradient(y_c, x_c)
            theta = np.arctan(dyc_dx)
            S_theta = np.sin(theta)
            C_theta = np.cos(theta)

            # Determine top thickness distribution
            thickness_error = 1
            while thickness_error > 1e-10:
                print(thickness_error)
                
                # Guess point on surface
                x_t_guess = x_c-t_t*S_theta
                y_t_guess = y_c+t_t*C_theta

                # Find true point
                y_t_true = self._get_closest_outline_point_vertically(x_t_guess, y_t_guess)

                # Estimate error
                correction = y_t_true-y_t_guess
                thickness_error = np.max(np.abs(correction))

                # Apply correction
                t_t += correction

            # Determine bottom thickness distribution
            thickness_error = 1
            while thickness_error > 1e-10:
                print(thickness_error)
                
                # Guess point on surface
                x_b_guess = x_c+t_b*S_theta
                y_b_guess = y_c-t_b*C_theta

                # Find true point
                y_b_true = self._get_closest_outline_point_vertically(x_b_guess, y_b_guess)

                # Estimate error
                correction = y_b_true-y_b_guess
                thickness_error = np.max(np.abs(correction))

                # Apply correction
                t_b -= correction

            # Calculate new camber line points and thickness
            x_c_new = (x_t_guess+x_b_guess)*0.5
            y_c_new = (y_t_guess+y_b_guess)*0.5
            t = np.sqrt((x_t_guess-x_b_guess)**2+(y_t_guess-y_b_guess)**2)
            t_t = t*0.5
            t_b = t*0.5

            # Update for next iteration
            x_diff = x_c-x_c_new
            y_diff = y_c-y_c_new
            camber_error = np.max(np.sqrt(x_diff*x_diff+y_diff*y_diff))
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
        plt.plot(x_t, y_t, 'r--', label='Top Fit')
        plt.plot(x_b, y_b, 'r--', label='Bottom Fit')
        plt.legend()
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    
    def _get_closest_outline_point_vertically(self, x, y):
        # Returns the y value on the outline that is closest to (x, y) and has the same x-coordinate
        y_t = self._top_surface(x)
        y_b = self._bottom_surface(x)

        return np.where(np.abs(y_t-y)<np.abs(y_b-y), y_t, y_b)