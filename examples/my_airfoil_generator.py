import airfoil_db as adb
import matplotlib.pyplot as plt
import numpy as np
import math as m
from mpl_toolkits.mplot3d import Axes3D

if __name__=="__main__":

    # Where the airfoil geometry is stored
    geometry_file = "examples/my_airfoil_geom.txt"

    # Declare input structure
    airfoil_input = {
        "type" : "database",
        "geometry" : {
            "outline_points" : geometry_file
        },
        "trailing_flap_type" : "parabolic"
    }

    # Initialize airfoil
    airfoil = adb.Airfoil("my_airfoil", airfoil_input, verbose=True)
    # Setting verbose to true will let us see how well the camber line/outline are detected and fit

    # Set up database degrees of freedom
    dofs = {
        "alpha" : {
            "range" : [m.radians(-5.0), m.radians(5.0)],
            "steps" : 11,
            "index" : 1
        },
        "Rey" : {
            "range" : [100000, 400000],
            "steps" : 2,
            "index" : 2
        },
        "Mach" : {
            "range" : [0.0, 0.4],
            "steps" : 3,
            "index" : 3
        },
        "trailing_flap_deflection" : {
            "range" : [m.radians(-20.0), m.radians(20.0)],
            "steps" : 7,
            "index" : 0
        },
        #"trailing_flap_fraction" : 0.25
        "trailing_flap_fraction" : {
            "range" : [0.0, 0.5],
            "steps" : 3,
            "index" : 3
        }
    }

    # Generate database using xfoil
    airfoil.generate_database(degrees_of_freedom=dofs, max_iter=100, show_xfoil_output=False)
    airfoil.export_database(filename="my_airfoil_database.txt")

    # If you just want a data database, stop here.
    # If you'd like to create a set of polynomial fits to your database, this next set of code will do that

    # Declare the order of the fit for each degree of freedom for each coefficient
    CL_fit_orders = {
        "alpha" : 3,
        "Rey" : 1,
        "trailing_flap" : 1,
        "Mach" : 2
    }

    CD_fit_orders = {
        "alpha" : 4,
        "Rey" : 1,
        "trailing_flap" : 4,
        "Mach" : 2
    }

    Cm_fit_orders = {
        "alpha" : 1,
        "Rey" : 1,
        "trailing_flap" : 1,
        "Mach" : 2
    }

    # Generate fits
    airfoil.generate_polynomial_fit(CL_degrees=CL_fit_orders, CD_degrees=CD_fit_orders, Cm_degrees=Cm_fit_orders)

    # Alternatively, you can have it automatically detect the order of fit that would be best for each
    #airfoil.generate_polynomial_fit(CL_degrees="auto", CD_degrees="auto", Cm_degrees="auto", max_order=6) 

    # Export fits
    airfoil.export_polynomial_fits(filename="my_airfoil_fits.json")

    # To use these database files in MachUpX, you'll want to declare the airfoil as follows
    #   "my_airfoil" : {
    #       "type" : "database",
    #       "input_file" : "my_airfoil_database.txt",
    #       "geometry" : {
    #           "outline_points" : "my_airfoil_geom.txt"
    #       }
    #   }
    # 
    # or
    #
    #   "my_airfoil" : {
    #       "type" : "poly_fit",
    #       "input_file" : "my_airfoil_fits.json",
    #       "geometry" : {
    #           "outline_points" : "my_airfoil_geom.txt"
    #       }
    #   }