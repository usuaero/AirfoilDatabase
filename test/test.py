import airfoil_db as adb
import matplotlib.pyplot as plt

geometry_file = "test/uCRM-9_wr0_xfoil.txt"
#geometry_file = "test/NACA_9412_geom.txt"
#geometry_file = "test/NACA_0012_geom.txt"
airfoil_input = {
    "type" : "database",
    "geometry" : {
        "outline_points" : geometry_file,
        "top_first" : True
        #"NACA" : "9412"
    },
    "trailing_flap" : {
        "type" : "linear",
        "x" : 0.7
    }
}

airfoil = adb.Airfoil("test_airfoil", airfoil_input, verbose=False)

degrees_of_freedom = {
    "alpha" : {
        "range" : [-5.0, 5.0],
        "steps" : 11,
        "index" : 0
    },
    "Rey" : {
        "range" : [100000, 500000],
        "steps" : 3,
        "index" : 1
    },
    "Mach" : {
        "range" : [0.0, 0.5],
        "steps" : 3,
        "index" : 2
    },
    "trailing_flap_deflection" : {
        "range" : [-20.0, 20.0],
        "steps" : 5,
        "index" : 3
    }
}

airfoil.generate_database(degrees_of_freedom=degrees_of_freedom, max_iter=50)