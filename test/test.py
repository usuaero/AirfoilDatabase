import airfoil_db as adb
import matplotlib.pyplot as plt
import numpy as np

geometry_file = "test/uCRM-9_wr0_xfoil.txt"
#geometry_file = "test/NACA_9412_geom.txt"
#geometry_file = "test/NACA_0012_geom.txt"
airfoil_input = {
    "type" : "linear",
    "geometry" : {
        "outline_points" : geometry_file,
        "top_first" : True
        #"NACA" : "9412"
    },
    "trailing_flap" : {
        "type" : "parabolic",
        "x" : 0.7
    }
}

airfoil = adb.Airfoil("test_airfoil", airfoil_input, verbose=False)
#airfoil.get_outline_points(trailing_flap_deflection=-10.0)

degrees_of_freedom = {
    "alpha" : {
        "range" : [-6.0, 6.0],
        "steps" : 5,
        "index" : 1
    },
    "Rey" : {
        "range" : [100000, 400000],
        "steps" : 4,
        "index" : 2
    },
    #"Mach" : {
    #    "range" : [0.0, 0.4],
    #    "steps" : 4,
    #    "index" : 3
    #},
    "trailing_flap" : {
        "range" : [-5.0, 5.0],
        "steps" : 5,
        "index" : 0
    }
}

#airfoil.generate_database(degrees_of_freedom=degrees_of_freedom)
#airfoil.export_database(filename="database.txt")
#airfoil.import_database(filename="database.txt")

CL = airfoil.get_CL(alpha=0.0, Rey=150000, Mach=0.1, trailing_flap=0.0)
print(CL)
alphas = np.linspace(0, np.pi/20, 10)
CL = airfoil.get_CL(alpha=alphas, Rey=150000, Mach=0.1, trailing_flap=np.pi*3.0/180)
print(CL)
CL = airfoil.get_CL(alpha=0.0, Rey=210000, Mach=0.1, trailing_flap=np.pi*5.0/180)
print(CL)