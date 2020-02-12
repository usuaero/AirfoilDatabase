import airfoil_db as adb
import matplotlib.pyplot as plt
import numpy as np
import math as m

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
        "type" : "parabolic",
        "x" : 0.7
    }
}

airfoil = adb.Airfoil("test_airfoil", airfoil_input, verbose=False)
#airfoil.get_outline_points(trailing_flap_deflection=-10.0)

degrees_of_freedom = {
    "alpha" : {
        "range" : [m.radians(-10.0), m.radians(10.0)],
        "steps" : 10,
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
        "range" : [m.radians(-15.0), m.radians(15.0)],
        "steps" : 5,
        "index" : 0
    }
}

#airfoil.generate_database(degrees_of_freedom=degrees_of_freedom, max_iter=100)
#airfoil.export_database(filename="database.txt")
airfoil.import_database(filename="database.txt")

# Check interpolation
alphas = np.radians(np.linspace(-10, 10, 10))
flaps = np.radians(np.linspace(10, -10, 10))
Re = np.linspace(0, 500000, 10)
CL = airfoil.get_CL(alpha=alphas, trailing_flap=flaps, Rey=Re, Mach=0.1)
CD = airfoil.get_CD(alpha=alphas, trailing_flap=flaps, Rey=Re, Mach=0.1)
Cm = airfoil.get_Cm(alpha=alphas, trailing_flap=flaps, Rey=Re, Mach=0.1)
CLa = airfoil.get_CLa(alpha=alphas, trailing_flap=flaps, Rey=Re, Mach=0.1)
aL0 = airfoil.get_aL0(trailing_flap=flaps, Rey=Re, Mach=0.1)

print(CL)
print(CD)
print(Cm)
print(CLa)
print(aL0)