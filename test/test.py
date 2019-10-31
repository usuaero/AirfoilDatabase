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
        "range" : [-20.0, 20.0],
        "steps" : 100,
        "index" : 3
    },
    "Rey" : {
        "range" : [100000, 500000],
        "steps" : 4,
        "index" : 0
    },
    "Mach" : {
        "range" : [0.0, 0.4],
        "steps" : 4,
        "index" : 1
    },
    "trailing_flap" : {
        "range" : [-20.0, 20.0],
        "steps" : 10,
        "index" : 2
    }
}

#airfoil.generate_database(degrees_of_freedom=degrees_of_freedom)
#airfoil.export_database(filename="database.txt")
airfoil.import_database(filename="database.txt")

CL = airfoil.get_CL(alpha=0.0, Rey=250000, Mach=0.1, trailing_flap=0.0)
print(CL)
CL = airfoil.get_CL(alpha=0.0, Rey=250000, Mach=0.1, trailing_flap=5.0)
print(CL)
CL = airfoil.get_CL(alpha=0.0, Rey=250000, Mach=0.1, trailing_flap=10.0)
print(CL)