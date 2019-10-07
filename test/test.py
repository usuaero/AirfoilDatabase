import airfoil_db as adb

#geometry_file = "test/uCRM-9_fr0_profile.txt"
geometry_file = "test/NACA_9412_geom.txt"
#geometry_file = "test/NACA_0012_geom.txt"
airfoil_input = {
    "type" : "database",
    "geometry" : {
        #"outline_points" : geometry_file
        "NACA" : "9412"
    }
}

airfoil = adb.Airfoil("test_airfoil", airfoil_input, verbose=True)

alphas = [i-5 for i in range(11)]
Reynolds = [10000.0, 50000.0, 100000.0]
Machs = [0.0, 0.3, 0.5]

CL, CD, Cm = airfoil.run_xfoil(alpha=alphas, Rey=Reynolds, Mach=Machs)