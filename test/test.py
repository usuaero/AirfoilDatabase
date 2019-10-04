import airfoil_db as adb

#geometry_file = "test/uCRM-9_fr0_profile.txt"
geometry_file = "test/NACA_9412_geom.txt"
#geometry_file = "test/NACA_0012_geom.txt"
airfoil_input = {
    "type" : "database",
    #"points" : geometry_file
    "NACA" : "9412"
}
airfoil = adb.Airfoil("test_airfoil", airfoil_input, verbose=True)
airfoil.set_verbosity(True)
airfoil.run_xfoil(alpha=[i-5 for i in range(11)], Rey=[10000, 50000, 100000], Mach=[0.0, 0.3, 0.5])
