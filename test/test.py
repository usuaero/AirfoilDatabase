import airfoil_db as adb

#geometry_file = "test/uCRM-9_fr0_profile.txt"
geometry_file = "test/NACA_9412_geom.txt"
#geometry_file = "test/NACA_0012_geom.txt"
airfoil = adb.Airfoil(geom_file=geometry_file, compare_to_NACA="9412", verbose=True)
airfoil.get_outline_points(export="outline_test.txt")