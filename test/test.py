import airfoil_db as adb

geometry_file = "test/NACA_9412_geom.txt"
airfoil = adb.Airfoil(geom_file=geometry_file, compare_to_NACA="9412", verbose=True)
airfoil.export_outline("outline_test.txt")