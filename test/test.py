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

airfoil = adb.Airfoil("test_airfoil", airfoil_input, verbose=True)

outline = airfoil.get_outline_points(200, trailing_flap_deflection=30.0, export="test.txt")
plt.figure()
plt.plot(outline[:,0], outline[:,1])
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

#alphas = [i-5 for i in range(11)]
#Reynolds = [10000.0, 50000.0, 100000.0]
#Machs = [0.0, 0.3, 0.5]
#dfs = [30.0]
#
#CL, CD, Cm = airfoil.run_xfoil(alpha=alphas, Rey=Reynolds, Mach=Machs, trailing_flap_deflection=dfs)
#print(CL)
#print(CD)
#print(Cm)
