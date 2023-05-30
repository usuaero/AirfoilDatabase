import airfoil_db as adb
import matplotlib.pyplot as plt
import numpy as np
import math as m
from mpl_toolkits.mplot3d import Axes3D

#geometry_file = "dev/uCRM-9_wr0_xfoil.txt"
#geometry_file = "dev/symmetric.dat"
#geometry_file = "dev/MMXX0700cTE4.muxf"
geometry_file = "dev/Eppler_335_profile.csv"
#geometry_file = "dev/64A204.txt"
#geometry_file = "dev/NACA_9412_geom.txt"
#geometry_file = "dev/NACA_0012_geom.txt"
airfoil_input = {
    "type" : "database",
    "geometry" : {
        #"outline_points" : geometry_file
        "NACA" : "0012"
        #"NACA_closed_te" : True
    },
    "trailing_flap_type" : "parabolic"
}

airfoil = adb.Airfoil("test_airfoil", airfoil_input, verbose=True, camber_termination_tol=1e-11)
#airfoil.get_outline_points(trailing_flap_fraction=0.3, trailing_flap_deflection=np.radians(0.0), export="dev/symmetric.dat")

dofs = {
    "alpha" : {
        "range" : [m.radians(-15.0), m.radians(15.0)],
        "steps" : 11,
        "index" : 1
    },
    "Rey" : {
        "range" : [1000000, 4000000],
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

# Generate or import database
#airfoil.generate_database(degrees_of_freedom=dofs, max_iter=100, show_xfoil_plots=True, show_xfoil_output=False, visc=True)
airfoil.run_xfoil(alpha=0.0, show_xfoil_output=True)
#airfoil.export_database(filename="database.txt")
#airfoil.import_database(filename="database.txt")
print("Imported database")
airfoil.generate_linear_model()

# Fit orders
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

# Generate or import fits
#airfoil.generate_polynomial_fit(CL_degrees=CL_fit_orders, CD_degrees=CD_fit_orders, Cm_degrees=Cm_fit_orders)
#airfoil.generate_polynomial_fit(CL_degrees="auto", CD_degrees="auto", Cm_degrees="auto", max_order=6)
#airfoil.export_polynomial_fits(filename="fits.json")
#airfoil.import_polynomial_fits(filename="database.json")

# Compare interpolation and fits
alphas = np.radians(np.linspace(-10, 10, 20))
flaps = np.radians(np.linspace(10, -10, 20))
Re = 2000000
c_tf = 0.2

fig = plt.figure(figsize=plt.figaspect(1.0))
ax = fig.gca(projection='3d')
for a in alphas:

    # Interpolation results
    airfoil.set_type("database")
    CL_int = airfoil.get_CL(alpha=a, trailing_flap_deflection=flaps, Rey=Re, trailing_flap_fraction=c_tf)
    CD_int = airfoil.get_CD(alpha=a, trailing_flap_deflection=flaps, Rey=Re, trailing_flap_fraction=c_tf)

    # Polynomial results
    #airfoil.set_type("poly_fit")
    #CL_fit = airfoil.get_CL(alpha=a, trailing_flap_deflection=flaps, Rey=Re, trailing_flap_fraction=c_tf)
    #CD_fit = airfoil.get_CD(alpha=a, trailing_flap_deflection=flaps, Rey=Re, trailing_flap_fraction=c_tf)

    # Plot
    alpha_graph = np.full(20, a)
    ax.plot(alpha_graph, flaps, CL_int, 'rx')
    #ax.plot(alpha_graph, flaps, CL_fit, 'b-')
    #ax.plot(alpha_graph, flaps, CD_int, 'x', color='orange')
    #ax.plot(alpha_graph, flaps, CD_fit, 'g-')

ax.set_xlabel("Angle of Attack")
ax.set_ylabel("Flap Deflection")
ax.set_zlabel("CL")
plt.show()