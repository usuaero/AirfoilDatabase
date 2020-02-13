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
        "range" : [m.radians(-20.0), m.radians(20.0)],
        "steps" : 5,
        "index" : 0
    }
}

#airfoil.generate_database(degrees_of_freedom=degrees_of_freedom, max_iter=100)
#airfoil.export_database(filename="database.txt")
airfoil.import_database(filename="database.txt")

# Check interpolation
alphas = np.radians(np.linspace(-10, 10, 100))
flaps = np.radians(np.linspace(10, -10, 100))
Re = np.linspace(0, 500000, 100)
CL_int = airfoil.get_CL(alpha=alphas, trailing_flap=flaps, Rey=Re)
CD_int = airfoil.get_CD(alpha=alphas, trailing_flap=flaps, Rey=Re)
Cm_int = airfoil.get_Cm(alpha=alphas, trailing_flap=flaps, Rey=Re)
CLa_int = airfoil.get_CLRe(alpha=alphas, trailing_flap=flaps, Rey=Re)
aL0_int = airfoil.get_aL0(trailing_flap=flaps, Rey=Re)

print("\n---Raw Interpolation Results---")
print(CL_int)
print(CD_int)
print(Cm_int)
print(CLa_int)
print(aL0_int)

CL_fit_orders = {
    "alpha" : 1,
    "Rey" : 1,
    "trailing_flap" : 1
}

CD_fit_orders = {
    "alpha" : 4,
    "Rey" : 1,
    "trailing_flap" : 4
}

Cm_fit_orders = {
    "alpha" : 1,
    "Rey" : 1,
    "trailing_flap" : 1
}

#airfoil.generate_polynomial_fit(CL_fit_orders, CD_fit_orders, Cm_fit_orders)
airfoil.import_polynomial_fits(filename="database.json")

CL_pol = airfoil.get_CL(alpha=alphas, trailing_flap=flaps, Rey=Re)
CD_pol = airfoil.get_CD(alpha=alphas, trailing_flap=flaps, Rey=Re)
Cm_pol = airfoil.get_Cm(alpha=alphas, trailing_flap=flaps, Rey=Re)
CLa_pol = airfoil.get_CLRe(alpha=alphas, trailing_flap=flaps, Rey=Re)
aL0_pol = airfoil.get_aL0(trailing_flap=flaps, Rey=Re)

print("\n---Polynomial Fit Results---")
print(CL_pol)
print(CD_pol)
print(Cm_pol)
print(CLa_pol)
print(aL0_pol)

plt.figure()
plt.plot(flaps, CL_int, 'xr')
plt.plot(flaps, CL_pol, 'b-')
plt.show()

plt.figure()
plt.plot(flaps, CD_int, 'xr')
plt.plot(flaps, CD_pol, 'b-')
plt.show()

plt.figure()
plt.plot(flaps, Cm_int, 'xr')
plt.plot(flaps, Cm_pol, 'b-')
plt.show()

plt.figure()
plt.plot(flaps, CLa_int, 'xr')
plt.plot(flaps, CLa_pol, 'b-')
plt.show()

plt.figure()
plt.plot(flaps, aL0_int, 'xr')
plt.plot(flaps, aL0_pol, 'b-')
plt.show()

airfoil.export_polynomial_fits(filename="database.json")