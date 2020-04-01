# Creating Input Files
As an initializer, the Airfoil class takes a JSON object or Python dictionary. This object describes the geometry of the airfoil, the trailing flap (if present), and how its aerodynamics are to be predicted, whether using linear predictions, a database, or polynomial fits. The main purpose of this input structure is to allow native interface with [MachUpX](https://www.github.com/usuaero/MachUpX). Generating or exporting a database cannot be specified through the input object. This is done using the Airfoil class methods (see [Airfoil Class](airfoil_class) for more information).

## Airfoil Object
The following are keys which can be specified in the airfoil JSON object or dictionary.

>**"type" : string**
>>The type of prediction to be used in determining aerodynamic properties. Can be "linear", "database", or "poly_fit". "linear" airfoils are defined by a linearization of airfoil parameters about the zero-lift angle of attack, with the exception of CD being a quadratic function of CL. All appropriate coefficients and derivatives should then be given. UNITS MAY NOT BE SPECIFIED BY THE USER FOR ANY AIRFOIL PARAMETERS. THESE VALUES MUST BE SPECIFIED IN THE UNITS GIVEN.
>>
>>"database" airfoils are defined by an array of data, that is aerodynamic coefficients, determined at various angles of attack, Reynolds numbers, Mach numbers, and flap deflections. Such a database must be generated using this software. More information, see generate_database(), export_database(), and import_database() in [Airfoil Class](airfoil_class).
>>
>>"poly_fit" airfoils are defined by polynomial functions of aerodynamic coefficients as a function of angle of attack, Reynolds number, Mach number, and flap deflection. These fits must be generated using this software. More information, see generate_polynomial_fits(), export_polynomial_fits(), and import_polynomial_fits() in [Airfoil Class](airfoil_class).
>
>**"input_file" : str, optional**
>>File to read database or fit data from. The file specified here should match the "type" of the airfoil specified (e.g. if "type" is "database", this should be the database file).
>
>**"geometry" : dict, optional**
>>Describes the geometry of the airfoil.
>>
>>**"outline_points" : str or array, optional**
>>>Path to a file containing airfoil outline points or array of outline points. The first column contains x-coordinates and the second column contains y-coordinates, where the x-axis originates at the leading edge and points back along the chord line and the y-axis points up. If a file, it should be comma-delimited or space-delimited. The airfoil will automatically determine the ordering of the points (whether top first or bottom first).
>>
>>**"NACA" : str, optional**
>>>NACA designation for the airfoil. If given, the airfoil will automatically generate outline points using the NACA equations. Can only be NACA 4-digit series. Cannot be specified along with "outline_points".
>
>**"trailing_flap_type" : str, optional**
>>The shape of the flap for this airfoil. Can be "linear" or "parabolic". Defaults to "linear".
>
>**"trailing_flap_hinge_height" : float, optional**
>>y location of the flap hinge, nondimensionalized by the chord length. Defaults to 0.0 (i.e. on the chord line).
>
>The following keys are pertinent to the "linear" type.
>
>**"aL0" : float, optional**
>>The zero-lift angle of attack in radians. Defaults to 0.0. Only for "linear".
>
>**"CLa" : float, optional**
>>The lift slope in radians^-1. Defaults to 2pi Only for "linear".
>
>**"CmL0" : float, optional**
>>The zero-lift moment coefficient. Defaults to 0.0. Only for "linear".
>
>**"Cma" : float, optional**
>>The moment slope in radians^-1. Defaults to 0.0. Only for "linear".
>
>**"CD0" : float, optional**
>>Constant coefficient in the quadratic fit of the CD/CL curve. Defaults to 0.0. Only for "linear".
>
>**"CD1" : float, optional**
>>Linear coefficient in the quadratic fit of the CD/CL curve. Defaults to 0.0. Only for "linear".
>
>**"CD2" : float, optional**
>>Quadratic coefficient in the quadratic fit of the CD/CL curve. Defaults to 0.0. Only for "linear".
>
>**"CL_max" : float, optional**
>>Maximum lift coefficient. Defaults to infinity. Only for "linear".