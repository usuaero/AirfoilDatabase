# Creating Input Files
The input to the airfoil database generator is a JSON file. A JSON file is structured using sets of key-value pairs, much like a Python dictionary. The structure of this file is described below.

## Input File Structure
The following are keys which can be specified in the scene JSON object. NOTE: all keys not marked as optional are required. Key names typed in all capitals between carats (e.g. <KEY_VALUE>) are to be deterimined by the user.

## Airfoil Object
The following are keys which can be specified in the airfoil JSON object.

>**"type" : string**
>>The type of information describing the airfoil. Can be "linear", "database", or "poly_fit". "linear" airfoils are defined by a linearization of airfoil parameters about the zero-lift angle of attack, with the exception of CD being a quadratic function of CL. All appropriate coefficients and derivatives should then be given. UNITS MAY NOT BE SPECIFIED BY THE USER FOR ANY AIRFOIL PARAMETERS. THESE VALUES MUST BE SPECIFIED IN THE UNITS GIVEN.
>>
>>"database" airfoils are defined by an array of data, that is aerodynamic coefficients, determined at various angles of attack, Reynolds numbers, Mach numbers, and flap deflections. Such a database must be generated using this software. More information, see generate_database(), export_database(), and import_database() in [Airfoil Class](airfoil_class).
>>
>>"poly_fit" airfoils are defined by polynomial functions of aerodynamic coefficients as a function of angle of attack, Reynolds number, Mach number, and flap deflection. These fits must be generated using this software. More information, see generate_polynomial_fits(), export_polynomial_fits(), and import_polynomial_fits() in [Airfoil Class](airfoil_class).
>
>**"input_file" : str, optional**
>>File to read database or fit data from. The file specified here should match the "type" of the airfoil specified.
>
>**"geometry" : dict, optional**
>>Describes the geometry of the airfoil.
>>
>>**"outline_points" : str or array, optional**
>>>Path to a file containing airfoil outline points or array of outline points. The first column contains x-coordinates and the second column contains y-coordinates, where the x-axis originates at the leading edge and points back along the chord line and the y-axis points up. If a file, it should be comma-delimited. AirfoilDatabase will automatically determine the ordering of the points (whether top first or bottom first).
>>
>>**"NACA" : str, optional**
>>>NACA designation for the airfoil. If given, MachUpX will automatically generate outline points using the NACA equations. Can only be NACA 4-digit series. Cannot be specified along with "outline_points".
>
>**"trailing_flap" : dict, optional**
>>Specifies a trailing-edge flap for the airfoil.
>>
>>**"type" : str, optional**
>>>The shape of the flap. Can be "linear" or "parabolic". Defaults to "linear".
>>
>>**"x" : float, optional**
>>>x location, nondimensionalized by the chord length, where the hinge is. Defaults to 1.0.
>>
>>**"y" : float, optional**
>>>y location, nondimensionalized by the chord length, where the hinge is. Defaults to 0.0.
>>
>>**"is_sealed" : bool, optional**
>>>Whether the hinge is sealed. Defaults to True.
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