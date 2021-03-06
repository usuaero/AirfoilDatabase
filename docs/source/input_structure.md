# Input Structure
As an initializer, the Airfoil class takes a JSON object (file) or Python dictionary. This object describes the geometry of the airfoil, the trailing flap (if present), and how its aerodynamics are to be predicted, whether using linear predictions, a database, or polynomial fits. The main purpose of this input structure is to allow native interface with [MachUpX](https://www.github.com/usuaero/MachUpX).

Generating or exporting a database cannot be specified through the input; this is done using the Airfoil class methods (see [Airfoil Class](airfoil_class) for more information). If you wish to use the Airfoil class to generate a database, all that is required to input are the geometry specifications.

The following are keys which can be specified in the airfoil JSON object or dictionary.

>**"type" : string**
>>The type of prediction to be used in determining aerodynamic properties. Can be "linear", "database", "poly_fit", or "functional". "linear" airfoils are defined by a linearization of airfoil parameters about the zero-lift angle of attack, with the exception of CD being a quadratic function of CL. All appropriate coefficients and derivatives should then be given. UNITS MAY NOT BE SPECIFIED BY THE USER FOR ANY AIRFOIL PARAMETERS. THESE VALUES MUST BE SPECIFIED IN THE UNITS GIVEN.
>>
>>For the linear class of airfoils, we make corrections to section data due to flap deflection based on Phillips' corrections (see Phillips, *Mechanics of Flight*, 2010, pp. 39-46).
>>
>>"database" airfoils are defined by an array of data, that is aerodynamic coefficients, determined at various angles of attack, Reynolds numbers, Mach numbers, and flap deflections. Such a database must be generated using this software. More information, see generate_database(), export_database(), and import_database() in [Airfoil Class](airfoil_class).
>>
>>"poly_fit" airfoils are defined by polynomial functions of aerodynamic coefficients as a function of angle of attack, Reynolds number, Mach number, and flap deflection. These fits must be generated using this software. More information, see generate_polynomial_fits(), export_polynomial_fits(), and import_polynomial_fits() in [Airfoil Class](airfoil_class).
>>
>>"functional" airfoils are defined by user-defined functions. This type can only be used if the airfoil input is declared as a dictionary within a Python script. The user defines three functions which give the coefficients of lift, drag, and moment as a function of various parameters.
>
>**"input_file" : str, optional**
>>File to read database or fit data from. The file specified here should match the "type" of the airfoil specified (e.g. if "type" is "database", this should be the database file).
>>
>>A note about database files: if the database was generated outside of AirfoilDatabase, be sure the angle of attack is stored in radians, not degrees.
>
>**"geometry" : dict, optional**
>>Describes the geometry of the airfoil.
>>
>>**"outline_points" : str or array, optional**
>>>Path to a file containing airfoil outline points or array of outline points. The first column contains x-coordinates and the second column contains y-coordinates, where the x-axis originates at the leading edge and points back along the chord line and the y-axis points up. If a file, it should be comma-delimited or space-delimited. The points should begin at the trailing edge, wrap around the leading edge, and then end at the trailing edge. The airfoil will automatically determine the ordering of the points (whether top first or bottom first). Cannot be specified along with "NACA".
>>>
>>>Experience has shown that airfoil outlines ending in a blunt tip at a closed trailing edge are poorly handled by AirfoilDatabase. We recommend having either an open trailing edge, or having the trailing edge taper to a smooth point.
>>
>>**"NACA" : str, optional**
>>>NACA designation for the airfoil. If given, the airfoil will automatically generate outline points using the NACA equations. Can only be NACA 4-digit series. Cannot be specified along with "outline_points".
>>
>>**"NACA_closed_te" : bool, optional**
>>>Whether to use the NACA 4-digit equations which force the trailing edge to be closed. Defaults to False.
>>
>>**"max_camber" : float, optional**
>>>Maximum camber of the airfoil as a fraction of the chord. Can be specified if "outline_points" and "NACA" are not specified. This and "max_thickness" affect the corrections for sweep applied in MachUpX. This has no effect on computations made within the Airfoil class. Defaults to 0.0.
>>
>>**"max_thickness" : float, optional**
>>>Maximum thickness of the airfoil as a fraction of the chord. Can be specified if "outline_points" and "NACA" are not specified. This and "max_camber" affect the corrections for sweep applied in MachUpX. This has no effect on computations made within the Airfoil class. Defaults to 0.0.
>
>**"trailing_flap_type" : str, optional**
>>The shape of the flap for this airfoil. Can be "linear" or "parabolic". Defaults to "linear".
>
>**"trailing_flap_hinge_height" : float, optional**
>>y location of the flap hinge, nondimensionalized by the chord length. Defaults to lie on the camber line.
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
>>The quarter-chord moment coefficient at the zero-lift angle of attack. Defaults to 0.0. Only for "linear".
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
>
>The following keys are pertinent to the "functional" type.
>
>**"CL" : func**
>>Function handle returning the coefficient of lift based on the given parameters. The function declaration should look like
>>
>>>def user_CL(**kwargs):
>>>    ...
>>>    return CL
>>
>>Kwargs may be given as floats or as 1-D numpy arrays. This function must be able to handle both (not difficult to implement). The following kwargs can be given to this function:
>>>
>>>"alpha" : float
>>>>Angle of attack in radians.
>>>
>>>"Rey" : float
>>>>Reynolds number.
>>>
>>>"Mach" : float
>>>>Mach number.
>>>
>>>"trailing_flap_deflection" : float
>>>>Trailing flap deflection in radians.
>>>
>>>"trailing_flap_fraction" : float
>>>>Trailing flap fraction of the chord length.
>>
>>This function should have defined defaults for each parameter in case they are not given.
>
>**"CD" : func**
>>Function handle returning the coefficient of drag based on the given parameters. Same as "CL".
>
>**"Cm" : func**
>>Function handle returning the moment coefficient based on the given parameters. Same as "CL".