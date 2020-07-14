# Error Handling

AirfoilDatabase has a few custom exceptions defined. These can all be accessed through the airfoil_db module. Information on each error and its attributes is provided below.

## DatabaseBoundsError

A **DatabaseBoundsError** will be raised if the inputs to the airfoil database fall outside the bounds of the provided data. AirfoilDatabase cannot perform extrapolation and so an error must be thrown if this occurs. **DatabaseBoundsError** has the following members:

* **airfoil : str** The name of the airfoil for which this exception occurred.
* **inputs_dict : dict** The arguments passed to the airfoil.
* **exception_indices : list** The indices at which the arguments fell outside the database bounds.