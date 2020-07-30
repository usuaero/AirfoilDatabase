"""Custom exceptions used in the airfoil class."""

class DatabaseBoundsError(Exception):
    """An exception thrown when the inputs to the airfoil database fall outside the database bounds.

    Attributes
    ----------

    airfoil : str
        The name of the airfoil for which this exception occurred.

    inputs_dict : dict
        The arguments passed to the airfoil.

    exception_indices : list
        The indices at which the arguments fell outside the database bounds.

    message : str
        A message about the error.
    """

    def __init__(self, airfoil, exception_indices, inputs_dict):
        self.airfoil = airfoil
        self.exception_indices = exception_indices
        self.inputs_dict = inputs_dict
        self.message = "The inputs to the airfoil database fell outside the bounds of available data."
        super().__init__(self.message)


    def __str__(self):
        return self.message+" Airfoil: {0}".format(self.airfoil)