"""Custom exceptions used in the airfoil class."""

class DatabaseBoundsError(Exception):
    """An exception thrown when the inputs to the airfoil database fall outside the database bounds."""

    def __init__(self, airfoil, exception_indices, inputs_dict, message="The inputs to the airfoil database fell outside the bounds of available data."):
        self.airfoil = airfoil
        self.exception_indices = exception_indices
        self.inputs_dict = inputs_dict
        self.message = message
        super().__init__(self.message)


    def __str__(self):
        return self.message+" Airfoil: {0}".format(self.airfoil)