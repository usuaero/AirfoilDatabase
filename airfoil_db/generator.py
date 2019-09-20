"""A class for collecting data on an airfoil from Xfoil and organizing it into a database."""

import numpy as np
from .airfoil import Airfoil
import json
import copy

class Generator(object):
    """A class that creates an airfoil coefficient database from the user's input.

    Parameters
    ----------
    settings : str or dict
        A set of parameters for generating the database. Can be specified either as a
        dictionary or a JSON object.
    """

    def __init__(self, settings):
        
        # Load settings
        if isinstance(settings, str):
            with open(settings, 'r') as input_handle:
                self._input_dict = json.load(input_handle)

        elif isinstance(settings, dict):
            self._input_dict = copy.deepcopy(settings)

        else:
            raise IOError("settings must be str or dict. Got {0}.".format(type(settings)))

        # Create airfoil object
        geometry_type = self._input_dict["geometry"]["type"]

        if geometry_type == "outline_points":
            outline_points = self._input_dict["geometry"]["data"]
            le_coords = self._input_dict["geometry"]["leading_edge_coords"]
            self.airfoil = Airfoil(geom_file=outline_points, leading_edge_coords=le_coords)
