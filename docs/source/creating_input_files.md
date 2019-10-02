# Creating Input Files
The input to the airfoil database generator is a JSON file. A JSON file is structured using sets of key-value pairs, much like a Python dictionary. The structure of this file is described below.

## Input File Structure
The following are keys which can be specified in the scene JSON object. NOTE: all keys not marked as optional are required. Key names typed in all capitals between carats (e.g. <KEY_VALUE>) are to be deterimined by the user.

## Airfoil Object
The following are keys which can be specified in the airfoil JSON object.
>**"geometry" : dict**
>>Contains information about the geometry of the airfoil.
>>
>>**"NACA" : str, optional
>>>NACA designation