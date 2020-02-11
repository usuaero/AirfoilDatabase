# Installation

## Prerequisites

### Getting Python

If you do not have Python installed on your machine, it can be downloaded from [https://www.anaconda.com/distribution/](https://www.anaconda.com/distribution/).

### Xfoil

AirfoilDatabase calculates section properties using a code called Xfoil, written by MIT. The Xfoil executable can be downloaded for Windows at [http://web.mit.edu/drela/Public/web/xfoil/](http://web.mit.edu/drela/Public/web/xfoil/). This is also available as a package on some Linux distributions and can be installed like any other package. For example on Mint:

    $ sudo apt-get install xfoil

For other Unix systems where a package is not available, compiling from the Fortran source code is required. Step by step instructions on how to do this can be found on the page [Compiling Xfoil on Unix](xfoil_installation).

**AirfoilDatabase must be able to send only the command ```xfoil``` to the system to run the Xfoil executable.** Please have this properly configured before trying to run AirfoilDatabase.

## Getting the Source Code

You can either download the source as a ZIP file and extract the contents, or clone the AirfoilDatabase repository using Git. If your system does not already have a version of Git installed, you will not be able to use this second option unless you first download and install Git. If you are unsure, you can check by typing `git --version` into a command prompt. We recommend the second option as this makes updating the code much easier.

### Downloading source as a ZIP file

1. Open a web browser and navigate to [https://github.com/usuaero/AirfoilDatabase](https://github.com/usuaero/AirfoilDatabase)
2. Make sure the Branch is set to `Master`
3. Click the `Clone or download` button
4. Select `Download ZIP`
5. Extract the downloaded ZIP file to a local directory on your machine

### Cloning the Github repository

1. From the command prompt navigate to the directory where AirfoilDatabase will be installed. Note: git will automatically create a folder within this directory called MachUpX. Keep this in mind if you do not want multiple nested folders called AirfoilDatabase.
2. Execute `git clone https://github.com/usuaero/AirfoilDatabase`

## Installing

Navigate to the root (AirfoilDatabase/) directory and execute

`pip install .`

Please note that any time you update the source code, the above command will need to be run.

## Testing the Installation

Once the installation is complete, run

`py.test test/`

to verify AirfoilDatabase is working properly on your machine.