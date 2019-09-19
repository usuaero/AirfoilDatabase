# Compiling Xfoil on Unix (Linux and Mac) Systems
The purpose of this document is to outline how to compile Xfoil on a Unix machine The system I am using to do this is Linux Mint 19.2 on an HP Z230 desktop. I cannot guarantee this process will work on your setup. It is possible to compile Xfoil using only the instructions in the READMEs (that's what I did), but hopefully this will be more clear.

## Compilers
I use the GCC (GNU Compiler Collection). On Debian systems, this can be installed with

    $ sudo apt-get install gcc

## Getting the Source Code
The source code for Xfoil can be downloaded from [http://web.mit.edu/drela/Public/web/xfoil/](http://web.mit.edu/drela/Public/web/xfoil/). Download the .tar.gz file which says it is for Unix. Open this with an archive manager and extract the contents into your home/ directory. Note that .tar.gz is a zip of a zip. Make sure you're extracting the contents of the .tar zip into your home/directory.

## Build Sequence
At this point, we will simply be following the instructions given in the READMEs. Open a terminal window and navigate to the Xfoil/ directory.

### Orr-Sommerfeld Database
Navigate to the orrs/ directory

    $ cd orrs

Determine the absolute path to this directory by typing

    $ pwd

Write this path down. Navigate to the src/ directory and open the osmap.f file using your favorite command-line text editor (mine's Vim).

    $ cd src
    $ vim osmap.f

We need to let the compiler know exactly where the osmap.dat file is. Around line 100 of osmap.f, edit the following line to contain the absolute path you got with the ```pwd``` command

    DATA OSFILE / '/var/local/codes/orrs/osmapDP.dat' /

When done it should look something like this

    DATA OSFILE / '/home/cory/Xfoil/orrs/osmapDP.dat' /

Around line 75, there are options for setting the precision (either double or single). We want to ensure double precision here, so make sure the following lines are not commented

        REAL RLSP, WLSP, HLSP,
        &    RINCR, WINCR, RL, WL, HL,
        &    A, AR, AW, AH, ARW, ARH, AWH, ARWH

Save and exit. Navigate to the Xfoil/bin directory and open the Makefile

    $ cd ../bin
    $ vim Makefile_DP

Lines 14-17 contain the compiler flags. You'll need to edit these to fit your compiler. I made changes to the first two lines to be

    FC = gfortran
    FLG = -0 -f-default-real-8

Save and exit Makefile_DP.

Now it's time to compile the database! Type

    $ make -f Makefile_DP osgen
    $ make -f Makefile_DP osmap.o
    $ cd ..
    $ bin/osgen osmaps_ns.lst

You may get a couple warnings; that's okay.

### Plot Library
Navigate to the plotlib/ directory

    $ cd ../../plotlib

Open Makefile

    $ vim Makefile

Edit line 72 to use the correct compiler

    FC = gfortran

You will also need to set double precision here. Uncomment line 77 and change it to be

    DP = -f-default-real-8

Save and exit. Do the same thing for config.make. Specify the compiler on line 54 and the double precision option on line 58. Save and exit. Now you can compile the plot library

    $ make libPlt.a

### Binaries
Navigate to the bin/ directory and open the Makefile

    $ cd ../bin
    $ vim Makefile

You will need to make changes to lines 102, 111-116 in the Makefile. Line 102 should be

    FC = gfortran

Lines 111-116 should be

    FFLAGS = -O -fdefault-real-8 -B
    FFLOPT = -O -fdefault-real-8 -B
    PLTOBJ = ../plotlib/libPlt.a

    FFLAGS = -O -fdefault-real-8 -ftrapv -fdec
    FFLOPT = -O -fdefault-real-8 -ftrapv -fdec

There are subtle differences between the above and what is originally in the file, so be careful. If you make mistakes, the compiler will usually let you know what it didn't understand. Save and exit.

Compile xfoil

    $ make xfoil

You now need to edit a bit of the source code. Open the pplot.f file

    $ vim ../src/pplot.f

And change line 57 to be

    IF(LERR.ne.0) THEN

Save and exit. Now you can compile the rest of the files

    $ make pplot
    $ make pxplot

### Running Xfoil
The Xfoil executable is located in the bin/ directory. Type

    $ bin/xfoil