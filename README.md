# Welcome to the TCLX wiki!

TCLX is a C++ Coulomb excitation code based on the Fortran CLX code (written by H. Ower and modified by J. Gerl), using similar syntax and the same methodology. It uses the ROOT and GSL libraries to perform the necessary integrations. 

# Compilation requirements

TCLX requires that the MathMore ROOT libraries are compiled on your machine. These include a number of GSL functions, which simplifies things considerably. Instructions on compilation of ROOT with the MathMore libraries can be found online, [here](https://root.cern.ch/drupal/content/installing-root-source).

You can check if the MathMore libraries are compiled with your version of ROOT by:

> root-config --features

This will list the features, which should include mathmore.

Additionally, in order to use the MathMore libraries, they must be loaded on startup:

> gSystem->Load("libMathMore.so");

> gSystem->Load("libMathCore.so");

Alternatively, you can just include the two lines above in your rootlogon.C file.

In order to compile the TCLX libraries, use the makefile. This has been shamelessly modified from the [GRSISort](https://github.com/GRIFFINCollaboration/GRSISort) makefile, and hence looks very similar.

Once the libTExperiment.so library has been compiled, one can simply load it in ROOT:
> .L libTExperiment.so 

Or alternatively in your rootlogon file.

It can also be included in a macro, for example in the Run.C macro included.

Instructions on running the code can be found in the [Wiki](https://github.com/jhenderson88/TCLX/wiki)
