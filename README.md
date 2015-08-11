# Welcome to the TCLX wiki!

TCLX is a C++ Coulomb excitation code based on the Fortran CLX code (written by H. Ower and modified by J. Gerl), using similar syntax and the same methodology. It uses the ROOT and GSL libraries to perform the necessary integrations. 

# Compilation requirements

TCLX requires that the MathMore ROOT libraries are compiled on your machine. These include a number of GSL functions, which simplifies things considerably. Instructions on compilation of ROOT with the MathMore libraries can be found online, [here](https://root.cern.ch/drupal/content/installing-root-source).

In order to compile the TCLX libraries, use the makefile. This has been shamelessly modified from the [GRSISort](https://github.com/GRIFFINCollaboration/GRSISort) makefile, and hence looks very similar.

Once the libTCLX.so library has been compiled, one can simply load it in ROOT:
> .L libTCLX.so 

Or alternatively in your rootlogon file.

# Running the TCLX code

TCLX requires a number of inputs from the user in order to be able to calculate Coulomb excitation cross sections. Required information:
Target and projectile A and Z values,
Beam energy (in MeV),
States of interest and associated energies (in MeV), J values and parities,
Transition matrix elements,
Theta range and the steps over which to integrate.

In order to pass this information to the code, there are three options.

## Using an input file

This is the standard method and is certainly the best option if you're likely to need a large number of states or will want to repeat the calculation with only slight tweaks. The format of the file is as follows:

Line 1: The title. This is ignored by the code, but is useful for keeping track.

Line 2: The number of states we're going to be using.

Line 3: The projectile Z and A

Line 4: The target Z and A

Line 5: The beam energy (in MeV)

Line 6: Theta minimum, maximum and the step size

Lines 7 - 7+number of states: State information in the following form:
State index (incremental integers starting from 1), State J, State E (MeV), State parity (1 = +, 0=-), State K
Note that parity and K are not yet properly implemented in the integration.

Remaining lines: Transition information, in the following form:
State index 1, State index 2, matrix element (NOT B(E2)), multipolarity

This is all of the information the code needs, and can be passed simply as:

> TCLX *myclx = new TCLX()
> myclx->GrabData("filename.txt","v")

Where "v" is the verbose option.

## Input from the terminal

This is not recommended for large calculations (the probability of a mistake tends to 1), but can be done, as follows:
> TCLX *myclx = new TCLX();
> myclx->SetupCLXReaction();

You will then be prompted to manually enter the various variables.

## Input using functions

It is also possible to change individual values using defined "Setter" functions within the code, for example:

>myclx->SetStateEnergy(1,2.5)

Will set the energy of state 1 (note that this is the SECOND state, due to the use of vectors to store data) to be 2.5 MeV.

# Performing the calculation

Now you've loaded all of the necessary information into TCLX, we're nearly ready to perform the calculation. In preparation though, we need to set up all of the matrices that are needed. This is not (yet) done automatically, to give you to the opportunity to check your input variables for errors.

> myclx->PrepareMatrices();

With the matrices now all in order, we can perform the calculation:

> myclx->IntegrateTheta();

Which will integrate across the entire range of theta, as defined in your inputs.

The code will then take a few seconds to calculate the excitation probabilities for all of your requested theta values, before storing them in vectors. They can be accessed individually (myclx->Probabilities() will return a vector of TVectorDs, which contain all of the probabilities).

In order to convert them into cross sections and print them, one can do:

> myclx->PrintCrossSections("outputfilename.txt")

Which will print them to file. Alternatively, if a graphical representation is all you're after:

> myclx->StateCrossSections(1)->Draw()

Will draw a TSpline3 of your cross section. Again note that the indices here are different from those in the input file.

# Defaults/assumptions/things not yet implemented

Currently, the default is for particle excitation. If you want to do target excitation, enter "t" as the second option for the GrabData function.

As of now, magnetic transitions are not included. They will be, once I get around to it.

All cross sections are given in b/sr.

# CoulEx information

For further information on the theory behind Coulomb excitation, see the early chapters of the [GOSIA manual](http://www-user.pas.rochester.edu/~gosia/mediawiki/index.php/Gosia_Manual), or _Electromagnetic excitation_, by K. Alder and A. Winther.
