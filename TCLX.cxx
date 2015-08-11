#define TCLX_cxx

#include "TCLX.h"

//ClassImp(TCLX);

//*****************************************************//
//	Functions herein are involved in the creation
//	of relevant matrices and the integration of the
//	amplitude derivatives, as defined in the GOSIA
//	manual, Eq. 3.18. 
//	Functions included in the TCLX.h file are
//	generally more related to functionalities, such
//	as grabbing the inputs, printing the data, and
//	evaluating the resultant cross sections
//*****************************************************//

//******************************************************** SET UP MATRICES BEFORE BEGINNING INTEGRATION **********************************************************//

//**********************************//
//	Master function prepares all
//	matrices for integration
//**********************************//


void TCLX::PrepareMatrices()
{

	CheckMatrix();
	ConstructLDNUMMatrix();

	ConstructMatrixElementPairsMatrices();

	SetEta();
	MakeXiMatrix();
	MakePsiMatrix();
	MakeMagneticSubstates();
	MakeMRange();
	MakeZeta();

}

//*********************************//
// 	Check matrix symmetry
// 	Identify non-zero matrices
// 	vs lambda
//*********************************//

//Double_t TCLX::Accuracy = 1e-8;

void TCLX::CheckMatrix()
{

	for(unsigned int i=0;i<TransitionMatrix.size();i++)
	{
		Double_t sum = 0.0;
		for(int j=0;j<N_States;j++)
			for(int k=0;k<N_States;k++)
				sum = sum + TransitionMatrix.at(i)[j][k];

		if(sum == 0.0)
			continue;

		if(!TransitionMatrix.at(i).IsSymmetric())
			continue;
		
		if(verbose){
			printf("Non-zero, symmetric matrix, lambda = %i",i+1);
			TransitionMatrix.at(i).Print();
		}
	}

}

//**********************************//
//	Construct a matrix 
//  (Lambda_Level) which 
//  contains information on the
//  number of transitions 
//  associated with each state
//	for each value of lambda
//**********************************//

void TCLX::ConstructLDNUMMatrix()
{

	Lambda_Level.ResizeTo(MaxLambda,N_States); // Resize the matrix to the appropriate size, columns correspond to states, rows to values of lambda
	
	for(unsigned int i=0;i<TransitionMatrix.size();i++) // Loop over the transition matrices for each value of lambda
	{	
		for(int j=0;j<N_States;j++)	// Loop over the states
		{
			Lambda_Level[i][j] = 0; // Set the value to zero - this is technically superfluous, but a safety check
			for(int k=0;k<N_States;k++) // ... Loop over the states again
			{	
				if(TransitionMatrix.at(i)[j][k]!=0) Lambda_Level[i][j]++; // For each transition matrix element associated with state j, increment the Lambda_Level matrix by one		
			}
		}
	}

	if(verbose){ // If we're being verbose, print out the matrix
		printf("LDNUM matrix:\n");
		Lambda_Level.Print();
	}

	MaxStateConnections = GetMaxMatrix(Lambda_Level); // Grab the maximum value from the matrix, this corresponds to the largest number of transition matrix elements associated with a given state

	if(MaxStateConnections>20) // 20 is a crazy-high number - probably would indicate an error in the input
		printf("More than 20 matrix elements connecting to a state");

}

//**********************************//
// 	Based on the Lambda_Level 
// 	matrix, we now construct a
//  matrix (MatrixElementPairs)
//	which tells us which states 
//  are linked, and another
//  (MatrixElementsByConnections) 
// 	which gives us the 
//	corresponding matrix elements
//**********************************//


void TCLX::ConstructMatrixElementPairsMatrices()
{

	Int_t connectioncounter=0; // Keep track of the number of connections

	// The matrices are stored in vectors to allow for one matrix for every value of lambda being used	
	MatrixElementPairs.resize(MaxLambda); // resize the std::vectors to allow for the values of lambda being used
	MatrixElementsByConnections.resize(MaxLambda);

	for(unsigned int i=0;i<TransitionMatrix.size();i++) // Loop over the transition matrix elements - i.e. loop over all values of lambda
	{	
		MatrixElementPairs.at(i).ResizeTo(N_States,MaxStateConnections); // Resize the two matrices
		MatrixElementsByConnections.at(i).ResizeTo(N_States,MaxStateConnections);

		for(int j=0;j<N_States;j++) // Loop over the number of states
		{
			for(int k=0;k<MaxStateConnections;k++) MatrixElementPairs.at(i)[j][k]=-1; // Initially set matrix connections to -1 - allows them to be ignored later if needed
			connectioncounter = 0; // Set the connection counter for the given state to 0

			for(int k=0;k<N_States;k++) // Now we loop over the states
			{	
				if(TransitionMatrix.at(i)[j][k]!=0){ // For a non-zero matrix element
					MatrixElementPairs.at(i)[j][connectioncounter]=k; // Set the MatrixElementsPairs matrix to give the state which state "j" is connected to ("k")
					MatrixElementsByConnections.at(i)[j][connectioncounter]=TransitionMatrix.at(i)[j][k]; // Set the MatrixElementsByConnections to contain the corresponding matrix elements
					connectioncounter++;
				}
			}
		}
		if(verbose && GetMaxMatrix(MatrixElementPairs.at(i))>0){ // Again, if we're being verbose (and the matrix is not empty) print the matrices
			printf("Rows = States, Columns = Connection counter");
			MatrixElementPairs.at(i).Print();
			printf("Rows = States, Columns = Matrix Elements");
			MatrixElementsByConnections.at(i).Print();
		}
	}	

}

//*********************************************************************************************************************//
// 	Create a matrix containing
//  the Eta values, where the Eta value between two states is 
// 	defined as:
//  (Z1 * Z2 * SQRT(A1) / 6.24977) - (1/(SQRT(Ebeam - s*StateEnergy1)) - (1/SQRT(Ebeam - s*StateEnergy2)))
//	here: s = (1 + (A1/A2))
//  the indices 1 & 2 depend on where we're doing target or projectile excitation
//  in this matrix we simply calculate the value for a single state, we'll subtract them from one another later
//*********************************************************************************************************************//

void TCLX::SetEta()
{

	Eta.ResizeTo(N_States); // Set up the TVector
	for(unsigned int i=0;i<level_E.size();i++){ // Loop over the number of states
		Eta[i] = reaction->EtaCalc(level_E.at(i)); // Calculate the Eta value from the TReactions package
	}
	
	if(verbose){ // If verbose, let's print the vector
		printf("Eta:");
		Eta.Print();
	}
}

//******************************************//
// 	Use the Eta matrix to calculate
// 	the differences as defined above
//	for every matrix element 
//  combination
//******************************************//

void TCLX::MakeXiMatrix()
{

	XiMatrix.ResizeTo(N_States,N_States); // Set up the matrix
	for(int i=0;i<MaxLambda;i++) // Loop over the values of lambda
	{
		for(int j=0;j<N_States;j++) // ... and the number of states
		{
			for(int k=0;k<Lambda_Level[i][j];k++) //  ... and the number of matrix element associated with that state
			{
				Int_t State = MatrixElementPairs.at(i)[j][k]; // Grab the state which "j" is connected to
				XiMatrix[j][State] = Eta[j] - Eta[State]; // Calculate the differences in Eta values and stores them in the Xi matrix
			}
		}
	}

	if(verbose){ // Again, if we're being verbose, print the matrix
		printf("Xi matrix");
		XiMatrix.Print();
	}

	XiMax = GetMaxMatrix(XiMatrix); // Get the maximum Xi value
	if(XiMax > 5) XiMax = 5; // It shouldn't be larger than 5

}

//******************************************//
//  Create the Psi matrix, which is used
//  to help determine coupling strengths
//  similar to Xi matrix, depends on 
// 	differences between Ebeam and Elevel
//******************************************//


void TCLX::MakePsiMatrix()
{

	bool IsM1 = false; // This isn't used yet - will be needed for magnetic transitions

	Double_t C; // Set up a variable for use later

	Double_t CPsi_Elec[6] = {5.169286,14.359366,56.982577,263.812653,1332.409500,7117.6915}; // These are variables, CPsi for E1, E2, E3, E4 and E5 transitions
	Double_t CPsi_Mag = reaction->beta * CPsi_Elec[0] / 95.0981942; // For M1 - not used yet

	PsiMatrix.resize(MaxLambda); // Set PsiMatrix for each value of lambda
	for(unsigned int i=0;i<PsiMatrix.size();i++) // Loop over lamba
	{

		Int_t lambda = i+1; // Set lambda to a real value

		Double_t AAZZ = (reaction->Proj_Z * reaction->Tar_Z) * ((Double_t)1 + reaction->Proj_A / reaction->Tar_A); 

		Double_t s = pow(AAZZ,lambda); // Calculate a variable, dependent only on the nuclear masses, charges and the lambda involved

		if(IsM1) // Again, not being used yet
			C = CPsi_Mag / s;
		else 
			C = CPsi_Elec[i] / s;

		PsiMatrix.at(i).ResizeTo(N_States,MaxStateConnections); // Resize the matrix as needed

		for(int j=0;j<N_States;j++) // Loop over states
		{

			Double_t PsiParticle1 = pow((reaction->Elab - (1 + reaction->Proj_A / reaction->Tar_A) * level_E.at(j)),0.25); // Calculate Ebeam - Elevel dependence for state j

			for(int k=0;k<Lambda_Level[i][j];k++)
			{
				Int_t State = MatrixElementPairs.at(i)[j][k]; // Grab connected state
				if(State == -1) // If the state == -1 (as we preset everything to earlier) we can ignore it
					continue;
				Double_t PsiParticle2 = pow((reaction->Elab - (1 + reaction->Proj_A / reaction->Tar_A) * level_E.at(State)),0.25);	// Calculate Ebeam - Elevel dependence for connected state

				Double_t Power = (2*((Double_t)i+1)-1); 

				Double_t PsiVal = C * reaction->Tar_Z * TMath::Sqrt(reaction->Proj_A) * pow((PsiParticle1 * PsiParticle2),Power) * MatrixElementsByConnections.at(i)[j][k];  // Calculate the psi value for the two states, defined as
				// (C * Z1 * SQRT(A1) / (s * Z1* Z2)^lambda) * ((Ebeam - s*Estate1)*(Ebeam-s*Estate2))^((2lambda-1)/4)
				// See GOSIA manual for easier-to-read definitions

				PsiMatrix.at(i)[j][k]=PsiVal; // Stick the value in the matrix
			}

		}
		if(verbose && GetMaxMatrix(PsiMatrix.at(i))>0){ // Again, if verbose, we'll print the matrix
			printf("Psi Matrix:");
			PsiMatrix.at(i).Print();
		}

	}

}

//**************************************************//
//	Now we want to work out how many 
// 	magnetic substates exist for each state.
//	We'll do this by making four matrices
//  one defining the magnetic substates and 
//  three others telling us which ones
// 	correspond to which states
//**************************************************//

void TCLX::MakeMagneticSubstates()
{

	MagneticSubstatesCatalogue.clear(); // Empty all of our vectors
	MagneticSubstatesStart.clear();
	MagneticSubstatesStop.clear();
	MagneticSubstatesHalt.clear();

	Int_t SubStateCounter=0; // Count substates - set to zero

	for(unsigned int i=0;i<level_J.size();i++) // Loop over the states
	{

		MagneticSubstatesStart.push_back(SubStateCounter); // Define the starte of the magnetic substates - i.e. for the ground state, magnetic substate 0 will correspond to the first associated substate of the ground state
		Double_t SpinQuantumNumber = 2*level_J.at(i); // Allows us to work with non-integer spins
		if(debug) // Quick debug statement
			printf("State: %f Spin Quantum Number: %f\n",level_J.at(i),SpinQuantumNumber);

		for(int j=-SpinQuantumNumber;j<=SpinQuantumNumber;j=j+2) // Loop over potential connected states, from -J to +J
		{
			MagneticSubstatesCatalogue.push_back(0.5*(Double_t)j); // Add these substates to the catalogue of substates
			SubStateCounter++; // Increment counter
		}
		MagneticSubstatesStop.push_back(SubStateCounter); // Tells us where the LAST substate associated with each state is in our catalogue of substates

	}

	if(level_J.at(0) == 0) // For even-even nuclei, lets use some symmetry properties
	{
		for(int i=0;i<N_States;i++)
		{
			MagneticSubstatesHalt.push_back(MagneticSubstatesStart.at(i) + level_J.at(i) + 1); // Only use states from -J to 0
		}
	}
	else // Otherwise, we'll use all substates
		MagneticSubstatesHalt = MagneticSubstatesStop;

	if(verbose){ // Standard verbose printing
		printf("\nTotal number of magnetic substates: %i\n\n",SubStateCounter);

		for(unsigned int i=0;i<MagneticSubstatesCatalogue.size();i++) printf("Magnetic substate %i: %lf\n",i,MagneticSubstatesCatalogue.at(i));

		for(unsigned int i=0;i<MagneticSubstatesStart.size();i++){
			printf("\nState %i:\n",i);
			for(int j=MagneticSubstatesStart.at(i);j<MagneticSubstatesHalt.at(i);j++) printf("Magnetic substate %i = %f\n",j,MagneticSubstatesCatalogue.at(j));
		}
		printf("\n");
	}

}

//****************************************************//
//	Not all transitions between substates
//  are allowed by angular momentum conservation
//  this matrix is used to define which ones are
//****************************************************//

void TCLX::MakeMRange()
{

	// Lets count the number of substates
	Int_t MR=0;	
	for(int i=0;i<MaxLambda;i++) // Loop over lambda
	{
		if(GetMaxMatrix(TransitionMatrix.at(i)) == 0) // If the transition matrix is empty, skip it
			continue;
		for(int j=0;j<N_States;j++)
		{
			MR = MR + Lambda_Level[i][j] * (Int_t)(MagneticSubstatesHalt.at(j) - MagneticSubstatesStart.at(j)); // Increment the counter by the number of matrix elements for a given state, multiplied by the number of substates of that state
		}
	}

	if(verbose) // Verbose condition
		printf("MRange elements: %i\n",MR);

	MRange.ResizeTo(MR,2); // Setup the matrix

	MR = 0; // Reset the counter

	ZetaCounter=0; // Set up the counter for the zeta matrix, which we'll construct shortly

	for(int i=0;i<MaxLambda;i++)
	{

		if(GetMaxMatrix(TransitionMatrix.at(i)) == 0) // If the transition matrix is empty, skip it
				continue;

		for(int j=0;j<N_States;j++) // Loop over states - this will be referred to as the "Final State"
		{

			if(GetMaxMatrix(TransitionMatrix.at(i)) == 0) // If the transition matrix element is zero, skip it
				continue;		
	
			for(int k=0;k<Lambda_Level[i][j];k++)
			{

				Int_t State = MatrixElementPairs.at(i)[j][k]; // Find a state connected to state j - "Initial State"

				if(State == -1) // If there's no state, skip it
					continue;

				Int_t Substate1 = MagneticSubstatesStart.at(State); // Grab the first substate of the initial state


				for(int l=MagneticSubstatesStart.at(j);l<MagneticSubstatesHalt.at(j);l++) // Loop over substates of the final state
				{

					// If Final substate - Initial substate (dSubstate) > Multipolarity, return 0
					// Else return (dSubstate) - Multipolarity
					Int_t MinSubstate = (Int_t)std::min((Double_t)0,MagneticSubstatesCatalogue.at(l)-MagneticSubstatesCatalogue.at(Substate1)-(i+1)); // Minimum difference between substates, taking into account the multipolarity

					// If Final substate - Initial substate (dSubstate) > Multipolarity, return (dSubstate) - Multipolarity
					// Else return 0
					Int_t PlusSubstate = (Int_t)std::max((Double_t)0,MagneticSubstatesCatalogue.at(l)-MagneticSubstatesCatalogue.at(Substate1)-(i+1));
					
					// Change substate so that it satisfies angular momentum conservation - if needed
					Int_t Substate2 = Substate1 + PlusSubstate;

					Int_t lambda = i+1; // Set lambda value

					MRange[MR][0]=std::min((Double_t)(2*(i+1) + 1 + MinSubstate), MagneticSubstatesStop.at(State)-(Double_t)Substate2); // whichever is smaller, the difference between the initial and final substate, or the multipolarity of the transition
					MRange[MR][1]=Substate2; // Final substate
					
					ZetaCounter = ZetaCounter + std::max((Double_t)0,MRange[MR][0]); // Increment zeta counter
				

					MR++;
				}
			}
		}
	}

	if(verbose)
		MRange.Print();
	if(verbose)
		printf("ZetaCounter: %i\n",ZetaCounter);

}

//***********************************************//
//	Now we've constructed our MRange matrix
//	we can use it to put together our zeta 
// 	matrix. This matrix defines the coupling
// 	between the intial and final (sub)states
// 	see Eq 3.21 in the GOSIA manual
//***********************************************//

void TCLX::MakeZeta()
{

	Zeta.ResizeTo(ZetaCounter); // Setup the vector

	Int_t MR = 0;
	Int_t NZ = 0;
	for(int i=0;i<MaxLambda;i++)
	{
		if(GetMaxMatrix(TransitionMatrix.at(i)) == 0) // Again, if the transition matrix for this value of lambda is empty, skip it
				continue;		

		Double_t LambdaSqrt = TMath::Sqrt(2*(Double_t)(i+1) + (Double_t)1); // Define a lambda dependent variable (first variable in Eq. 3.21 in GOSIA manual)
		
		for(int j=0;j<N_States;j++) // Loop over states - this will be referred to as the "Final State"
		{

			for(int k=0;k<Lambda_Level[i][j];k++)
			{

				Int_t State = MatrixElementPairs.at(i)[j][k]; // Find a state connected to state j

				if(State == -1) // Non-connected state - skip it
					continue;

				Double_t Psi = PsiMatrix.at(i)[j][k] * LambdaSqrt; // Grab the psi value between those two states

				for(int l=MagneticSubstatesStart.at(j);l<MagneticSubstatesHalt.at(j);l++) // For each substate of the "Final State"
				{

					for(int m=0;m<MRange[MR][0];m++) // Loop over the connected states
					{

						Int_t InitialSubstate = MRange[MR][1] + m; // Grab the substate

						Int_t IIEX = level_J.at(State) - MagneticSubstatesCatalogue.at(InitialSubstate); 
						Int_t Phase = pow((-1),IIEX); // Second variable in Eq. 3.21

						// Define the inputs for the Wigner 3J calculator - all inputs are of 2*J form
						Int_t InitialSpin = 2*level_J.at(State); 
						Int_t FinalSpin = 2*level_J.at(j);
						Int_t InitialSubstateValue = 2*MagneticSubstatesCatalogue.at(InitialSubstate);
						Int_t FinalSubstate = 2*MagneticSubstatesCatalogue.at(l);
						Int_t Lambda = 2*(i+1);
						Int_t Mu = 2*(MagneticSubstatesCatalogue.at(InitialSubstate) - MagneticSubstatesCatalogue.at(l));
						
						// Calculate the Zeta value
						Zeta[NZ] = (Double_t)Phase * Psi * reaction->ThreeJ(InitialSpin,Lambda,FinalSpin,-InitialSubstateValue,Mu,FinalSubstate);

						NZ++;
					}

					MR++;
				}

			}

		}


	}
	
	if(verbose){ // Standard verbose
		printf("Zeta:");
		Zeta.Print();
	}

}

//******************************************************** INTEGRATION SECTION BEGINS **********************************************************//

//*****************************************************//
//	Right, now we've set up all of the matrices
//	we need to perform the integration, let's 
//  start integrating. We need to perform this 
// 	integration for every value of theta we want
//	to consider
//*****************************************************//

void TCLX::IntegrateTheta()
{

	TVectorD temp; // We set up a TVectorD, this will be a temporary storage for the population probabilities so needs to be resized to allow a probability for each state (not substate)
	temp.ResizeTo(N_States);

	stepfile.open("CoulExStepChanges.txt"); // Open an output file, in here we'll put all of the details of the calculations

	printf("Integrating between Theta = %f and Theta = %f, in steps of %f\n",Theta_Min,Theta_Max,Theta_Step); // Output the theta range we're working with, just to make sure

	for(Double_t i=Theta_Min; i<=Theta_Max; i=i+Theta_Step){ // Loop over the defined theta range, in the defined number of steps

		temp = Integration(i); // Perform the integration for the theta value given
		Probabilities.push_back(temp); // Push the resulting TVectorD into a std::vector
		stepfile << "\n";
		if((int)i % 1 == 0) std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(3) << 100 * (i-Theta_Min)/(Theta_Max-Theta_Min) << " % complete" << "\r" << std::flush; // Progress bar - important to keep distracted

	}

	stepfile.close(); // Close our output file, now we're done

	printf("Integration complete!!!\n"); // SUCCESS!
}

//****************************************************//
//	This is where we're going to actually do the
//	integration.
//
//	Firstly, we need to redefine our coordinate
//	system into something more sensible than x,y,z
//	
//	Then we can begin our integration, using some
//	numerical differentiation methods. We begin 
//	with the more accurate Runge-Kutta method,
//	this gives us the starting point.
//	Once we have that, we move to the faster
//	Adams-Boulton method. If we notice that our 
// 	values lie outside a defined acceptable
//	accuracy, we go back to the Runge-Kutta
//  and redefine our incrementation to be
// 	smaller (or larger)
//
//	Once we've completed our integration, we can
// 	calculate probabilities and send them back
//	to the IntegrateTheta() function.
//****************************************************//


TVectorD TCLX::Integration(Double_t theta)
{

	// Set up a print flag. We will print the details of the calculation to the terminal only for the first and last integration steps. The rest of the details will be printed elsewhere, don't worry!
	bool PrintFlag=false;
	if(theta + Theta_Step > Theta_Max) // Define the good values of PrintFlag
		PrintFlag=true;
	if(theta == Theta_Min)
		PrintFlag=true;

	// Here we define ourselves std::vectors of TMatrices and TVectors, note that we define vector.at(0) to be the real part and vector.at(1) to be the imaginary part
	std::vector<TMatrixD> Amplitude;
	std::vector<TMatrixD> AmplitudeP; // Used for Adams-Boulton
	std::vector<TMatrixD> AmplitudeDerivative;
	std::vector<TMatrixD> RealAmpF; // These two matrices (RealAmpF and ImagAmpF) give us the previous four steps for our numerical differentiation - differ from the aforementioned real/imaginary vector convention
	std::vector<TMatrixD> ImagAmpF;
	std::vector<TMatrixD> Q1_Matrix;

	// These are the vectors which we will eventually dump our probabilities into, here we define and resize them to appropriate values
	TVectorD SubStateProbabilities; // Substate excitation probabilities
	TVectorD StateProbabilities; // State excitation probabilities
	SubStateProbabilities.ResizeTo(MagneticSubstatesCatalogue.size());
	StateProbabilities.ResizeTo(N_States);

	// Define an initial TotalProbability - soon to be overwritten
	Double_t TotalProbability = 0;
	
	// Resize the vectors for their purpose (2 for most - real, imaginary, and 4 for RealAmpF and ImagAmpF, for the numerical differentiation)
	Amplitude.resize(2);
	AmplitudeP.resize(2);
	AmplitudeDerivative.resize(2);
	RealAmpF.resize(4);
	ImagAmpF.resize(4);
	Q1_Matrix.resize(2);

	// Resize the various matrices to the appropriate sizes
	Amplitude.at(0).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	Amplitude.at(1).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	AmplitudeP.at(0).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	AmplitudeP.at(1).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	AmplitudeDerivative.at(0).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	AmplitudeDerivative.at(1).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	RealAmpF.at(0).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	RealAmpF.at(1).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	RealAmpF.at(2).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	RealAmpF.at(3).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	ImagAmpF.at(0).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	ImagAmpF.at(1).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	ImagAmpF.at(2).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	ImagAmpF.at(3).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	Q1_Matrix.at(0).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	Q1_Matrix.at(1).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	
	// Now we're onto the serious stuff, here we define the accuracy and the range and steps of the integration we're going to perform. We also redefine everything into a better coordinate system.
	// For more details on the coordinate system especially, see GOSIA manual section 3.3 - but note that omega = 0 corresponds to the closest approach
	Double_t Accuracy = 1e-8; // This is our accuracy. We could go more accurate, but it will slow down the calculation.
	Double_t Accuracy_50 = Accuracy / 50.; // Used to define whether to increase the integration steps
	Double_t Epsilon = 1 / (TMath::Sin(theta * TMath::DegToRad() / 2.0)); // This is a redefinition of theta into units of eliptical eccentricity, more useful for what we're doing
	Double_t RootEpsilon = TMath::Sqrt(pow(Epsilon,2) - 1);
	Double_t Distance = reaction->ClosestApproach() * (1 + Epsilon) / 2.0; // Defines the closest approach given our value of theta
	Double_t Up = TMath::Log(1 / (Epsilon * TMath::Sqrt(Accuracy))); // This is our start/end "distance", as defined by omega in the GOSIA manual
	Double_t ABW = 0.0;
	Double_t dOmega = (Double_t) 40.0 * (pow(Accuracy,0.2)) / (10.0 + 48.0 * XiMax + 16.0 * XiMax * Epsilon); // Our initial step value is defined here - for reasoning, see GOSIA manual
	Double_t Steps = Up / (8.0*dOmega) + 1;
	Steps = Steps*8.0;
	dOmega = Up / (Steps - 0.25);
	Up = dOmega * Steps;
	Double_t Omega = -Up; // Set our initial time/distance to be -the end time/distance

	if(PrintFlag){ // Print these variables, assuming the printflag is good
		printf("\nUp: %f\n",Up);
		printf("Xi Max: %f\n",XiMax);
		printf("Epsilon: %f\n",Epsilon);
		printf("Accuracy: %E\n",Accuracy);
		printf("Omega (start): %f\n", Omega);
		printf("dOmega: %lf\n",dOmega);
		printf("Steps: %f\n",Steps);
	}

	// Output everything so far to the output file for any required debugging
	stepfile << "Theta: " << theta << "\n";
	stepfile << "\nHyperbola eccentricity: " << Epsilon << ", Closest approach: " << Distance << " fm\n"; 
	stepfile << "Integrate between Omega = " << Omega << " and " << Up << ", initial stepwidth = " << 2*dOmega << "\n\n";

	// To begin with, the only non-zero amplitude should be the ground state, which should have a real amplitude of 1. Everything else should be zero.
	for(int i=0;i<LMax;i++) Amplitude.at(0)[0][i] = 1;
	
	if(theta == Theta_Min && debug) // Just some debugging checks
		Amplitude.at(0).Print();
	if(theta == Theta_Min && debug)
		Amplitude.at(1).Print();

	// Now we get into the integration process itself, this is all based on numerical methods, no real physics going on

	while(Omega <= Up) // Loop over Omega until it hits the end point
	{
		if(theta == Theta_Min && debug) // Debug stuff
		{
			printf("Real Amp Matrix:\n");
			Amplitude.at(0).Print();
			printf("Imag Amp Matrix:\n");
			Amplitude.at(1).Print();
		}

		// First we're going to use the accurate Runge-Kutta method to get four intitial points: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
		// I'm not going to describe the method in too much detail

		bool dOmegaFlag = true; // This flag tells us that dOmega is initially "good", that is, it fits with our required accuracy
		AmplitudeDerivative = ComputeAmpDerivativeMatrices(Amplitude, Epsilon, Omega); // Determine the derivatives of the amplitude matrices
		RealAmpF.at(0) = AmplitudeDerivative.at(0); // Stick them into the Real and Imaginary matrices as value 0
		ImagAmpF.at(0) = AmplitudeDerivative.at(1);

		for(int n=1;n<4;n++) // Now we loop over the next three points to get our four start points
		{

			Q1_Matrix.at(0) = dOmega * AmplitudeDerivative.at(0);
			Q1_Matrix.at(1) = dOmega * AmplitudeDerivative.at(1); 
			Amplitude.at(0) = Amplitude.at(0) + Q1_Matrix.at(0); 
			Amplitude.at(1) = Amplitude.at(1) + Q1_Matrix.at(1);

			if( level_J.at(0) == 0 ) // Use the symmetry of even-even nuclei to speed things up a bit, just integrate substates between -j and 0
			{
				for(int i=0;i<LMax;i++)
				{
					for(int j=0;j<N_States;j++)
					{
						for(int k=MagneticSubstatesStart.at(j);k<MagneticSubstatesHalt.at(j);k++)
						{
							Int_t IR = k - 2*MagneticSubstatesCatalogue.at(k);
							Amplitude.at(0)[IR][i] = (Double_t)IFAC.at(j) * Amplitude.at(0)[k][i];
							Amplitude.at(1)[IR][i] = (Double_t)IFAC.at(j) * Amplitude.at(1)[k][i];
						}
					}
				}
			}
			Omega = Omega + dOmega; // Increment Omega

			AmplitudeDerivative = ComputeAmpDerivativeMatrices(Amplitude, Epsilon, Omega); 

			Amplitude.at(0) = Amplitude.at(0) + 0.5857864 * (dOmega*AmplitudeDerivative.at(0) - Q1_Matrix.at(0)); 
			Amplitude.at(1) = Amplitude.at(1) + 0.5857864 * (dOmega*AmplitudeDerivative.at(1) - Q1_Matrix.at(1)); 
			Q1_Matrix.at(0) = 0.5857864 * dOmega * AmplitudeDerivative.at(0) + 0.1213204 * Q1_Matrix.at(0); 
			Q1_Matrix.at(1) = 0.5857864 * dOmega * AmplitudeDerivative.at(1) + 0.1213204 * Q1_Matrix.at(1); 

			if( level_J.at(0) == 0 ) // Use the symmetry of even-even nuclei to speed things up a bit, just integrate substates between -j and 0
			{
				for(int i=0;i<LMax;i++)
				{
					for(int j=0;j<N_States;j++)
					{
						for(int k=MagneticSubstatesStart.at(j);k<MagneticSubstatesHalt.at(j);k++)
						{
							Int_t IR = k - 2*MagneticSubstatesCatalogue.at(k);
							Amplitude.at(0)[IR][i] = (Double_t)IFAC.at(j) * Amplitude.at(0)[k][i];
							Amplitude.at(1)[IR][i] = (Double_t)IFAC.at(j) * Amplitude.at(1)[k][i];
						}
					}
				}
			}
			AmplitudeDerivative = ComputeAmpDerivativeMatrices(Amplitude, Epsilon, Omega); 
			
			Amplitude.at(0) = Amplitude.at(0) + 3.414214 * (dOmega * AmplitudeDerivative.at(0) - Q1_Matrix.at(0)); 
			Amplitude.at(1) = Amplitude.at(1) + 3.414214 * (dOmega * AmplitudeDerivative.at(1) - Q1_Matrix.at(1)); 
			Q1_Matrix.at(0) = 3.414214 * dOmega * AmplitudeDerivative.at(0) - 4.1213204 * Q1_Matrix.at(0); 
			Q1_Matrix.at(1) = 3.414214 * dOmega * AmplitudeDerivative.at(1) - 4.1213204 * Q1_Matrix.at(1); 

			if( level_J.at(0) == 0 ) // Use the symmetry of even-even nuclei to speed things up a bit, just integrate substates between -j and 0
			{
				for(int i=0;i<LMax;i++)
				{
					for(int j=0;j<N_States;j++)
					{
						for(int k=MagneticSubstatesStart.at(j);k<MagneticSubstatesHalt.at(j);k++)
						{
							Int_t IR = k - 2*MagneticSubstatesCatalogue.at(k);
							Amplitude.at(0)[IR][i] = (Double_t)IFAC.at(j) * Amplitude.at(0)[k][i];
							Amplitude.at(1)[IR][i] = (Double_t)IFAC.at(j) * Amplitude.at(1)[k][i];
						}
					}
				}
			}
			Omega = Omega + dOmega;

			AmplitudeDerivative = ComputeAmpDerivativeMatrices(Amplitude, Epsilon, Omega);

			Amplitude.at(0) = Amplitude.at(0) + (1./3.) * dOmega * AmplitudeDerivative.at(0) - (2./3.) * Q1_Matrix.at(0); 
			Amplitude.at(1) = Amplitude.at(1) + (1./3.) * dOmega * AmplitudeDerivative.at(1) - (2./3.) * Q1_Matrix.at(1);  

			if( level_J.at(0) == 0 ) // Use the symmetry of even-even nuclei to speed things up a bit, just integrate substates between -j and 0
			{
				for(int i=0;i<LMax;i++)
				{
					for(int j=0;j<N_States;j++)
					{
						for(int k=MagneticSubstatesStart.at(j);k<MagneticSubstatesHalt.at(j);k++)
						{
							Int_t IR = k - 2.0*MagneticSubstatesCatalogue.at(k);
							Amplitude.at(0)[IR][i] = (Double_t)IFAC.at(j) * Amplitude.at(0)[k][i];
							Amplitude.at(1)[IR][i] = (Double_t)IFAC.at(j) * Amplitude.at(1)[k][i];
						}
					}
				}
			}
			AmplitudeDerivative = ComputeAmpDerivativeMatrices(Amplitude, Epsilon, Omega);

			RealAmpF.at(n) = AmplitudeDerivative.at(0);
			ImagAmpF.at(n) = AmplitudeDerivative.at(1);

		} 

		// Now that the Runge-Kutta has given us our initial points, we can use the fast Adams-Moulton method to get the rest done: https://en.wikipedia.org/wiki/Linear_multistep_method#Adams.E2.80.93Moulton_methods

		while( Omega < Up && dOmegaFlag) // As long as the dOmegaFlag condition is still true, that is
		{

			AmplitudeP.at(0) = Amplitude.at(0) + (dOmega/12.) * (55.0 * RealAmpF.at(3) - 59.0 * RealAmpF.at(2) + 37.0 * RealAmpF.at(1) - 9.0 * RealAmpF.at(0));
			AmplitudeP.at(1) = Amplitude.at(1) + (dOmega/12.) * (55.0 * ImagAmpF.at(3) - 59.0 * ImagAmpF.at(2) + 37.0 * ImagAmpF.at(1) - 9.0 * ImagAmpF.at(0)); 

			if( level_J.at(0) == 0 ) // Use the symmetry of even-even nuclei to speed things up a bit, just integrate substates between -j and 0
			{
				for(int i=0;i<LMax;i++)
				{
					for(int j=0;j<N_States;j++)
					{
						for(int k=MagneticSubstatesStart.at(j);k<MagneticSubstatesHalt.at(j);k++)
						{
							Int_t IR = k - 2*MagneticSubstatesCatalogue.at(k);
							AmplitudeP.at(0)[IR][i] = (Double_t)IFAC.at(j) * AmplitudeP.at(0)[k][i];
							AmplitudeP.at(1)[IR][i] = (Double_t)IFAC.at(j) * AmplitudeP.at(1)[k][i];
						}
					}
				}
			}
			Omega = Omega + dOmega + dOmega;

			AmplitudeDerivative = ComputeAmpDerivativeMatrices(AmplitudeP, Epsilon, Omega);

			Amplitude.at(0) = Amplitude.at(0) + (dOmega/12.0) * (9.0 * AmplitudeDerivative.at(0) + 19.0 * RealAmpF.at(3) - 5.0 * RealAmpF.at(2) + RealAmpF.at(1));
			Amplitude.at(1) = Amplitude.at(1) + (dOmega/12.0) * (9.0 * AmplitudeDerivative.at(1) + 19.0 * ImagAmpF.at(3) - 5.0 * ImagAmpF.at(2) + ImagAmpF.at(1));
	
			if( level_J.at(0) == 0 ) // Use the symmetry of even-even nuclei to speed things up a bit, just integrate substates between -j and 0
			{
				for(int i=0;i<LMax;i++)
				{
					for(int j=0;j<N_States;j++)
					{
						for(int k=MagneticSubstatesStart.at(j);k<MagneticSubstatesHalt.at(j);k++)
						{
							Int_t IR = k - 2.*MagneticSubstatesCatalogue.at(k);
							Amplitude.at(0)[IR][i] = (Double_t)IFAC.at(j) * Amplitude.at(0)[k][i];
							Amplitude.at(1)[IR][i] = (Double_t)IFAC.at(j) * Amplitude.at(1)[k][i];
						}
					}
				}
			}			
		
			AmplitudeDerivative = ComputeAmpDerivativeMatrices(Amplitude, Epsilon, Omega);

			// Shift all of the differential points along by one - Adams-Moulton only requires three points
			RealAmpF.at(0) = RealAmpF.at(1);
			ImagAmpF.at(0) = ImagAmpF.at(1);
			RealAmpF.at(1) = RealAmpF.at(2);
			ImagAmpF.at(1) = ImagAmpF.at(2);	
			RealAmpF.at(2) = RealAmpF.at(3);
			ImagAmpF.at(2) = ImagAmpF.at(3);	
			RealAmpF.at(3) = AmplitudeDerivative.at(0);
			ImagAmpF.at(3) = AmplitudeDerivative.at(1);	

			// Check amplitudes to see if the stepsize needs changing
			if(Omega + dOmega <= Up)
			{
				Double_t FF=0;
				for(int i=0;i<LMax;i++)
				{
					for(unsigned int j=0;j<MagneticSubstatesCatalogue.size();j++)
					{
						Double_t FZR = AmplitudeP.at(0)[j][i] - Amplitude.at(0)[j][i];
						Double_t FZI = AmplitudeP.at(1)[j][i] - Amplitude.at(1)[j][i];
						Double_t FZ = TMath::Sqrt(pow(FZR,2) + pow(FZI,2)) / 14.;
						if(FZ > FF) FF = FZ;

					}
				}
				if(FF <= Accuracy_50)
				{
					dOmegaFlag = false;
					dOmega = 2*dOmega;
					stepfile << "At Omega: " << Omega << " Stepwidth doubled to: " << 2*dOmega << "\n";
					if(verbose && PrintFlag)
						printf("At Omega: %f, Stepwidth doubled to: %f\n",Omega,2*dOmega);
				}
				if(FF > Accuracy)
				{
					dOmegaFlag = false;
					dOmega = 0.5*dOmega;
					stepfile << "At Omega: " << Omega << " Stepwidth halved to: " << 2*dOmega << "\n";
					if(verbose && PrintFlag)
						printf("At Omega: %f, Stepwidth halved to: %f\n",Omega,2*dOmega);
				}				
			}

			// Check excitation probabilities
			SubStateProbabilities.Clear();
			SubStateProbabilities.ResizeTo(MagneticSubstatesCatalogue.size());
			for(unsigned int i=0;i<MagneticSubstatesCatalogue.size();i++)
			{
				for(int j=0;j<LMax;j++)
				{
					SubStateProbabilities[i] = SubStateProbabilities[i] + 2.0*(pow(Amplitude.at(0)[i][j],2) + pow(Amplitude.at(1)[i][j],2)) / NormalisationFactor;
				}
				if(IntegerSpin == 0)
				{
					SubStateProbabilities[i] = SubStateProbabilities[i] - (pow(Amplitude.at(0)[i][LMax-1],2) + pow(Amplitude.at(1)[i][LMax-1],2)) / NormalisationFactor;
				}
			}
			StateProbabilities.Clear();
			StateProbabilities.ResizeTo(N_States);
			for(int i=0;i<N_States;i++)
			{
				for(int j=MagneticSubstatesStart.at(i);j<MagneticSubstatesStop.at(i);j++)
				{
					StateProbabilities[i] = StateProbabilities[i] + SubStateProbabilities[j];
				}
			}

			TotalProbability = 0;
			for(int i=0;i<N_States;i++)
				TotalProbability = TotalProbability + StateProbabilities[i];

			if(abs(TotalProbability - 1) > abs(ABW))
				ABW = TotalProbability - 1.0;

		}

	}

	// The below elements just write all of the calculated values to our output file, allowing for comparison

	stepfile << "Integration Complete!!!\n\n";

	std::vector< std::complex<Double_t> > Q_Matrix = CollisionFunction(2,Epsilon,Omega);
	stepfile << "Up: " << Omega << ", Epsilon: " << Epsilon << ", Q1: " << std::real(Q_Matrix.at(0)) << ", Q2: " << std::imag(Q_Matrix.at(1)) << ", Q3: " << std::real(Q_Matrix.at(2)) << "\n\n";

	stepfile << "Largest deviation from 1: " << ABW << ", Final deviation from 1: " << TotalProbability - 1 << "\n";
	
	stepfile << "\n" << std::setw(6) << "Spin:\t" << std::setw(3) << "M:\t\t" << std::setw(15) << "Real Amplitude:\t\t" << std::setw(15) << "Imag. Amplitude:\t" << std::setw(15) << "Population:\n";
	for(int i=0;i<N_States;i++)
	{
		for(int j=MagneticSubstatesStart.at(i);j<MagneticSubstatesStop.at(i);j++)
		{
			for(int k=0;k<LMax;k++)
			{
				Double_t pop = pow(Amplitude.at(0)[j][k],2) + pow(Amplitude.at(1)[j][k],2);

				stepfile << std::setw(6) << level_J.at(i) << "\t" << std::setw(3) << MagneticSubstatesCatalogue.at(j) << "\t\t" << std::setw(15) << Amplitude.at(0)[j][k] << "\t\t" << std::setw(15) << Amplitude.at(1)[j][k] << "\t\t" << std::setw(15) << pop << "\n";

			}
		}
	}

	stepfile <<  std::setw(20) << "Scattering angle in CM system: " << theta << " degree(s) \n";
	stepfile <<  std::setw(12) << "Level index:" << "\t " <<  std::setw(23) << "Excitation probability: " << "\t" <<  std::setw(24) << "CM cross-section (b/sr):\n";
	for(int i=0;i<N_States;i++)
	{
		Double_t CS = reaction->EvalRutherfordLevel(theta,level_E.at(i)) * StateProbabilities[i];
		stepfile << std::setw(12) << (i+1) << "\t" << std::setw(23) << StateProbabilities[i] << "\t" << std::setw(24) << CS << "\n";
	}

	// Push the probabilities and return the stateprobabilities

	SubStateProbabilitiesTheta.push_back(SubStateProbabilities);

	FinalRealAmplitude.push_back(Amplitude.at(0));
	FinalImagAmplitude.push_back(Amplitude.at(1));

	return StateProbabilities;

}

//******************************************************//
//	This function evaluates the Amplitude Derivatives
// 	at a given point in the integration process.
// 	The amplitude derivative function is defined in
//	the GOSIA manual, Eq. 3.18
//******************************************************//

std::vector<TMatrixD> TCLX::ComputeAmpDerivativeMatrices(std::vector<TMatrixD> AmplitudeIn, Double_t Epsilon, Double_t Omega)
{

	Double_t RAlfa = Epsilon * TMath::SinH(Omega) + Omega; // Here we calculate the exponent from Eq. 3.18 in the GOSIA manual

	std::vector<TMatrixD> AmpDot; // Prepare the vector of TMatrices which we will pass back to the integration process
	AmpDot.resize(2);
	AmpDot.at(0).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	AmpDot.at(1).ResizeTo(MagneticSubstatesCatalogue.size(),LMax);

	Int_t MR = 0, NZ = 0;
	for(int i=0;i<MaxLambda;i++) // Loop over multipolarities
	{

		if(GetMaxMatrix(TransitionMatrix.at(i)) == 0) // If the transition matrix element is zero, skip
			continue;

		std::vector< std::complex<Double_t> > Q_Matrix = CollisionFunction((i+1),Epsilon,Omega); // Grab the collision functions for this lambda, omega, epsilon combination

		for(int j=0;j<N_States;j++) // Loop over all states
		{
			for(int k=0;k<Lambda_Level[i][j];k++) // For this (final) state, loop over all connections
			{
				Int_t State = MatrixElementPairs.at(i)[j][k]; // Grab connected state (initial)

				if(State == -1) // If there is no state, we can skip
					continue;

				Double_t Xi = Eta[j] - Eta[State]; // Determine Xi parameter between these states

				Double_t RealEx = TMath::Cos(Xi * RAlfa); // Here we calculate the real part of the exponent from Eq. 3.18 of the GOSIA manual
				Double_t ImagEx = TMath::Sin(Xi * RAlfa); // ... and the imaginary part
	
				for(int l=MagneticSubstatesStart.at(j);l<MagneticSubstatesHalt.at(j);l++) // Loop over the substates of final state
				{

					for(int m=0;m<MRange[MR][0];m++) // For this substate determine which are connected
					{

						Int_t InitialSubstate = MRange[MR][1] + m; // Grab initial substate (valid within multipolarity conditions)
						
						Int_t MuAbs = (Int_t)(abs(MagneticSubstatesCatalogue.at(InitialSubstate) - MagneticSubstatesCatalogue.at(l))); // Determine ang. mom. difference between initial and final substates

						Double_t Real_RC = (std::real(Q_Matrix.at(MuAbs))*RealEx - std::imag(Q_Matrix.at(MuAbs))*ImagEx) * Zeta[NZ]; // Calculate the real initial part
						Double_t Imag_RC = (std::real(Q_Matrix.at(MuAbs))*ImagEx + std::imag(Q_Matrix.at(MuAbs))*RealEx) * Zeta[NZ]; // ... and the imaginary initial part
						NZ++;

						for(int n=0;n<LMax;n++)
						{
							Double_t RRC = Real_RC * AmplitudeIn.at(0)[InitialSubstate][n] - Imag_RC * AmplitudeIn.at(1)[InitialSubstate][n]; // Determine the differences
							Double_t IRC = Real_RC * AmplitudeIn.at(1)[InitialSubstate][n] + Imag_RC * AmplitudeIn.at(0)[InitialSubstate][n];
							AmpDot.at(0)[l][n] = AmpDot.at(0)[l][n] + IRC; // ... and push them back into the vector
							AmpDot.at(1)[l][n] = AmpDot.at(1)[l][n] - RRC;

						}		
					}
					MR++;
				}
			}
		}
	}

	if(NZ != ZetaCounter)
		printf("NZ not equal to ZetaCounter\n");

	return AmpDot; // Push back the amplitude differential

}

//********************************************************//
// 	Calculates collision functions, dependent upon
//	lambda, epsilon and omega for the possible values
//	of mu. This is returned as a vector of complex
//	numbers. Even mu values correspond to real 
//	collision functions, odd mu values correspond
//	to imaginary collision functions.
//********************************************************//

std::vector< std::complex<Double_t> > TCLX::CollisionFunction(Int_t Lambda, Double_t Epsilon, Double_t Omega)
{

	std::vector< std::complex<Double_t> > output; // The vector of complex numbers which we will be filling

	// Define the Epsilon and Omega dependent variables, as described in Table 2.1a of the GOSIA manual
	Double_t a = TMath::CosH(Omega) + Epsilon; 
	Double_t b = Epsilon * TMath::CosH(Omega) + 1.;
	Double_t c = TMath::Sqrt(pow(Epsilon,2) - 1.) * TMath::SinH(Omega);

	// Temporary real and imaginary values
	Double_t Real; 
	Double_t Imag;

	// Polarisation factor to account for dipole polarisation contributions
	Double_t PolFactor = (1. - E1PolFactor/b);


	// E1 Collision Functions:
	if(Lambda == 1)
	{
		
		output.resize(2);
		Real = 0.5 * (a / pow(b,2) );
		Imag = 0.0;		
		output.at(0) = std::complex<Double_t>(Real,Imag);

		Real = 0.0;
		Imag = -(1.0 / (2.0 * TMath::Sqrt(2) ) ) * c / pow(b,2);
		output.at(1.0) = std::complex<Double_t>(Real,Imag);
		

	}

	// E2 Collision Functions:
	if(Lambda == 2)
	{
		
		output.resize(3);
		Real = (3./4.) * (2.*pow(a,2) - pow(c,2))/pow(b,4) * PolFactor;
		Imag = 0.;
		output.at(0) = std::complex<Double_t>(Real,Imag);

		Real = 0.;
		Imag = -((3. * TMath::Sqrt(3))/(2. * TMath::Sqrt(2) )) * (a * c)/pow(b,4) * PolFactor;
		output.at(1) = std::complex<Double_t>(Real,Imag);

		Real = -((3. * TMath::Sqrt(3)) / (4. * TMath::Sqrt(2) )) * pow(c,2)/pow(b,4) * PolFactor;
		Imag = 0.;		
		output.at(2) = std::complex<Double_t>(Real,Imag);
		

	}

	// E3 Collision Functions:
	if(Lambda == 3)
	{
		
		output.resize(4);
		Real = (15./8.) * a*(2*pow(a,2) - 3*pow(c,2))/pow(b,6);
		Imag = 0.;		
		output.at(0) = std::complex<Double_t>(Real,Imag);

		Real = 0.;
		Imag = -((15. * TMath::Sqrt(3))/16.) * c * (4 * pow(a,2) - pow(c,2))/pow(b,6);
		output.at(1) = std::complex<Double_t>(Real,Imag);

		Real = -((15. * TMath::Sqrt(15)) / (8 * TMath::Sqrt(2) )) * a*pow(c,2)/pow(b,6);
		Imag = 0.;		
		output.at(2) = std::complex<Double_t>(Real,Imag);

		Real = 0.;
		Imag = ((15. * TMath::Sqrt(5))/16.) * pow(c,3) / pow(b,6);
		output.at(3) = std::complex<Double_t>(Real,Imag);
		

	}	

	// E4 Collision Functions:
	if(Lambda == 4)
	{
		
		output.resize(5);
		Real = (35/32) * (8*pow(a,4) - 24*pow(a,2)*pow(c,2) + 3*pow(c,4)) / pow(b,8);
		Imag = 0.;		
		output.at(0) = std::complex<Double_t>(Real,Imag);

		Real = 0.;
		Imag = -((35 * TMath::Sqrt(5))/16) * a*c * (4 * pow(a,2) - 3*pow(c,2))/pow(b,8);
		output.at(1) = std::complex<Double_t>(Real,Imag);

		Real = -((35 * TMath::Sqrt(5)) / (16 * TMath::Sqrt(2) )) * pow(c,2) * (6*pow(a,2) - pow(c,2)) / pow(b,8);
		Imag = 0.;		
		output.at(2) = std::complex<Double_t>(Real,Imag);

		Real = 0.;
		Imag = ((35 * TMath::Sqrt(35))/16) * a * pow(c,3) / pow(b,8);
		output.at(3) = std::complex<Double_t>(Real,Imag);

		Real = (35 * TMath::Sqrt(35))/(32 * TMath::Sqrt(2)) * pow(c,4)/pow(b,8);
		Imag = 0.;		
		output.at(4) = std::complex<Double_t>(Real,Imag);
		

	}	

	// This will be extended. In addition, magnetic transitions will eventually be included.
	if(Lambda > 4)
		printf("Error: Cannot process multipolarities greater than E4");

	return output;

}
//*******************************************************//
//	Calculate Electric Dipole Normalization value 
// 	this accounts for E1 strength. See GOSIA manual
// 	Eqs. 3.25 and 3.26 for description.
//*******************************************************//

Double_t TCLX::ElectricDipoleNormalization()
{

	Double_t factor;
	
	factor = TReactions::dipole * (reaction->Elab * reaction->Tar_A) / (pow(reaction->Proj_Z,2) * (1 + reaction->Proj_A/reaction->Tar_A));

	return factor;

}
