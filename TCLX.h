#ifndef TCLX_h
#define TCLX_h

#include <TSpline.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <complex>
#include <TH1.h>
#include <TMatrix.h>
#include <TMatrixD.h>
#include <TVector.h>
#include <TVectorD.h>
#include <TSystem.h>
#include <Math/SpecFunc.h>
#include <TCanvas.h>

#include "TCLX_DataInput.h" // The TCLX_DataInput class grabs the data from the input file, in order to be processed by the TCLX class.
#include "TReactions.h" // The TReactions class includes a load of useful reactions (duh) information, such as Rutherford cross sections, as well as a load of information defining the nuclear properties


//**********************************************//
//	The TCLX class enables the integration of
//	Coulomb excitation amplitudes, based on 
//	user inputs of nuclear information, such as
//	matrix elements. It is based (and in many
//	parts, directly ported) on the clx fortran
//	code produced by H. Ower et al.
//
//**********************************************//

class TCLX {

	private :

	public :
	
		TCLX(){;} 
		virtual ~TCLX(){;}

		ofstream stepfile; // This output file is where we are going to print all of the steps of the integration, allowing for you to check things more carefully than just grabbing the cross-sections

		std::vector<TVectorD> SubStateProbabilitiesTheta; // A std::vector, one entry for each step in theta, which contains the excitation probabilities for each individual substate

		bool verbose, debug; // Flags to be defined by input options - at the moment, debug is not used

		void WriteMatrix(std::ofstream& outfile, TMatrixD Matrix); // A function which prints a matrix in ASCII format to a predefined outfile
	
		std::vector<TMatrixD> TransitionMatrix; // Vector of Transition TMatrices, for lambda = 1 -> max lambda
		Int_t N_States; // Number of states
		Int_t MaxLambda; // Maximum lambda value

		Double_t Theta_Min, Theta_Max, Theta_Step; // These are the angles we are going to integrate between, and the steps we are going to increment in
		
		void GrabData(const char* s, Option_t *opt="", Option_t *opt_bt="b"); // Grab data from file via TCLX_DataInput class

		void GetMatrices(TCLX_DataInput *input); // Load the vector of TransitionMatrices from input
		void GetNuclei(TCLX_DataInput *input, Option_t *opt="", Option_t *opt_bt="b");   // Load the nuclei into TReactions from input
		void GetKinematics(TCLX_DataInput *input); // Load the kinematics into TReactions from input
		void GetOtherVariables(TCLX_DataInput *input); // Load other useful variables

		bool BeamExcitation; // Are we dealing with BeamExcitation?

		TReactions *reaction; // Reaction information, including nuclei and kinematics

		Double_t GetMaxMatrix(TMatrixD max); // Grabs the maximum value of a TMatrixD

		void CheckMatrix();	// Check's symmetry of TransitionMatrix

		void PrintCalculationDetails(); // Print the details of the calculation input file in some semblence of English, to check things stack up

		Double_t ElectricDipoleNormalization(); // Calculates E1PolFactor
		Double_t E1PolFactor; // Factor which accounts for dipole polarisation

		TMatrixD Lambda_Level; // Matrix showing the number of matrix elements of a given multipolarity populating each level
		void ConstructLDNUMMatrix(); // Fill the above matrix
		Int_t MaxStateConnections; // Maximum number of matrix elements connected to a single state

		std::vector<TMatrixD> MatrixElementPairs; // Matrices showing which other states a given state is "connected" to
		std::vector<TMatrixD> MatrixElementsByConnections; // Matrices containing the matrix elements corresponding to the "connections" in MatrixElementPairs
		void ConstructMatrixElementPairsMatrices(); // Construct MatrixElementPairs and MatrixElementsByConnections

		TVectorD Eta; // Vector containing Eta (similar to GOSIA manual, Eq. 3.19) values for each state (combined to give difference)
		void SetEta(); // Creates the above vector

		Int_t LMax; // Maximum L - based on ground-state spin (i.e. J=0 gs, LMax = 1)

		// Nuclear level input properties, energy, J and parity.
		std::vector<Double_t> level_E;
		std::vector<Double_t> level_J;
		std::vector<Int_t> level_P;
		std::vector<Int_t> IFAC;

		TMatrixD XiMatrix; // Xi matrix gives the full Eta values (i.e. diff. between TVectorD Eta) for any two states
		void MakeXiMatrix(); // Make the above matrix
		Double_t XiMax; // The maximum value of the Xi matrix

		std::vector<TMatrixD> PsiMatrix; // Psi matrix (GOSIA manual, Eq. 3.21)
		void MakePsiMatrix(); // Make Psi matrix
		Double_t PsiMax; // Maximum value of Psi

		// These vectors contain the magnetic substate information
		std::vector<Double_t> MagneticSubstatesCatalogue; // Contains all of the magnetic substates
		std::vector<Int_t> MagneticSubstatesStart; // Gives the FIRST substate catalogue value for any given STATE
		std::vector<Int_t> MagneticSubstatesStop; // Gives the LAST substate catalogue value for any given state
		std::vector<Int_t> MagneticSubstatesHalt; // Used to speed things up with symmetry arguments
		void MakeMagneticSubstates(); // Fill the above vectors

		TMatrixD MRange; // Matrix defining allowed combinations of substates for the given transition multipolarity
		void MakeMRange(); // Make the above

		Int_t ZetaCounter; // Count Zeta values
		TVectorD Zeta; // Zeta contains all of the coupling information for the allowed combinations of magnetic substates as defined in MRange
		void MakeZeta(); // Make the Zeta matrix

		void IntegrateTheta(); // Function which integrates over theta by calling...
		TVectorD Integration(Double_t Theta); // Integrates over time (in the form of omega, see GOSIA manual) for a given value of theta

		std::vector< std::complex<Double_t> > CollisionFunction(Int_t Lambda, Double_t Epsilon, Double_t Omega); // Determine the collision functions (see GOSIA manual, Table 2.1a)

		std::vector<TMatrixD> ComputeAmpDerivativeMatrices(std::vector<TMatrixD> Amplitude, Double_t Epsilon, Double_t Omega); // Eq. 3.18 of the GOSIA manual, used to perform the integration with numerical differentiation methods

		Double_t NormalisationFactor; // Factor for normalisation
		Double_t IntegerSpin; // Are we using integer spin?

		std::vector<TVectorD> Probabilities; // These are the excitation probabilities

		std::vector<TMatrixD> FinalRealAmplitude; // Final real and imaginary amplitudes at the end of the integration process
		std::vector<TMatrixD> FinalImagAmplitude;

		//****************************************//
		//	Master function which prepares all
		//	matrices.
		//****************************************//
		
		void PrepareMatrices();

		//****************************************//
		//	The below functions take the outputs
		// 	from our integration and print/process
		//	them into a useful format, as defined
		// 	by the user.
		//****************************************//

		TSpline3* StateProbability(Int_t State); // Creates a TSpline3 containing the probabilities vs theta for a given state

		TGraph* StateCrossSection(Int_t State);	// Creates a TSpline3 containing the cross-sections (b/sr) vs theta for a given state

		TSpline3* Rutherford(); // Grabs the Rutherford from TReactions and creates a TSpline3

		TReactions* GetReaction() { return reaction; } // Grab the TReaction being used in the TCLX class - this is particularly useful if you want to use it externally

		void PrintCrossSections(const char *filename); // Print our resultant cross sections (and probabilities) for all states

		//****************************************//
		//	The below functions are intended to 
		//	set up the variables without reading
		//	in an input file, by taking inputs
		//	from the terminal.
		//
		//	Note that this may be slightly less
		//	reliable, not tested properly yet.
		//****************************************//
		
		void SetupCLXReaction();
		void SetNStates(Int_t n) { N_States = n ;}
		void SetStateVariables();
		void SetMatrixElements();
		void SetReactionInformation();
		void SetThetaRange();

		//****************************************//
		//	The below functions are intended to 
		//	set up the variables without reading
		//	in an input file. These functions are
		//	more suited to being used as part of
		//	a macro, for example.
		//
		//	Note that this may be slightly less
		//	reliable, not tested properly yet.
		//****************************************//

		void SetNStatesMaxLambda(Int_t n, Int_t maxlambda); // This will clear and reset ALL information, if you just want to reset a single value, use one of the following functions ALONE
		void SetStateEnergy(Int_t i, Double_t E) { level_E.at(i) = E; }
		void SetStateParity(Int_t i, Int_t P) { level_P.at(i) = P; }
		void SetStateJ(Int_t i, Double_t J) { level_J.at(i) = J; }
		void SetMatrixElement(Int_t lambda, Int_t i, Int_t j, Double_t ME) { TransitionMatrix.at(lambda)[i][j] = ME; }
		void SetThetaRange(Double_t thetamin, Double_t thetamax, Double_t thetastep)
		{
			Theta_Min = thetamin;
			Theta_Max = thetamax;
			Theta_Step = thetastep;
		}
		void SetupCLXReaction(Int_t PA, Int_t PZ, Int_t TA, Int_t TZ, Double_t LabE)
		{
			reaction->SetNuclei(PZ,PA,TZ,TA);
			reaction->SetElab(LabE);
		}
		void SetupVariables(); // This needs to be run after setting up the information and before performing the integration

		//ClassDef(TCLX,1);

};

#endif

#ifdef TCLX_cxx

//******************************//
//	Grabs the data from a
//	define file. Takes options
// 	which allow for verbose/debug 
//	mode.
//******************************//

void TCLX::GrabData(const char* filename, Option_t *opt, Option_t *opt_bt)
{

	// Check options and change flags accordingly
	verbose = false;
	debug = false;
	if(strncmp(opt,"v",1)==0)
		verbose = true;
	if(strncmp(opt,"d",1)==0){
		debug = true;
		verbose = true;
	}

	// Make TCLX_DataInput and read the file
	TCLX_DataInput *inputdata = new TCLX_DataInput();
	inputdata->ReadFile(filename,opt);

	if(verbose)
		printf("Data grabbed\n");

	// Grab the matrices from the defined TCLX_DataInput
	GetMatrices(inputdata);

	// Initialise our TReactions, now we can start creating nuclei
	reaction = new TReactions();

	// Grab the nuclei, the kinematics and other variables from the TCLX_DataInput
	GetNuclei(inputdata,opt,opt_bt);
	GetKinematics(inputdata);
	GetOtherVariables(inputdata);

	// Make the dipole polarization value
	E1PolFactor = ElectricDipoleNormalization();

	if(verbose)
		printf("Got data\n");


}

//*****************************//
// Grab assorted other variables 
// from the input file - of 
// varying importance
//*****************************//

void TCLX::GetOtherVariables(TCLX_DataInput *input)
{

	// Grab theta values from input
	Theta_Min = input->Theta_Start;
	Theta_Max = input->Theta_End;
	Theta_Step = input->Theta_Step;

	// Grab nuclear level information from input
	level_E = input->level_E;
	level_J = input->level_J;
	level_P = input->level_P;

	// Fill IFAC vector - defines how symmetry short-cuts are performed
	IFAC.resize(level_J.size());
	for(unsigned int i=0;i<level_J.size();i++)
	{
		Int_t IDPAR = 0;
		if(level_P.at(i) != level_P.at(0)) IDPAR = 1;
		IFAC.at(i) = pow((-1),(IDPAR + Int_t(level_J.at(0)-level_J.at(i)))); 

	}

	// Check whether we're dealing with a J=0 ground state, if not we have to consider other substates too
	LMax = level_J.at(0) + 1;

	// Define a few factors
	NormalisationFactor = 2.0 * level_J.at(0) + 1.0;
	IntegerSpin = (Int_t)(2*level_J.at(0)) % 2;

}

//*************************//
//	Grabs the matrices from
//	the input file. Nothing
//	fancy.
//*************************//

void TCLX::GetMatrices(TCLX_DataInput *input)
{
	
	MaxLambda = input->GetNLambda();
	N_States = input->GetNStates();
	TransitionMatrix.resize(MaxLambda);
	for(unsigned int i=0;i<TransitionMatrix.size();i++) TransitionMatrix.at(i).ResizeTo(N_States,N_States);
	TransitionMatrix = input->GetTransitionMatrix();
	
}

//******************************//
//	Grabs the nuclei information
//	from the input file and puts
//	it in our newly initialised
//	TReactions for use later.
//******************************//

void TCLX::GetNuclei(TCLX_DataInput *input, Option_t *opt, Option_t *opt_bt)
{

	verbose = false;
	if(strncmp(opt,"v",1)==0)
		verbose = true;

	if(strncmp(opt_bt,"b",1)==0)
		BeamExcitation = true;
	else if(strncmp(opt_bt,"t",1)==0)
		BeamExcitation = false;
	else
		printf("Not valid beam/target excitation option\n");	
	
	if(!reaction)
		reaction = new TReactions();

	Int_t temp_ZP;
	Int_t temp_AP;
	Int_t temp_ZT;
	Int_t temp_AT;

	if(BeamExcitation)
	{
		temp_ZP = input->ZP;
		temp_AP = input->AP;
		temp_ZT = input->ZT;
		temp_AT = input->AT;
	}
	else if(!BeamExcitation)
	{
		temp_ZT = input->ZP;
		temp_AP = input->AP;
		temp_ZP = input->ZT;
		temp_AT = input->AT;
	}
	reaction->SetNuclei(temp_ZP,temp_AP,temp_ZT,temp_AT);	

	if(verbose)
		printf("Projectile Z: %i, Projectile A: %i, Target Z: %i, Target A: %i\n",reaction->GetNuclei(3),reaction->GetNuclei(2),reaction->GetNuclei(1),reaction->GetNuclei(0));

}

//******************************//
//	Pretty simple right now, pass
//	kinematics information to the
//	TReactions class.
//******************************//


void TCLX::GetKinematics(TCLX_DataInput *input)
{

	if(!reaction)
		reaction = new TReactions();

	reaction->SetElab(input->EP);

}

//******************************************//
//	Create a text output file (which will
//	be overwritten) and write the various
//	calculation details into it. This will
//	hopefully be done in a reasonably
//	sensible manner, so you can easily
//	identify any mistakes made in the input
//******************************************//

void TCLX::PrintCalculationDetails()
{

	ofstream outfile;
	outfile.open("TCLX_CalculationDetails.txt");

	outfile << "Proj. Z:\tProj. A:\tTar. Z:\t\tTar. A:\t\tE Lab (MeV):\tEcm (MeV):\n";

	outfile << reaction->GetNuclei(3) << "\t\t\t" << reaction->GetNuclei(2) << "\t\t\t" << reaction->GetNuclei(1) << "\t\t\t" << reaction->GetNuclei(0) << "\t\t\t" << reaction->GetElab() << "\t\t\t\t" << reaction->GetEcm() << "\n";

	outfile << "\nTheta min: " << (Double_t)Theta_Min << " (degrees), Theta max: " << (Double_t)Theta_Max << " (degrees), in steps of " << (Double_t)Theta_Step << " degrees\n";

	outfile << "\nClosest separation distance (fm): " << reaction->ClosestApproach() << "\n";

	outfile << "\nSommerfield Parameter (g.s.): " << reaction->EtaCalc(0) << "\n";

	outfile << "\nElectric dipole polarization factor: " << ElectricDipoleNormalization() << "\n";

	outfile << "\nMatrix Elements (empty matrices not shown):\n";

	outfile << "\nEta Matrix:\n";

	for(int i=0;i<Eta.GetNrows();i++)
		outfile << i << "\t" << Eta[i] << "\n";

	outfile << "\nLargest Xi: " << XiMax << "\n";

	outfile << "\nXi matrix:\n";
	WriteMatrix(outfile,XiMatrix);
	
	for(int l=0;l<MaxLambda;l++)
	{

		if(GetMaxMatrix(TransitionMatrix.at(l))==0)
			continue;
		outfile << "\nLambda: " << l+1 << "\n";

		WriteMatrix(outfile,TransitionMatrix.at(l));
		
	}

	outfile << "\nLDNUM matrix (Transition elements connected to each state):\n";
	WriteMatrix(outfile,Lambda_Level);

	outfile << "\nConnections matrices (rows indicate states, elements indicate other states that are connected to this state):\n";
	
	for(int l=0;l<MaxLambda;l++)
	{

		if(GetMaxMatrix(MatrixElementPairs.at(l))==0)
			continue;
		outfile << "\nLambda: " << l+1 << "\n";

		WriteMatrix(outfile,MatrixElementPairs.at(l));
		
	}

	outfile << "\nConnections matrix elements (rows indicate states, elements contain the matrix element for the connection shown in connections matrices):\n";
	
	for(int l=0;l<MaxLambda;l++)
	{

		if(GetMaxMatrix(MatrixElementsByConnections.at(l))==0)
			continue;
		outfile << "\nLambda: " << l+1 << "\n";

		WriteMatrix(outfile,MatrixElementsByConnections.at(l));
		
	}


	outfile << "\nPsi Matrix:\n";
	for(int l=0;l<MaxLambda;l++)
	{

		if(GetMaxMatrix(PsiMatrix.at(l))==0)
			continue;
		outfile << "\nLambda: " << l+1 << "\n";

		WriteMatrix(outfile,PsiMatrix.at(l));
		
	}

	outfile << "\nZeta Matrix:\n";
	for(int i=0;i<Zeta.GetNrows();i++)
		outfile << Zeta[i] << "\n";

	outfile << "\n\nMagnetic substates catalogue:\n" << std::endl;
	for(int i=0;i<N_States;i++)
	{
		outfile << "State " << i << ", J=" << level_J.at(i) << ":\t";
		for(int j=MagneticSubstatesStart.at(i);j<MagneticSubstatesStop.at(i);j++)
		{
			outfile << MagneticSubstatesCatalogue.at(j) << ",\t";
		}
		outfile << "\n";
	}

	outfile << "Substate range matrix:\n";
	WriteMatrix(outfile,MRange);
	outfile << "\n";

	TMatrixD realamp,imagamp;
	realamp.ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	imagamp.ResizeTo(MagneticSubstatesCatalogue.size(),LMax);
	for(unsigned int i=0;i<FinalRealAmplitude.size();i++)
	{
	
		outfile << "Theta: " << i*Theta_Step + Theta_Min << "\n";
		outfile << std::setw(12) << "Substate:\t\t" << std::setw(12) << "RealAmp:\t\t" << std::setw(12) << "ImagAmp:\n";

		for(unsigned int j=0;j<MagneticSubstatesCatalogue.size();j++)
		{


			for(int k=0;k<LMax;k++)
			{
				outfile << std::setw(12) << MagneticSubstatesCatalogue.at(j) << "\t\t";

				realamp[j][k] = realamp[j][k] + FinalRealAmplitude.at(i)[j][k];
				imagamp[j][k] = imagamp[j][k] + FinalImagAmplitude.at(i)[j][k];

				outfile << std::setw(12) << FinalRealAmplitude.at(i)[j][k] << "\t\t" << std::setw(12) << FinalImagAmplitude.at(i)[j][k] << "\n";

			}
		}
	}
	for(unsigned int j=0;j<MagneticSubstatesCatalogue.size();j++)
	{
		for(int k=0;k<LMax;k++)
		{
			realamp[j][k] = realamp[j][k] / (Double_t)FinalRealAmplitude.size();
			imagamp[j][k] = imagamp[j][k] / (Double_t)FinalRealAmplitude.size();
		}
	}

	outfile << "\nFinal amplitudes:\n\n";
	outfile << std::setw(6) << "SPIN:\t" << std::setw(6) << "M:\t" << std::setw(15) << "REAL AMP.:\t" << std::setw(15) << "IMAG AMP.:\n";
	for(int i=0;i<N_States;i++)
	{
		for(int j=MagneticSubstatesStart.at(i);j<MagneticSubstatesStop.at(i);j++)
		{			
			outfile << std::setw(6) << level_J.at(i) << "\t" << std::setw(6) << MagneticSubstatesCatalogue.at(j) << "\t" << std::setw(15) << realamp[j][0] << "\t" << std::setw(15) << imagamp[j][0] << "\n";
		}
	}

	Double_t Amplitudes[SubStateProbabilitiesTheta.at(0).GetNrows()];
	for(unsigned int i=0;i<SubStateProbabilitiesTheta.size();i++)
	{
		for(int j=0;j<SubStateProbabilitiesTheta.at(i).GetNrows();j++)
		{
			Amplitudes[j] = Amplitudes[j] + SubStateProbabilitiesTheta.at(i)[j];
		}
	}	

	outfile << "\nSubstate probabilities:\n\n";
	for(int i=0;i<N_States;i++)
	{
		for(int j=MagneticSubstatesStart.at(i);j<MagneticSubstatesStop.at(i);j++)
		{			
			outfile << std::setw(6) << level_J.at(i) << "\t" << std::setw(6) << MagneticSubstatesCatalogue.at(j) << "\t" << std::setw(15) << Amplitudes[j] << "\n";
		}
	}

	outfile.close();

}

//****************************************//
//	A quick function which writes a matrix
//	to a define outfile.
//****************************************//


void TCLX::WriteMatrix(std::ofstream& outfile, TMatrixD Matrix)
{

	if(!outfile.is_open())
		return;

	char line[2048]= "-------";
	for(int i=0;i<Matrix.GetNcols();i++) strcat(line,"-----------------");

	outfile << "\nState:\t";
	for(int i=0;i<Matrix.GetNcols();i++) outfile << "|" << std::setw(12) << i+1 << "\t";
	outfile << "|\n";

	outfile << line << "\n";

	for(int i=0;i<Matrix.GetNrows();i++)
	{
		outfile << "\t" << i+1 << "\t";
		for(int j=0;j<Matrix.GetNcols();j++)
		{
			outfile << "|" << std::setw(12) << Matrix[i][j] << "\t";

		}
		outfile << "|\n";

	}

}

//****************************************//
//	Grabs the maximum value of a matrix
//****************************************//

Double_t TCLX::GetMaxMatrix(TMatrixD mat){

	Double_t max = 0;
	for(int i=0;i<mat.GetNcols();i++){
		for(int j=0;j<mat.GetNrows();j++){
			if(mat[j][i] > max) max = mat[j][i];
		}
	}
	return max;

}

//****************************************//
//	Create a TSpline of probabilities
//	for a given state
//****************************************//

TSpline3* TCLX::StateProbability(Int_t State)
{

	if(State >= N_States || State < 0){
		printf("Error, state doesn't exist\n");
	}

	Double_t Probs[Probabilities.size()];
	Double_t Theta[Probabilities.size()];

	Double_t step = (Theta_Max - Theta_Min)/Probabilities.size();

	TVectorD Prob;	
	Prob.ResizeTo(N_States);

	for(unsigned int i=0;i<Probabilities.size();i++)
	{
		Prob = Probabilities.at(i);
		Probs[i] = (Double_t)Prob[State];
		Theta[i] = (Double_t)Theta_Min + i*step;
	}

	char sname[32];
	sprintf(sname,"Probability_%i",State);
	TSpline3 *ProbGraph = new TSpline3(sname,Theta,Probs,Probabilities.size());

	return ProbGraph;
	
}

//****************************************//
//	Create a TSpline of the Rutherford 
//	cross section
//****************************************//

TSpline3* TCLX::Rutherford()
{

	TSpline3 *ruth = reaction->Rutherford(Theta_Min, Theta_Max, Theta_Step);
	
	return ruth;

}

//*****************************************//
//	Create a TSpline of the cross-section
//	for a given state
//*****************************************//

TGraph* TCLX::StateCrossSection(Int_t State)
{

	if(State >= N_States || State < 0){
		printf("Error, state doesn't exist\n");
	}

	TVectorD Prob;	
	Prob.ResizeTo(N_States);

	Double_t CS[Probabilities.size()];
	Double_t Theta[Probabilities.size()];

	Double_t step = (Theta_Max - Theta_Min)/Probabilities.size();


	for(unsigned int i=0;i<Probabilities.size();i++)
	{
		Prob = Probabilities.at(i);		

		for(int j=0;j<N_States;j++)
		{
			Double_t ruth = reaction->EvalRutherfordLevel((Theta_Min+i*step),level_E.at(j));
			CS[i] = Prob[State] * ruth;
			Theta[i] = Theta_Min + i*step;
		}
	}

	char sname[32];
	sprintf(sname,"CrossSection_%i",State);
	TGraph *CSSpline = new TGraph(Probabilities.size(),Theta,CS);

	return CSSpline;
	
}

//*****************************************//
//	Print the state cross-sections and 
//	probabilities to file for all values
//	of theta.
//*****************************************//

void TCLX::PrintCrossSections(const char *filename)
{

	TVectorD Prob;	
	Prob.ResizeTo(N_States);

	Double_t CS[Probabilities.size()];
	Double_t Theta[Probabilities.size()];

	Double_t step = (Theta_Max - Theta_Min)/(Probabilities.size()-1);

	printf("Steps: %f\n",step);
	printf("Prob size: %i\n",Probabilities.size());

	ofstream outfile;
	outfile.open(filename);	

	outfile << std::setw(6) << "Theta\t";
	for(int i=0;i<N_States;i++)
		outfile << std::setw(12) << "State_" << i << "_prob.\t" << std::setw(12) << "State_" << i << std::setw(6) << "(b/sr)\t";
	outfile << "\n";

	for(unsigned int i=0;i<Probabilities.size();i++)
	{
		Prob = Probabilities.at(i);		

		outfile << Theta_Min + i*step << "\t\t";

		for(int j=0;j<N_States;j++)
		{

			Double_t ruth = reaction->EvalRutherfordLevel((Theta_Min+i*step),level_E.at(j));
			outfile << std::setw(12) << Prob[j] << "\t" << std::setw(12) << Prob[j] * ruth << "\t";

		}
		outfile << "\n";
	}

}

//****************************************//
//	The below functions are intended to 
//	set up the variables without reading
//	in an input file.
//
//	Note that this may be slightly less
//	reliable, not tested properly yet.
//****************************************//

void TCLX::SetupCLXReaction()
{

	Int_t n;
	std::cout << "Number of states: " << std::endl;
	std::cin >> n;

	SetNStates(n);

	SetStateVariables();
	SetMatrixElements();
	SetReactionInformation();
	SetThetaRange();

	SetupVariables();

	// Make the dipole polarization value
	E1PolFactor = ElectricDipoleNormalization();

}

void TCLX::SetStateVariables()
{

	level_E.clear();
	level_J.clear();
	level_P.clear();

	Double_t temp_E, temp_J;
	Int_t temp_P;

	for(int i=0;i<N_States;i++)
	{
		std::cout << "Input state " << (i+1) << " energy (MeV):\n";
		std::cin >> temp_E;
		std::cout << "Input state " << (i+1) << " J:\n";
		std::cin >> temp_J;
		std::cout << "Input state " << (i+1) << " parity (1 = +, 0 = -)\n";
		std::cin >> temp_P;

		level_E.push_back(temp_E);
		level_J.push_back(temp_J);
		level_P.push_back(temp_P);
	}

}

void TCLX::SetMatrixElements()
{

	MaxLambda = 0;

	Double_t temp_ME;
	
	TransitionMatrix.clear();
	std::cout << "Maximum multipolarity needed:" << std::endl;
	std::cin >> MaxLambda;

	TransitionMatrix.resize(MaxLambda);

	for(int i=0;i<MaxLambda;i++)
	{
		TransitionMatrix.at(i).ResizeTo(N_States,N_States);
		for(int j=0;j<N_States;j++)
		{
			for(int k=j;k<N_States;k++)
			{
				std::cout << "Matrix element between state: " << j+1 << " and state " << k+1 << " for multipolarity " << (i+1) << std::endl;
				std::cin >> temp_ME;
				TransitionMatrix.at(i)[j][k] = temp_ME;
				TransitionMatrix.at(i)[k][j] = temp_ME;
			}
		}
	}

}

void TCLX::SetReactionInformation()
{

	if(!reaction)
		reaction = new TReactions();

	Double_t Ebeam;
	std::cout << "Beam energy: " << std::endl;
	std::cin >> Ebeam;
	
	reaction->SetElab(Ebeam);

	Int_t TA,TZ,PA,PZ;

	std::cout << "Target A:" << std::endl;
	std::cin >> TA;
	std::cout << "Target Z:" << std::endl;
	std::cin >> TZ;
	std::cout << "Projectile A" << std::endl;
	std::cin >> PA;
	std::cout << "Projectile Z" << std::endl;
	std::cin >> PZ;

	reaction->SetNuclei(PZ,PA,TZ,TA);

}

void TCLX::SetThetaRange()
{

	Double_t mintheta,maxtheta,steptheta;
	std::cout << "Theta minimum: " << std::endl;
	std::cin >> mintheta;

	std::cout << "Theta maximum: " << std::endl;
	std::cin >> maxtheta;

	std::cout << "Theta in steps of: " << std::endl;
	std::cin >> steptheta;

	Theta_Min = mintheta;
	Theta_Max = maxtheta;
	Theta_Step = steptheta;

}

//****************************************//
//	The below functions are intended to 
//	set up the variables without reading
//	in an input file. These functions are
//	more suited to being used as part of
//	a macro, for example.
//
//	Note that this may be slightly less
//	reliable, not tested properly yet.
//****************************************//

void TCLX::SetNStatesMaxLambda(Int_t n, Int_t maxlambda)
{

	N_States = n;
	MaxLambda = maxlambda;

	level_E.clear();
	level_J.clear();
	level_P.clear();
	level_E.resize(N_States);
	level_J.resize(N_States);
	level_P.resize(N_States);

	TransitionMatrix.clear();
	TransitionMatrix.resize(MaxLambda);
	for(int i=0;i<MaxLambda;i++)
		TransitionMatrix.at(i).ResizeTo(N_States,N_States);

	if(!reaction)
		reaction = new TReactions();

}

// Sets up the other various variables required for the integration
void TCLX::SetupVariables()
{

	// Make the dipole polarization value
	E1PolFactor = ElectricDipoleNormalization();

	// Fill IFAC vector - defines how symmetry short-cuts are performed
	IFAC.resize(level_J.size());
	for(unsigned int i=0;i<level_J.size();i++)
	{
		Int_t IDPAR = 0;
		if(level_P.at(i) != level_P.at(0)) IDPAR = 1;
		IFAC.at(i) = pow((-1),(IDPAR + Int_t(level_J.at(0)-level_J.at(i)))); 

	}

	// Check whether we're dealing with a J=0 ground state, if not we have to consider other substates too
	LMax = level_J.at(0) + 1;

	// Define a few factors
	NormalisationFactor = 2.0 * level_J.at(0) + 1.0;
	IntegerSpin = (Int_t)(2*level_J.at(0)) % 2;

}


#endif
