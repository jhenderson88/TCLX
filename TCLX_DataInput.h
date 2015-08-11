#ifndef TCLX_DataInput_h
#define TCLX_DataInput_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <complex>
#include <TMatrix.h>
#include <TMatrixD.h>
#include <TVector.h>
#include <TVectorD.h>
#include <TSystem.h>

class TCLX_DataInput {

	private :

	public :

		TCLX_DataInput(){ 
			if(level_index.size() || trans_ME.size())
				ClearAll();
		}
		virtual ~TCLX_DataInput(){;}

		void ClearAll();

		void ReadFile(const char* s, Option_t *opt="");

		std::string title; // Experiment title

		Int_t MaxLambda; // Maximum lambda
		Int_t N_States; // Number of states
		Int_t ZP, AP, ZT, AT; // Bulk nucleus information
		Double_t EP; // Beam energy
		Double_t Theta_Start, Theta_End, Theta_Step; // Range of theta and step size

		std::vector<Int_t> level_index, level_P; // Level information
		std::vector<Double_t> level_J, level_E, level_K; // Level information

		std::vector<Int_t> trans_N, trans_M, trans_Lambda; // Transition information (initial state, final state, multipolarity)
		std::vector<Double_t> trans_ME; // Transition matrix elements
		
		std::vector<TMatrixD> TransitionMatrix;
		
		std::vector<TMatrixD> GetTransitionMatrix() { return TransitionMatrix; }
		
		Int_t GetNStates() { return N_States; }
		Int_t GetNLambda() { return MaxLambda; }
			
		//ClassDef(TCLX_DataInput,1);

};

#endif

#ifdef TCLX_DataInput_cxx

void TCLX_DataInput::ClearAll(){ 
	
	// Level information
	level_index.clear();
	level_P.clear();
	level_J.clear();
	level_E.clear();
	level_K.clear();
	
	// Transition information
	trans_N.clear();
	trans_M.clear();
	trans_Lambda.clear();
	trans_ME.clear();

}

void TCLX_DataInput::ReadFile(const char* filename, Option_t *opt){
	
	bool verbose = false;
	if(strncmp(opt,"v",1)==0)
		verbose = true;
	
	std::ifstream infile;
	infile.open(filename);
	
	if(!infile.is_open()){
		printf("Error: File not opened\n");
		return;
	}
	
	if(verbose)
		printf("Clearing vectors\n");

	ClearAll();

	Int_t i=0;

	Int_t temp_index, temp_P;
	Double_t temp_J, temp_E, temp_K;

	Int_t temp_N, temp_M, temp_Lambda;
	Double_t temp_ME;
	
	if(verbose)
		printf("Reading input file:\n");

	std::string line;
	while(getline(infile,line))
	{

		std::stringstream linestream(line);

		if(line.length() == 0 || (line[0] == '/' && line[1] == '/'))
			continue;

		if(verbose)
			std::cout << line << std::endl;

		if(i==0) title = line;

		if(i==1) linestream >> N_States;
		if(i==2) linestream >> ZP >> AP;
		if(i==3) linestream >> ZT >> AT;
		if(i==4) linestream >> EP;
		if(i==5) linestream >> Theta_Start >> Theta_End >> Theta_Step;
	
		if(verbose && i==5)
			std::cout << Theta_Start << " " << Theta_End << " " << Theta_Step << std::endl;

		if(i>5 && i <= 5+N_States)
		{
			linestream >> temp_index >> temp_J >> temp_E >> temp_P >> temp_K;
			level_index.push_back(temp_index);
			level_J.push_back(temp_J);
			level_E.push_back(temp_E);
			level_P.push_back(temp_P);
			level_K.push_back(temp_K);
			
		}
		if(i>5+N_States)
		{
			linestream >> temp_N >> temp_M >> temp_ME >> temp_Lambda;
			if(temp_N < temp_M)
			{
				trans_N.push_back(temp_N);
				trans_M.push_back(temp_M);
			}
			else
			{
				trans_N.push_back(temp_M);
				trans_M.push_back(temp_N);
			}
			trans_ME.push_back(temp_ME);
			if(temp_Lambda>=0) trans_Lambda.push_back(temp_Lambda);
			else trans_Lambda.push_back(5);
		}

		i++;

	}

	if(verbose)
		printf("Input file read successfully\n");
	
	if(Theta_End == 0) Theta_End = Theta_Start;
	if(Theta_Step == 0) Theta_Step = Theta_End - Theta_Start;
	
	if(level_E.at(0) != 0){
		printf("ERROR: Non-zero ground-state energy!\n"); 
		return;
	}
	if((int)level_E.size() != N_States){
		printf("ERROR: Inconsistent NMAX and number of states!\n"); 
		return;
	}
	if(level_J.at(0) > 7){ 
		printf("ERROR: Ground-state spin is greater than 7!\n"); 
		return;
	}

	if(verbose){
		printf("//=====================================//\n");
		printf("\t\tEXPERIMENT TITLE:\n");
		printf("//=====================================//\n\n");
		std::cout << title << std::endl;
		printf("\n");

		printf("//===========================================//\n");
		printf("\t\tINPUT STATES:\n");
		printf("//===========================================//\n");

		printf("Index:\tJ:\t\tEnergy:\t\tParity:\n");	

		for(unsigned int i=0;i<level_J.size();i++){
			printf("%li\t%lf\t%lf\t%li\n",(long int)level_index.at(i),level_J.at(i),level_E.at(i),(long int)level_P.at(i));
		}
		printf("\n");

		printf("//===========================================//\n");
		printf("\t\tINPUT MATRIX ELEMENTS::\n");
		printf("//===========================================//\n");

		printf("Init:\tFinal:\tMatrix Element:\t\tL:\n");	
		
		for(unsigned int i=0;i<trans_ME.size();i++){
			printf("%li\t%li\t%lf\t\t%li\n",(long int)trans_N.at(i),(long int)trans_M.at(i),trans_ME.at(i),(long int)trans_Lambda.at(i));
		}
		printf("\n");
	}
	
	if(verbose)
		printf("Creating transition matrices...\n");
	
	MaxLambda=0;
	for (unsigned int i = 0; i < trans_Lambda.size(); i++)
	  if (trans_Lambda.at(i) > MaxLambda)
	    MaxLambda = trans_Lambda.at(i);
	
	if(verbose)
		printf("Maximum multipolarity: %i\n",MaxLambda);
	
	TransitionMatrix.resize(MaxLambda);
	
	
	for(int i=0;i<MaxLambda;i++){
		TMatrixD tempmatrix(N_States,N_States);
		for(unsigned int j=0;j<trans_ME.size();j++){
			if(trans_Lambda.at(j) == (i+1)){
				tempmatrix[trans_N.at(j)-1][trans_M.at(j)-1]=trans_ME.at(j);
				tempmatrix[trans_M.at(j)-1][trans_N.at(j)-1]=trans_ME.at(j);
			}
		}
		if(verbose)
			printf("Matrix for multipolarity = %i:\n",i+1);
		TransitionMatrix.at(i).ResizeTo(N_States,N_States);
		TransitionMatrix.at(i) = tempmatrix;
		if(verbose)
			TransitionMatrix.at(i).Print();
	}
	if(verbose)
		printf("Transition matrices built\n");

}


#endif
