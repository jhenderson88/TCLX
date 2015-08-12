#ifndef TReactions_h
#define TReactions_h

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

//*****************************************//
//	TReactions class includes a load
//	of nuclear definitions, such as
//	A and Z, as well as some kinematic
//	properties.
//
//	Also included are the Wigner 3j 
//	symbols, and the Rutherford
//	cross section
//*****************************************//

class TReactions {

	
	public :

		TReactions(){ ;}
		virtual ~TReactions(){ ;}

		Double_t Adi_Par; // Adiabacity parameter

		Double_t beta; // v/c

		Double_t Ecm; // Center of mass energy

		Double_t Elab; // Lab energy

		Double_t Epsilon; // Orbit eccentricity


		// Stored as floats to simplify arithmetic and stop (reduce) casting errors
		Double_t Proj_A; // Projectile A
		Double_t Tar_A; // Target A
		Double_t Proj_Z; // Projectile Z
		Double_t Tar_Z; // Target Z

		void SetNuclei(Int_t Z1, Int_t A1, Int_t Z2, Int_t A2); // Set Nuclei, index 2 = target, index 1 = projectile

		//*********************************//
		// Get Nuclei indexing:
		// 0 = Target A
		// 1 = Target Z
		// 2 = Projectile A
		// 3 = Projectile Z
		//*********************************//

		Int_t GetNuclei(Int_t i) {
			if(i==0) return Tar_A;
			if(i==1) return Tar_Z;
			if(i==2) return Proj_A;
			if(i==3) return Proj_Z;
		}

		Double_t GetElab() { return Elab; } // Grab lab energy
		Double_t GetEcm()  { return Ecm;  } // Grab CoM energy
		
		static Double_t ThreeJ(Int_t I1, Int_t I2, Int_t I3, Int_t M1, Int_t M2, Int_t M3); // Three j value - inputs in 2*J format

		TSpline3* Rutherford(Double_t theta_start, Double_t theta_end, Double_t stepsize); // TSpline of Rutherford cross-section
		Double_t EvalRutherford(Double_t theta); // Rutherford for a given theta
		Double_t EvalRutherfordLevel(Double_t theta, Double_t E_level); // Total cross section for a given theta, accounting for relative velocities

		static Double_t hbarc; // MeV fm
		static Double_t finestruc; // Fine structure constant 
		static Double_t nuclearmagneton; // Nuclear magneton e.fm
		static Double_t dipole; // Normalisation factor for electric dipole polarization - nominally 0.005

		Double_t ClosestApproach(Int_t PA, Int_t TA, Int_t PZ, Int_t TZ, Double_t E_cm); // Define distance of closest approach
		Double_t ClosestApproach();	// As above, but uses pre-defined values

		Double_t RelativeVelocity(Double_t State_E);

		void SetElab(Double_t x) { // Set Elab, sets Ecm and beta too, if possible
			Elab = x; 
			if(Proj_A && Tar_A)
			{
				Ecm = Elab/(1+(Proj_A/Tar_A));
				beta = 0.046337*TMath::Sqrt(Elab / Proj_A);
			}
		}

		void SetEcm(Double_t x) { // Set Ecm, sets Elab and beta too, if possible
			Ecm = x; 
			if(Proj_A && Tar_A)
			{
				Elab = Ecm*(1+(Proj_A/Tar_A));
				beta = 0.046337*TMath::Sqrt(Elab / Proj_A);
			}
		}

		Double_t EtaCalc(Double_t v); // Calculates Eta, using relative velocities
	
		TSpline3* PlotSafeCoulEx(Double_t Sep, Double_t Theta_Min, Double_t Theta_Max, Double_t StepSize); // Plots the safe Coulex angles for a given distance of separation in fm (i.e. 5 fm is the "traditional" safe Coulex limit)

		//ClassDef(TReactions,1);

};

#endif


#ifdef TReactions_cxx

void TReactions::SetNuclei(Int_t Z1, Int_t A1, Int_t Z2, Int_t A2){

	Tar_A=A2; 
	Tar_Z=Z2; 
	Proj_A=A1, 
	Proj_Z=Z1;

	return;

}

Double_t TReactions::RelativeVelocity(Double_t State_E){

	Double_t v = 0.046337*TMath::Sqrt((Elab - State_E*(1 + ((Double_t)Proj_A/(Double_t)Tar_A) ) )/Proj_A);

	return v;

}

Double_t TReactions::ClosestApproach() { return ClosestApproach(Proj_A,Tar_A,Proj_Z,Tar_Z,Elab); }

Double_t TReactions::ClosestApproach(Int_t PA, Int_t TA, Int_t PZ, Int_t TZ, Double_t Elab){

	Double_t A = (2. * 0.71999 * (1. + ((Double_t)PA)/((Double_t)TA)) * ((Double_t)PZ) * ((Double_t)TZ)) /Elab;

	return A;

}

Double_t TReactions::EtaCalc(Double_t E){

	Double_t v = RelativeVelocity(E);

	Double_t s = finestruc * Tar_Z * Proj_Z / v;

	return s;

}

TSpline3* TReactions::PlotSafeCoulEx(Double_t Sep, Double_t Theta_Min, Double_t Theta_Max, Double_t StepSize)
{

	Double_t Contact = 1.25 * (pow(Proj_A,1/3) + pow(Tar_A,1/3));

	Int_t Steps = (Int_t) (Theta_Max-Theta_Min)/StepSize;

	Double_t Theta[Steps];
	Double_t SafeEnergy[Steps];

	Double_t Distance = Contact + Sep;

	for(int i=0;i<Steps;i++)
	{		

		Theta[i] = i * StepSize + Theta_Min;
		SafeEnergy[i] = 0.72 * (Tar_Z*Proj_Z/Distance) * ((Tar_A + Proj_A)/(Tar_A)) * (1 + (1 / (TMath::Sin((TMath::DegToRad() * Theta[i])/2))));

	}

	TSpline3 *spline = new TSpline3("SafeCoulExEnergy",Theta,SafeEnergy,Steps);
	
	return spline;

}


#endif
