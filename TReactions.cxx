#define TReactions_cxx

#include "TReactions.h"
#include "Math/SpecFunc.h"

// Set up some useful nuclear constants
Double_t TReactions::hbarc = 197.326;
Double_t TReactions::finestruc = 0.007297352;
Double_t TReactions::nuclearmagneton = 0.105155;
Double_t TReactions::dipole=0.005;

Double_t TReactions::ThreeJ(Int_t I1, Int_t I2, Int_t I3, Int_t M1, Int_t M2, Int_t M3){

	Double_t threej = (Double_t)ROOT::Math::wigner_3j(I1,I2,I3,M1,M2,M3); // Uses the GNU science libraries

	return threej;

}


TSpline3* TReactions::Rutherford(Double_t theta_start, Double_t theta_end, Double_t stepsize){

	Double_t a = pow(((Tar_Z * Proj_Z / Ecm)*(finestruc * hbarc)),2);
	
	Int_t steps = (theta_end-theta_start)/stepsize;

	Double_t Ruth[steps];
	Double_t Theta[steps];

	for(int i=0;i<steps;i++){
		Ruth[i] = a*(1/(4*pow(TMath::Sin(TMath::DegToRad()*(theta_start + i*stepsize)/2),4)));
		Theta[i] = theta_start + i*stepsize;
	}

	TSpline3 *ruth_cs = new TSpline3("Rutherford",Theta,Ruth,steps);

	return ruth_cs;

}

Double_t TReactions::EvalRutherford(Double_t theta)
{
	Double_t a = pow(((Tar_Z * Proj_Z / Ecm)*(finestruc * hbarc)),2);
	Double_t ruth = a*(1/(4*pow(TMath::Sin(TMath::DegToRad()*theta/2),4)));

	return ruth;

}


Double_t TReactions::EvalRutherfordLevel(Double_t theta, Double_t E_level)
{
	Double_t Epsilon = 1 / (TMath::Sin(theta * TMath::DegToRad() / 2.0));
	
	Double_t Separation = ClosestApproach();

	Double_t CS = 0.000625*TMath::Sqrt(Elab / (Elab - (1 + Proj_A/Tar_A) * E_level)) * pow(Separation,2) * pow(Epsilon,4);

	return CS;
}
