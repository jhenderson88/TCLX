#ifndef TExperiment_h
#define TExperiment_h

#include "TCLX.h" // The TCLX class performs Coulomb excitation calculations
#include "TReactions.h" // The TReactions class includes a load of useful reactions (duh) information, such as Rutherford cross sections, as well as a load of information defining the nuclear properties

//*****************************************************//
//	The TExperiment class is intended to grab
//	experimental classes, such as TCLX, and 
//	turn them into some realistic experimental
//	values, such as cross-sections in the lab
//	frames, and total yield values.
//*****************************************************//


class TExperiment {

	public :
		
		TExperiment(){;} 
		virtual ~TExperiment(){;}

		//***********************************//
		//	These are the variables which
		//	define our experiment.
		//***********************************//
		Int_t BeamIntensity;
		Double_t TargetDensity;
		Double_t ExperimentLength;

		//***********************************//
		//	Experimental variables set
		//	with these functions.
		//***********************************//
		void SetBeamIntensity(Int_t intensity) { BeamIntensity = intensity; }
		void SetTargetDensity(Double_t density) { TargetDensity = density; }
		void SetExperimentLength(Double_t length) {ExperimentLength = length; }


		//*************************************//
		//	Defined and ready to be grabbed
		//*************************************//
		TReactions *reaction;
		TCLX *clx;

		//*************************************//
		//	Grab a TCLX which has already
		//	been calculated, ready for
		//	turning into some real values.
		//*************************************// 
		void GrabCLX(TCLX *clx_experiment) 
		{
			clx = clx_experiment; 
			reaction = clx->GetReaction();
		}

		//*************************************//
		//	Plot differential cross-sections.
		//	These are the same plots which are
		//	produced in the TCLX code
		//*************************************//
		TSpline3* PlotDiffCrossSection(Int_t State);
		TSpline3* PlotDiffCrossSectionLab(Int_t State);


		//*************************************//
		//	Turn those differential cross-
		//	sections into cross-sections.
		//*************************************//
		TSpline3* PlotCrossSection(Int_t State);
		TSpline3* PlotCrossSectionLab(Int_t State);

		//*************************************//
		//	Plot yields over the course of 
		//	a defined experiment vs theta.
		//	These are the PARTICLE yields, 
		//	and take no account of gammma-
		//	ray detection efficiency.
		//*************************************//
		TSpline3* PlotYieldLab(Int_t State); 

};
#endif

#ifdef TExperiment_cxx

TSpline3* TExperiment::PlotDiffCrossSection(Int_t State)
{

	if(State >= clx->N_States || State < 0){
		printf("Error, state doesn't exist\n");
	}

	TVectorD Prob;	
	Prob.ResizeTo(clx->N_States);

	Double_t CS[clx->Probabilities.size()];
	Double_t Theta[clx->Probabilities.size()];

	Double_t step = (clx->Theta_Max - clx->Theta_Min)/clx->Probabilities.size();


	for(unsigned int i=0;i<clx->Probabilities.size();i++)
	{
		Prob = clx->Probabilities.at(i);		

		Double_t ruth = reaction->EvalRutherfordLevel((clx->Theta_Min+i*step),clx->level_E.at(State));
		CS[i] = Prob[State] * ruth;
		Theta[i] = clx->Theta_Min + i*step;
	}

	char sname[64];
	sprintf(sname,"DiffCrossSection_CoM_%i",State);
	TSpline3 *CSSpline = new TSpline3(sname,Theta,CS,clx->Probabilities.size());

	return CSSpline;
	
}

TSpline3* TExperiment::PlotDiffCrossSectionLab(Int_t State)
{

	if(State >= clx->N_States || State < 0){
		printf("Error, state doesn't exist\n");
	}

	TVectorD Prob;	
	Prob.ResizeTo(clx->N_States);

	Double_t CS[clx->Probabilities.size()];
	Double_t Theta[clx->Probabilities.size()];

	Double_t step = (clx->Theta_Max - clx->Theta_Min)/clx->Probabilities.size();

	Double_t Tau = reaction->Proj_A / reaction->Tar_A;

	for(unsigned int i=0;i<clx->Probabilities.size();i++)
	{
		Prob = clx->Probabilities.at(i);		

		Double_t ruth = reaction->EvalRutherfordLevel((clx->Theta_Min+i*step),clx->level_E.at(State));
		CS[i] = Prob[State] * ruth;
		Double_t thetacm = clx->Theta_Min + i*step;
		Double_t thetalab = TMath::RadToDeg() * TMath::ATan( (TMath::Sin(TMath::DegToRad() * thetacm) / ( ( TMath::Cos(TMath::DegToRad() * thetacm) + Tau ) )));//thetacm
		if(thetalab > 0)
			Theta[i] = thetalab;
		else
			Theta[i] = thetalab + 180;
	}

	char sname[64];
	sprintf(sname,"DiffCrossSection_Lab_State_%i",State);
	TSpline3 *CSSpline = new TSpline3(sname,Theta,CS,clx->Probabilities.size());

	return CSSpline;
	
}

TSpline3* TExperiment::PlotCrossSection(Int_t State)
{

	if(State >= clx->N_States || State < 0){
		printf("Error, state doesn't exist\n");
	}

	TVectorD Prob;	
	Prob.ResizeTo(clx->N_States);

	Double_t CS[clx->Probabilities.size()];
	Double_t Theta[clx->Probabilities.size()];

	Double_t step = (clx->Theta_Max - clx->Theta_Min)/clx->Probabilities.size();


	for(unsigned int i=0;i<clx->Probabilities.size();i++)
	{
		Prob = clx->Probabilities.at(i);		

		Double_t ruth = reaction->EvalRutherfordLevel((clx->Theta_Min+i*step),clx->level_E.at(State));
		Double_t tempCS = Prob[State] * ruth;
		CS[i] = tempCS * 2 * TMath::Pi() * TMath::Sin(TMath::DegToRad() * (clx->Theta_Min + i*step)) * (TMath::DegToRad());
		Theta[i] = clx->Theta_Min + i*step;

	}

	char sname[64];
	sprintf(sname,"DiffCrossSection_CoM_%i",State);
	TSpline3 *CSSpline = new TSpline3(sname,Theta,CS,clx->Probabilities.size());

	return CSSpline;
	
}

TSpline3* TExperiment::PlotCrossSectionLab(Int_t State)
{

	if(State >= clx->N_States || State < 0){
		printf("Error, state doesn't exist\n");
	}

	TVectorD Prob;	
	Prob.ResizeTo(clx->N_States);

	Double_t CS[clx->Probabilities.size()];
	Double_t Theta[clx->Probabilities.size()];

	Double_t step = (clx->Theta_Max - clx->Theta_Min)/clx->Probabilities.size();

	Double_t Tau = reaction->Proj_A / reaction->Tar_A;

	for(unsigned int i=0;i<clx->Probabilities.size();i++)
	{
		Prob = clx->Probabilities.at(i);		

		Double_t ruth = reaction->EvalRutherfordLevel((clx->Theta_Min+i*step),clx->level_E.at(State));
		Double_t thetacm = clx->Theta_Min + i*step;
		Double_t thetalab = TMath::RadToDeg() * TMath::ATan( (TMath::Sin(TMath::DegToRad() * thetacm) / ( ( TMath::Cos(TMath::DegToRad() * thetacm) + Tau ) )));//thetacm
		if(thetalab > 0)
			Theta[i] = thetalab;
		else
			Theta[i] = thetalab + 180;
		Double_t tempCS = Prob[State] * ruth;
		CS[i] = tempCS * 2 * TMath::Pi() * TMath::Sin(TMath::DegToRad() * Theta[i]) * (TMath::DegToRad());

	}

	char sname[64];
	sprintf(sname,"DiffCrossSection_Lab_State_%i",State);
	TSpline3 *CSSpline = new TSpline3(sname,Theta,CS,clx->Probabilities.size());

	return CSSpline;
	
}

TSpline3* TExperiment::PlotYieldLab(Int_t State)
{

	TVectorD Prob;	
	Prob.ResizeTo(clx->N_States);

	Double_t CS[clx->Probabilities.size()];
	Double_t Theta[clx->Probabilities.size()];
	Double_t Yield[clx->Probabilities.size()];

	if(!TargetDensity || !BeamIntensity || !ExperimentLength)
	{
		printf("Experiment details not set!\n");
	}
	else
	{

		Double_t TargetNucleons = 0.0006022 * TargetDensity / reaction->Tar_A; // Target nucleon density (Avg. const * size of a barn * target density / target mass)

		Double_t step = (clx->Theta_Max - clx->Theta_Min)/clx->Probabilities.size();

		Double_t Tau = reaction->Proj_A / reaction->Tar_A;

		for(unsigned int i=0;i<clx->Probabilities.size();i++)
		{
			Prob = clx->Probabilities.at(i);		

			Double_t ruth = reaction->EvalRutherfordLevel((clx->Theta_Min+i*step),clx->level_E.at(State));
			Double_t thetacm = clx->Theta_Min + i*step;
			Double_t thetalab = TMath::RadToDeg() * TMath::ATan( (TMath::Sin(TMath::DegToRad() * thetacm) / ( ( TMath::Cos(TMath::DegToRad() * thetacm) + Tau ) )));//thetacm
			if(thetalab > 0)
				Theta[i] = thetalab;
			else
				Theta[i] = thetalab + 180;

			Double_t tempCS = Prob[State] * ruth;
			CS[i] = tempCS * 2 * TMath::Pi() * TMath::Sin(TMath::DegToRad() * Theta[i]) * (TMath::DegToRad());
			Double_t tempyield = CS[i] * TargetNucleons * BeamIntensity * ExperimentLength * 60 * 60 * 24;
			Yield[i] = tempyield;

		}
	}
	
	TSpline3 *yield = new TSpline3("Yield vs theta",Theta,Yield,clx->Probabilities.size());

	return yield;

}

#endif
