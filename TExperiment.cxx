#define TExperiment_cxx

#include "TExperiment.h"

TGraph* TExperiment::PlotDiffCrossSection(Int_t State)
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
	TGraph *CSSpline = new TGraph(clx->Probabilities.size(),Theta,CS);
  	CSSpline->SetTitle(Form("Differential cross section; #theta_{CoM} [deg]; d#sigma/d#omega [b/sr];"));
	CSSpline->GetYaxis()->SetTitleOffset(1.2);
	CSSpline->GetYaxis()->SetTitleSize(0.04);
	CSSpline->GetXaxis()->SetTitleSize(0.04);
	CSSpline->GetYaxis()->CenterTitle();
	CSSpline->GetXaxis()->CenterTitle();

	return CSSpline;
	
}

TGraph* TExperiment::PlotDiffCrossSectionLab(Int_t State)
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
	TGraph *CSSpline = new TGraph(clx->Probabilities.size(),Theta,CS);
  	CSSpline->SetTitle(Form("Differential cross section lab; #theta_{lab} [deg]; d#sigma/d#omega [b/sr];"));
	CSSpline->GetYaxis()->SetTitleOffset(1.2);
	CSSpline->GetYaxis()->SetTitleSize(0.04);
	CSSpline->GetXaxis()->SetTitleSize(0.04);
	CSSpline->GetYaxis()->CenterTitle();
	CSSpline->GetXaxis()->CenterTitle();

	return CSSpline;
	
}

TGraph* TExperiment::PlotCrossSection(Int_t State)
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
	TGraph *CSSpline = new TGraph(clx->Probabilities.size(),Theta,CS);
  	CSSpline->SetTitle(Form("Cross section; #theta_{CoM} [deg]; #sigma [b];"));
	CSSpline->GetYaxis()->SetTitleOffset(1.2);
	CSSpline->GetYaxis()->SetTitleSize(0.04);
	CSSpline->GetXaxis()->SetTitleSize(0.04);
	CSSpline->GetYaxis()->CenterTitle();
	CSSpline->GetXaxis()->CenterTitle();

	return CSSpline;
	
}

TGraph* TExperiment::PlotCrossSectionLab(Int_t State)
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
	TGraph *CSSpline = new TGraph(clx->Probabilities.size(),Theta,CS);
  	CSSpline->SetTitle(Form("Cross section lab; #theta_{lab} [deg]; #sigma [b];"));
	CSSpline->GetYaxis()->SetTitleOffset(1.2);
	CSSpline->GetYaxis()->SetTitleSize(0.04);
	CSSpline->GetXaxis()->SetTitleSize(0.04);
	CSSpline->GetYaxis()->CenterTitle();
	CSSpline->GetXaxis()->CenterTitle();

	return CSSpline;
	
}

TGraph* TExperiment::PlotYieldLab(Int_t State)
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
	
	TGraph *yield = new TGraph(clx->Probabilities.size(),Theta,Yield);
  	yield->SetTitle(Form("Experiment yield; #theta_{lab} [deg]; Yield;"));
	yield->GetYaxis()->SetTitleOffset(1.2);
	yield->GetYaxis()->SetTitleSize(0.04);
	yield->GetXaxis()->SetTitleSize(0.04);
	yield->GetYaxis()->CenterTitle();
	yield->GetXaxis()->CenterTitle();
	
	return yield;

}

TGraph* TExperiment::PlotYieldLabGamma(Int_t State, Int_t N_Detectors)
{

	TVectorD Prob;	
	Prob.ResizeTo(clx->N_States);

	Double_t CS[clx->Probabilities.size()];
	Double_t Theta[clx->Probabilities.size()];
	Double_t Yield[clx->Probabilities.size()];

	Double_t GammaEff =	TigressEfficiency(N_Detectors)->Eval(clx->level_E.at(State)) / 100;

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
			Yield[i] = GammaEff * tempyield;

		}
	}
	
	TGraph *yield = new TGraph(clx->Probabilities.size(),Theta,Yield);
  	yield->SetTitle(Form("Experiment yield; #theta_{lab} [deg]; Yield;"));
	yield->GetYaxis()->SetTitleOffset(1.2);
	yield->GetYaxis()->SetTitleSize(0.04);
	yield->GetXaxis()->SetTitleSize(0.04);
	yield->GetYaxis()->CenterTitle();
	yield->GetXaxis()->CenterTitle();

	return yield;

}
