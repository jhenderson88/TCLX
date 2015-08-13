#define TExperiment_cxx

#include "TExperiment.h"

TGraph* TExperiment::PlotDiffCrossSection(Int_t State)
{

	if(State >= clx->N_States || State < 0){
		printf("Error, state doesn't exist\n");
	}

	std::vector<TVectorD> Probabilities;
	if(clx->E_Loss_inc) 
		Probabilities = clx->ELossProbabilities.at(0);
	else
		Probabilities = clx->Probabilities;

	TVectorD Prob;	
	Prob.ResizeTo(clx->N_States);

	Double_t CS[Probabilities.size()];
	Double_t Theta[Probabilities.size()];

	Double_t step = (clx->Theta_Max - clx->Theta_Min)/Probabilities.size();


	for(unsigned int i=0;i<Probabilities.size();i++)
	{
		Prob = Probabilities.at(i);		

		Double_t ruth = reaction->EvalRutherfordLevel((clx->Theta_Min+i*step),clx->level_E.at(State));
		CS[i] = Prob[State] * ruth;
		Theta[i] = clx->Theta_Min + i*step;
	}

	char sname[64];
	sprintf(sname,"DiffCrossSection_CoM_%i",State);
	TGraph *CSSpline = new TGraph(Probabilities.size(),Theta,CS);
  	CSSpline->SetTitle(Form("Differential cross section - Incident E_{Beam}; #theta_{CoM} [deg]; d#sigma/d#omega [b/sr];"));
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

	std::vector<TVectorD> Probabilities;
	if(clx->E_Loss_inc) 
		Probabilities = clx->ELossProbabilities.at(0);
	else
		Probabilities = clx->Probabilities;

	TVectorD Prob;	
	Prob.ResizeTo(clx->N_States);

	Double_t CS[Probabilities.size()];
	Double_t Theta[Probabilities.size()];

	Double_t step = (clx->Theta_Max - clx->Theta_Min)/Probabilities.size();

	Double_t Tau = reaction->Proj_A / reaction->Tar_A;

	for(unsigned int i=0;i<Probabilities.size();i++)
	{
		Prob = Probabilities.at(i);		

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
	TGraph *CSSpline = new TGraph(Probabilities.size(),Theta,CS);
  	CSSpline->SetTitle(Form("Differential cross section lab - Incident E_{Beam}; #theta_{lab} [deg]; d#sigma/d#omega [b/sr];"));
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

	std::vector<TVectorD> Probabilities;
	if(clx->E_Loss_inc) 
		Probabilities = clx->ELossProbabilities.at(0);
	else
		Probabilities = clx->Probabilities;

	TVectorD Prob;	
	Prob.ResizeTo(clx->N_States);

	Double_t CS[Probabilities.size()];
	Double_t Theta[Probabilities.size()];

	Double_t step = (clx->Theta_Max - clx->Theta_Min)/Probabilities.size();


	for(unsigned int i=0;i<Probabilities.size();i++)
	{
		Prob = Probabilities.at(i);		

		Double_t ruth = reaction->EvalRutherfordLevel((clx->Theta_Min+i*step),clx->level_E.at(State));
		Double_t tempCS = Prob[State] * ruth;
		CS[i] = tempCS * 2 * TMath::Pi() * TMath::Sin(TMath::DegToRad() * (clx->Theta_Min + i*step)) * (TMath::DegToRad());
		Theta[i] = clx->Theta_Min + i*step;

	}

	char sname[64];
	sprintf(sname,"DiffCrossSection_CoM_%i",State);
	TGraph *CSSpline = new TGraph(Probabilities.size(),Theta,CS);
  	CSSpline->SetTitle(Form("Cross section - Incident E_{Beam}; #theta_{CoM} [deg]; #sigma [b];"));
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

	std::vector<TVectorD> Probabilities;
	if(clx->E_Loss_inc) 
		Probabilities = clx->ELossProbabilities.at(0);
	else
		Probabilities = clx->Probabilities;

	if(clx->E_Loss_inc)
		printf("Printing cross section for incident energy!\n");

	
	TVectorD Prob;	
	Prob.ResizeTo(clx->N_States);

	Double_t CS[Probabilities.size()];
	Double_t Theta[Probabilities.size()];

	printf("Probabilities size: %i\n",(Int_t)Probabilities.size());

	Double_t step = (clx->Theta_Max - clx->Theta_Min)/Probabilities.size();

	Double_t Tau = reaction->Proj_A / reaction->Tar_A;

	for(unsigned int i=0;i<Probabilities.size();i++)
	{
		Prob = Probabilities.at(i);		

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
	TGraph *CSSpline = new TGraph(Probabilities.size(),Theta,CS);
  	CSSpline->SetTitle(Form("Cross section lab - Incident E_{Beam}; #theta_{lab} [deg]; #sigma [b];"));
	CSSpline->GetYaxis()->SetTitleOffset(1.2);
	CSSpline->GetYaxis()->SetTitleSize(0.04);
	CSSpline->GetXaxis()->SetTitleSize(0.04);
	CSSpline->GetYaxis()->CenterTitle();
	CSSpline->GetXaxis()->CenterTitle();

	return CSSpline;
	
}

TGraph* TExperiment::PlotYieldLab(Int_t State)
{

	TGraph *yield;

	std::vector<TVectorD> Probabilities;
	if(clx->E_Loss_inc) 
	{

		Double_t Theta[clx->ELossProbabilities.at(0).size()];
		Double_t Yield[clx->ELossProbabilities.at(0).size()];

		for(int i=0;i<clx->ELossProbabilities.at(0).size();i++)
			Yield[i] = 0;
	
		Double_t TargetNucleons = 0.0006022 * TargetDensity / reaction->Tar_A; // Target nucleon density (Avg. const * size of a barn * target density / target mass)

		TargetNucleons = TargetNucleons / (Double_t)clx->ELossProbabilities.size();

		for(int de = 0; de<clx->ELossProbabilities.size(); de++)
		{
			Probabilities = clx->ELossProbabilities.at(de);
			TVectorD Prob;	
			Prob.ResizeTo(clx->N_States);

			reaction->SetElab(clx->TargetEnergies.at(de));

			if(!TargetDensity || !BeamIntensity || !ExperimentLength)
			{
				printf("Experiment details not set!\n");
			}
			else
			{


				Double_t step = (clx->Theta_Max - clx->Theta_Min)/Probabilities.size();

				Double_t Tau = reaction->Proj_A / reaction->Tar_A;

				for(unsigned int i=0;i<Probabilities.size();i++)
				{
					Prob = Probabilities.at(i);		

					Double_t ruth = reaction->EvalRutherfordLevel((clx->Theta_Min+i*step),clx->level_E.at(State));
					Double_t thetacm = clx->Theta_Min + i*step;
					Double_t thetalab = TMath::RadToDeg() * TMath::ATan( (TMath::Sin(TMath::DegToRad() * thetacm) / ( ( TMath::Cos(TMath::DegToRad() * thetacm) + Tau ) )));//thetacm
					if(thetalab > 0)
						Theta[i] = thetalab;
					else
						Theta[i] = thetalab + 180;

					Double_t tempCS = Prob[State] * ruth * 2 * TMath::Pi() * TMath::Sin(TMath::DegToRad() * Theta[i]) * (TMath::DegToRad());
					Double_t tempyield = tempCS * TargetNucleons * BeamIntensity * ExperimentLength * 60 * 60 * 24;
					Yield[i] += (Double_t)(tempyield);

				}
			}
		}
	
		yield = new TGraph(Probabilities.size(),Theta,Yield);
	  	yield->SetTitle(Form("Experiment yield - Incident E_{Beam}; #theta_{lab} [deg]; Yield;"));
		yield->GetYaxis()->SetTitleOffset(1.2);
		yield->GetYaxis()->SetTitleSize(0.04);
		yield->GetXaxis()->SetTitleSize(0.04);
		yield->GetYaxis()->CenterTitle();
		yield->GetXaxis()->CenterTitle();

		return yield;

	}
	else if(!clx->E_Loss_inc)
	{
		Probabilities = clx->Probabilities;

		TVectorD Prob;	
		Prob.ResizeTo(clx->N_States);

		Double_t CS[Probabilities.size()];
		Double_t Theta[Probabilities.size()];
		Double_t Yield[Probabilities.size()];

		for(int i=0;i<Probabilities.size();i++)
			Yield[i]=0;

		if(!TargetDensity || !BeamIntensity || !ExperimentLength)
		{
			printf("Experiment details not set!\n");
		}
		else
		{

			Double_t TargetNucleons = 0.0006022 * TargetDensity / reaction->Tar_A; // Target nucleon density (Avg. const * size of a barn * target density / target mass)

			Double_t step = (clx->Theta_Max - clx->Theta_Min)/Probabilities.size();

			Double_t Tau = reaction->Proj_A / reaction->Tar_A;

			for(unsigned int i=0;i<Probabilities.size();i++)
			{
				Prob = Probabilities.at(i);		


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
	
		yield = new TGraph(Probabilities.size(),Theta,Yield);
	  	yield->SetTitle(Form("Experiment yield - Incident E_{Beam}; #theta_{lab} [deg]; Yield;"));
		yield->GetYaxis()->SetTitleOffset(1.2);
		yield->GetYaxis()->SetTitleSize(0.04);
		yield->GetXaxis()->SetTitleSize(0.04);
		yield->GetYaxis()->CenterTitle();
		yield->GetXaxis()->CenterTitle();

		return yield;

	}

	return yield;

}

TGraph* TExperiment::PlotYieldLabGamma(Int_t State, Int_t N_Detectors)
{

	TGraph *yield;

	std::vector<TVectorD> Probabilities;
	if(clx->E_Loss_inc) 
	{

		Double_t Theta[clx->ELossProbabilities.at(0).size()];
		Double_t Yield[clx->ELossProbabilities.at(0).size()];

		for(int i=0;i<clx->ELossProbabilities.at(0).size();i++)
			Yield[i] = 0;
	
		Double_t TargetNucleons = 0.0006022 * TargetDensity / reaction->Tar_A; // Target nucleon density (Avg. const * size of a barn * target density / target mass)

		TargetNucleons = TargetNucleons / (Double_t)clx->ELossProbabilities.size();

		for(int de = 0; de<clx->ELossProbabilities.size(); de++)
		{
			Probabilities = clx->ELossProbabilities.at(de);
			TVectorD Prob;	
			Prob.ResizeTo(clx->N_States);

			reaction->SetElab(clx->TargetEnergies.at(de));

			Double_t GammaEff =	TigressEfficiency(N_Detectors)->Eval(clx->level_E.at(State)) / 100;

			if(!TargetDensity || !BeamIntensity || !ExperimentLength)
			{
				printf("Experiment details not set!\n");
			}
			else
			{


				Double_t step = (clx->Theta_Max - clx->Theta_Min)/Probabilities.size();

				Double_t Tau = reaction->Proj_A / reaction->Tar_A;

				for(unsigned int i=0;i<Probabilities.size();i++)
				{
					Prob = Probabilities.at(i);		

					Double_t ruth = reaction->EvalRutherfordLevel((clx->Theta_Min+i*step),clx->level_E.at(State));
					Double_t thetacm = clx->Theta_Min + i*step;
					Double_t thetalab = TMath::RadToDeg() * TMath::ATan( (TMath::Sin(TMath::DegToRad() * thetacm) / ( ( TMath::Cos(TMath::DegToRad() * thetacm) + Tau ) )));//thetacm
					if(thetalab > 0)
						Theta[i] = thetalab;
					else
						Theta[i] = thetalab + 180;

					Double_t tempCS = Prob[State] * ruth * 2 * TMath::Pi() * TMath::Sin(TMath::DegToRad() * Theta[i]) * (TMath::DegToRad());
					Double_t tempyield = tempCS * TargetNucleons * BeamIntensity * ExperimentLength * 60 * 60 * 24;
					Yield[i] += (Double_t)(GammaEff * tempyield);

				}
			}
		}
	
		yield = new TGraph(Probabilities.size(),Theta,Yield);
	  	yield->SetTitle(Form("Experiment yield - Incident E_{Beam}; #theta_{lab} [deg]; Yield;"));
		yield->GetYaxis()->SetTitleOffset(1.2);
		yield->GetYaxis()->SetTitleSize(0.04);
		yield->GetXaxis()->SetTitleSize(0.04);
		yield->GetYaxis()->CenterTitle();
		yield->GetXaxis()->CenterTitle();

		return yield;

	}
	else if(!clx->E_Loss_inc)
	{
		Probabilities = clx->Probabilities;

		TVectorD Prob;	
		Prob.ResizeTo(clx->N_States);

		Double_t CS[Probabilities.size()];
		Double_t Theta[Probabilities.size()];
		Double_t Yield[Probabilities.size()];

		for(int i=0;i<Probabilities.size();i++)
			Yield[i]=0;

		Double_t GammaEff =	TigressEfficiency(N_Detectors)->Eval(clx->level_E.at(State)) / 100;

		if(!TargetDensity || !BeamIntensity || !ExperimentLength)
		{
			printf("Experiment details not set!\n");
		}
		else
		{

			Double_t TargetNucleons = 0.0006022 * TargetDensity / reaction->Tar_A; // Target nucleon density (Avg. const * size of a barn * target density / target mass)

			Double_t step = (clx->Theta_Max - clx->Theta_Min)/Probabilities.size();

			Double_t Tau = reaction->Proj_A / reaction->Tar_A;

			for(unsigned int i=0;i<Probabilities.size();i++)
			{
				Prob = Probabilities.at(i);		


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
	
		yield = new TGraph(Probabilities.size(),Theta,Yield);
	  	yield->SetTitle(Form("Experiment yield - Incident E_{Beam}; #theta_{lab} [deg]; Yield;"));
		yield->GetYaxis()->SetTitleOffset(1.2);
		yield->GetYaxis()->SetTitleSize(0.04);
		yield->GetXaxis()->SetTitleSize(0.04);
		yield->GetYaxis()->CenterTitle();
		yield->GetXaxis()->CenterTitle();

		return yield;

	}

	return yield;

}

TGraph* TExperiment::PlotCrossSectionLabdE(Int_t State, Int_t x)
{

	if(State >= clx->N_States || State < 0){
		printf("Error, state doesn't exist\n");
	}

	std::vector<TVectorD> Probabilities;
	if(clx->E_Loss_inc && x < clx->ELossProbabilities.size()) 
		Probabilities = clx->ELossProbabilities.at(x);
	else if(!clx->E_Loss_inc)
		Probabilities = clx->Probabilities;
	else if(clx->E_Loss_inc && x >= clx->ELossProbabilities.size())
	{
		printf("X beyond target range. Printing cross section for incident energy!\n");
		Probabilities = clx->ELossProbabilities.at(0);
	}
	
	TVectorD Prob;	
	Prob.ResizeTo(clx->N_States);

	Double_t CS[Probabilities.size()];
	Double_t Theta[Probabilities.size()];

	printf("Probabilities size: %i\n",(Int_t)Probabilities.size());

	Double_t step = (clx->Theta_Max - clx->Theta_Min)/Probabilities.size();

	Double_t Tau = reaction->Proj_A / reaction->Tar_A;

	for(unsigned int i=0;i<Probabilities.size();i++)
	{
		Prob = Probabilities.at(i);		

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
	TGraph *CSSpline = new TGraph(Probabilities.size(),Theta,CS);
  	CSSpline->SetTitle(Form("Cross section lab - Incident E_{Beam}; #theta_{lab} [deg]; #sigma [b];"));
	CSSpline->GetYaxis()->SetTitleOffset(1.2);
	CSSpline->GetYaxis()->SetTitleSize(0.04);
	CSSpline->GetXaxis()->SetTitleSize(0.04);
	CSSpline->GetYaxis()->CenterTitle();
	CSSpline->GetXaxis()->CenterTitle();

	return CSSpline;
	
}
