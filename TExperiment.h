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
		TGraph* PlotDiffCrossSection(Int_t State);
		TGraph* PlotDiffCrossSectionLab(Int_t State);


		//*************************************//
		//	Turn those differential cross-
		//	sections into cross-sections.
		//*************************************//
		TGraph* PlotCrossSection(Int_t State);
		TGraph* PlotCrossSectionLab(Int_t State);

		//*************************************//
		//	Plot yields over the course of 
		//	a defined experiment vs theta.
		//	These are the PARTICLE yields, 
		//	and take no account of gammma-
		//	ray detection efficiency.
		//*************************************//
		TGraph* PlotYieldLab(Int_t State); 

		//*************************************//
		//	Plot the particle yield, assuming 
		//	gamma-ray detection requirement.
		//	i.e. take into account gamma-
		//	ray detection efficiency.
		//*************************************//
		TGraph* PlotYieldLabGamma(Int_t State, Int_t N_Detectors=12);

		//********************************************************************//
		// 	Graph of the TIGRESS efficiency based on the GEANT simulations 
		//	in the paper: C. E. Svensson et al., TIGRESS: TRIUMF-ISAC 
		//	gamma-ray escape-suppressed spectrometer,
		//	J. Phys. G: Nucl. Part. Phys. 31 (2005) S1663
		//
		//	These simulations were for a 12 detector configuration in
		//	an unsuppressed mode. Can be normalised to a different number
		//	of detectors (and eventually the suppressed mode).
		//********************************************************************//
		TGraph* TigressEfficiency(Int_t N_Detectors=12);


		//********************************************************************//
		//	Here we are going to define the detector coverage. At some
		//	point I would like to autmatically include SHARC, BAMBINO 
		//	and TIP detector geometries in here so they can be included
		//	easily, with the only variable being the distance from the
		//	target.
		//********************************************************************//


		std::vector<Double_t> det_cov_theta_min;
		std::vector<Double_t> det_cov_theta_max;
		void DefineDetectorCoverageTheta(Double_t theta_min, Double_t theta_max)
		{
			det_cov_theta_min.push_back(theta_min);
			det_cov_theta_max.push_back(theta_max);
		}

		void RemoveParticleDetectors() 
		{
			det_cov_theta_min.clear(); 
			det_cov_theta_max.clear(); 
		}
		
		//********************************************//
		//	This function will print the yields,
		//	separated by the included particle
		//	detectors.
		//********************************************//
		void PrintYieldsByDetector(Int_t N_Detectors=12);

};
#endif

#ifdef TExperiment_cxx

TGraph* TExperiment::TigressEfficiency(Int_t N_Detectors)
{
	Double_t Tig_Eff[12] = {0, 36.57, 40.21, 40.32, 35.97, 27.72, 23.35, 18.42, 12.64, 8.00, 5.83, 3.66};
	Double_t Tig_E[12] = {0, 40.92, 61.74, 103.72, 204.86, 411.31, 612.66, 1017.62, 2026.93, 4076.33, 6077.89, 10200.48};
	for(int i=0;i<12;i++)
		Tig_Eff[i] = Tig_Eff[i] * N_Detectors/12;
	for(int i=0;i<12;i++)
		Tig_E[i] = Tig_E[i] / 1000;
	TGraph *TigEff = new TGraph(11,Tig_E,Tig_Eff);

	return TigEff;
}

void TExperiment::PrintYieldsByDetector(Int_t N_Detectors)
{

	TVectorD Prob;
	Prob.ResizeTo(clx->N_States);

	std::ofstream outfile;
	outfile.open("TCLX_Detection_Yields.txt");

	Double_t step = (clx->Theta_Max - clx->Theta_Min)/clx->Probabilities.size();
	Double_t Tau = reaction->Proj_A / reaction->Tar_A;
	Double_t TargetNucleons = 0.0006022 * TargetDensity / reaction->Tar_A; // Target nucleon density (Avg. const * size of a barn * target density / target mass)

	for(int i=0;i<det_cov_theta_min.size();i++)
	{
		Double_t MinTheta = det_cov_theta_min.at(i);
		Double_t MaxTheta = det_cov_theta_max.at(i);

		outfile << "Detector " << (i+1) << " coverage theta between " << MinTheta << " and " << MaxTheta << "\n\n";


		for(int k=0;k<clx->N_States;k++)
		{

			outfile << "State " << (k+1) << " J: " << std::setw(3) << clx->level_J.at(k) << " | State E: " << std::setw(8) << clx->level_E.at(k) << " MeV |";

			Double_t tempCS = 0;
			Double_t CS = 0;
			Double_t Yield = 0;
			Double_t GammaEff =	TigressEfficiency(N_Detectors)->Eval(clx->level_E.at(k)) / 100;

			for(int j=0;j<clx->Probabilities.size();j++)
			{
				Prob = clx->Probabilities.at(j);		
				Double_t ruth = reaction->EvalRutherfordLevel((clx->Theta_Min+j*step),clx->level_E.at(k));
				Double_t thetacm = clx->Theta_Min + j*step;
				Double_t thetalab = TMath::RadToDeg() * TMath::ATan( (TMath::Sin(TMath::DegToRad() * thetacm) / ( ( TMath::Cos(TMath::DegToRad() * thetacm) + Tau ) )));//thetacm

				if(thetalab > 0)
					thetalab = thetalab;
				else
					thetalab = thetalab + 180;

				if(thetalab > MinTheta && thetalab < MaxTheta){
					tempCS = Prob[k] * ruth * 2 * TMath::Pi() * TMath::Sin(TMath::DegToRad() * thetalab) * (TMath::DegToRad());
					CS = CS + tempCS;
				}
			}

			Yield = CS * TargetNucleons * BeamIntensity * ExperimentLength * 60 * 60 * 24;
	
			outfile << " Particles detected: " << std::setw(12) << Yield << " | ";
			
			Yield = Yield * GammaEff;

			outfile << "Coincident gamma rays: " << std::setw(8) << Yield << " |\n";

		}

		outfile << "\n\n";

	}

	outfile.close();

}

#endif
