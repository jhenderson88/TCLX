#include "TCLX.cxx"
#include "TCLX_DataInput.cxx"
#include "TReactions.cxx"
#include <Math/SpecFunc.h>
#include "linkdef.h"
#include <TFile.h>

void Run(const char* filename, Int_t State = 1, Option_t *opt=""){

	TCLX *myclx = new TCLX();
	
	myclx->GrabData(filename,opt);

	myclx->PrepareMatrices();

	myclx->IntegrateTheta();

	myclx->StateCrossSection(State)->Draw();

	TGraph *graphs[myclx->N_States];
	for(int i=0;i<myclx->N_States;i++)
		graphs[i] = myclx->StateCrossSection(i);

	TFile *myfile = new TFile("TCLX_Plots.root","RECREATE");
	for(int i=0;i<myclx->N_States;i++)
		graphs[i]->Write();
	myfile->Close();

	myclx->PrintCalculationDetails();

	myclx->PrintCrossSections("22Mg_CS.txt");

}
