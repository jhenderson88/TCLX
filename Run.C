//*************************************************************************//
//	This macro does the basic jobs required of the code. It loads
//	the required libraries, then creates a TCLX. This has it's matrices
//	populated and the integration is performed.
//
//	A TExperiment is then created, in which we defined our experimental
//	variables. We then grab our previously created TCLX, which we
//	pass to the TExperiment. We then plot the yields vs theta.
//*************************************************************************//

{
	
	// You might find that you also need to include the following two commented lines
	//gSystem->Load("libMatrix.so");
	//gSystem->Load("libHist.so");
    	gSystem->Load("libMathMore.so");
    	gSystem->Load("libMathCore.so");
	gSystem->Load("libTExperiment.so");

	TCLX *myclx = new TCLX();
	myclx->GrabData("22Mg_104Pd_example.txt");
	myclx->PrepareMatrices();
	myclx->IntegrateTheta();

	TExperiment *exp = new TExperiment();
	exp->SetBeamIntensity(1e5);
	exp->SetTargetDensity(2);
	exp->SetExperimentLength(5);
	exp->GrabCLX(myclx);
	exp->PlotYieldLab(1)->Draw();

}
