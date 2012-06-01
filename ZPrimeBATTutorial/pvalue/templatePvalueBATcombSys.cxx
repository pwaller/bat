// ***************************************************************
// This file was created using the ./CreateProject.sh script
// for project BCMultiTemplateFitter
// ./CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>
#include <BAT/BCSummaryPriorModel.h>
#include <BAT/BCModelOutput.h>
#include <BAT/BCH1D.h>
#include <BAT/BCH2D.h>

#include <BCMTFAnalysisFacility.h>
#include <BCMultiTemplateFitter.h>
#include <BCChannel.h>

#include <TROOT.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TLegend.h>
#include <TColor.h>
#include <TArrow.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <stdlib.h>
 
#include "plot2Dscan.h"
#include "plot_pvalue.h"

using std::cerr;
using std::cout;
using std::endl;
using std::vector;


int main(int argc, char ** argv){
  
  if(argc<11){
    cerr<<"usage:  templatePvalueBATcombSys.exe iRun doEnsemble nEnsemble iMass doPlots doMCMC doSYS doTheory doElectron doMuon [model]\n";
    exit(1);
  }
  /* Note: Typically run PE before Data.

     Examples:
     ./templatePvalueBATcombSys 2 1 5 0 1 1 0 0 0 1    # PE,   muon channel
     ./templatePvalueBATcombSys 3 0 5 0 1 1 0 0 0 1    # Data, muon channel
     ./templatePvalueBATcombSys 4 1 5 0 1 1 0 0 1 0    # PE,   electron channel
     ./templatePvalueBATcombSys 5 0 5 0 1 1 0 0 1 0    # Data, electron channel
  */


  cerr << "++ Running: " << argv[0] << endl;
  //cout << "++ with iMass: " << argv[1] << " iRun "<< argv[2] << endl;
	
  //Some flags to set (could come from cmd line)
  int iRun = atoi(argv[1]);   		// iRun number (for parallel running)
  bool doEnsemble=atoi(argv[2]);//true; // true=doPEs, false=data analysis
  int nEnsemble=atoi(argv[3]);//150;	// Number of pseudo-experiments (PEs) if doEnsemble is true
  int iMass = atoi(argv[4]);            // template number
  bool doPlots=atoi(argv[5]);//true;	// Make plots of posterior PDF etc, for data analysis
  bool doMCMC=atoi(argv[6]);//true;	// only relevant for PEs (use MarkovChain MC for marginalization or Minuit profiling)
  bool doSys =atoi(argv[7]);//true;
  bool doTheoryUncertainty =atoi(argv[8]);// false= Z' theory uncertainty outside of likelihood , true= Z' theory uncertainty included in likelihood function 
  bool doElectron=atoi(argv[9]);
  bool doMuon=atoi(argv[10]);

  int model=0; //0 - SSM Z' 1 - Graviton
  if(argc==12){
    model=atoi(argv[11]);
  }
	
	
  bool doXsecLimit=true;		// true=set xsec limits [pb] false=set N_zprime limits [counts]
  double lumiRatio=1.0;			// only relevant if limit is set in units of N_zprime
  double maxNsig=100.0;			// max parameter range of signal (in counts)
  bool activateEE=doElectron;		// enable EE channel
  bool activateMM=doMuon;		// enable MM channel


  cerr<<"doEnsemble: "<<doEnsemble<<"\n";
  cerr<<"nEnsemble: "<<nEnsemble<<"\n";
  cerr<<"iMass: "<< iMass <<"\n";
  cerr<<"doPlots: "<< doPlots <<"\n";
  cerr<<"doMCMC: "<< doMCMC <<"\n";
  cerr<<"doSys: "<< doSys <<"\n";
  cerr<<"doTheoryUncertainty: "<< doTheoryUncertainty <<"\n";
  cerr<<"doElectron: "<<doElectron<<"\n";
  cerr<<"doMuon: "<<doMuon<<"\n";
  cerr<<"lumiRatio: "<<lumiRatio<<"\n";
  cerr<<"maxNsig: "<<maxNsig<<"\n";
  cerr<<"model: "<<model<<"\n\n\n";


  // Define label based on settings
  std::string labelTxt = ""; 
  if (activateEE && !activateMM)
    labelTxt = "_ee";
  if (activateMM && !activateEE)
    labelTxt = "_mm";
  if (activateMM && activateEE)
    labelTxt = "_comb";
	
  if (doSys)
    labelTxt +="Sys";  

  if (doTheoryUncertainty)
    labelTxt +="_inclTheo";  

  // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();

  // remember old directory
  TDirectory* f = gDirectory;

  // ----------------------------------------------------
  // Get electron templates
  // ----------------------------------------------------

  // open ele template file
  TFile * file_ee;
  if(model==0){
    file_ee= new TFile("../input/Template2DObjee_0710_130.root", "READ");
  } 
  else if(model==1){
    file_ee= new TFile("../input/Gee_Template.root", "READ");
  }
  // check if file is open
  if (!file_ee->IsOpen()) {
    cout << "Could not open file. Exit." << endl;
    return 1;
  }

  TObjArray* obj = new TObjArray();
  obj->Read("template");

  // Get ele background 
  TFile * file_ee_sm = new TFile("../input/Bkg_Template_0710_130.root","READ");
  TH1D h_ee_bkg = *((TH1D*) file_ee_sm->Get("FullBackgroundTeV"));

  // Get ele data
  TH1D h_ee_data = *((TH1D*) file_ee_sm->Get("DataTeV"));

  //skip for now
  //ele resolution systematics
  // TFile * file_sys = new TFile("../input/Template2DObjee_res_smearup.root","READ");
  // file_sys->cd();
  // TObjArray* ores = new TObjArray();
  // ores->Read("res");
  // TH1D* h_ee_sys_res =((TH1D*)(TObjArray*)ores->At(iMass));

  //ele QCD background systematics
  //TFile * file_bkgsys = new TFile("../input/bkg_sys_ee_2011_0710_130.root","READ");
  //TH1D* h_sys_bkg =((TH1D*) file_bkgsys->Get("FullBackgroundTeV"));

  // close file
  file_ee->Close();

  // ----------------------------------------------------
  // Get ele signal template
  // ----------------------------------------------------
  // mass = 0.04 * iMass + 0.13 [GeV/c^2];
  TH1D* h_ee_sgn = ((TH1D*)(TObjArray*)obj->At(iMass));

  delete file_ee;


  // ----------------------------------------------------
  // Get muon templates
  // ----------------------------------------------------

  // open muon template file
  TFile * file_mm;
  if(model==0){
    file_mm= new TFile("../input/Template2DObjmm_130.root","READ");//Template2DObjmm_130.root", "READ");
  }
  else if(model==1){
    file_mm= new TFile("../input/GravTemplatesMM_130.root", "READ");
  }
  // check if file is open
  if (!file_mm->IsOpen()) {
    cout << "Could not open file. Exit." << endl;
    return 1;
  }

  TObjArray* obj_mm = new TObjArray();
  obj_mm->Read("template");

  // Get total background 
  TFile * file_mm_sm = new TFile("../input/totalmc_130.root", "READ");
  TH1D h_mm_bkg = *((TH1D*) file_mm_sm->Get("TotalMC"));

  // Get data
  TFile * file_mm_data = new TFile("../input/muonZprime_130.root", "READ");
  TH1D h_mm_data = *((TH1D*) file_mm_data->Get("Data"));

/*// Neglect resolution systematics for now

  // open resolution systematic for each test mass
  TFile * file_mm_sys;
  if(model==0){
    file_mm_sys = new TFile("../input/resSyst2_130.root", "READ");
  }
  else if(model==1){
    file_mm_sys= new TFile("../input/Gmm_resSyst_130.root", "READ"); 
  }
  TObjArray* ores_mm = new TObjArray();
  ores_mm->Read("res");
  TH1D* h_mm_sys_res = ((TH1D*)(TObjArray*)ores_mm->At(iMass));
*/
/*// Offset systematics are negligible

  // open offset systematic for each test mass
  TFile * file_offsys = new TFile("../input/resOffsetSyst2_130.root", "READ");
  TObjArray* ooff = new TObjArray();
  ooff->Read("res");
  TH1D* h_sys_off = ((TH1D*)(TObjArray*)ooff->At(iMass));

  // open offset systematic for background
  TFile * file_offbkgsys = new TFile("../input/resOffsetBkgSyst2_130.root", "READ");
  TH1D h_sys_bkgoff = *((TH1D*) file_offbkgsys->Get("BkgOffset"));
*/

  // close files
  file_mm->Close();
  file_mm_sm->Close();
//  file_mm_sys->Close();
//  file_offsys->Close();
//  file_offbkgsys->Close();


  // ----------------------------------------------------
  // Get signal template
  // ----------------------------------------------------
  // mass = 0.04 * iMass + 0.13 [GeV/c^2];
  TH1D* h_mm_sgn = ((TH1D*)(TObjArray*)obj_mm->At(iMass));

  // go back to old directory for memory handling
  f->cd();

  // ----------------------------------------------------
  // configure BAT
  // ----------------------------------------------------
  // set nice style for drawing than the ROOT default
  BCAux::SetStyle();
  std::ostringstream iMassTxt; iMassTxt << iMass <<"_run"<<iRun;
  // open log file
  BCLog::OpenLog(("log_mass"+iMassTxt.str()+labelTxt+".txt").c_str());
  BCLog::SetLogLevel(BCLog::summary);

  // ----------------------------------------------------
  // Normalization
  // ----------------------------------------------------
  double Neebkg = h_ee_bkg.Integral();
  double Nmmbkg = h_mm_bkg.Integral();
	
  // ----------------------------------------------------
  // create new BCMultiTemplateFitter object
  // ----------------------------------------------------
  BCMultiTemplateFitter * m = new BCMultiTemplateFitter();

  m->MCMCSetPrecision(BCEngineMCMC::kMedium);
//  m->MCMCSetPrecision(BCEngineMCMC::kHigh);

  BCLog::OutSummary("Test model created");

  // create a new summary tool object
  BCSummaryTool * summary = new BCSummaryTool(m);

  // ----------------------------------------------------
  // Define required input
  // ----------------------------------------------------
  double mass = 0.04 * iMass + 0.13;
  if(model==1){
    double temp_mass[10]={0.3, 0.5, 0.7, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25};
    mass=temp_mass[iMass];
  }
  
  TF1* accEE=0;
  if (model==0){
    accEE= new TF1("accEE","pol6",0.07,3.07);
    accEE->SetParameters(0.230274,1.26172,-1.49026,0.844042,-0.215362,0.013045,0.00194791);       // for cut at 130 GeV
  } 
  else if (model==1){
    accEE= new TF1("accEE","pol1",0.07,3.07);
    accEE->SetParameters(0.701244,-0.0107504);
  }
  
  TF1* accMM=0;
  if (model==0){
    accMM= new TF1("accMM","pol6",0.07,3.07);
    accMM->SetParameters(0.248482,0.613077,-0.960227,0.774202,-0.339321,0.0756178,-0.00670578);   // for cut at 130 GeV
  } 
  else if (model==1){
    accMM= new TF1("accMM","pol1",0.07,3.07);
    accMM->SetParameters(0.455201,-0.0114039);
  }
  
  
  double zxsec = 0.98*989; // NNLO x-sec for m>70 GeV
  double NZee = 2933.3;    // found Z->ee's - other backgrounds
  double NZmm = 2696.5;    // found Z->mm's - other backgrounds
  double Aee = 0.0034008;  // Acceptance above 130 GeV from Zee inclusive sample
  double Amm = 0.0026469;  // Acceptance above 130 GeV from Zmm inclusive sample

  double maxSigmaSig = maxNsig * zxsec * Aee / ( NZee * accEE->Eval(mass) );
  cout<<"max "<<maxSigmaSig<<endl;

  double maxSig = maxNsig;
  if (doXsecLimit) maxSig = maxSigmaSig;


  m->AddChannel("ee");
  m->AddChannel("mm");

  m->AddProcess("eeBkg", Neebkg, Neebkg);
  m->AddProcess("mmBkg", Nmmbkg, Nmmbkg);
  if (doXsecLimit)	
    m->AddProcess("signal",   0., maxSigmaSig);
  else
    m->AddProcess("signal",   0., maxNsig);
	
  m->SetData("ee", h_ee_data); 
  m->SetData("mm", h_mm_data); 
        	
  m->SetTemplate("ee", "eeBkg", h_ee_bkg, 1.0); 
  m->SetTemplate("mm", "mmBkg", h_mm_bkg, 1.0); 
	
  if (doXsecLimit)
    {
      m->SetTemplate("ee", "signal", obj, ( accEE->Eval(mass)  * NZee * lumiRatio ) / (Aee * zxsec) ); 
      m->SetTemplate("mm", "signal", obj_mm , accMM->Eval(mass) * NZmm / (Amm * zxsec )); 
    }
  else
    { 
      m->SetTemplate("ee", "signal", obj, 1.0 * lumiRatio); 
      m->SetTemplate("mm", "signal", obj_mm, 1.0); 
    }  

  // ----------------------------------------------------
  // Specify active channels and priors
  // ----------------------------------------------------
  m->GetChannel(0)->SetFlagChannelActive(activateEE);
  m->GetChannel(1)->SetFlagChannelActive(activateMM);

  // set priors
  m->SetPriorGauss("eeBkg", Neebkg, 0.1);
  m->SetPriorGauss("mmBkg", Nmmbkg, 0.1);
  m->SetPriorConstant("signal");
	
  //increase number of marginalization bins (higher precision)
  m->SetNbins("signal",1000);
  m->SetNbins("Mass",150);


  if (doSys)
    {
      // ----------------------------------------------------
      // Specify systematic uncertainties
      // ----------------------------------------------------

      m->AddSystematic("EFF",-5,5);m->SetPriorGauss("EFF", 0., 1.);
      m->AddSystematic("PDF",-5,5);m->SetPriorGauss("PDF", 0., 1.);
      m->AddSystematic("KEWK",-5,5);m->SetPriorGauss("KEWK", 0., 1.);
      //m->AddSystematic("KQCD",-5,5);m->SetPriorGauss("KQCD", 0., 1.);
      
      //m->AddSystematic("EERES",-5,5);m->SetPriorGauss("EERES", 0., 1.);
      //m->AddSystematic("MMRES",-5,5);m->SetPriorGauss("MMRES", 0., 1.);
      //m->AddSystematic("SIGOFF",-5,5);m->SetPriorGauss("SIGOFF", 0., 1.);
      //m->AddSystematic("BKGOFF",-5,5);m->SetPriorGauss("BKGOFF", 0., 1.);
      m->AddSystematic("ZXSEC",-5,5);m->SetPriorGauss("ZXSEC", 0., 1.);
      //m->AddSystematic("MODELQCD",-5,5);m->SetPriorGauss("MODELQCD", 0., 1.);
      //m->AddSystematic("EEISOEFF",-5,5);m->SetPriorGauss("EEISOEFF", 0., 1.);

      // ----------------------------------------------------
      // Make up histograms (BAT wants TH1Ds)
      // ----------------------------------------------------
      TH1D h_sys_eff  = TH1D(*(h_mm_sgn));	  // const slope efficiency uncertainty
      TH1D h_sys_pdfZ = TH1D(*(h_ee_sgn));	  // PDF systematic
      TH1D h_sys_pdfG = TH1D(*(h_ee_sgn));	  // PDF systematic
      TH1D h_sys_kEWK = TH1D(*(h_ee_sgn));	  // EWK K-factor
      TH1D h_sys_kQCD = TH1D(*(h_ee_sgn));	  // QCD K-factor
      TH1D h_sys_zXsec = TH1D(*(h_ee_sgn));	  // Z boson theory xsec uncertainty from ratio method
      TH1D h_sys_eeIsoEff = TH1D(*(h_ee_sgn));	  // isolation in ee channel

      //Functions to map out linear increase of systematic vs mass (relative uncertainty)
      TF1 *fsys_eff  = new TF1("fsys_eff" ,"0.03*x",0.0,3.0);
      TF1 *fsys_pdfG = new TF1("fsys_pdfG","0.0641*x+0.0320",0.0,3.0);
      TF1 *fsys_kEWK = new TF1("fsys_kEWK","0.03*x",0.0,3.0);
      TF1 *fsys_eeIsoEff = new TF1("fsys_fsys_eeIsoEff","0.01*x",0.0,3.0);
	  
      //PDF, alpha_s and kQCD (=scale) error TGraph, input from T. Nunnemann
      double Vmass[39] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.091, 0.1, 0.125, 0.15, 0.175, 
      			  0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 
			  2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5};
      double VPDF[39] = {0.161273, 0.076909, 0.0583181, 0.0517107, 0.0485283, 0.047392, 0.0457165, 0.0449222, 
      			   0.0440114, 0.0434281, 0.0421782, 0.0410122, 0.0399375, 0.0399249, 0.0397618, 0.041, 
			   0.0428719, 0.0454423, 0.0488774, 0.0524023, 0.0572276, 0.0628967, 0.0695557, 0.0887525, 
			   0.115953, 0.151611, 0.19758, 0.25315, 0.313861, 0.377535, 0.442274, 0.507875, 0.576834, 
			   0.648193, 0.726607, 0.818732, 0.945838, 1.20089, 2.00758};
      TGraph *graph_pdfZ = new TGraph(39, Vmass, VPDF);


      //Fill histograms with relative uncertainties (in %)
      for (int i=1; i<=h_sys_eff.GetNbinsX();i++){
	double mass = h_sys_eff.GetBinCenter(i);
	// all uncertainties cancel to first order in norm region (70,110GeV)
	if (h_sys_eff.GetBinCenter(i) < 0.11){
	  h_sys_zXsec.SetBinContent(i,0.0);
	  h_sys_eff.SetBinContent(i,0.0);
	  h_sys_pdfZ.SetBinContent(i,0.0);
	  h_sys_pdfG.SetBinContent(i,0.0);
	  h_sys_kEWK.SetBinContent(i,0.0);
	  //h_sys_kQCD.SetBinContent(i,0.0);  //Now included in pdfZ
	  h_sys_eeIsoEff.SetBinContent(i,0.0);
	}
	// fill mass dependent uncertainties
	else{
	  h_sys_zXsec.SetBinContent(i,0.05); // flat 5%
	  h_sys_eff.SetBinContent(i,fsys_eff->Eval(mass));
          h_sys_pdfZ.SetBinContent(i,graph_pdfZ->Eval(mass, 0, "S"));
	  h_sys_pdfG.SetBinContent(i,fsys_pdfG->Eval(mass));
	  h_sys_kEWK.SetBinContent(i,fsys_kEWK->Eval(mass));
	  //h_sys_kQCD.SetBinContent(i,fsys_kQCD->Eval(mass));  //Now included in pdfZ
	  h_sys_eeIsoEff.SetBinContent(i,fsys_eeIsoEff->Eval(mass));
	}
      }
	  
      //EFF (only mm)
      m->SetSystematicVariation("mm", "mmBkg", "EFF", h_sys_eff,h_sys_eff);
      m->SetSystematicVariation("mm", "signal", "EFF", h_sys_eff,h_sys_eff);
	  
      //Isolation efficiency (only ee)
      //m->SetSystematicVariation("ee", "signal", "EEISOEFF", h_sys_eff,h_sys_eff);

      //PDF
      m->SetSystematicVariation("ee", "eeBkg", "PDF", h_sys_pdfZ,h_sys_pdfZ);
      if (doTheoryUncertainty){
	if (model==0) m->SetSystematicVariation("ee", "signal", "PDF", h_sys_pdfZ,h_sys_pdfZ);
	if (model==1) m->SetSystematicVariation("ee", "signal", "PDF", h_sys_pdfG,h_sys_pdfG);
      }
      m->SetSystematicVariation("mm", "mmBkg", "PDF", h_sys_pdfZ,h_sys_pdfZ);
      if (doTheoryUncertainty){
	if (model==0) m->SetSystematicVariation("mm", "signal", "PDF", h_sys_pdfZ,h_sys_pdfZ);
	if (model==1) m->SetSystematicVariation("mm", "signal", "PDF", h_sys_pdfG,h_sys_pdfG);
      }

/*    // Now included in pdfZ
      //K-factor (QCD)
      m->SetSystematicVariation("ee", "eeBkg", "KQCD", h_sys_kQCD,h_sys_kQCD);
      if (doTheoryUncertainty && model==0)
	m->SetSystematicVariation("ee", "signal", "KQCD", h_sys_kQCD,h_sys_kQCD);
	    
      m->SetSystematicVariation("mm", "mmBkg", "KQCD", h_sys_kQCD,h_sys_kQCD);
      if (doTheoryUncertainty && model==0)
	m->SetSystematicVariation("mm", "signal", "KQCD", h_sys_kQCD,h_sys_kQCD);
*/
      //K-factor (EWK)
      m->SetSystematicVariation("ee", "eeBkg", "KEWK", h_sys_kEWK,h_sys_kEWK);
      m->SetSystematicVariation("mm", "mmBkg", "KEWK", h_sys_kEWK,h_sys_kEWK);

      //Resolution smearing (ee, mm uncorrelated)
      //m->SetSystematicVariation("ee", "signal", "EERES", *h_ee_sys_res,*h_ee_sys_res);  //off -under 3%
      //m->SetSystematicVariation("mm", "signal", "MMRES", *h_mm_sys_res,*h_mm_sys_res);  //off -under 3%

      //Momentum offset (mm)
      //m->SetSystematicVariation("mm", "signal", "SIGOFF", *h_sys_off,   *h_sys_off);    //off -under 3%
      //m->SetSystematicVariation("mm", "mmBkg",  "BKGOFF",  h_sys_bkgoff, h_sys_bkgoff); //off -under 3%
          
      //Flat DY xsec uncertainty
      m->SetSystematicVariation("ee", "signal", "ZXSEC", h_sys_zXsec,h_sys_zXsec);
      m->SetSystematicVariation("mm", "signal", "ZXSEC", h_sys_zXsec,h_sys_zXsec);

      //QCD background systematic
      //m->SetSystematicVariation("ee", "eeBkg", "MODELQCD", *h_sys_bkg,*h_sys_bkg);

    }
		
		

  // ----------------------------------------------------
  // perform analysis
  // ----------------------------------------------------
  if (!doEnsemble){


    // run MCMC
    m->SetData("ee", h_ee_data); 
    m->SetData("mm", h_mm_data); 
    m->MarginalizeAll();

    // find global mode
    //m->FindMode();
    m->FindMode( m->GetBestFitParameters() );
	   
    // print a summary
    m->PrintSummary();
    //double mass = 0.04 * iMass + 0.13;

    // ----------------------------------------------------
    // Print results a) 2D scan and result of test statistic
    // ----------------------------------------------------
    // S+B hypothesis
    vector<double> test_parameters;
    test_parameters = m->GetBestFitParameters();
    cout<<"LogLike "<<m->LogLikelihood(test_parameters)<<endl;
    cout<<"S+B "<<exp(m->LogLikelihood(test_parameters))<<endl;

    cout << "Test parameters:" << endl;
    cout << "  NmmBkg = " << test_parameters[0] << endl;
    cout << "  NeeBkg = " << test_parameters[1] << endl;
    if (doXsecLimit)	
      cout << "  SignalXsec = " << test_parameters[2] << " pb" << endl;
    else
      cout << "  NSignal = " << test_parameters[2] << endl;
    cout << "  At Mass = " << test_parameters[3] * 40 + 130 << " GeV" << endl;


    if (doPlots){
      m->PrintAllMarginalized(("model_marginalized_mass"+iMassTxt.str()+labelTxt+".eps").c_str()); 
      // print all summary plots
      summary->PrintParameterPlot(("BCMultiTemplateFitter_parameters"+iMassTxt.str()+labelTxt+".eps").c_str());
      summary->PrintCorrelationPlot(("BCMultiTemplateFitter_correlation"+iMassTxt.str()+labelTxt+".eps").c_str());
      summary->PrintKnowledgeUpdatePlots(("BCMultiTemplateFitter_update"+iMassTxt.str()+labelTxt+".ps").c_str());

      // print results of the analysis into a text file
      m->PrintResults(("BCMultiTemplateFitter_results"+iMassTxt.str()+labelTxt+".txt").c_str());
    }

    TH2D* signalVSmass = new TH2D();
    TFile* fout = new TFile(("usefulPlots_pvalue"+labelTxt+".root").c_str(),"RECREATE");
    signalVSmass = m->GetMarginalized("signal","Mass")->GetHistogram();
    fout->cd();
    signalVSmass->Write();
    fout->Close();

    // B only (Null) hypothesis
    m->SetParameterRange(2,0.0,0.0); // set signal to 0
    // run MCMC
    m->MarginalizeAll();	 
    // find global mode
    m->FindMode();
    vector<double> null_parameters; 
    null_parameters = m->GetBestFitParameters();
    for (unsigned int i = 0;i<null_parameters.size();i++)
    {
       cout<<i<<" null parameter "<<null_parameters[i]<<endl;
    }
    //null_parameters.push_back(Nbkg);// bkg norm
    cout<<"LogLike "<<m->LogLikelihood(null_parameters)<<endl;
    cout<<"B "<<exp(m->LogLikelihood(null_parameters))<<endl;

    double TSobs = -2*(m->LogLikelihood(test_parameters) - m->LogLikelihood(null_parameters));
    cout << "TS " << TSobs << endl;

    if (doPlots && activateEE){
      double pvalue = plot_pvalue(TSobs, 1, 0);
      plot2Dscan (TSobs, pvalue, 1, 0);
    }    
    if (doPlots && activateMM){
      double pvalue = plot_pvalue(TSobs, 0, 1);
      plot2Dscan (TSobs, pvalue, 0, 1);       
    }
    
    return 0;
  }

  // ----------------------------------------------------
  // create prior model for pseudo-experiments
  // ----------------------------------------------------

  // create new prior model
  BCSummaryPriorModel* pm = new BCSummaryPriorModel();
  // set model (and make adjustment suitable for background only pseudo-experiments (PE))
  m->SetParameterRange(0,Neebkg,Neebkg);
  m->SetParameterRange(1,Nmmbkg,Nmmbkg);
  m->SetParameterRange(2,0.0,0.0);//set signal to zero here (flat prior would incl signal in PE)
  m->SetParameterRange(3,0.0,0.0);//set signal mass to zero
  pm->SetModel(m);

  // ----------------------------------------------------
  // create output object
  // ----------------------------------------------------
  BCModelOutput* pmout = new BCModelOutput(pm, ("prior"+iMassTxt.str()+labelTxt+".root").c_str());

  // switch writing of Markov Chains on
  pmout->WriteMarkovChain(true);

  // set precision
  pm->MCMCSetPrecision(BCEngineMCMC::kMedium);
  //pm->MCMCSetPrecision(BCEngineMCMC::kHigh);

  // perform marginalization
  pm->MarginalizeAll(); 

  // get tree
  TTree* priortree = (TTree*) pmout->GetFile()->Get("MarkovChainTree_0");

  //undo PE modifications
  m->SetParameterRange(0,Neebkg,Neebkg);
  m->SetParameterRange(1,Nmmbkg,Nmmbkg);
  if (doXsecLimit)	
    m->SetParameterRange(2,0.0,maxSigmaSig);
  else
    m->SetParameterRange(2,0.0,maxNsig);

  m->SetParameterRange(3,0.0,73);


  // ----------------------------------------------------------------
  // Perform ensemble test
  // ----------------------------------------------------------------
  // create new analysis facility 
  BCMTFAnalysisFacility* facility = new BCMTFAnalysisFacility(m); 
  facility->SetFlagMCMC(doMCMC);

  // create ensembles
  //	TTree* tree = facility->BuildEnsembles( m->GetBestFitParameters(), 10000 );
  TTree* tree = facility->BuildEnsembles( priortree, 10000 );

  // run ensemble test
  TTree* tree_out = facility->PerformEnsembleTest(tree, nEnsemble, maxSig);
	
  // open new file
  TFile *file = new TFile(("zprime_ensembles_pvalue_"+iMassTxt.str()+labelTxt+".root").c_str(), "RECREATE");
  file->cd(); 

  // write trees into file
  tree->Write();
  tree_out->Write(); 

  // close file
  file->Close();
	
  // free memory
  delete file;

  // -----------

  // close log file
  BCLog::CloseLog();

  // close output file
  pmout->Close();

  // free memory
  delete pm;
  delete pmout;
  delete facility;
  delete m;
  delete summary;

  BCLog::OutSummary("Test program ran successfully");
  BCLog::OutSummary("Exiting");

  // no error
  return 0;

}

