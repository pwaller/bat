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

#include <BCMTFAnalysisFacility.h>
#include <BCMultiTemplateFitter.h>
#include <BCChannel.h>

#include <TROOT.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <stdlib.h>
 
using std::cerr;
using std::cout;
using std::endl;
using std::vector;


int main(int argc, char ** argv){
  
  if(argc<11){
    cerr<<"usage:  templateLimitBATcombSys.exe iRun doEnsemble nEnsemble iMass doPlots doMCMC doSYS doTheory doElectron doMuon [model]\n";
    exit(1);
  }
      

  std::cerr << "++ Running: " << argv[0] << std::endl;
  //std::cout << "++ with iMass: " << argv[1] << " iRun "<< argv[2] << std::endl;
	
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

  int model=1; //0 - SSM Z' 1 - Graviton
  if(argc==12){
    model=atoi(argv[11]);
  }
	
	
  bool doXsecLimit=true;		// true=set xsec limits [pb] false=set N_zprime limits [counts]
  double lumiRatio=1.0;			// only relevant if limit is set in units of N_zprime
  double maxNsig=300.0;			// max parameter range of signal (in counts)
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

  cerr<<"ee template file: "<<file_ee->GetName()<<"\n";

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
  TFile * file_bkgsys = new TFile("../input/bkg_sys_ee_2011_0710_130.root","READ");
  TH1D* h_sys_bkg =((TH1D*) file_bkgsys->Get("FullBackgroundTeV"));

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
    file_mm= new TFile("../input/Template2DObjmm_130.root","READ");
  }
  else if(model==1){
    file_mm= new TFile("../input/GravTemplatesMM_130.root", "READ");
  }
  // check if file is open
  if (!file_mm->IsOpen()) {
    cout << "Could not open file. Exit." << endl;
    return 1;
  }

  cerr<<"mm template file: "<<file_mm->GetName()<<"\n";

  TObjArray* obj_mm = new TObjArray();
  obj_mm->Read("template");

  // Get total background 
  TFile * file_mm_sm = new TFile("../input/totalmc_130.root", "READ");
  TH1D h_mm_bkg = *((TH1D*) file_mm_sm->Get("TotalMC"));

  // Get data
  TFile * file_mm_data = new TFile("../input/muonZprime_130.root", "READ");
  TH1D h_mm_data = *((TH1D*) file_mm_data->Get("Data"));

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

/*// Offset systematics are negligible

  // open offset systematic for each test mass
  TFile * file_mm_offsys;
  if(model==0){
    file_mm_offsys = new TFile("../input/resOffsetSyst2_130.root", "READ");
  }
  else if(model==1){
    file_mm_offsys= new TFile("../input/Gmm_resOffsetSyst_130.root", "READ");
  }
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
  file_mm_sys->Close();
//  file_mm_offsys->Close();
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
  double Aee = 0.0034008;  // Acceptance from Zee inclusive sample
  double Amm = 0.0026469;  // Acceptance from Zmm inclusive sample

  double maxSigmaSig = maxNsig * zxsec * Aee / ( NZee * accEE->Eval(mass) );
  std::cerr<<"mass: "<<mass<<"\n";
  std::cerr<<"ACCee: "<<accEE->Eval(mass)<<"\n";
  std::cerr<<"ACCmm: "<<accMM->Eval(mass)<<"\n";
  
  std::cerr<<"maxSigmaSig "<<maxSigmaSig<<std::endl;
  if(mass>0.79){
    maxSigmaSig=0.1;
  }
  if(mass>0.99){
    maxSigmaSig=0.05;
  }
  std::cerr<<"after lowering maxSigmaSig: "<<maxSigmaSig<<std::endl;

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
      m->SetTemplate("ee", "signal", *h_ee_sgn, ( accEE->Eval(mass)  * NZee * lumiRatio ) / (Aee * zxsec) ); 
      m->SetTemplate("mm", "signal", *h_mm_sgn, accMM->Eval(mass) * NZmm / (Amm * zxsec )); 
    }
  else
    { 
      m->SetTemplate("ee", "signal", *h_ee_sgn, 1.0 * lumiRatio); 
      m->SetTemplate("mm", "signal", *h_mm_sgn, 1.0); 
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
  m->SetNbins("signal",2000);
  //m->SetNbins("eeBkg",100);
  //m->SetNbins("mmBkg",100);


  if (doSys)
    {
      // ----------------------------------------------------
      // Specify systematic uncertainties
      // ----------------------------------------------------

      m->AddSystematic("EFF",-5,5);m->SetPriorGauss("EFF", 0., 1.);
      m->AddSystematic("PDF",-5,5);m->SetPriorGauss("PDF", 0., 1.);
      m->AddSystematic("KEWK",-5,5);m->SetPriorGauss("KEWK", 0., 1.);
      m->AddSystematic("KQCD",-5,5);m->SetPriorGauss("KQCD", 0., 1.);
      //m->AddSystematic("EERES",-5,5);m->SetPriorGauss("EERES", 0., 1.);
      //m->AddSystematic("MMRES",-5,5);m->SetPriorGauss("MMRES", 0., 1.);
      //m->AddSystematic("SIGOFF",-5,5);m->SetPriorGauss("SIGOFF", 0., 1.);
      //m->AddSystematic("BKGOFF",-5,5);m->SetPriorGauss("BKGOFF", 0., 1.);
      m->AddSystematic("ZXSEC",-5,5);m->SetPriorGauss("ZXSEC", 0., 1.);
      m->AddSystematic("MODELQCD",-5,5);m->SetPriorGauss("MODELQCD", 0., 1.);
      //m->AddSystematic("EEISOEFF",-5,5);m->SetPriorGauss("EEISOEFF", 0., 1.);

      // ----------------------------------------------------
      // Make up histograms (BAT wants TH1Ds)
      // ----------------------------------------------------
      TH1D h_sys_eff  = TH1D(*(h_ee_sgn));	  // const slope efficiency uncertainty
      TH1D h_sys_pdfZ = TH1D(*(h_ee_sgn));	  // PDF systematic
      TH1D h_sys_pdfG = TH1D(*(h_ee_sgn));	  // PDF systematic
      TH1D h_sys_kEWK = TH1D(*(h_ee_sgn));	  // EWK K-factor
      TH1D h_sys_kQCD = TH1D(*(h_ee_sgn));	  // QCD K-factor
      TH1D h_sys_zXsec = TH1D(*(h_ee_sgn));	  // Zboson theory xsec uncertainty from ratio method
      TH1D h_sys_eeIsoEff= TH1D(*(h_ee_sgn)); // isolation in ee channel

      //Functions to map out linear increase of systematic vs mass (relative uncertainty)
      TF1 *fsys_eff  = new TF1("fsys_eff" ,"0.03*x",0.0,3.0);
      TF1 *fsys_pdfZ = new TF1("fsys_pdfZ","0.0592*x+0.0176",0.0,3.0);
      TF1 *fsys_pdfG = new TF1("fsys_pdfG","0.0641*x+0.0320",0.0,3.0);
      TF1 *fsys_kEWK = new TF1("fsys_kEWK","0.03*x",0.0,3.0);
      TF1 *fsys_kQCD = new TF1("fsys_kQCD","0.02*x",0.0,3.0);
      TF1 *fsys_eeIsoEff = new TF1("fsys_fsys_eeIsoEff","0.01*x",0.0,3.0);
	  
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
	  h_sys_kQCD.SetBinContent(i,0.0);
	  h_sys_eeIsoEff.SetBinContent(i,0.0);
	}
	// fill mass dependent uncertainties
	else{
	  h_sys_zXsec.SetBinContent(i,0.05); // flat 5%
	  h_sys_eff.SetBinContent(i,fsys_eff->Eval(mass));
	  h_sys_pdfZ.SetBinContent(i,fsys_pdfZ->Eval(mass));
	  h_sys_pdfG.SetBinContent(i,fsys_pdfG->Eval(mass));
	  h_sys_kEWK.SetBinContent(i,fsys_kEWK->Eval(mass));
	  h_sys_kQCD.SetBinContent(i,fsys_kQCD->Eval(mass));
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

      //K-factor (QCD)
      m->SetSystematicVariation("ee", "eeBkg", "KQCD", h_sys_kQCD,h_sys_kQCD);
      if (doTheoryUncertainty && model==0)
	m->SetSystematicVariation("ee", "signal", "KQCD", h_sys_kQCD,h_sys_kQCD);
	    
      m->SetSystematicVariation("mm", "mmBkg", "KQCD", h_sys_kQCD,h_sys_kQCD);
      if (doTheoryUncertainty && model==0)
	m->SetSystematicVariation("mm", "signal", "KQCD", h_sys_kQCD,h_sys_kQCD);

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
      m->SetSystematicVariation("ee", "eeBkg", "MODELQCD", *h_sys_bkg,*h_sys_bkg);

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
	   
    // print templates and stacks
    //for (int i = 0; i < m->GetNChannels(); ++i) {
    //	BCChannel* channel = m->GetChannel(i);
    //	channel->PrintTemplates(Form("%s_templates.ps", channel->GetName().c_str())); 
    //	m->PrintStack(i, m->GetBestFitParameters(), Form("%s_stack.eps", channel->GetName().c_str()), "logxlogy"); 
    //   } 

    // print a summary
    m->PrintSummary();
    //double mass = 0.04 * iMass + 0.13;
    std::cout<<"limit "<<mass<<" "<<m->GetMarginalized("signal")->GetLimit(0.95)<<std::endl;

    // create summary tool
    //BCSummaryTool* st = new BCSummaryTool(model); 
	
    if (doPlots){
      m->PrintAllMarginalized(("model_marginalized_mass"+iMassTxt.str()+labelTxt+".eps").c_str()); 
      // print all summary plots
      summary->PrintParameterPlot(("BCMultiTemplateFitter_parameters"+iMassTxt.str()+labelTxt+".eps").c_str());
      summary->PrintCorrelationPlot(("BCMultiTemplateFitter_correlation"+iMassTxt.str()+labelTxt+".eps").c_str());
      summary->PrintKnowledgeUpdatePlots(("BCMultiTemplateFitter_update"+iMassTxt.str()+labelTxt+".ps").c_str());

      // print results of the analysis into a text file
      m->PrintResults(("BCMultiTemplateFitter_results"+iMassTxt.str()+labelTxt+".txt").c_str());
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
  pm->SetModel(m);

  // ----------------------------------------------------
  // create output object
  // ----------------------------------------------------
  BCModelOutput* pmout = new BCModelOutput(pm, ("prior"+iMassTxt.str()+labelTxt+".root").c_str());

  // switch writing of Markov Chains on
  pmout->WriteMarkovChain(true);

  // set precision
  pm->MCMCSetPrecision(BCEngineMCMC::kMedium);

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
  TTree* tree_out = facility->PerformEnsembleTest(tree, nEnsemble);
	
  // open new file
  TFile *file = new TFile(("zprime_ensembles_mass"+iMassTxt.str()+labelTxt+".root").c_str(), "RECREATE");
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

