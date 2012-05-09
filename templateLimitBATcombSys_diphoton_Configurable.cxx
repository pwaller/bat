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

int main(int argc, char ** argv){
  
  if(argc<12){
    cerr<<"usage:  templateLimitBATcombSys.exe iRun doEnsemble nEnsemble iMass doPlots doMCMC doSYS doTheory doElectron doMuon doPhoton [model] outputdir\n";
    exit(1);
  }
      

  std::cerr << "++ Running: " << argv[0] << std::endl;
  //std::cout << "++ with iMass: " << argv[1] << " iRun "<< argv[2] << std::endl;
	
  //Some flags to set (could come from cmd line)
  int iRun = atoi(argv[1]);   		// iRun number (for parallel running)
  bool doEnsemble=atoi(argv[2]);//true;			// true=doPEs, false=data analysis
  int nEnsemble=atoi(argv[3]);//150;	// Number of pseudo-experiments (PEs) if doEnsemble is true
  int iMass = atoi(argv[4]);              // template number
  bool doPlots=atoi(argv[5]);//true;	// Make plots of posterior PDF etc, for data analysis
  bool doMCMC=atoi(argv[6]);//true;	// only relevant for PEs (use MarkovChain MC for marginalization or Minuit profiling)
  bool doSys =atoi(argv[7]);//true;
  bool doTheoryUncertainty =atoi(argv[8]);// false= Z' theory uncertainty outside of likelihood , true= Z' theory uncertainty included in likelihood function 
  bool doElectron=atoi(argv[9]);
  bool doMuon=atoi(argv[10]);
  bool doPhotons=atoi(argv[11]);
  TString outputdir=argv[12];
  TString configuration=argv[13];
  

  int model=1; 
	
	
  bool doXsecLimit=true;			// true=set xsec limits [pb] false=set N_zprime limits [counts]
  double lumiRatio=1.0;			// only relevant if limit is set in units of N_zprime
  double maxNsig=100.0;			// max parameter range of signal (in counts)
  bool activateGG=doPhotons;			// enable MM channel

  int lumiGG = 1075.0;       //B_H4
  if(configuration.SubString("B_K6") =="B_K6")
    lumiGG= 2116;
  
  lumiGG = 4910; // Full 2011
  
  cerr<<"doEnsemble: "<<doEnsemble<<"\n";
  cerr<<"nEnsemble: "<<nEnsemble<<"\n";
  cerr<<"iMass: "<< iMass <<"\n";
  cerr<<"doPlots: "<< doPlots <<"\n";
  cerr<<"doMCMC: "<< doMCMC <<"\n";
  cerr<<"doSys: "<< doSys <<"\n";
  cerr<<"doTheoryUncertainty: "<< doTheoryUncertainty <<"\n";
  cerr<<"doElectron: "<<doElectron<<"\n";
  cerr<<"doMuon: "<<doMuon<<"\n";
  cerr<<"doPhotons: "<<doPhotons<<"\n";
  cerr<<"lumiRatio: "<<lumiRatio<<"\n";
  cerr<<"maxNsig: "<<maxNsig<<"\n";
  cerr<<"model: "<<model<<"\n\n\n";
  cerr<<"outputdir: "<<outputdir<<"\n\n\n";
  cerr<<"Configuration: "<<configuration<<"\n\n\n";

  std::string labelTxt = ""; 
  labelTxt = "_gg";
  
  if (doSys)
    labelTxt +="Sys";  

  if (doTheoryUncertainty)
    labelTxt +="_inclTheo";  

  // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();

  // remember old directory
  TDirectory* f = gDirectory;

  //################################################################################################################################
  // ----------------------------------------------------
  // Get photon templates
  // ----------------------------------------------------
  
  // Open Signal Templates#####################################
  TFile * file_gg;
  if(configuration.SubString("IsoTemplate") =="IsoTemplate"){
    if(configuration.SubString("c0.1") =="c0.1"){
      if(configuration.SubString("B_K6") =="B_K6")
	file_gg = new TFile("Ggg_Templates_IsoTemplates_c0.1_B_K6.root", "READ");   //Graviton Model
      else
	file_gg = new TFile("Ggg_Templates_IsoTemplates_c0.1.root", "READ");   //Graviton Model
    }
    else if(configuration.SubString("c0.05") =="c0.05"){
      file_gg = new TFile("Ggg_Templates_IsoTemplates_c0.05.root", "READ");   //Graviton Model
    }
    else if(configuration.SubString("c0.01") =="c0.01"){
      file_gg = new TFile("Ggg_Templates_IsoTemplates_c0.01.root", "READ");   //Graviton Model
    }
  }
  else if(configuration.SubString("Method1") =="Method1"){
    if(configuration.SubString("c0.1") =="c0.1"){
      file_gg = new TFile("Ggg_Templates_Method1_c0.1.root", "READ");   //Graviton Model
    }
  }

  //  check if file is open
  if (!file_gg->IsOpen()) {
    std::cout << "Could not open file -puto-. Exit." << std::endl;
    return 1;
  }
  
  TObjArray* obj_gg = new TObjArray();
  obj_gg->Read("template");

  // close file
  file_gg->Close();

  // Get gg Background template#################################
  TFile * file_gg_data;

  if(configuration.SubString("IsoTemplate") =="IsoTemplate"){
    if(configuration.SubString("B_K6") =="B_K6")
      file_gg_data = new TFile("Bkg_gg_Templates_IsoTemplates_B_K6_newCaloC.root","READ");
    else
      file_gg_data = new TFile("Bkg_gg_Templates_IsoTemplates_newCaloC.root","READ");
  }
  else if(configuration.SubString("Method1") =="Method1"){
    file_gg_data = new TFile("Bkg_gg_Templates_Method1_newCaloC.root","READ");
  }

  TH1D h_gg_bkg = *((TH1D*) file_gg_data->Get("bkg_total_gg"));
  TH1D h_gg_data = *((TH1D*) file_gg_data->Get("hh_data"));
  
  //Photon QCD background systematics
  TH1D hh_bkg_syst_gg   = *((TH1D*) file_gg_data->Get("bkg_total_syst_gg"));


  // ----------------------------------------------------
  // Get photon signal template
  // ----------------------------------------------------
  TH1D* h_gg_sgn = ((TH1D*)(TObjArray*)obj_gg->At(iMass));
  delete file_gg;
  
  
  //################################################################################################################################
  // ----------------------------------------------------
  // configure BAT
  // ----------------------------------------------------
  // set nice style for drawing than the ROOT default
  BCAux::SetStyle();
  std::ostringstream iMassTxt; iMassTxt << iMass <<"_run"<<iRun;
  // open log file
  BCLog::OpenLog(outputdir+("/log_mass"+iMassTxt.str()+labelTxt+".txt").c_str());
  BCLog::SetLogLevel(BCLog::summary);

  // ----------------------------------------------------
  // Normalization
  // ----------------------------------------------------
  double Nggbkg   = h_gg_bkg.Integral();

  if(configuration.SubString("bkgScale") =="bkgScale")
    Nggbkg*=1.3;

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
  double mass = 0.;
  //double temp_mass[47]={ 0.2 , 0.25 , 0.3 , 0.35 , 0.4 , 0.45 , 0.5 , 0.55 , 0.6 , 0.65 , 0.7 , 0.75 , 0.8 , 0.85 , 0.9 , 0.95 , 1 , 1.05 , 1.1 , 1.15 , 1.2 , 1.25 , 1.3 , 1.35 , 1.4 , 1.45 , 1.5 , 1.55 , 1.6 , 1.65 , 1.7 , 1.75 , 1.8 , 1.85 , 1.9 , 1.95 , 2 , 2.05 , 2.1 , 2.15 , 2.2 , 2.25 , 2.3 , 2.35 , 2.4 , 2.45 , 2.5 };
  double temp_mass[71];
  //  for(int ii=0;ii<71;ii++)   temp_mass[ii]= 0.45 + (ii*0.025); //New templates from Evan 
  //mass=temp_mass[iMass];
  mass =  0.45+ (iMass*0.025);

  TF1* accGG=0;
  //Aceptance vs mass parameterization for phootns
  accGG= new TF1("accEE","pol1",0.07,3);
  accGG->SetParameters(0.527305 , 1.56236e-05 ); //new parameterization Note_v1  
  
  double maxSigmaSig = 1;

  
  m->AddChannel("gg");
  m->AddProcess("ggBkg", Nggbkg, Nggbkg);
  
  cerr<<"Predicted ngg bkg "<< Nggbkg <<"\n";
 
  if (doXsecLimit)	
    m->AddProcess("signal",   0., maxSigmaSig);
  else
    m->AddProcess("signal",   0., maxNsig);
	
  m->SetData("gg", h_gg_data); 
  m->SetTemplate("gg", "ggBkg"  , h_gg_bkg, 1.0); 
	
  
  if (doXsecLimit)
    {
      m->SetTemplate("gg", "signal", *h_gg_sgn, ( accGG->Eval(mass * 1000) *lumiGG *lumiRatio ) ); 
    }
  else
    { 
      m->SetTemplate("gg", "signal", *h_gg_sgn, 1.0); 
    }  

  // ----------------------------------------------------
  // Specify active channels and priors
  // ----------------------------------------------------
  m->GetChannel(0)->SetFlagChannelActive(activateGG);

  // set priors
  m->SetPriorGauss("ggBkg", Nggbkg, Nggbkg*0.2);
  m->SetPriorConstant("signal");
	
  //increase number of marginalization bins (higher precision)
  m->SetNbins("signal",2000);
  m->SetNbins("ggBkg_reducible",2000);
  m->SetNbins("ggBkg_irreducible",2000);


  if (doSys &&  (configuration.SubString("NoSyst")!="NoSyst"))
    {
      // ----------------------------------------------------
      // Specify systematic uncertainties
      // ----------------------------------------------------
      m->AddSystematic("gg_EFF",-5,5);   m->SetPriorGauss("gg_EFF", 0., 1.);           // Efficiency : including PileUp / FF / Calo / Distorted Material
      m->AddSystematic("gg_LUMI",-5,5);  m->SetPriorGauss("gg_LUMI", 0., 1.);          // Luminosity 
      m->AddSystematic("gg_BKG",-5,5);  m->SetPriorGauss("gg_BKG", 0., 1.);          // Irreducible/reducible and purity

     
      // ----------------------------------------------------
      // Make up histograms (BAT wants TH1Ds)
      // ----------------------------------------------------
      TH1D h_sys_gg_eff     = TH1D(*(h_gg_sgn));        // const slope efficiency uncertainty
      TH1D h_sys_gg_Lumi    = TH1D(*(h_gg_sgn));        // lumi for gg channel
      TH1D h_sys_gg_Stat    = TH1D(*(h_gg_sgn));        // lumi for gg channel
      TH1D h_sys_gg_Super   = TH1D(*(h_gg_sgn));        // lumi for gg channel
    
      //For the moment mass independent values
      for (int i=1;i<=h_sys_gg_Lumi.GetNbinsX();i++){
	h_sys_gg_Lumi .SetBinContent(i, 0.037);    //Given by Collaboration
	h_sys_gg_eff  .SetBinContent(i, 0.05);     // 
	h_sys_gg_Stat .SetBinContent(i, 0.02);     // 
	h_sys_gg_Super.SetBinContent(i, 1.);     // 
      }
      if(configuration.SubString("FullSyst") =="FullSyst"){
	m->SetSystematicVariation("gg","ggBkg","gg_BKG",hh_bkg_syst_gg,hh_bkg_syst_gg);
	m->SetSystematicVariation("gg","signal","gg_EFF" ,h_sys_gg_eff,h_sys_gg_eff);
	m->SetSystematicVariation("gg","signal","gg_LUMI",h_sys_gg_Lumi,h_sys_gg_Lumi);
      }
      else  if(configuration.SubString("StatSyst") =="StatSyst"){
	m->SetSystematicVariation("gg","ggBkg","gg_BKG",h_sys_gg_Stat,h_sys_gg_Stat);
	m->SetSystematicVariation("gg","signal","gg_EFF" ,h_sys_gg_eff,h_sys_gg_eff);
	m->SetSystematicVariation("gg","signal","gg_LUMI",h_sys_gg_Lumi,h_sys_gg_Lumi);
      }
      else  if(configuration.SubString("SuperSyst") =="SuperSyst"){
	m->SetSystematicVariation("gg","ggBkg","gg_BKG",h_sys_gg_Super,h_sys_gg_Super);
	m->SetSystematicVariation("gg","signal","gg_EFF" ,h_sys_gg_eff,h_sys_gg_eff);
	m->SetSystematicVariation("gg","signal","gg_LUMI",h_sys_gg_Lumi,h_sys_gg_Lumi);
      }

    }

  // ----------------------------------------------------
  // perform analysis
  // ----------------------------------------------------
  if (!doEnsemble){
    // run MCMC
    m->SetData("gg", h_gg_data); 
    m->MarginalizeAll();

    // find global mode
    m->FindMode( m->GetBestFitParameters() );
	   
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
  m->SetParameterRange(0,Nggbkg,Nggbkg);
  m->SetParameterRange(1,0.0,0.0);//set signal to zero here (flat prior would incl signal in PE)
  pm->SetModel(m);

  // ----------------------------------------------------
  // create output object
  // ----------------------------------------------------
  BCModelOutput* pmout = new BCModelOutput(pm, outputdir+("/prior"+iMassTxt.str()+labelTxt+".root").c_str());

  // switch writing of Markov Chains on
  pmout->WriteMarkovChain(true);

  // set precision
  pm->MCMCSetPrecision(BCEngineMCMC::kMedium);

  // perform marginalization
  pm->MarginalizeAll(); 

  // get tree
  TTree* priortree = (TTree*) pmout->GetFile()->Get("MarkovChainTree_0");

  //undo PE modifications
  m->SetParameterRange(0,Nggbkg,Nggbkg);
  if (doXsecLimit)	
    m->SetParameterRange(1,0.0,maxSigmaSig);
  else
    m->SetParameterRange(1,0.0,maxNsig);

  // ----------------------------------------------------------------
  // Perform ensemble test
  // ----------------------------------------------------------------
  // create new analysis facility 
  BCMTFAnalysisFacility* facility = new BCMTFAnalysisFacility(m); 
  facility->SetFlagMCMC(doMCMC);

  TFile *file = new TFile(outputdir+("/zprime_ensembles_mass"+iMassTxt.str()+labelTxt+".root").c_str(), "RECREATE");
  file->cd(); 


  TTree* tree = facility->BuildEnsembles( priortree, nEnsemble);

  // run ensemble test
  TTree* tree_out = facility->PerformEnsembleTest(tree, nEnsemble);
	
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

