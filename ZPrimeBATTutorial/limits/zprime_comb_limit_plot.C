#include "TString.h"
#include "TSystem.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TFile.h"
#include "TLine.h"
#include "TH1.h"
#include "TH1F.h"
#include "TF1.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "AtlasStyle.h"
#include "TStyle.h"

using std::string;
using std::cerr;

Double_t kFactor(Double_t mass)  // NNLO/2008LO
{
  double truemass = mass;
  double kF = 1.;
  // QCD
  if(truemass > 400) kF *= 1.34815e+00 + 8.58823e-05*mass - 2.11201e-07*mass*mass + 1.11821e-10*mass*mass*mass - 1.76464e-14*mass*mass*mass*mass;
  else kF *= 1.21330e+00 + 8.17233e-04*mass - 1.17503e-06*mass*mass;

  // EW not applied to signal

  return kF;
}

Double_t findIntersection ( TGraph* graphOfData,  TGraph* graphOfTheory, Double_t xMin, Double_t xMax) {
  //gDirectory->GetList()->ls();
  //gDirectory->GetList()->Delete();
  double x = xMin;
  double intersection = -1;
  const int nSteps = 15000;
  double deltaX = (xMax-xMin)/nSteps;
  double deltaOld = 10;

  for (int i=0;i<nSteps;i++) {
    x += deltaX;
    double data = graphOfData->Eval(x);
    double theory = graphOfTheory->Eval(x);
    double delta = fabs(data-theory)/fabs(data);
    if (delta<deltaOld) {
      intersection = x;
      deltaOld = delta;
    }
  }

  return intersection;
}


void zprime_comb_limit_plot(bool ratio, bool xsec, bool doLogmass, int spin,int channel){
/*
  bool ratio=true;
  bool xsec=true;
  bool doLogmass=true;
  int spin=1;
  channel: 0-ee, 1-mm, 2-combo
*/
  double xsec_scalefactor = 1;
  
  double g = 0.2;  // Coupling
  std::string gstr = "02";
  //double g = 0.7;  // Coupling
  //std::string gstr = "07";

  bool sneutrino = (spin==0);
  bool zprime = (spin==1);
  bool graviton = (spin==2);
  bool techni = (spin==3);
  bool torsion = (spin==4);


  float ymin = 1;
  float ymax = 1000;
  if (ratio) {ymin = 0.000001; ymax = 0.0007;}
  if (xsec) {
    ymin=0.0001; ymax=0.4;

    if (graviton) {ymin=0.002; ymax=0.6;}
    if (sneutrino) {ymin=0.002; ymax=0.6;}
  } 

  std::string ch_names[3];
  ch_names[0]="ee";
  ch_names[1]="mm";
  ch_names[2]="combo";
  
  

  TGraph* graph = new TGraph(); 
  TGraph* loggraph = new TGraph(); 
  TGraph* datagraph = new TGraph();
  TGraph* logdatagraph = new TGraph();


  //PDF, alpha_s and kQCD (=scale) error TGraph, input from T. Nunnemann
  double Vmass[39] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.091, 0.1, 0.125, 0.15, 0.175, 
  		      0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 
    		      2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5};
  double VPDF[39] = {0.161273, 0.076909, 0.0583181, 0.0517107, 0.0485283, 0.047392, 0.0457165, 0.0449222, 
  		       0.0440114, 0.0434281, 0.0421782, 0.0410122, 0.0399375, 0.0399249, 0.0397618, 0.041, 
    		       0.0428719, 0.0454423, 0.0488774, 0.0524023, 0.0572276, 0.0628967, 0.0695557, 0.0887525, 
    		       0.115953, 0.151611, 0.19758, 0.25315, 0.313861, 0.377535, 0.442274, 0.507875, 0.576834, 
    		       0.648193, 0.726607, 0.818732, 0.945838, 1.20089, 2.00758};
  TGraph *theory_uncertainty_pdfZ = new TGraph(39, Vmass, VPDF);

  TF1 *theory_uncertainty_pdfG = new TF1("theory_uncertainty_pdfG","pol2",0.07,3.07);
  theory_uncertainty_pdfG->SetParameters(0.0460,0.0054,0.0306);

  TF1 *theory_uncertainty_pdfTS = new TF1("theory_uncertainty_pdfTS","pol2",0.1,3.0);
  if(g==0.2) theory_uncertainty_pdfTS->SetParameters(0.0237639,0.024988,0.014378);        //eta =  0.2
  if(g==0.5) theory_uncertainty_pdfTS->SetParameters(0.024426,0.0341878,-0.00500503);     //eta = 0.5
  if(g==0.7) theory_uncertainty_pdfTS->SetParameters(0.0248354,0.02603,-0.00557187);      //eta = 0.7

 
  int nZoomPointsMass = 73;
  if(graviton){
    nZoomPointsMass = 10;
  }

  std::vector<double> vecLimit[nZoomPointsMass];
  std::vector<double> Mass;   
  std::vector<double> vec68pos, vec95pos;
  std::vector<double> vec68neg, vec95neg;
  //double pelimit = 0;
  //double mass = 0;
  double datalimit = 0;
  double data = 0;
  double datamass = 0;
  

  for (int j=0;j<nZoomPointsMass;j++) {
    char buffer [52];
    //int k = j+1;
    int k = j;
    sprintf (buffer, "%d",k);
    string iMassTxt = buffer;
    string inFile;

    if(zprime || techni){
      if(channel==0){
	inFile = ("PE_out_ee_Zp/zprime_ensembles_mass"+iMassTxt+"_run0_ee_Sys.root").c_str();
      }
      if(channel==1){
	inFile = ("PE_out_mm_Zp/zprime_ensembles_mass"+iMassTxt+"_run0_mm_Sys.root").c_str();
      }
      if(channel==2){
	inFile = ("PE_out_combo_Zp/zprime_ensembles_mass"+iMassTxt+"_run0_comb_Sys.root").c_str();
      }
    }
    if(graviton){
      if(channel==0){
	inFile = ("PE_out_ee_G/zprime_ensembles_mass"+iMassTxt+"_run0_ee_Sys.root").c_str();
      }
      if(channel==1){
	inFile = ("PE_out_mm_G/zprime_ensembles_mass"+iMassTxt+"_run0_mm_Sys.root").c_str();
      }
      if(channel==2){
	inFile = ("PE_out_combo_G/zprime_ensembles_mass"+iMassTxt+"_run0_comb_Sys.root").c_str();
      }
    }
    if(torsion){
      if(channel==0){
	inFile = (Form("PE_out_eeTorsion%s/zprime_ensembles_mass",gstr.c_str())+iMassTxt+"_run0_ee_Sys.root").c_str();
      }
      if(channel==1){
	inFile = (Form("PE_out_mmTorsion%s/zprime_ensembles_mass",gstr.c_str())+iMassTxt+"_run0_mm_Sys.root").c_str();
      }
      if(channel==2){
	inFile = (Form("PE_out_combTorsion%s/zprime_ensembles_mass",gstr.c_str())+iMassTxt+"_run0_comb_Sys.root").c_str();
      }
    }


    
    TFile* myfile = new TFile(inFile.c_str(),"READ");

    TTree* tree = (TTree*)myfile->Get("ensemble_test");

    //cerr<<"looking for branch\n";
//    string branchname="95quantile_marginalized_2";
    string branchname="95quantile_marginalized_3";   // since implementation of the MuonLoose option

    TBranch* b_marg=tree->GetBranch(branchname.c_str());
    TObjArray* objarr_leaves=b_marg->GetListOfLeaves();
    //cerr<<"num leaves: "<<objarr_leaves->GetEntries()<<"\n";
    //cerr<<"first leaf name: "<<objarr_leaves->At(0)->GetName()<<"\n";
    TLeaf* l_marg=b_marg->GetLeaf(objarr_leaves->At(0)->GetName());


    for (int i=0;i<tree->GetEntries();i++) {
      tree->GetEntry(i);
      double pelimit=l_marg->GetValue();
      vecLimit[j].push_back(pelimit); 
    }

    cout << j << "  ";

    //cerr<<"num entries: "<<tree->GetEntries()<<"\n";
    double mass = double(j)*0.04 + 0.13;
    if(graviton){
      double temp_mass[10]={0.3,0.5,0.7,0.8,1.0,1.25,1.5, 1.75, 2.0, 2.25};
      mass=temp_mass[j];
    }
    //cout<<"mass "<<mass<<" limit "<<pelimit<<endl;
    Mass.push_back(mass); 
  }
  
  cout << endl;


  string data_limit_filename;
  if(zprime || techni){
    if(channel==0){
      data_limit_filename="./zprime_Logmassdata_ee.txt";
    }
    if(channel==1){
      data_limit_filename="./zprime_Logmassdata_mm.txt";
    }
    if(channel==2){
      data_limit_filename="./zprime_Logmassdata_combo.txt";
    }
  }
  if(graviton){
    if(channel==0){
      data_limit_filename="./GravitonEE_ObservedLimit.txt";
    }
    if(channel==1){
      data_limit_filename="./GravitonMM_ObservedLimit.txt";
    }
    if(channel==2){
      data_limit_filename="./GravitonCombo_ObservedLimit.txt";
    }
  }
  if(torsion){
    if(channel==0){
      data_limit_filename= Form("./Obs_out_eeTorsion$s/Electron_ObservedLimit.txt",gstr.c_str());
    }
    if(channel==1){
      data_limit_filename= Form("./Obs_out_mmTorsion%s/Muon_ObservedLimit.txt",gstr.c_str());
    }
    if(channel==2){
      data_limit_filename= Form("./Obs_out_combTorsion%s/Combo_ObservedLimit.txt",gstr.c_str());
    }
  }

  ifstream datafile (data_limit_filename.c_str());


  double zxsec = 0.98*989; // NNLO x-sec
 
 
  for (int j=0;j<nZoomPointsMass;j++) {

    datafile >> datamass >> datalimit; 
    
    std::sort(vecLimit[j].begin(), vecLimit[j].end()); 
    unsigned neg95 = unsigned(vecLimit[j].size() * (1 - .9545) / 2);
    unsigned neg68 = unsigned(vecLimit[j].size() * (1 - .6827) / 2);     
    unsigned central = unsigned(vecLimit[j].size() * .5);
    unsigned pos68 = unsigned(vecLimit[j].size() * (1 + .6827) / 2);
    unsigned pos95 = unsigned(vecLimit[j].size() * (1 + .9545) / 2);   

    std::cout << datamass << "  " << vecLimit[j][neg95] << "  " << vecLimit[j][neg68] << "  " << central << "  " << vecLimit[j][pos68] << "  " << vecLimit[j][pos95] << "  " << std::endl;
    
    vec68pos.push_back(vecLimit[j][pos68]);
    vec68neg.push_back(vecLimit[j][neg68]);
    vec95pos.push_back(vecLimit[j][pos95]);
    vec95neg.push_back(vecLimit[j][neg95]); 
    
    if((zprime || techni || torsion) && j==0)  continue; //skip the point at 130GeV
    
    
    int ipoint=graph->GetN();
    datagraph->SetPoint(ipoint, datamass, datalimit);   
    logdatagraph->SetPoint(ipoint, datamass, log(datalimit));   

    graph->SetPoint(ipoint, Mass[j], vecLimit[j][central]);     
    loggraph->SetPoint(ipoint, Mass[j], log(vecLimit[j][central]));     
    //std::cout << Mass[j] << " " << vecLimit[j][central] << std::endl;

  }

  datafile.close();

  int numTheories = 7;
  if(graviton)    numTheories=4;
  if(techni)      numTheories=1;
  if(torsion && g< 0.4)     numTheories=3;
  if(torsion && g>=0.4)     numTheories=5;

  const int numTheoriesConst = numTheories;
  
  ifstream smTheoryfile;
  
  if(zprime)   smTheoryfile.open("./zprime_masstheory.txt");
  if(graviton) smTheoryfile.open("./Graviton_Xsec.txt");
  if(techni)   smTheoryfile.open("./Technicolor_Xsec.txt");
  if(torsion && g< 0.4)  smTheoryfile.open("./TSXsecMassNoTitle_01-03.txt");
  if(torsion && g>=0.4)  smTheoryfile.open("./TSXsecMassNoTitle_04-08.txt");

  
  TGraph* smTheorygraph[numTheoriesConst]; 
  TGraph* smTheorygraphHigh[numTheoriesConst];
  TGraph* smTheorygraphLow[numTheoriesConst];   
  TGraph* logsmTheorygraph[numTheoriesConst]; 
  TGraph* logsmTheorygraphLow[numTheoriesConst];   
  for (int j=0;j<numTheories;j++) {
    smTheorygraph[j] = new TGraph(); 
    smTheorygraphHigh[j] = new TGraph();
    smTheorygraphLow[j] = new TGraph();     
    logsmTheorygraph[j] = new TGraph();
    logsmTheorygraphLow[j] = new TGraph();       
  }
  TGraph* smTheorygraphFill = new TGraph();
  smTheorygraphFill->SetLineWidth(0);
  smTheorygraphFill->SetLineStyle(0);
  smTheorygraphFill->SetLineColor(kWhite);
  smTheorygraphFill->SetFillColor(12);

  double theory[numTheoriesConst];
  int nTheoryPoints = 59;
  if(graviton)
      nTheoryPoints = 93;
  if(techni)
      nTheoryPoints = 11;
  if(torsion)
      nTheoryPoints = 59;
  
  for (int j=0;j<nTheoryPoints;j++) {
    smTheoryfile >> datamass;
    if(!graviton){
      datamass *= 0.001;   // GeV -> TeV
    }
    for (int k=0;k<numTheories;k++) {
      smTheoryfile >> theory[k];
      if(zprime){
	theory[k] *= 1000000000.; // input mb, plot pb
      } 
      else if(graviton) {
	theory[k] *= 0.001; //fb->pb
        theory[k] *= xsec_scalefactor; // Factor for ee -> ee+mm+gg = 4
      }
      else if(techni) {
	theory[k] *= 0.001; //fb->pb
      }
      else if(torsion) {
	theory[k] *= 1.; // pb
	//cout << datamass << " " << theory[k] << endl;
      }
      if (ratio && !xsec) theory[k] /= zxsec;
      if (zprime && k!=7) theory[k] *= kFactor(datamass*1000);
      if (sneutrino) theory[k] *= 2.0;
    }
    for (int k=0;k<numTheories;k++) {
      smTheorygraph[k]->SetPoint(j, datamass, theory[k]);
//      double current_theo_uncert=TMath::Sqrt(theory_uncertainty_pdf->Eval(datamass)*theory_uncertainty_pdf->Eval(datamass)+
//					     theory_uncertainty_kQCD->Eval(datamass)*theory_uncertainty_kQCD->Eval(datamass));
      double current_theo_uncert = 0;
      if(zprime)   current_theo_uncert = theory_uncertainty_pdfZ->Eval(datamass,0,"S");
      if(graviton) current_theo_uncert = theory_uncertainty_pdfG->Eval(datamass);
      if(techni)   current_theo_uncert = theory_uncertainty_pdfZ->Eval(datamass,0,"S");
      if(torsion)  current_theo_uncert = theory_uncertainty_pdfTS->Eval(datamass);
      //if(torsion)  current_theo_uncert = theory_uncertainty_pdfZ->Eval(datamass,0,"S");
      
      smTheorygraphHigh[k]->SetPoint(j, datamass, (1.+current_theo_uncert)*theory[k]); 
      smTheorygraphLow[k]->SetPoint(j, datamass, (1.-current_theo_uncert)*theory[k]);     
      if((k==0 && zprime) || (k==3 && graviton) || (k==0 && techni) || (k==0 && torsion)){ 
        smTheorygraphFill->SetPoint(j, datamass, (1.+current_theo_uncert)*theory[k]); 
        smTheorygraphFill->SetPoint(2*nTheoryPoints-j-1, datamass, (1.-current_theo_uncert)*theory[k]);
      }         
      logsmTheorygraph[k]->SetPoint(j, datamass, log(theory[k])); 
      logsmTheorygraphLow[k]->SetPoint(j, datamass, log((1.-current_theo_uncert)*theory[k]));         
    }
  }
  smTheoryfile.close();
 
  if (ratio && xsec){
    float lowScanRange = 0.2;
    float highScanRange = 3.0;
    for (int k=0;k<numTheories;k++) {
      double xIntersectionExp = findIntersection (loggraph, logsmTheorygraph[k], lowScanRange, highScanRange);
      double xIntersection = findIntersection (logdatagraph, logsmTheorygraph[k], lowScanRange, highScanRange);     
      std::cout << "Expected Mass Limit (TeV) = " << xIntersectionExp << " Observed Mass Limit (TeV) = " << xIntersection << std::endl;   
      std::cout << "Expected Cross Section (pb) = " <<  smTheorygraph[k]->Eval(xIntersectionExp) << " Observed cross section (pb) = " << smTheorygraph[k]->Eval(xIntersection) << std::endl;                 
    }
  }

 //-----Make the pretty plot-----//

   gDirectory->GetList()->Delete();

   gROOT->ProcessLine(".x ./AtlasStyle.C");
   gROOT->SetStyle("ATLAS");
   gROOT->ForceStyle();

   TCanvas *c1 = new TCanvas("c1","C1",800,600);
   TGraph* hist95 = new TGraph();
   
   
  TH1F* AxisRanges = new TH1F("AxisRanges","",1000,0.1,3.1); // draw this first to define axes ranges
  AxisRanges->SetLineWidth(0);      
  AxisRanges->SetLineStyle(0);      
  AxisRanges->SetLineColor(kWhite); 
  AxisRanges->SetFillColor(kBlack);
  AxisRanges->GetYaxis()->SetTitleSize(0.045);     
  AxisRanges->GetYaxis()->SetLabelSize(0.045);     
  AxisRanges->GetXaxis()->SetTitleSize(0.05);      
  AxisRanges->GetXaxis()->SetLabelSize(0.045);     
  AxisRanges->GetXaxis()->SetTitleOffset(0.95);    
  AxisRanges->GetYaxis()->SetTitleOffset(1.1);  

  //AxisRanges->SetMinimum(ymin);                    
  //AxisRanges->SetMaximum(ymax);
  AxisRanges->SetMinimum(0.001);                    
  AxisRanges->SetMaximum(4.);
  AxisRanges->SetTitle("");
  AxisRanges->Draw("");
  
  AxisRanges->GetXaxis()->SetTitle("m [TeV]");
  AxisRanges->GetYaxis()->SetTitle("#sigma B [pb]");
  AxisRanges->Draw("");
  
  
  hist95->SetLineWidth(0);
  hist95->SetLineStyle(0);
  hist95->SetLineColor(kWhite);
  hist95->SetFillColor(kYellow);
  
   for (unsigned ibin = 0; ibin < vec95pos.size(); ++ibin)
     {
       if((zprime || techni || torsion) && ibin==0)  continue; //skip the point at 130GeV
       
       int ctr=hist95->GetN();
       hist95->SetPoint(ctr, Mass[ibin], vec95pos[ibin]);
     }
   unsigned counter95 = hist95->GetN();
   for (int bin = vec95neg.size() - 1; bin >= 0; --bin)
     {
       if((zprime || techni || torsion) && bin==0)  continue; //skip the point at 130GeV
       
       hist95->SetPoint(counter95++, Mass[bin], vec95neg[bin]);
     }
   
   hist95->Draw("sameF"); 
   // hist95->SetMaximum(ymax);
   // hist95->SetMinimum(ymin); 
   // hist95->SetTitle("");
/*
   if (graviton) {
     if(ratio && !xsec)  hist95->GetYaxis()->SetTitle("95% C.L. Limits on #sigma B(G*#rightarrow l^{+} l^{-})/#sigma B(Z#rightarrow l^{+} l^{-})");
     if(ratio && xsec)   hist95->GetYaxis()->SetTitle("95% C.L. Limits on #sigma B(G*#rightarrow l^{+} l^{-}) [pb]");
   }
   else if (zprime) {
     if(!ratio && !xsec) hist95->GetYaxis()->SetTitle("95% C.L. Limits on {N}_{Z'}");
     if(ratio && !xsec)  hist95->GetYaxis()->SetTitle("95% C.L. Limits on #sigma B(Z'#rightarrow l^{+} l^{-})/#sigma B(Z#rightarrow l^{+} l^{-})");
     if(ratio && !xsec)  hist95->GetYaxis()->SetTitle("R_{#sigma}");
     if(ratio && xsec)   hist95->GetYaxis()->SetTitle("#sigma B [pb]");
   }
   else if (sneutrino) {
     if(ratio && !xsec)  hist95->GetYaxis()->SetTitle("95% C.L. Limits on #sigma B(#tilde{#nu}#rightarrow l^{+} l^{-})/#sigma B(Z#rightarrow l^{+} l^{-})");
     if(ratio && xsec)   hist95->GetYaxis()->SetTitle("95% C.L. Limits on #sigma B(#tilde{#nu}#rightarrow l^{+} l^{-}) [pb]");
   }
*/
   TGraph* hist68 = new TGraph();
   hist68->SetLineWidth(0);
   hist68->SetLineStyle(0);
   hist68->SetLineColor(kWhite);
   hist68->SetFillColor(kGreen);
   for (unsigned gbin = 0; gbin < vec68pos.size(); ++gbin)
   {
     if((zprime || techni || torsion) && gbin==0)  continue; //skip the point at 130GeV
     
     int ctr=hist68->GetN();
     hist68->SetPoint(ctr, Mass[gbin], vec68pos[gbin]);
   }
   unsigned counter68 = hist68->GetN();
   for (int fbin = vec68neg.size() - 1; fbin >= 0; --fbin)
   {
     if((zprime || techni || torsion) && fbin==0)  continue; //skip the point at 130GeV
     
     hist68->SetPoint(counter68++, Mass[fbin], vec68neg[fbin]);
   }
  

   int theoryColor[8]  = {12,6,11,4,65,40,8,95};
   if(graviton){
     theoryColor[3]=12;
     theoryColor[0]=4;
   }
   int theoryMarker[8] = {24,25,22,26,23,27,21,28};


   for (int k=0;k<numTheories;k++) {
     smTheorygraph[k]->SetMarkerStyle(theoryMarker[k]);
     smTheorygraph[k]->SetMarkerSize(.7);
     //smTheorygraph[k]->SetMarkerSize(0);
     smTheorygraph[k]->SetMarkerColor(theoryColor[k]);
     smTheorygraph[k]->SetLineColor(theoryColor[k]);
     smTheorygraph[k]->SetLineWidth(2);
     if((k==0 && zprime) || (k==3 && graviton)) smTheorygraph[k]->SetLineWidth(6);  //TURN OFF THEORY BAND
     smTheorygraph[k]->SetFillColor(theoryColor[k]);
     //if(k==0 && !graviton) smTheorygraph[k]->SetLineStyle(1);    // 4 for dash-dotted     
     smTheorygraphHigh[k]->SetMarkerSize(0);
     smTheorygraphHigh[k]->SetMarkerColor(theoryColor[k]);
     smTheorygraphHigh[k]->SetLineColor(theoryColor[k]);
     smTheorygraphHigh[k]->SetFillColor(theoryColor[k]); 
     smTheorygraphHigh[k]->SetLineWidth(2);
     smTheorygraphLow[k]->SetMarkerSize(0);
     smTheorygraphLow[k]->SetMarkerColor(theoryColor[k]);
     smTheorygraphLow[k]->SetLineColor(theoryColor[k]);
     smTheorygraphLow[k]->SetFillColor(theoryColor[k]); 
     smTheorygraphLow[k]->SetLineWidth(2);         
     //smTheorygrapherrors[k]->SetFillStyle(3002);   
     
     //smTheorygrapherrors[k]->SetFillColor(kYellow);
     //if (xsec && ratio) smTheorygrapherrors[k]->Draw("3 same");     
     //if (ratio && (k==3 || k==4 || k==7)) smTheorygraph[k]->Draw("samel");
     //if (xsec && ratio && (k==0)) smTheorygraphHigh[k]->Draw("samel");
     //if (xsec && ratio && (k==0)) smTheorygraphLow[k]->Draw("samel");       
   }   

   graph->SetLineStyle(2);
   graph->SetLineWidth(2);
   graph->SetMarkerStyle(28);
   graph->SetMarkerSize(.7);
   graph->Draw("sameL");
   datagraph->SetMarkerStyle(20);
   datagraph->SetMarkerSize(.7);
   datagraph->SetMarkerColor(2);
   datagraph->SetLineColor(2);
   datagraph->SetLineWidth(2);
   datagraph->SetFillColor(2);
   datagraph->Draw("samel");   // UNCOMMENT this line to include datagraph !!!

   float legendXmin, legendXmax, legendYmin, legendYmax;
   legendXmin = 0.7; legendXmax = 0.88; legendYmin = 0.5; legendYmax = 0.9;
   if(ratio && !xsec) legendYmin = 0.55;
   TLegend* l = new TLegend(legendXmin, legendYmin, legendXmax, legendYmax);
   l->SetFillColor(10);
   l->SetLineColor(0);   
   l->SetTextFont(42);     
   l->AddEntry(graph,"Expected limit","l");
   l->AddEntry(hist68,"Expected #pm 1#sigma","f");
   l->AddEntry(hist95,"Expected #pm 2#sigma","f");
   l->AddEntry(datagraph,"Observed limit","l");  // Change to "Observed limit" or "Pseudo-Data"
   l->SetLineColor(0);
   l->SetFillColor(10);
   l->SetBorderSize(0); 
   l->SetTextSize(0.04);
       

   for(int k=numTheories-1;k>=0;k--){
     if(ratio){
       if (graviton) {
	 if (k==0) l->AddEntry(smTheorygraph[k],"k/#bar{M}_{Pl} = 0.01","l");
	 if (k==1) l->AddEntry(smTheorygraph[k],"k/#bar{M}_{Pl} = 0.03","l");
	 if (k==2) l->AddEntry(smTheorygraph[k],"k/#bar{M}_{Pl} = 0.05","l");
	 if (k==3) l->AddEntry(smTheorygraph[k],"k/#bar{M}_{Pl} = 0.1","l");
       }
     }
   }
   for (int k=0;k<numTheories;k++) {
     if (ratio) {
       if (zprime) {
         if (k==0 && numTheories>7) l->AddEntry(smTheorygraph[7],"Z*","l");
         if (k==0) l->AddEntry(smTheorygraph[k],"Z'_{SSM}","l");
         if (k==6){
           l->AddEntry(smTheorygraph[4],"Z'_{#chi}","l");
           //l->AddEntry(smTheorygraph[1],"Z'_{S}","l");
           //l->AddEntry(smTheorygraph[6],"Z'_{I}","l");
           //l->AddEntry(smTheorygraph[5],"Z'_{#eta}","l");
           //l->AddEntry(smTheorygraph[2],"Z'_{N}","l");
           l->AddEntry(smTheorygraph[3],"Z'_{#psi}","l");
         }
       }
       else if (sneutrino) {
         if (k==0) l->AddEntry(smTheorygraph[k],"#lambda^{2}BR = 0.0001","l");
         if (k==1) l->AddEntry(smTheorygraph[k],"#lambda^{2}BR = 0.0002","l");
         if (k==2) l->AddEntry(smTheorygraph[k],"#lambda^{2}BR = 0.0005","l");
         if (k==3) l->AddEntry(smTheorygraph[k],"#lambda^{2}BR = 0.001","l");
         if (k==4) l->AddEntry(smTheorygraph[k],"#lambda^{2}BR = 0.002","l");
         if (k==5) l->AddEntry(smTheorygraph[k],"#lambda^{2}BR = 0.005","l");
         if (k==6) l->AddEntry(smTheorygraph[k],"#lambda^{2}BR = 0.01","l");
       }
       else if (techni) {
         if (k==0) l->AddEntry(smTheorygraph[k],"#rho_{T}/#omega_{T}, no a_{T}","l");
       }
       else if (torsion && g< 0.4) {
         if (k==0) l->AddEntry(smTheorygraph[k],"#eta = 0.1","l");
         if (k==1) l->AddEntry(smTheorygraph[k],"#eta = 0.2","l");
         if (k==2) l->AddEntry(smTheorygraph[k],"#eta = 0.3","l");
       }
       else if (torsion && g>=0.4) {
         if (k==0) l->AddEntry(smTheorygraph[k],"#eta = 0.4","l");
         if (k==1) l->AddEntry(smTheorygraph[k],"#eta = 0.5","l");
         if (k==2) l->AddEntry(smTheorygraph[k],"#eta = 0.6","l");
         if (k==3) l->AddEntry(smTheorygraph[k],"#eta = 0.7","l");
         if (k==4) l->AddEntry(smTheorygraph[k],"#eta = 0.8","l");
       }
     }
   }
   l->SetFillColor(10);
   l->SetLineColor(0);
   l->Draw();   

   if(xsec && ratio) {
     char writetext1[100];
     char writetext2[100];
     char writetext3[100];
     char writetext4[100];
     char writetext5[100];
     char writetext6[100];
     TLatex *t = new TLatex();
     t->SetNDC(1);
     t->SetTextAlign(13);
     t->SetTextColor(kBlack);
     
     //sprintf(writetext1,"#font[72]{ATLAS}");               // ONLY AFTER PAPER SECOND ATLAS READING!
     //sprintf(writetext1,"#font[72]{ATLAS} Preliminary");   // ONLY AFTER PAPER FIRST ATLAS READING!
     //sprintf(writetext1,"#font[72]{ATLAS} For Approval");
     sprintf(writetext1,"");
      
     if(zprime){ 
       if(channel==0)
         sprintf(writetext2,"Z' #rightarrow ee");
       if(channel==1)
        sprintf(writetext2,"Z' #rightarrow #mu#mu");
       if(channel==2)
         sprintf(writetext2,"Z' #rightarrow ll");
     }
     else if(graviton){
       if(channel==0)
         sprintf(writetext2,"G* #rightarrow ee");
       if(channel==1)
         sprintf(writetext2,"G* #rightarrow #mu#mu");
       if(channel==2)
         sprintf(writetext2,"G* #rightarrow ll");
     }
     else if(techni){
       if(channel==0)
         sprintf(writetext2,"#rho_{T} #rightarrow ee");
       if(channel==1)
         sprintf(writetext2,"#rho_{T} #rightarrow #mu#mu");
       if(channel==2)
         sprintf(writetext2,"#rho_{T} #rightarrow ll");
     }
     else if(torsion){
       if(channel==0)
         sprintf(writetext2,"TS #rightarrow ee");
       if(channel==1)
         sprintf(writetext2,"TS #rightarrow #mu#mu");
       if(channel==2)
         sprintf(writetext2,"TS #rightarrow ll");
     }
     else sprintf(writetext2,"");

     sprintf(writetext3,"#sqrt{s} = 7 TeV");
     sprintf(writetext4,"ee: #int L dt = 4.9 fb^{-1}");
     sprintf(writetext5,"#mu#mu: #int L dt = 5.0 fb^{-1}");

     if(torsion) sprintf(writetext6,"Templates with #eta = %g",g);
     else sprintf(writetext6,"");

     double xtext = 0.45; 
     double ytext = 0.85; 
     t->SetTextSize(0.045);
     t->DrawLatex(xtext-0.07,ytext+0.06,writetext1);
     //t->DrawLatex(xtext+0.05,ytext+0.06,writetext1);
     t->DrawLatex(xtext+0.05,ytext-0.06,writetext2);
     t->DrawLatex(xtext+0.05,ytext     ,writetext3);
     t->SetTextSize(0.035);
     if(channel==0)
       t->DrawLatex(0.18,0.31,writetext4);
     if(channel==1)
       t->DrawLatex(0.18,0.31,writetext5);
     if(channel==2){
       t->DrawLatex(0.18,0.23,writetext5);
       t->DrawLatex(0.18,0.31,writetext4);
     }
     t->SetTextSize(0.030);
     t->DrawLatex(xtext+0.02,ytext-0.12,writetext6);    
   }
   

//   c1->SetGridx();
//   c1->SetGridy();


   c1->SetLogy(1); 
   hist95->Draw("sameF");    
   hist68->Draw("sameF");   
   graph->Draw("sameL");
   datagraph->Draw("samel");
   if(zprime){
     //smTheorygraph[0]->Draw("samel");
     smTheorygraphFill->Draw("sameF");  // comment to turn off theory band
     smTheorygraph[3]->Draw("samel"); 
     smTheorygraph[4]->Draw("samel");    
     //hist95->GetXaxis()->SetRangeUser(0.10,1.60);
   } 
   else if(graviton) {
     smTheorygraph[0]->Draw("samel");
     smTheorygraph[1]->Draw("samel");
     smTheorygraph[2]->Draw("samel");
     //smTheorygraph[3]->Draw("samel");
     smTheorygraphFill->Draw("sameF");  // comment to turn off theory band
   }
   else if(techni) {
     smTheorygraphFill->Draw("sameF"); 
   }   
   else if(torsion) {
     smTheorygraphFill->Draw("sameF"); // comment to turn off theory band
     //smTheorygraph[0]->Draw("samel");
     smTheorygraph[1]->Draw("samel");
     smTheorygraph[2]->Draw("samel");
     if(g>=0.4) {
       smTheorygraph[3]->Draw("samel");
       smTheorygraph[4]->Draw("samel");
     }
   }   
   
   c1->RedrawAxis();
   c1->Update();   
   if (graviton) { 
     if (doLogmass) {
       if(!ratio && !xsec) {
         c1->Print(Form("Logmasslimit_ngraviton_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Logmasslimit_ngraviton_%s.png",ch_names[channel].c_str()));
       }   
       if (ratio && !xsec) {
         c1->Print(Form("Logmasslimit_gravitonratio_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Logmasslimit_gravitonratio_%s.png",ch_names[channel].c_str()));
       }   
       if (xsec && ratio) {
	 c1->Print(Form("Logmasslimit_gravitonxsec_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Logmasslimit_gravitonxsec_%s.png",ch_names[channel].c_str())); 
       }       
     }
     else {
       if(!ratio && !xsec) {
         c1->Print(Form("Limit_ngraviton_vsMass_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Limit_ngraviton_vsMass_%s.png",ch_names[channel].c_str()));
       }   
       if (ratio && !xsec) {
         c1->Print(Form("Limit_gravitonratio_vsMass_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Limit_gravitonratio_vsMass_%s.png",ch_names[channel].c_str()));
       }   
       if (xsec && ratio) {
	 c1->Print(Form("Limit_gravitonxsec_vsMass_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Limit_gravitonxsec_vsMass_%s.png",ch_names[channel].c_str())); 
       }
     }
   }
   else if (zprime) {
     if (doLogmass) {
       if(!ratio && !xsec) {
         c1->Print(Form("Logmasslimit_nzprime_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Logmasslimit_nzprime_%s.png",ch_names[channel].c_str()));
       }   
       if (ratio && !xsec) {
         c1->Print(Form("Logmasslimit_zprimeratio_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Logmasslimit_zprimeratio_%s.png",ch_names[channel].c_str()));
       }   
       if (xsec && ratio) {
	 c1->Print(Form("Logmasslimit_zprimexsec_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Logmasslimit_zprimexsec_%s.png",ch_names[channel].c_str())); 
       }       
     }
     else {
       if(!ratio && !xsec) {
         c1->Print(Form("Limit_nzprime_vsMass_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Limit_nzprime_vsMass_%s.png",ch_names[channel].c_str()));
       }   
       if (ratio && !xsec) {
         c1->Print(Form("Limit_zprimeratio_vsMass_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Limit_zprimeratio_vsMass_%s.png",ch_names[channel].c_str()));
       }   
       if (xsec && ratio) {
	 c1->Print(Form("Limit_zprimexsec_vsMass_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Limit_zprimexsec_vsMass_%s.png",ch_names[channel].c_str())); 
       }
     }
   }
   else if (sneutrino) {
     if (doLogmass) {
       if(!ratio && !xsec) {
         c1->Print(Form("Logmasslimit_nsneutrino_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Logmasslimit_nsneutrino_%s.png",ch_names[channel].c_str()));
       }   
       if (ratio && !xsec) {
         c1->Print(Form("Logmasslimit_sneutrinoratio_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Logmasslimit_sneutrinoratio_%s.png",ch_names[channel].c_str()));
       }   
       if (xsec && ratio) {
	 c1->Print(Form("Logmasslimit_sneutrinoxsec_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Logmasslimit_sneutrinoxsec_%s.png",ch_names[channel].c_str())); 
       }       
     }
     else {
       if(!ratio && !xsec) {
         c1->Print(Form("Limit_nsneutrino_vsMass_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Limit_nsneutrino_vsMass_%s.png",ch_names[channel].c_str()));
       }   
       if (ratio && !xsec) {
         c1->Print(Form("Limit_sneutrinoratio_vsMass_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Limit_sneutrinoratio_vsMass_%s.png",ch_names[channel].c_str()));
       }   
       if (xsec && ratio) {
	 c1->Print(Form("Limit_sneutrinoxsec_vsMass_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Limit_sneutrinoxsec_vsMass_%s.png",ch_names[channel].c_str())); 
       }
     }
   }
   else if (techni) {
     if (doLogmass) {
       if(!ratio && !xsec) {
         c1->Print(Form("Logmasslimit_ntechnicolor_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Logmasslimit_ntechnicolor_%s.png",ch_names[channel].c_str()));
       }   
       if (ratio && !xsec) {
         c1->Print(Form("Logmasslimit_technicolorratio_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Logmasslimit_technicolorratio_%s.png",ch_names[channel].c_str()));
       }   
       if (xsec && ratio) {
	 c1->Print(Form("Logmasslimit_technicolorxsec_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Logmasslimit_technicolorxsec_%s.png",ch_names[channel].c_str())); 
       }       
     }
     else {
       if(!ratio && !xsec) {
         c1->Print(Form("Limit_ntechnicolor_vsMass_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Limit_ntechnicolor_vsMass_%s.png",ch_names[channel].c_str()));
       }   
       if (ratio && !xsec) {
         c1->Print(Form("Limit_technicolorratio_vsMass_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Limit_technicolorratio_vsMass_%s.png",ch_names[channel].c_str()));
       }   
       if (xsec && ratio) {
	 c1->Print(Form("Limit_technicolorxsec_vsMass_%s.eps",ch_names[channel].c_str())); 
	 c1->Print(Form("Limit_technicolorxsec_vsMass_%s.png",ch_names[channel].c_str())); 
       }
     }
   }
   else if (torsion) {
     if (doLogmass) {
       if(!ratio && !xsec) {
         c1->Print(Form("Logmasslimit_ntorsion_%s_%s.eps",gstr.c_str(),ch_names[channel].c_str())); 
	 c1->Print(Form("Logmasslimit_ntorsion_%s_%s.png",gstr.c_str(),ch_names[channel].c_str()));
       }   
       if (ratio && !xsec) {
         c1->Print(Form("Logmasslimit_torsionratio_%s_%s.eps",gstr.c_str(),ch_names[channel].c_str())); 
	 c1->Print(Form("Logmasslimit_torsionratio_%s_%s.png",gstr.c_str(),ch_names[channel].c_str()));
       }   
       if (xsec && ratio) {
	 c1->Print(Form("Logmasslimit_torsionxsec_%s_%s.eps",gstr.c_str(),ch_names[channel].c_str())); 
	 c1->Print(Form("Logmasslimit_torsionxsec_%s_%s.png",gstr.c_str(),ch_names[channel].c_str())); 
       }       
     }
     else {
       if(!ratio && !xsec) {
         c1->Print(Form("Limit_ntorsion_vsMass_%s_%s.eps",gstr.c_str(),ch_names[channel].c_str())); 
	 c1->Print(Form("Limit_ntorsion_vsMass_%s_%s.png",gstr.c_str(),ch_names[channel].c_str()));
       }   
       if (ratio && !xsec) {
         c1->Print(Form("Limit_torsionratio_vsMass_%s_%s.eps",gstr.c_str(),ch_names[channel].c_str())); 
	 c1->Print(Form("Limit_torsionratio_vsMass_%s_%s.png",gstr.c_str(),ch_names[channel].c_str()));
       }   
       if (xsec && ratio) {
	 c1->Print(Form("Limit_torsionxsec_vsMass_%s_%s.eps",gstr.c_str(),ch_names[channel].c_str())); 
	 c1->Print(Form("Limit_torsionxsec_vsMass_%s_%s.png",gstr.c_str(),ch_names[channel].c_str())); 
       }
     }
   }
   //gROOT->ProcessLine(".q");
}
