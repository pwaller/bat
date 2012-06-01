
int plot2Dscan(double TSobs, double pvalue, bool doEE, bool doMM)
{

  gROOT->ProcessLine(".x SetStyle.C");

  TCanvas *c1 = new TCanvas("c1", "c1",600,500);

    const Int_t NRGBs = 5;
    const Int_t NCont = 99;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);


  TFile * file;
  if (doEE)
    //file = new TFile("usefulPlots_pvalue_eeSys_inclTheo.root", "READ");
    file = new TFile("usefulPlots_pvalue_ee.root", "READ");
  else 
    //file = new TFile("usefulPlots_pvalue_mmSys_inclTheo.root", "READ");
    file = new TFile("usefulPlots_pvalue_mm.root", "READ");
  

  TH2D* tmp = (TH2D*)(file->Get("hist_model_signal_Mass"));
  tmp->SetContour(100);


  double norm  = fabs(TSobs) / tmp->GetBinContent(tmp->GetMaximumBin());
  tmp->Scale(norm);


  tmp->Draw("colz");
  tmp->GetXaxis()->SetTitle("#sigma_{Z'} [pb]");
  tmp->SetTitle("");
  tmp->GetYaxis()->SetTitle("M_{Z'} [Tev]");
  //tmp->GetYaxis()->SetTitle("M_{Z'} [Tev/c^{2}]");
  tmp->GetYaxis()->SetTitleOffset(1.2);
  tmp->GetYaxis()->SetLabelFont(42);
  tmp->GetXaxis()->SetLabelFont(42);
  tmp->GetYaxis()->SetTitleFont(42);
  tmp->GetXaxis()->SetTitleFont(42);

  TLatex *t = new TLatex();
  t->SetNDC(1);
  t->SetTextAlign(13);
  t->SetTextColor(kBlack);

  if (doEE)
     t->DrawLatex(0.5,0.52,"Z' #rightarrow ee");
  else
     t->DrawLatex(0.5,0.52,"Z' #rightarrow #mu#mu");


  char writetext1[100];
  char writetext2[100];
  char writetext3[100];
  char writetext4[100];
  //	 sprintf(writetext1,"#font[72]{ATLAS} #font[42]{For Approval}");
       sprintf(writetext1," ");
  sprintf(writetext2,"p = %g", pvalue);
  sprintf(writetext3,"#sqrt{s} = 7 TeV");
/*
  if (doEE)
    sprintf(writetext4,"#int L dt = 1.08 fb^{-1}");
  else
    sprintf(writetext4,"#int L dt = 2.49 fb^{-1}");
*/

  double xtext = 0.4; 
  double ytext = 0.73; 
  t->SetTextSize(0.05);
  t->DrawLatex(xtext     ,ytext+0.10,writetext1);
  t->SetTextSize(0.045);
  t->DrawLatex(xtext+0.20,ytext-0.10,writetext2);
  t->DrawLatex(xtext	 ,ytext-0.10,writetext3);
  t->DrawLatex(xtext	 ,ytext     ,writetext4);
  

   tmp->GetYaxis()->SetLimits(0.13,3.05);
  
   tmp->GetYaxis()->SetRangeUser(0.13,2.1);
   tmp->GetXaxis()->SetRangeUser(0.0,0.16);
   

   double modex=0.01;
   double modey=0.85;

   TH1D* fHistogram = (TH1D*)(tmp->Clone());
   int maximumbin = fHistogram->GetMaximumBin();
  
   int binx = maximumbin % (fHistogram->GetNbinsX() + 2);
   int biny = maximumbin / (fHistogram->GetNbinsX() + 2);
  
   modex = fHistogram->GetXaxis()->GetBinCenter(binx);
   modey = fHistogram->GetYaxis()->GetBinCenter(biny);
  
   fHistogram->SetLineWidth(1);
   fHistogram->Draw("CONT3same");
  

   // set contours
   int nz = 100;
  
   double zmax = fHistogram->GetMaximum();
   double dz   = zmax / double(nz);
  
   double nx = fHistogram->GetNbinsX();
   double ny = fHistogram->GetNbinsY();

   TH1D* fIntegratedHistogram = new TH1D("", "", nz, 0.0, zmax);
   fIntegratedHistogram->SetXTitle("z");
   fIntegratedHistogram->SetYTitle("Integrated probability");
   fIntegratedHistogram->SetStats(kFALSE);
  
   // loop over histogram
   for (int ix = 1; ix <= nx; ix++) {
      for (int iy = 1; iy <= ny; iy++) {
	 int binmin = int(fHistogram->GetBinContent(ix, iy) / dz);
	 for (int i = binmin; i <= nz; i++) {
	    fIntegratedHistogram->SetBinContent(i,
						fIntegratedHistogram->GetBinContent(i) +
						fHistogram->GetBinContent(ix, iy));
	 }
      }
   }


  double levels[2];

  double quantiles[1];
  double probsum[1];
  probsum[0] = 0.32;

  fIntegratedHistogram->GetQuantiles( 1, quantiles, probsum);
  	  
  levels[0] = 0.;
  levels[1] = quantiles[0];

  fHistogram->SetContour(2, levels);


  TMarker * marker0 = new TMarker(modex, modey, 8);
  marker0->SetMarkerColor(0);
  marker0->SetMarkerSize(.7);
  marker0->Draw();
  TMarker * marker1 = new TMarker(modex, modey, 4);
  marker1->SetMarkerColor(1);
  marker1->SetMarkerSize(.7);
  marker1->Draw();

  TLegend* legend = new TLegend(0.47, 0.33, 0.87, 0.46);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetFillColor(kWhite);
  legend->AddEntry(marker1, "Signal Scan, Best Fit", "P");
  legend->AddEntry(fHistogram, "68% Contour", "L");
  legend->Draw();


  c1->Update();
  if (doEE) {
    c1->Print("scan2d_ee.eps");
    c1->Print("scan2d_ee.gif");
  }
  else {
    c1->Print("scan2d_mm.eps");
    c1->Print("scan2d_mm.gif");
  }

  return 0;

}
