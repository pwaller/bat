
double plot_pvalue(double TSobs, bool do_ee, bool do_mm)
{
  gROOT->ProcessLine(".x SetStyle.C");
  TCanvas *c1 = new TCanvas("c1", "c1",600,550);
  c1->SetLogy();
  gROOT->ForceStyle();


  TFile * file;
  if (do_ee)
    file = new TFile("zprime_ensembles_pvalue_0_run4_ee.root", "READ");
  else 
    file = new TFile("zprime_ensembles_pvalue_0_run2_mm.root", "READ");

  
  TTree* tree = ((TTree*) file->Get("ensemble_test"));
  
  tree->Draw("-2*log(Lsb/Lb)>>tmp(40,-19,1)");
  tree->Draw("-2*log(Lsb/Lb)>>tmp2(10000,-19,1)");
  TH1F *tmp  = (TH1F*)gDirectory->Get("tmp");
  TH1F *tmp2 = (TH1F*)gDirectory->Get("tmp2");

  double num = tmp2->Integral(0,tmp2->FindBin(TSobs));
  double den =  tmp2->Integral();
  double pvalue = num / den;
  std::cout << "Pvalue = " << pvalue << " = " << num << "/" << den << std::endl;
  
  tmp->SetFillColor(kYellow);  
  tmp->SetMaximum(tmp->GetMaximum()*8);
  tmp->Draw();
  tmp->GetXaxis()->SetTitle("LLR");
  tmp->SetTitle("");
  tmp->GetYaxis()->SetTitle("Pseudo-Experiments");
  tmp->GetYaxis()->SetLabelFont(42);
  tmp->GetXaxis()->SetLabelFont(42);
  tmp->GetYaxis()->SetTitleFont(42);
  tmp->GetXaxis()->SetTitleFont(42);

  TArrow* tl = new TArrow(TSobs, tmp->GetMaximum()*0.005, TSobs, 0);
  tl->SetLineWidth(2);
  tl->SetLineColor(4);
  tl->Draw();

  TLatex *t = new TLatex();
  t->SetNDC(1);
  t->SetTextAlign(13);
  t->SetTextColor(kBlack);
  char writetext2[100];
  sprintf(writetext2,"p = %g", pvalue);
  t->SetTextSize(0.045);
  t->DrawLatex(0.20,0.65,writetext2);

  TLegend* l = new TLegend(0.14, 0.72, 0.87, 0.86);
  l->SetTextFont(42);
  l->AddEntry(tmp,"Pseudo-Experiments","lfp");
  if (do_ee)
    l->AddEntry(tl,"Observed value in e+e- Data","l");
  else
    l->AddEntry(tl,"Observed value in #mu+#mu- Data","l");
  l->SetLineColor(0);
  l->SetFillColor(10);
  l->SetBorderSize(0);     
  l->Draw("same");
 

  c1->Update();
  if (do_ee) {
    c1->Print("pvalue_ee.eps");
    c1->Print("pvalue_ee.gif");
  }
  else {
    c1->Print("pvalue_mm.eps");
    c1->Print("pvalue_mm.gif");
  }

  return pvalue;
  
}
