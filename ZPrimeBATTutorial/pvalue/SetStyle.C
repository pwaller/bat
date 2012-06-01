{

  //****************************STYLE***********************************

  const char* modified = "Borrowed and adapted from paus et al";

  TStyle *RootStyle = new TStyle("Root-Style","The Perfect Style for Plots ;-)");
  
#ifdef __CINT__
  TStyle *GloStyle;
  GloStyle = gStyle;                          // save the global style reference
  gStyle = RootStyle;
#endif

// otherwise you need to call TROOT::SetStyle("Root-Style")

  // Paper size

  RootStyle->SetPaperSize(TStyle::kUSLetter);
  //RootStyle->SetHatchesSpacing(0.05);

  // Canvas

  RootStyle->SetCanvasColor     (0);
  RootStyle->SetCanvasBorderSize(10);
  RootStyle->SetCanvasBorderMode(0);
  RootStyle->SetCanvasDefH      (600);
  RootStyle->SetCanvasDefW      (600);
  RootStyle->SetCanvasDefX      (10);
  RootStyle->SetCanvasDefY      (10);

  // Pads

  RootStyle->SetPadColor       (0);
  RootStyle->SetPadBorderSize  (2);
  RootStyle->SetPadBorderMode  (0);
  RootStyle->SetPadBottomMargin(0.10);
  RootStyle->SetPadTopMargin   (0.10);
  RootStyle->SetPadLeftMargin  (0.12);
  RootStyle->SetPadRightMargin (0.12);
  RootStyle->SetPadGridX       (0);
  RootStyle->SetPadGridY       (0);
  RootStyle->SetPadTickX       (1);
  RootStyle->SetPadTickY       (1);

  // Frames

  RootStyle->SetFrameFillStyle ( 0);
  RootStyle->SetFrameFillColor ( 0);
  RootStyle->SetFrameLineColor ( 1);
  RootStyle->SetFrameLineStyle ( 0);
  RootStyle->SetFrameLineWidth ( 2);
  RootStyle->SetFrameBorderSize(10);
  RootStyle->SetFrameBorderMode( 0);


  // Histograms

  // RootStyle->SetHistFillColor(2);
//   RootStyle->SetHistFillStyle(1);
  //RootStyle->SetHistLineColor(1);
  //RootStyle->SetHistLineStyle(0);
  RootStyle->SetHistLineWidth(1);

  // Functions

  RootStyle->SetFuncColor(1);
  RootStyle->SetFuncStyle(0);
  RootStyle->SetFuncWidth(2);

  //Legends 

  RootStyle->SetStatBorderSize(1);
  RootStyle->SetStatFont      (42);
//  RootStyle->SetOptStat       (1111);
   RootStyle->SetOptStat       (0);
  RootStyle->SetStatColor     (0);
//  RootStyle->SetStatX         (1.2);
//  RootStyle->SetStatY         (1.2);
   RootStyle->SetStatW         (0.25);
   RootStyle->SetStatH         (0.20);

  // Labels, Ticks, and Titles

  RootStyle->SetTickLength ( 0.015,"X");
  RootStyle->SetTitleSize  ( 0.045,"X");
  RootStyle->SetTitleColor ( 1    ,"X");
//  RootStyle->SetTitleOffset( 1.100,"X");
  RootStyle->SetTitleOffset( 1.000,"X");
  RootStyle->SetLabelOffset( 0.015,"X");
  RootStyle->SetLabelSize  ( 0.045,"X");
  RootStyle->SetLabelFont  ( 42   ,"X");
  RootStyle->SetNdivisions ( 505 ,"X");

  RootStyle->SetTickLength ( 0.015,"Y");
  RootStyle->SetTitleSize  ( 0.045,"Y");
  RootStyle->SetTitleOffset( 1.100,"Y");
  RootStyle->SetLabelOffset( 0.015,"Y");
  RootStyle->SetLabelSize  ( 0.045,"Y");
  RootStyle->SetLabelFont  ( 42   ,"Y");
  RootStyle->SetTitleFont  ( 42   ,"Y");
  RootStyle->SetNdivisions ( 505   ,"Y");

  RootStyle->SetTickLength ( 0.015,"Z");
  RootStyle->SetTitleSize  ( 0.060,"Z");
  RootStyle->SetTitleOffset( 1.100,"Z");
  RootStyle->SetLabelOffset( 0.015,"Z");
  RootStyle->SetLabelSize  ( 0.050,"Z");
  RootStyle->SetLabelFont  ( 42   ,"Z");
  RootStyle->SetTitleFont  ( 42   ,"Z");
  RootStyle->SetNdivisions ( 707   ,"Z");


  RootStyle->SetTitleBorderSize  (0);
  RootStyle->SetTitleFillColor  (0);  
  RootStyle->SetTitleFont  (42);
  RootStyle->SetTitleColor  (1);

  // Options

  RootStyle->SetOptFit     (11);
  RootStyle->SetOptStat    (0);
  RootStyle->SetStatFormat("8.5g");
//RootStyle->SetOptFit     (0);
// RootStyle->SetOptStat    (0);
 //  RootStyle->SetMarkerStyle(20);
  RootStyle->SetMarkerSize(0.9);

  RootStyle->SetPalette(42, NULL);

}
