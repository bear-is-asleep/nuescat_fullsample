#include "TStyle.h"
#include "TColor.h"

void selectionstyle()
{
  //  TStyle *henrySBNDStyle= new TStyle("henrySBND","My SBND style for final plots");
  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameLineWidth(2);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  //gStyle->SetFillColor(0);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendTextSize(0.04);
  

  // set the paper & margin sizes
  gStyle->SetPaperSize(30,39);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.2);

  // use large Times-Roman fonts
  gStyle->SetTextFont(52);
  gStyle->SetTextSize(0.08);
  gStyle->SetLabelFont(52,"x");
  gStyle->SetLabelFont(52,"y");
  gStyle->SetLabelFont(52,"z");
  gStyle->SetLabelSize(0.06,"x");
  gStyle->SetTitleSize(0.06,"x");
  gStyle->SetLabelSize(0.06,"y");
  gStyle->SetTitleSize(0.06,"y");
  gStyle->SetLabelSize(0.06,"z");
  gStyle->SetTitleSize(0.06,"z");
  gStyle->SetLabelFont(52,"t");
  gStyle->SetTitleFont(52,"x");
  gStyle->SetTitleFont(52,"y");
  gStyle->SetTitleFont(52,"z");
  gStyle->SetTitleFont(52,"t"); 
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetTitleOffset(1.4,"y");
  gStyle->SetTitleOffset(1,"z");
  gStyle->SetTitleOffset(1.4,"pad");
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetTitleFont(52,"pad");
  gStyle->SetTitleBorderSize(0);

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetHistLineWidth(2);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  //gStyle->SetErrorX(0.001);



  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);


  // Add a greyscale palette for 2D plots
  int ncol=50;
  double dcol = 1./float(ncol);
  double gray = 1;
  TColor **theCols = new TColor*[ncol];
  for (int i=0;i<ncol;i++) theCols[i] = new TColor(999-i,0.0,0.7,0.7);
  for (int j = 0; j < ncol; j++) {
    theCols[j]->SetRGB(gray,gray,gray);
    gray -= dcol;
  }
  int ColJul[100];
  for  (int i=0; i<100; i++) ColJul[i]=999-i;
  // gStyle->SetPalette(ncol,ColJul);


 
  // Define a nicer color palette (red->blue)
  // Uncomment these lines for a color palette (default is B&W)
  gStyle->SetPalette(1,0);  // use the nice red->blue palette
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
				   NCont);
  gStyle->SetNumberContours(NCont); 
}