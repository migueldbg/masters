#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TApplication.h>  
#include <TROOT.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TProfile.h>
#include <TRandom2.h>
#include <TArrow.h>
#include <TEllipse.h>
#include <TGraphAsymmErrors.h>

//--------------
//Root Init 
//--------------
//TApplication theApp("App",NULL,NULL);
TLatex l;

int ww=500;
int wh=500;

int wtopx=15;
int wtopy=15;

int canvas=0;
const int MaxCanvas=50;
TCanvas *cc[MaxCanvas];

void SetStyle()
{

  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0000000000);
 
  gStyle->SetMarkerStyle(20);
  gStyle->SetErrorX(0.);

  gStyle->SetTitleOffset(1.7,"Y"); 
  gStyle->SetLabelSize(0.035,"X");  gStyle->SetLabelSize(0.035,"Y");
  gStyle->SetTitleSize(0.045,"X");  gStyle->SetTitleSize(0.045,"Y");
  gStyle->SetPalette(1);

  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15); 

  gStyle->SetNdivisions(505,"X");

  gStyle->SetTextSize(0.035);
  
  l.SetNDC(true); 
  l.SetTextSize(0.04);

  
}
