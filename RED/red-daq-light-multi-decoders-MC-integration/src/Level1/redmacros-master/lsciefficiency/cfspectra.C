#include "TTree.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVirtualPad.h"
#include "TChain.h"
#include "TMath.h"
#include <iostream>

using namespace std;

Double_t K(Double_t t)
{
  if (t<56.8)
    return 0.49078+0.0099946*t+0.0008073*t*t -2.6375e-5*t*t*t + 2.1063e-7*t*t*t*t;
  else
    return 0.973+0.2563*TMath::Exp(-0.02902*t);      
}

Double_t CubFunction(Double_t *x, Double_t *par)
{
  Double_t t =x[0];
  Double_t a = 8.345;
  Double_t b = 3681.04;
  Double_t tsq=t*t;


  Double_t f = TMath::Sqrt(b*b*b/TMath::Pi())*4./(tsq*tsq)*(1+(5*a)/(2*tsq))*
    TMath::Exp(-(b/tsq)*(1+(a/tsq)))*K(t);
  cout << t << " " << f << " " << K(t) << " " << par[0]*f << endl;

  return (par[0]*f);
}


int cfspectra()
{
  const Double_t mN = 939.5653; //MeV
  const Double_t c_light = 30.0; //cm/ns

  TF1* cffun = new TF1("Cf252-watt",
		       "[0]*TMath::Exp(-0.88*x)*TMath::SinH(TMath::Sqrt(2.*x))",
		       0,12);
  cffun->SetParameter(0,1);
  Double_t int1 = cffun->Integral(0.,12.);
  cffun->SetParameter(0,1./int1);

  TF1* maxfun = new TF1("Cf252-watt",
			"[0]*TMath::Sqrt(x)*TMath::Power([1],1.5)*TMath::Exp(-x/[1])",
			0,12);
  
  maxfun->SetParameter(1,1.406);
  maxfun->SetParameter(0,1);
  int1 = maxfun->Integral(0.,12.);
  maxfun->SetParameter(0,1./int1);

  
  
  TF1 *cub= new TF1("Cf252-Cub",CubFunction,10,200,1);
  cub->SetParameter(0,1);
  TH1D* h1 = new TH1D("Cf252-Cub-E","Energy spectrum",1000,0,12);
  Int_t nloops = 10000000;
  for (Int_t i=0;i<nloops;i++)
    {
      Double_t t = cub->GetRandom();
      //t is time on a 100 cm baseline
      Double_t beta = (100./t)/c_light; //
      Double_t gamma = 1./TMath::Sqrt(1-beta*beta);
      Double_t ekin = mN*(gamma-1);
      h1->Fill(ekin);
    }
  int1 = h1->Integral(1,h1->GetNbinsX(),"width");
  h1->Scale(1./int1);
  cffun->Draw();
  maxfun->SetLineColor(kBlue);
  maxfun->Draw("");
  cffun->SetLineColor(kRed);
  cffun->Draw("same");
  h1->Draw("same");



  TCanvas* c1 = new TCanvas();
  c1->cd();
  cub->Draw();

  

  return 1;


}
