#include "TROOT.h"
#include "TH1D.h"
#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"
#include <iostream>
#include "TCanvas.h"
#include "TLegend.h"

using namespace std;

/* 
   This macro is meant to compare the SiMon and the SiTel spectra 
   taken in different runs, to check for possible drifts in the energy 
   scale. The runs are supposed to be taken either 
   in "flower" or in "triple" mode (i.e. the Si detectors are channels 16 & ff 
   instead of 0 & ff, as in "banana" mode). 

   The 1D SiMon spectra are normalized such to have the same area.
   The 2D E/DeltaE spectra are displayed as a color-plot and a 
   countour-plot. 

   The inputs are the two run numbers. ROOT files are expected to be found 
   in the current directory.
*/

int checkSiDetectors(Int_t run1,Int_t run2=193)
{
  TFile f1(Form("run_%d.root",run1));
  //TFile f1("run_199.root");
  TFile f2(Form("run_%d.root",run2));
  
  TTree* r1 = (TTree*) f1.Get("reco");
  TTree* r2 = (TTree*) f2.Get("reco");
  
  TH1D* h1 = new TH1D("h1",Form("Run %d, Monitor",run1),
		      1000,0,12000);
  TH1D* h2 = new TH1D("h2",Form("Run %d, Monitor",run2),1000,0,12000);

  TH2D* h3 = new TH2D("h3",Form("Run %d, Telescope",run1),
		      400,200,10000,400,0,8000);
  TH2D* h4 = new TH2D("h4",Form("Run %d, Monitor",run2),
		      400,200,10000,400,0,8000);

  h2->SetLineColor(kRed);
  
  r1->Project("h1","ymax[20]-baseline_mean[20]","ymax[19]-baseline_mean[19]<100 && ymax[18]-baseline_mean[18]<100");
  r2->Project("h2","ymax[20]-baseline_mean[20]","ymax[19]-baseline_mean[19]<100 && ymax[18]-baseline_mean[18]<100");
  r1->Project("h3","ymax[18]-baseline_mean[18]:ymax[19]-baseline_mean[19]","ymax[20]-baseline_mean[20]<100");
  r2->Project("h4","ymax[18]-baseline_mean[18]:ymax[19]-baseline_mean[19]","ymax[20]-baseline_mean[20]<100");

  cout << "SiMon, Run " << run1 <<" : " << h1->GetEntries() << "/" << r1->GetEntries() << " " << 
    h1->GetEntries()/(Double_t)r1->GetEntries() << endl;
  cout << "SiMon, Run " << run2 << " : " << h2->GetEntries() << "/" << r2->GetEntries() << " " <<
    h2->GetEntries()/(Double_t)r2->GetEntries() << endl;
  
  h1->Sumw2();
  h2->Sumw2();
  h1->Scale(1./h1->GetEntries());
  h2->Scale(1./h2->GetEntries());
  h2->SetLineColor(kRed);

  h1->DrawCopy("HIST");
  h2->DrawCopy("HISTSAME");
  
  TLegend* l1 = new TLegend(0.1,0.7,0.48,0.9); 
  l1->AddEntry(h1,Form("Run %d",run1),"l");
  l1->AddEntry(h2,Form("Run %d",run2),"l");
  l1->Draw("same");

  TCanvas* c2 = new TCanvas();
  //c2->SetLogz();
  c2->cd();
  h3->SetFillColor(1);
  h4->SetLineColor(kRed);
  h4->SetContour(5);
  h3->DrawCopy("col");
  h4->DrawCopy("cont3 same");

  return 0;
}
