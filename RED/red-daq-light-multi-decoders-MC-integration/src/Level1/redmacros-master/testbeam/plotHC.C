#include "TH2Poly.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCut.h"
#include "TCutG.h"
#include "TChain.h"
#include "TLatex.h"
#include <iostream>
#include "red-daq/EvRec0.hh"

/*
  This macro plots the monitor-Si spectrum for the four runs 
  which used the 353 ug/cm2 target.

  The C peaks are normalized, such to show the depletion of the 
  H component.

*/




int plotHC(Bool_t useAu=false)
{
  
 
  gStyle->SetOptStat(0);
  vector<TH1D*> *spectra = new vector<TH1D*>;
  vector<TH1D*> *spectraAu = new vector<TH1D*>;
  Int_t index = (useAu) ? 4 : 20;
  vector<TString> fnames;
  if (useAu)
    fnames =  
    {
      "run_99",
      "run_105",
      "run_108",
      "run_111"
    };
  else
    fnames = {
      "run_102",
      "run_104",
      "run_109",
      "run_112"
    };
  
 

  vector<TString> fnamesAu = {
    "run_100",
    "run_106",
    "run_107","run_110"};
  vector<Int_t> rawCounter(fnames.size(),0);
  for (size_t i=0;i<fnames.size();i++)
    {
      spectra->push_back(new TH1D(Form("h%d",i),fnames.at(i),
				  1000,0,5200));
      spectraAu->push_back(new TH1D(Form("Au%d",i),fnamesAu.at(i),
				  1000,0,5200));
    }

  for (size_t i=0;i<fnamesAu.size();i++)
    {  
      TFile* f = new TFile(fnamesAu.at(i)+".root");
      TTree* reco = (TTree*) f->Get("reco");
      EvRec0* evReco = new EvRec0();
      reco->SetBranchAddress("recoevent",&evReco);
      
      for (Int_t iloop=0;iloop<reco->GetEntries();iloop++)
	{
	  reco->GetEntry(iloop);
	  vector<double> ymax = evReco->GetYmax();
	  vector<double> basemean = evReco->GetBaseMean();
	  if (ymax.at(4)-basemean.at(4) > 100)
	    spectraAu->at(i)->Fill(ymax.at(4)-basemean.at(4));
	}
      f->Close();
    }


  for (size_t i=0;i<fnames.size();i++)
    {  
      TFile* f = new TFile(fnames.at(i)+".root");
      TTree* reco = (TTree*) f->Get("reco");
      EvRec0* evReco = new EvRec0();
      reco->SetBranchAddress("recoevent",&evReco);
      
      for (Int_t iloop=0;iloop<reco->GetEntries();iloop++)
	{
	  reco->GetEntry(iloop);
	  vector<double> ymax = evReco->GetYmax();
	  vector<double> basemean = evReco->GetBaseMean();
	  if (ymax.at(index)-basemean.at(index) > 100)
	    spectra->at(i)->Fill(ymax.at(index)-basemean.at(index));
	  //if (ymax.at(20)-basemean.at(20) > 4700)
	  // (rawCounter.at(i))++;
	}
      f->Close();
    }
  
  cout << spectra->at(0)->GetEntries() << endl;
  
  TF1* gaus = new TF1("gaus","gaus(0)");
  vector<Double_t> mu;
  vector<Double_t> N;
  for (size_t i=0;i<spectra->size();i++)
    {
      if (useAu && fnamesAu.at(i) == "run_106")
	spectraAu->at(i)->Scale(92./82.); //thinner Au
      if (useAu)
	spectraAu->at(i)->Fit(gaus,"N","",3350,3500);
      else
	spectra->at(i)->Fit(gaus,"N","",4700,4850);
      mu.push_back(gaus->GetParameter(1));
      N.push_back(TMath::Sqrt(2.*TMath::Pi())*gaus->GetParameter(0)*
		   gaus->GetParameter(2));
    }

   for (size_t i=0;i<spectra->size();i++)
    {
      cout << "Mu: " << mu.at(i) << " N: " << N.at(i) << " " << 
	rawCounter.at(i) << " " << N.at(i)/N.at(0) << " " << endl;
      spectra->at(i)->Sumw2();
      spectra->at(i)->Scale(N.at(0)/N.at(i));
      if (useAu)
	{
	  spectraAu->at(i)->Sumw2();
	  spectraAu->at(i)->Scale(N.at(0)/N.at(i));
	}
    }

   TCanvas* c1 = new TCanvas();
   if (useAu)
     {
       c1->Divide(1,2);
       c1->cd(1);
     }
   else 
     c1->cd();
   for (size_t i=0;i<fnames.size();i++)
     {
       spectra->at(i)->SetLineColor(kBlack+i);
       if (!i)
	 {
	   spectra->at(i)->SetTitle("CH_{2}");	   
	   spectra->at(i)->GetXaxis()->SetTitle("Amplitude (a.u.)");
	   spectra->at(i)->Draw("HIST");
	 }
       else
	 spectra->at(i)->Draw("sameHIST");
    }
  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  for (size_t i=0;i<fnames.size();i++)
    legend->AddEntry(spectra->at(i),fnames.at(i),"lep");
  legend->Draw("same");

  if (useAu)
    {      
      TVirtualPad* p1 = c1->cd(2);
      p1->SetLogy();
      for (size_t i=0;i<fnamesAu.size();i++)
	{
	  spectraAu->at(i)->SetLineColor(kBlack+i);
	  if (!i)
	    {
	      spectraAu->at(i)->SetTitle("Au");	  
	      spectraAu->at(i)->GetXaxis()->SetTitle("Amplitude (a.u.)");
	      spectraAu->at(i)->Draw("HIST");
	    }
	  else
	spectraAu->at(i)->Draw("sameHIST");
	}
      auto legend2 = new TLegend(0.1,0.7,0.48,0.9);
      for (size_t i=0;i<fnamesAu.size();i++)
	legend2->AddEntry(spectraAu->at(i),fnamesAu.at(i),"lep");
      legend2->Draw("same");
    }      


  return 0;
  
}
