#include "TTree.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TVirtualPad.h"
#include "TChain.h"
#include "TMath.h"
#include <iostream>

/*
  This macro reads a "light" file produced by lscicalibrator(), 
  builds the tof spectrum for gamma events (f90 < 0.13 on both near 
  and far detector) and fit it with the two-peak function DT by 
  Marco, to derive baseline and offset. The fit is run channel-by-channel.

  The "high resolution" version of the input file is requires, i.e. in 
  which the reconstruction was performed with the interpolation of 
  the start_time.

  TOF spectra are written in a root file.
*/

//Fitting function
Double_t DT(Double_t *x, Double_t *par);

using namespace std;

int timing(TString filename)
{
  Double_t threshold = 100; //threshold on far detectors, keV
  //Constants and parameters

  vector<TH1D*> tofs; //vector of tof spectra
  for (size_t i=16;i<24;i++)
    {
      tofs.push_back(new TH1D(Form("tof%d",i),
			      Form("Tof, far chan %d",i),
			      400,-20,20));    
    }  

  //Read data from file!
  TChain* tree = new TChain("t1");
  
  if (filename!=" ")
    tree->AddFile(filename);
  else
    {
      tree->AddFile("run_75_int.root.light");
      tree->AddFile("run_77_int.root.light");
      tree->AddFile("run_78_int.root.light");
      tree->AddFile("run_79_int.root.light");
    }

  cout << "Opening data files. Total entries: " << tree->GetEntries() << endl;

  Double_t tof = 0;
  Double_t chargeFar = 0;
  Double_t chargeNear = 0;
  Double_t f90Far = 0;
  Double_t f90Near = 0;
  Int_t multiplicity = 0;
  Int_t nchan = 0;
  Int_t eventNb = 0;
  Int_t issaturated=0;
  tree->SetBranchAddress("eventNb",&eventNb);
  tree->SetBranchAddress("tof",&tof);
  tree->SetBranchAddress("chargeNear",&chargeNear);
  tree->SetBranchAddress("chargeFar",&chargeFar);
  tree->SetBranchAddress("f90Near",&f90Near);
  tree->SetBranchAddress("f90Far",&f90Far);
  tree->SetBranchAddress("nchan",&nchan);
  tree->SetBranchAddress("multiplicity",&multiplicity);
  tree->SetBranchAddress("issaturated",&issaturated);

  // Event loop
  for (Int_t iev=0;iev<tree->GetEntries();iev++)
    {
      tree->GetEntry(iev);

      if (f90Near < 0.13 && multiplicity==1 && f90Far<0.13 && chargeFar>threshold)
	{
	  // Fill tof spectrum
	  tofs.at(nchan-16)->Fill(tof); 	  
	}
    }
  
  //Done: data read.
  
  //Setup fitting function
  TF1 * d = new TF1("DT",DT,-20.,20.,8);

  d->SetParName(0,"Baseline");
  d->SetParName(1,"Offset");
  d->SetParName(2,"#Cf Photon");
  d->SetParName(3,"sigma1");
  d->SetParName(4,"sigma2");
  d->SetParName(5,"fraction");
  d->SetParName(6,"#Bk Photon");
  d->SetParName(7,"#Random Bk");
  d->SetParLimits(5,0,1);

  for (size_t i=0;i<tofs.size();i++)
    {
      if (tofs.at(i)->GetEntries())
	{
	  Double_t maxpos = tofs.at(i)->GetBinCenter(tofs.at(i)->GetMaximumBin());
	  //cout << i << " " << tofs.at(i)->GetBinCenter(tofs.at(i)->GetMaximumBin()) << endl;
	  //d->SetParameters(100.,1.2,1000.,0.5,1,0.4,10,0.5);
	  d->SetParameters(100.,maxpos-3.,1000.,0.5,1,0.4,1,0.5);
	  //Marco's line
	  //tofs.at(i)->Fit(d,"LQ","",-20,7.5);
	  tofs.at(i)->Fit(d,"LQ","",-18.,maxpos+4);
	  cout << "Ch" << i << "--> " << 
	    "Baseline: (" << d->GetParameter(0) << " +/- " << d->GetParError(0) << 
	    ") Offset: (" << d->GetParameter(1) << " +/- " << d->GetParError(1) << ")" << endl; 
	}
    }

  TCanvas* c1 = new TCanvas();
  c1->Divide(3,3);
  for (size_t i=0;i<tofs.size();i++)
    {
      if (tofs.at(i)->GetEntries())
	{
	  TVirtualPad* pad = c1->cd(i+1);
	  pad->SetLogy();
	  tofs.at(i)->DrawCopy();
	}
    }
   
  TCanvas* c2 = new TCanvas();
  c2->cd();
  tofs.at(0)->DrawCopy();

  //Save everything on file!
  TFile* output = new TFile("timing.root","RECREATE");
  output->cd();
  for (size_t i=0;i<tofs.size();i++)
    tofs.at(i)->Write(tofs.at(i)->GetName());
  output->Write();
  output->Close();
    

  return 0;
 
}

Double_t DT(Double_t *x, Double_t *par)
{

  Double_t DTCf = par[0]/29.98+par[1];
  Double_t NphCf = par[2];
  Double_t sigmaCf = par[3];
  Double_t sigmaCf2 = par[4];
  Double_t f1 = par[5];


  //  Double_t phCf = NphCf/TMath::Sqrt(2*3.1415)/sigmaCf*TMath::Exp(-0.5*(TMath::Power((x[0]-DTCf)/sigmaCf,2)));
  Double_t phCf = NphCf*(f1*1./TMath::Sqrt(2*3.1415)/sigmaCf*TMath::Exp(-0.5*(TMath::Power((x[0]-DTCf)/sigmaCf,2))) + 
			 (1-f1)*1./TMath::Sqrt(2*3.1415)/sigmaCf2*TMath::Exp(-0.5*(TMath::Power((x[0]-DTCf)/sigmaCf2,2))));


  Double_t NphBk = par[6];
  Double_t sigmaBk = par[7];
  Double_t DTBk = par[1]-par[0]/29.98;

  Double_t phBk = NphBk*(f1*1./TMath::Sqrt(2*3.1415)/sigmaCf*TMath::Exp(-0.5*(TMath::Power((x[0]-DTBk)/sigmaCf,2))) + 
			 (1-f1)*1./TMath::Sqrt(2*3.1415)/sigmaCf2*TMath::Exp(-0.5*(TMath::Power((x[0]-DTBk)/sigmaCf2,2))));

  //  Double_t phBk = NphBk/TMath::Sqrt(2*3.1415)/sigmaBk*TMath::Exp(-0.5*(TMath::Power((x[0]-DTBk)/sigmaBk,2)));

  //  Double_t NphCfBs = par[6];
  //  Double_t sigmaCfBs = par[7];
  //  Double_t DTCfBs = par[8];
  //
  //Double_t phCfBs = NphCfBs/TMath::Sqrt(2*3.1415)/sigmaCfBs*TMath::Exp(-0.5*(TMath::Power((x[0]-DTCfBs)/sigmaCfBs,2)));

  Double_t costBk=par[7];

 
  return phCf+phBk+costBk;


}
