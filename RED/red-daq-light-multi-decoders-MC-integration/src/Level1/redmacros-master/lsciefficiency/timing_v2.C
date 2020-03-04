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
  This macro reads a PAIR (Cf/background) of "light" files produced by 
  lscicalibrator(). It builds the tof spectrum for gamma events (f90 < 0.13 on both near 
  and far detector), by applying a +40 ns shift to the Cf data, and fits it with the four-peak 
  function DT2d by Marco, to derive baselines (near/far) and offset. The fit is run 
  channel-by-channel.

  The "high resolution" version of the input file is requires, i.e. in 
  which the reconstruction was performed with the interpolation of 
  the start_time.

  TOF spectra are written in a root file.
*/

//Fitting functions
Double_t DT2d(Double_t *x, Double_t *par);

using namespace std;

int timing_v2(TString fileBck,TString fileCf)
{
  Double_t threshold = 100; //threshold on far detectors, keV
  //Constants and parameters

  vector<TH1D*> tofs; //vector of tof spectra
  for (size_t i=16;i<24;i++)
    {
      tofs.push_back(new TH1D(Form("tof%d",i),
			      Form("Tof, far chan %d",i),
			      800,-20,60));    
    }  

  //Read data from file!
  TChain* tree1 = new TChain("t1");
  TChain* tree2 = new TChain("t1");

  //Background 
  if (fileBck)
    tree1->AddFile(fileBck);

  //Cf
  if (fileCf== " ")
    {
      tree2->AddFile("run_75_int.root.light");
      tree2->AddFile("run_77_int.root.light");
      tree2->AddFile("run_78_int.root.light");
      tree2->AddFile("run_79_int.root.light");
    }
  else
    tree2->AddFile(fileCf);

  cout << "Opening data files. Total entries: " << tree1->GetEntries() << " and " << 
    tree2->GetEntries() << endl;

  Double_t tof = 0;
  Double_t chargeFar = 0;
  Double_t chargeNear = 0;
  Double_t f90Far = 0;
  Double_t f90Near = 0;
  Int_t multiplicity = 0;
  Int_t nchan = 0;
  Int_t eventNb = 0;
  Int_t issaturated=0;
  tree1->SetBranchAddress("eventNb",&eventNb);
  tree1->SetBranchAddress("tof",&tof);
  tree1->SetBranchAddress("chargeNear",&chargeNear);
  tree1->SetBranchAddress("chargeFar",&chargeFar);
  tree1->SetBranchAddress("f90Near",&f90Near);
  tree1->SetBranchAddress("f90Far",&f90Far);
  tree1->SetBranchAddress("nchan",&nchan);
  tree1->SetBranchAddress("multiplicity",&multiplicity);
  tree1->SetBranchAddress("issaturated",&issaturated);
  tree2->SetBranchAddress("eventNb",&eventNb);
  tree2->SetBranchAddress("tof",&tof);
  tree2->SetBranchAddress("chargeNear",&chargeNear);
  tree2->SetBranchAddress("chargeFar",&chargeFar);
  tree2->SetBranchAddress("f90Near",&f90Near);
  tree2->SetBranchAddress("f90Far",&f90Far);
  tree2->SetBranchAddress("nchan",&nchan);
  tree2->SetBranchAddress("multiplicity",&multiplicity);
  tree2->SetBranchAddress("issaturated",&issaturated);

  // Event loop: first file
  for (Int_t iev=0;iev<tree1->GetEntries();iev++)
    {
      tree1->GetEntry(iev);

      if (f90Near < 0.13 && multiplicity==1 && f90Far<0.13 && chargeFar>threshold)
	{
	  if (tof > -20 && tof < 20)
	  // Fill tof spectrum
	    tofs.at(nchan-16)->Fill(tof); 	  
	}
    }

  for (Int_t iev=0;iev<tree2->GetEntries();iev++)
    {
      tree2->GetEntry(iev);

      if (f90Near < 0.13 && multiplicity==1 && f90Far<0.13 && chargeFar>threshold)
	{
	  if (tof > -20 && tof < 20)
	    // Fill tof spectrum
	    tofs.at(nchan-16)->Fill(tof+40.); 	  
	}
    }
  //Done: data read.
  cout << "All data read! " << endl;

  //Setup fitting function
  TF1 * d2d = new TF1("d2d",DT2d,-20.,60.,12);
  d2d->SetNpx(1000);
  Double_t mypar[12]={5.,100.,1.2,1000.,0.5,1,0.8,10,0.5,10,10,0.1};
  d2d->SetParameters(mypar);
  d2d->SetParName(0,"L_{NEAR}");
  d2d->SetParName(1,"L_{FAR}");
  d2d->SetParName(2,"Offset");
  d2d->SetParName(3,"N ph Cf");
  d2d->SetParName(4,"sigma1");
  d2d->SetParName(5,"sigma2");
  d2d->SetParName(6,"fraction 1");
  d2d->SetParName(7,"N ph Cf Bk");
  d2d->SetParName(8,"Random  Cf");
  d2d->SetParName(9,"N ph bkg 1");
  d2d->SetParName(10,"N ph bkg 2");
  d2d->SetParName(11,"Random Bkg");

  d2d->SetParLimits(0,0,20);
  d2d->SetParLimits(1,80,120);
  //d2d->SetParLimits(6,0,1);

  Int_t validChannels=0;
  for (size_t i=0;i<tofs.size();i++)
    {
      if (tofs.at(i)->GetEntries())
	{
	  validChannels++;
	  Double_t maxpos = tofs.at(i)->GetBinCenter(tofs.at(i)->GetMaximumBin());
	  //cout << i << " " << tofs.at(i)->GetBinCenter(tofs.at(i)->GetMaximumBin()) << endl;
	  //d->SetParameters(100.,1.2,1000.,0.5,1,0.4,10,0.5);
	  //d->SetParameters(100.,maxpos-3.,1000.,0.5,1,0.4,1,0.5);
	  //Marco's line
	  //tofs.at(i)->Fit(d,"LQ","",-20,7.5);
	  tofs.at(i)->Fit(d2d,"LQ","",-20.,maxpos+2.5);
	  //tofs.at(i)->Fit(d2d,"LQ","",-20.,46.);
	  cout << "Ch" << i << "--> " << 
	    "L_NEAR: (" << d2d->GetParameter(0) << " +/- " << d2d->GetParError(0) << 
	    ") L_FAR: (" << d2d->GetParameter(1) << " +/- " << d2d->GetParError(1) << 
	    ") Offset: (" << d2d->GetParameter(2) << " +/- " << d2d->GetParError(2) << 
	    ")" << endl; 

	}
    }

  TCanvas* c1 = new TCanvas();
  
  Int_t pixel = ((Int_t) std::sqrt(validChannels))+1;  
  c1->Divide(pixel,pixel);
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
  TFile* output = new TFile("timingv2.root","RECREATE");
  output->cd();
  for (size_t i=0;i<tofs.size();i++)
    tofs.at(i)->Write(tofs.at(i)->GetName());
  output->Write();
  output->Close();
    

  return 0;
 
}


Double_t DT2d(Double_t *x, Double_t *par)
{
  Double_t xx[1];
  Double_t fun;

  Double_t Lnear=par[0];
  Double_t BaseLine=par[1];
  Double_t Offset=par[2];

  if (x[0]>20) {
    xx[0]=x[0]-40.;
    
    Double_t DTCf = (BaseLine-Lnear)/29.98+Offset;  // taking into account in the CF252 fit the flight path from the source
    //    Double_t DTCf = (par[1]-2*Lnear)/29.98+par[2];  // taking into account in the CF252 fit the flight path from the source
    Double_t NphCf = par[3];
    Double_t sigmaCf = par[4];
    Double_t sigmaCf2 = par[5];
    Double_t f1 = par[6];
    Double_t phCf = NphCf*(f1*1./TMath::Sqrt(2*3.1415)/sigmaCf*TMath::Exp(-0.5*(TMath::Power((xx[0]-DTCf)/sigmaCf,2))) + (1-f1)*1./TMath::Sqrt(2*3.1415)/sigmaCf2*TMath::Exp(-0.5*(TMath::Power((xx[0]-DTCf)/sigmaCf2,2))));
    
    Double_t NphBkCf = par[7];
    Double_t DTBkCf = Offset-(BaseLine+Lnear)/29.98;
    //    Double_t DTBkCf = par[2]-par[1]/29.98;
    
    Double_t phBkCf = NphBkCf*(f1*1./TMath::Sqrt(2*3.1415)/sigmaCf*TMath::Exp(-0.5*(TMath::Power((xx[0]-DTBkCf)/sigmaCf,2))) + (1-f1)*1./TMath::Sqrt(2*3.1415)/sigmaCf2*TMath::Exp(-0.5*(TMath::Power((xx[0]-DTBkCf)/sigmaCf2,2))));

    // costant random background
    Double_t costBk=par[8];
    
    fun = phCf+phBkCf+costBk;
    
  } else {
    xx[0]=x[0];

    Double_t DTBk1 = (BaseLine+Lnear)/29.98+Offset; 
    //    Double_t DTBk1 = par[1]/29.98+par[2]; 
    Double_t NphBk1 = par[9];
    Double_t sigmaCf = par[4];
    Double_t sigmaCf2 = par[5];
    Double_t f1 = par[6];

    // Two gauss resolution function
    Double_t phBk1 = NphBk1*(f1*1./TMath::Sqrt(2*3.1415)/sigmaCf*TMath::Exp(-0.5*(TMath::Power((xx[0]-DTBk1)/sigmaCf,2))) + (1-f1)*1./TMath::Sqrt(2*3.1415)/sigmaCf2*TMath::Exp(-0.5*(TMath::Power((xx[0]-DTBk1)/sigmaCf2,2))));
    

    Double_t NphBk2 = par[10];
    Double_t DTBk2 = Offset-(BaseLine+Lnear)/29.98; 
    //    Double_t DTBk2 = par[2]-par[1]/29.98;

    Double_t phBk2 = NphBk2*(f1*1./TMath::Sqrt(2*3.1415)/sigmaCf*TMath::Exp(-0.5*(TMath::Power((xx[0]-DTBk2)/sigmaCf,2))) + (1-f1)*1./TMath::Sqrt(2*3.1415)/sigmaCf2*TMath::Exp(-0.5*(TMath::Power((xx[0]-DTBk2)/sigmaCf2,2))));

    // costant random background
    Double_t costBk=par[11];

    // correlated
    // Double_t NphCorBk = par[8];
    // Double_t sigmaCorBk = par[9]; // use same resol function
    // Double_t phCorBk = NphCorBk/TMath::Sqrt(2*3.1415)/sigmaCorBk*TMath::Exp(-0.5*(TMath::Power((x[0]-0.0)/sigmaCorBk,2)));

    fun= phBk1+phBk2+costBk;//+phCorBk;
    
  }
  
  return fun;

}

