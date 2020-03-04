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

using namespace std;

/* 
   This macro reads "light" files from the Cf252 campaign of February 2019 and produced 
   a set of distribution of the PSD parameter in different ranges of S1. 

   The energy range between 0 and 3000 keV is divided in 20 slices (50 keVee each). For 
   each slice, the PSD distribution is built, which is fitted with the sum of two gaussians. 
   Since the high-energy part of the spectrum is usually empty, only the first 16 slices are 
   considered in the procedure (and in some cases the fits fail due to the lack of 
   statistics). 

   The upper edge can be changed by
      Double_t emax = 3000.; //keV
   
   The input parameters are the filename and the number of slices, e.g.
    scanf90("run_78_int.root.light",20)
    scanf90("run_48_int.root.light",20)

  The part for the near detector is presented but commented, FoM distribution 
  are less clear because of the presence of mixed n-gamma events.

*/

int scanf90(TString filename, Int_t nslices=20)
{
  TFile* ifile = new TFile(filename);
  TTree* tree = (TTree*) ifile->Get("t1");
  
  Double_t tof = 0;
  Double_t chargeFar = 0;
  Double_t chargeNear = 0;
  Double_t f90Far = 0;
  Double_t f90Near = 0;
  Int_t multiplicity = 0;
  Int_t nchan = 0;
  Int_t eventNb = 0;
  tree->SetBranchAddress("eventNb",&eventNb);
  tree->SetBranchAddress("tof",&tof);
  tree->SetBranchAddress("chargeNear",&chargeNear);
  tree->SetBranchAddress("chargeFar",&chargeFar);
  tree->SetBranchAddress("f90Near",&f90Near);
  tree->SetBranchAddress("f90Far",&f90Far);
  tree->SetBranchAddress("nchan",&nchan);
  tree->SetBranchAddress("multiplicity",&multiplicity);

  vector<TH1D*> f90s;
  vector<TH2D*> f90vsE;

  Double_t emax = 3000.; //keV

  for (size_t i=16;i<24;i++)
    {
      f90s.push_back(new TH1D(Form("psd%d",i),
			      Form("PSD, far chan %d",i),
			      200,0.,0.5));      

      f90vsE.push_back(new TH2D(Form("psdvsE%d",i),
				Form("PSD vs. energy, chan %d",i),
				nslices,0,emax,
				100,0,0.5));
    }

  
  //Read data
  for (Int_t iev=0;iev<tree->GetEntries();iev++)
    {
      tree->GetEntry(iev);
      if (multiplicity==1)
	{
	  f90s.at(nchan-16)->Fill(f90Far);
	  f90vsE.at(nchan-16)->Fill(chargeFar,f90Far);
	}
    }

  //Find minimum of PSD spectra:
  //Full PSD spectra (over all energies)
  TF1* gaus = new TF1("gaus","gaus(0)",0.,0.3);
  TF1* fun = new TF1("fun","gaus(0)+gaus(3)",0.,0.3);
  Double_t pars[6];
  Double_t xlow = 0.055;
  Double_t xmed = 0.12;
  Double_t xhigh = 0.23;    
  Int_t filledSpectra = 0;
  for (size_t i=0;i<f90s.size();i++)
    {
      if (f90s.at(i)->GetEntries() < 10)
	continue;
      filledSpectra++;
      f90s.at(i)->Fit(gaus,"Q","",xlow,xmed);
      gaus->GetParameters(&pars[0]);
      f90s.at(i)->Fit(gaus,"Q","",xmed,xhigh);
      gaus->GetParameters(&pars[3]);
      //for (size_t j=0;j<6;j++)
      //cout << "par[" << j << "] = " << pars[j] << endl;
      fun->SetParameters(pars);
      fun->SetLineColor(kBlue);
      f90s.at(i)->Fit(fun,"Q","",xlow,xhigh);
      //Write results
      /*
      cout << "Chan: " << i+16 << " gamma:   " << fun->GetParameter(1) 
	   << " +/-" << fun->GetParameter(2) << 
	endl; 
      cout << "         neutron: " << fun->GetParameter(4) 
	   << " +/-" << fun->GetParameter(5) << 
	endl; 
      */
    }
  cout << "I have found: " << filledSpectra << " valid channels" << endl;

  TCanvas* c0 = new TCanvas();
  Int_t siz = (Int_t) TMath::Sqrt(filledSpectra) +1;
  c0->Divide(siz,siz);
  for (size_t i=0;i<filledSpectra;i++)
    {
      c0->cd(i+1);
      f90s.at(i)->DrawCopy();
    }
  
  //Slices in energy. All channels
  TCanvas* c1 = new TCanvas();
  c1->Divide(siz,siz);
  for (size_t i=0;i<filledSpectra;i++)
    {
      TVirtualPad* pad = c1->cd(i+1);
      pad->SetLogz();
      f90vsE.at(i)->DrawCopy("COLZ");
    }

  vector<TCanvas*> theSlices(filledSpectra,nullptr);
  vector<TGraphErrors*> FoM;
  for (size_t idet=0;idet<filledSpectra;idet++)
    {
      theSlices.at(idet) = new TCanvas();
      theSlices.at(idet)->Divide(4,4);
      xlow = 0.035;
      xmed = 0.14;
      xhigh = 0.28;   
      FoM.push_back(new TGraphErrors());
      for (int i=1;i<=std::min(16,nslices);i++)
	{
	  theSlices.at(idet)->cd(i);
	  cout << "Processing slice: " << i << endl;
	  TH1D* htemp = f90vsE.at(idet)->ProjectionY("",i,i+1);
	  htemp->Fit(gaus,"Q","",xlow,xmed);
	  gaus->GetParameters(&pars[0]);
	  htemp->Fit(gaus,"Q","",xmed,xhigh);
	  gaus->GetParameters(&pars[3]);
	  //for (size_t j=0;j<6;j++)
	  //cout << "par[" << j << "] = " << pars[j] << endl;
	  fun->SetParameters(pars);
	  fun->SetLineColor(kBlue);
	  htemp->Fit(fun,"QL","",xlow,xhigh);
	  htemp->DrawCopy();
	  cout << "Slice: " << i << " gamma:   " << fun->GetParameter(1) 
	       << " +/-" << fun->GetParameter(2) << 
	    endl; 
	  cout << "         neutron: " << fun->GetParameter(4) 
	       << " +/-" << fun->GetParameter(5) << endl;
	  double fom = (fun->GetParameter(4)-fun->GetParameter(1))
	    /TMath::Sqrt(fun->GetParameter(2)*fun->GetParameter(2)+fun->GetParameter(5)*fun->GetParameter(5));
	  FoM.at(idet)->SetPoint(i-1,(i-0.5)*emax/nslices,fom);
	  FoM.at(idet)->SetPointError(i-1,0.5*emax/nslices,0.);
	}
    }
  cout << "I have created " << FoM.size() << " FOM plots" << endl;

  if (FoM.size())
    {
      TCanvas* c2 = new TCanvas();
      c2->cd();
      FoM.at(0)->Draw("AZP");
      FoM.at(0)->SetMarkerStyle(20);
      for (size_t i=1;i<FoM.size();i++)
	{
	  FoM.at(i)->SetMarkerStyle(20+i);	  
	  FoM.at(i)->Draw("ZP");
	}

    }

  /*

  //Now take the near detector! The 2D distribution is already available 
  //in the root file
  TH2D* f90vsEnear = (TH2D*) ifile->Get("h2");
  
  //Slices in energy. 
  TCanvas* c3 = new TCanvas();
  TVirtualPad* pad2 = c3->cd();
  pad2->SetLogz();
  f90vsEnear->DrawCopy("COLZ");
  
  Int_t nbinsX = f90vsEnear->GetNbinsX();

  cout << "The distribution has " << nbinsX << " bins in x" << endl;
  Int_t span = nbinsX*(emax/5000.)/nslices;
  Double_t binsize = 5000./nbinsX ;
  cout << "Spanning in " << span << ". Total bin size: " << binsize << endl;

  xlow = 0.035;
  xmed = 0.14;
  xhigh = 0.28;   
  TGraphErrors* FoMnear = new TGraphErrors();
  FoMnear->SetMarkerStyle(20);
  TCanvas* c4 = new TCanvas();
  c4->Divide(siz,siz);
  for (int i=1;i<=std::min(16,nslices);i++)
   {
     c4->cd(i);
     TH1D* htemp = f90vsEnear->ProjectionY("",1+(i-1)*span,i*span);
     cout << 1+(i-1)*span << " to" << i*span << endl;
     htemp->Fit(gaus,"Q","",xlow,xmed);
     gaus->GetParameters(&pars[0]);
     htemp->Fit(gaus,"Q","",xmed,xhigh);
     gaus->GetParameters(&pars[3]);
     //for (size_t j=0;j<6;j++)
     //cout << "par[" << j << "] = " << pars[j] << endl;
     fun->SetParameters(pars);
     fun->SetLineColor(kBlue);
     htemp->Fit(fun,"QL","",xlow,xhigh);
     htemp->DrawCopy();
     cout << "Slice: " << i << " gamma:   " << fun->GetParameter(1) 
	  << " +/-" << fun->GetParameter(2) << 
       endl; 
     cout << "         neutron: " << fun->GetParameter(4) 
	  << " +/-" << fun->GetParameter(5) << endl;
     double fom = (fun->GetParameter(4)-fun->GetParameter(1))
       /TMath::Sqrt(fun->GetParameter(2)*fun->GetParameter(2)+fun->GetParameter(5)*fun->GetParameter(5));
     FoMnear->SetPoint(i-1,i,fom);
     FoMnear->SetPointError(i-1,0.,0.);
   }

  TCanvas* c5 = new TCanvas();
  c5->cd();
  FoMnear->Draw("AZP");

  */

  ifile->Close();
  

  return 0;
}
