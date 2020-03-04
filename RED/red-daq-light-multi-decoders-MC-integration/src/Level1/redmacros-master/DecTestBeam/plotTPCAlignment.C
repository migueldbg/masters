#include "TH2Poly.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCut.h"
#include "TCutG.h"
#include "TLatex.h"
#include <iostream>
#include "red-daq/EvRec0.hh"
#include "red-daq/RDCluster.hh"

/*
  This macro reads the output reco files from the LNS Test Beam 
  of Dec 2019. It is assumed that the "full" configuration is used 
  (i.e. trigger on Si only and TPC in slave). 

  The argument of the macro is the run number. 

  The macro select coincidence events in the TPC requesting:
  - Be events in the Si
  - valid cluster found within 100 ns with respect to DeltaE-fast
  - PSD compatible with neutron (fprompt > 0.4)

  The Be7 is selected using a graphical cut read from BeCut_Dec19.root. 
  Different cuts are used for different runs.
  If the file is not available, a square box is used in DeltaE/E.

  The macro produces: 2D histo with the full bananas (DeltaE/E) 
  superimposed with the events in coincidence.

*/



int plotTPCAlignment(Int_t runNb, Bool_t useLiveTime = false)
{ 
  gStyle->SetOptStat(0);
  gROOT->SetStyle("Plain");
  
 			
  TString filename;
  filename.Form("run_%d.root",runNb);

  TH2D* fBanana = new TH2D(Form("h0_%d",runNb),
			   Form("Full banana, run %d",runNb),
			   400,1500,7000,400,200,4000);
  TH2D* fCBanana = new TH2D(Form("h1_%d",runNb),
			    "Coincident banana",400,1500,7000,400,200,4000);
  fCBanana->SetMarkerStyle(1);
  fCBanana->SetMarkerColor(kBlack);

  TH1D* h90 = new TH1D("h90","fprompt for Be events",100,0.,1.);
  TH1D* dT = new TH1D("dT","DeltaT",400,-200,200);
  TH1D* hS1 = new TH1D("hS1","S1 for coincident events",200,0,1000); 

  //
  TString cutname = "be";

  //Open cut file
  TFile* fcut = new TFile("BeCut_Dec19-v2.root");
  //TFile* fcut = new TFile("Cuts_Dec19.root");
  TCutG* cut = nullptr;
  if (fcut->IsOpen())    
    { 
      cut = static_cast<TCutG*>(fcut->Get(cutname));
      fcut->Close();
    }


  //Open data file
  TFile* f = new TFile(filename);
  if (!(f->IsOpen()))
    {
      cout << "could not open file: " << filename << endl;
      return 0;
    }

  vector<TString> *chanIDs = new vector<TString>;
  //Read metadata
  TTree* metaevent = (TTree*) f->Get("metaevent");
  metaevent->SetBranchAddress("chanID",&chanIDs);
  metaevent->GetEntry(0);


  TTree* reco = (TTree*) f->Get("reco");
  EvRec0* evReco = new EvRec0();
  reco->SetBranchAddress("recoevent",&evReco);
  //cout << reco->GetEntries() << " entries" << endl;
  
  Int_t counter = 0;
  Int_t counterBck = 0;
  Int_t becounts = 0;

  //Event loop
  for (Int_t iloop=0;iloop<reco->GetEntries();iloop++)
    {            
      //if (!(iloop%1000000))
      //	cout << "Processing event #" << iloop << "/" << 
      //	  reco->GetEntries() << endl; 
      reco->GetEntry(iloop);
      vector<double> f90 = evReco->GetF90();
      vector<double> start_time = evReco->GetStartTime();
      vector<double> charge = evReco->GetCharge();
      vector<double> ymin = evReco->GetYmin();
      vector<double> basemean = evReco->GetBaseMean();
      vector<RDCluster*> clusters = evReco->GetClusters();
      Double_t totCharge = evReco->GetChargeTot();
      int nclu =  evReco->GetNClusters();
      Double_t f90_fixed = evReco->GetF90Tot();
      

      //This is the info of the Si detectors
      Double_t deltaE = basemean.at(30)-ymin.at(30);
      Double_t E = basemean.at(31)-ymin.at(31);
      fBanana->Fill(E,deltaE);
      

      //Apply cuts
      Bool_t isBe = false;
      if (cut) //The cut is found
	isBe = cut->IsInside(E,deltaE);
      else
	{	
	  if (deltaE >1300 && deltaE<1600 && E>3500 && E<4400)
	    isBe = true;
	}

      if (!isBe)
	continue;

      becounts++; //counts inside the Be window
      
      //FILL HISTOS
      if (nclu == 0)
	continue;
      

      Double_t fprompt = 0;
      Double_t time = clusters.at(0)->cdf_time;
      Double_t S1 = -1000; 
      Double_t deltaT = (time-start_time.at(30))*2.; //ns
      Int_t startT = clusters.at(0)->start_time;
      
     
      //Check for the presence of S1, in time
      for (size_t icl=0;icl<clusters.size();icl++)
	{
	  Double_t atime = clusters.at(icl)->cdf_time;
	  Double_t adeltaT = (atime-start_time.at(30))*2.; //ns
	 
	  
	  if (TMath::Abs(adeltaT)<200)
	    {
	      fprompt = clusters.at(icl)->f90;
	      if (clusters.at(icl)->f90>0.4) 
		//found a good coupling: timing + neutron-like
		{		  
		  S1 = clusters.at(icl)->charge;
		  time = clusters.at(icl)->cdf_time;		 
		  deltaT = (time-start_time.at(30))*2.; //ns
		  startT = clusters.at(icl)->start_time;
		  //cout << "Trovato: " << deltaT << " " << S1 << " " << fprompt << " " << icl << endl;
		}
	      h90->Fill(fprompt);
	    }
	}
      /*
      if (totCharge > 20)
	S1 = totCharge;
      fprompt = f90_fixed;
      Double_t deltaT = 0; //just for check
      */
      //No good matching found
      if (S1 < 0)
	continue;

     

      dT->Fill(deltaT);
      hS1->Fill(S1);
      
      if (TMath::Abs(deltaT)>50 && 
	  TMath::Abs(deltaT)<200)
	{
	  counterBck++;
	}
      
      if (TMath::Abs(deltaT)<50)
	{
	  //	  cout << "Sono qui: " << fprompt << " " << deltaT << " " << E << " " << deltaE << endl;
	  fCBanana->Fill(E,deltaE);	      
	  counter++;
	}
	     
    }


  //Now subtract background and fill histos
  Double_t counts = counter - counterBck*0.25; //background: 300 ns, signal: 100 ns
 
  Int_t totEvents = reco->GetEntries();
  f->Close();

  cout << "Total events: " << totEvents << endl;
  cout << "Total events in Be band: " << becounts << endl;
  cout << "Fraction of Be events with n-signal in TPC: " << (double)counts/becounts << " (" << 
    counts << "/" << becounts << ")" << endl;
  /*
  TCanvas* c1 = new TCanvas();
  c1->cd();
  c1->SetLogz();
  fBanana->Draw("COLZ");
  fCBanana->Draw("Psame");
  //fCBanana->Draw("P");
  if (cut)
    cut->Draw("same");
  
  TCanvas* c2 = new TCanvas();
  h90->Draw();
  
  
  TCanvas* c3 = new TCanvas();
  dT->Draw();
  
  TCanvas* c4 = new TCanvas();
  hS1->Draw();
  */

  TCanvas* c1 = new TCanvas();
  c1->Divide(2,2);
  c1->cd(1);
  c1->cd(1)->SetLogz();
  fBanana->Draw("COLZ");
  fCBanana->Draw("Psame");
  //fCBanana->Draw("P");
  if (cut)
    cut->Draw("same");
  c1->cd(2);
  h90->Draw();
  
  
  c1->cd(3);
  dT->Draw();
  
  c1->cd(4);
  hS1->Draw();


  //of.cd();
  //fBanana->Write(fBanana->GetName());
  //fCBanana->Write(fCBanana->GetName());
  /*
    for (size_t ich=0;ich<timeDiff.size();ich++)
    timeDiff.at(ich)->Write(timeDiff.at(ich)->GetName());
  */
  //of.Close();

  return 0;
  
}
