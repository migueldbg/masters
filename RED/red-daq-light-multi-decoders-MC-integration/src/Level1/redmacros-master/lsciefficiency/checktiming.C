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
  This macro reads two output files produced by efficiency() and 
  draws the (superimposed) tof spectra for each channel. Two kind 
  of zooms are available, depending on the flag isgamma: one for the 
  gamma region [-20,20] ns and one for the neutrons [-20,200] ns.

  In the case of gamma spectra, it fits the position of the main gamma 
  peak with a Gaussian and reports the distance (=t*c) and the width of the 
  peak (rms). 
 */

void checktiming(TString file1,TString file2,bool isgamma=true)
{
  /*
  TFile f("run48out.root");
  TFile g("run58out.root");
  */
  /*
  TFile f("runs75-79.out.root");
  TFile g("runs74.out.root");
  */
  vector<TH1D*> h1;
  vector<TH1D*> h2;

  TFile f(file1);
  TFile g(file2); 

  Double_t xmax = (isgamma) ? 20 : 200;
  for (size_t i=16;i<24;i++)
    {
      if (isgamma)
	{
	  h1.push_back((TH1D*) f.Get(Form("tofg%d",i)));
	  h2.push_back((TH1D*) g.Get(Form("tofg%d",i)));
	}
      else
	{
	  h1.push_back((TH1D*) f.Get(Form("tof%d",i)));
	  h2.push_back((TH1D*) g.Get(Form("tof%d",i)));
	}
    }

  for (size_t i=0;i<h1.size();i++)
    {
      if (h1.at(i)->GetEntries() && h2.at(i)->GetEntries())
	{
	  Int_t b1 = h1.at(i)->FindBin(-20);
	  Int_t b2 = h1.at(i)->FindBin(xmax);
	  Double_t x1 = h1.at(i)->Integral(b1,b2);
	  Double_t x2 = h2.at(i)->Integral(b1,b2);      
	  h1.at(i)->Sumw2();
	  h2.at(i)->Sumw2();
	  h1.at(i)->Scale(x2/x1);
	  //h1.at(i)->Rebin(2);
	  //h2.at(i)->Rebin(2);
	}
    }
  //cout << h1.size() << " " << h2.size() << endl;
  TCanvas* c1 = new TCanvas(); 
  c1->Divide(3,3);
  //c1->Divide(2,2);
  for (size_t i=0;i<h1.size();i++)
    {
      if (h1.at(i)->GetEntries())
	{		
	  TVirtualPad* p1 = c1->cd(i+1);
	  p1->SetLogy();
	  h1.at(i)->GetXaxis()->SetRangeUser(-20,xmax);
	  if (isgamma)
	    {
	      h1.at(i)->Fit("gaus","Q","",1,6);
	      TF1* gg = (TF1*) h1.at(i)->GetFunction("gaus");
	      cout << "ch" << i << "(" << file1 << ") --> (" << 
		gg->GetParameter(1)*30. << " +/- " << gg->GetParError(1)*30 << ") cm, (" <<
		gg->GetParameter(2)<< " +/- " << gg->GetParError(2) << ") ns" << endl;
	      if (h2.at(i)->GetEntries())
		{
		  h2.at(i)->Fit("gaus","Q","",1,6);
		  TF1* gg = (TF1*) h2.at(i)->GetFunction("gaus");
		  cout << "ch" << i << "(" << file2 << ") --> (" << 
		    gg->GetParameter(1)*30. << " +/- " << gg->GetParError(1)*30. << ") cm, (" <<
		    gg->GetParameter(2)<< " +/- " << gg->GetParError(2) << ") ns" << endl;		  
		}
	      
	    }
	  h1.at(i)->DrawCopy("HIST");
	  if (h2.at(i)->GetEntries())
	    {
	      h2.at(i)->SetLineColor(kRed);
	      h2.at(i)->DrawCopy("HISTsameE");
	    }
	}            
    }
}
