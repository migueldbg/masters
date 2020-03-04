#include <iostream>
#include <cstdlib>
#include <cstring>
#include <algorithm> 
#include <vector> 
#include <sstream> 

#include <TFile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

#include <TStyle.h>
#include <TCut.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TMath.h>
#include <TStyle.h>
#include "red-daq/EvRec0.hh"

//This is a square convolved with a Gaussian
double FuncDrift(double *v, double *par)
{
  double t=v[0];
  double width=par[0];
  double sigma=par[1];
  
  double A[2]={ t/(TMath::Sqrt2()*sigma), (t-width)/(TMath::Sqrt2()*sigma) };
  double val=  (TMath::Erf(A[0])- TMath::Erf(A[1]) )/2.0*width;
  return par[2]*val;
}

/*
  Luciano's version of the macro
*/

using namespace std;

int Kr(Int_t runsignal, Int_t runbck, bool isKr=true, Double_t xlow=430,
       Double_t xhigh=530)
{
  gROOT->SetStyle("Plain");
  gStyle->SetPaintTextFormat("2.2f");
  cout << "Warning: specularly swapping pos_y" << endl;
    
  //Signal file
  TFile *f = new TFile(Form("run_%d.root",runsignal), "read"); //sorgente
  TTree *sigT = (TTree*)f->Get("reco"); //sorgente

  TTree* bckT = nullptr;

  //Background file (Kr only)
  if (isKr)
    {
      TFile *g = new TFile(Form("run_%d.root",runbck), "read"); //bkg
      if (g)
	bckT = (TTree*)g->Get("reco"); //bkg
    }
  if (!bckT)
    cout << "Using no background subtraction " << endl;

  Int_t nbins= 150;
  //Book histograms
  //Energy spectra
  TH1D *hsg = new TH1D("hsg","S1 spectrum, signal",nbins,0,4000); //sorgente
  TH1D *hbg = new TH1D("hbg","S1 spectrum, background",nbins,0,4000); //bkg
  TH1D *hsubt =  new TH1D("hsubt","S1 spectrum, background-subtracted",
			  nbins,0,4000); //bkg

  TH2D *xysource = new TH2D("xysource","XY distribution, source",4,0,5,6,0,5);
  TH2D *xybck = new TH2D("xybck","XY distribution, background",4,0,5,6,0,5);

  Int_t zbins = (isKr) ? 40 : 200;
  TH1D *zsource = new TH1D("zsource","Tdrift, source",zbins,0,40000/500);
  TH1D *zbck = new TH1D("zbck","Tdrift, background",zbins,0,40000/500);

  Double_t s2s1max=100;
  TH1D *s2s1source = new TH1D("s1s2source","S2/S1, source",zbins,0,s2s1max);
  TH1D *s2s1bck = new TH1D("s1s2bck","S2/S1, background",zbins,0,s2s1max);

  //Fill everything
  //Signal
  EvRec0* evReco = new EvRec0();
  sigT->SetBranchAddress("recoevent",&evReco);
  cout << "Signal entries: " << sigT->GetEntries() << endl;
  //Loop!
  for (int i=0;i<sigT->GetEntries();i++)
    {      
      sigT->GetEntry(i); 
      vector<RDCluster *> clusters = evReco->GetClusters();
      int nclusters = evReco->GetNClusters();
      
      if (nclusters==2) //Select only events with two clusters
	{      
	  if (clusters.at(0)->f90 > 0.2 && clusters.at(1)->f90<0.2 &&
	      clusters.at(0)->rep == 1 && clusters.at(1)->rep == 1)
	    {
	      Double_t S1 = clusters.at(0)->charge; //charge;
	      Double_t S2 = clusters.at(1)->charge; //tot_charge_top;
	      hsg->Fill(S1);
	      hsubt->Fill(S1);
	      if (S1 > xlow && S1 < xhigh) //select peak
		{
		  Double_t tdrift = 
		    (clusters.at(1)->cdf_time-clusters.at(0)->cdf_time)*0.002;
		  xysource->Fill(clusters.at(1)->pos_x,
				 //Use specular reflection
				 5.-clusters.at(1)->pos_y);
		  zsource->Fill(tdrift);
		  s2s1source->Fill(S2/S1);
		}
	    }
	}
    }
  cout << "Signal file read " << endl;
  if (bckT)
    {
      bckT->SetBranchAddress("recoevent",&evReco);
      cout << "Background entries: " << bckT->GetEntries() << endl;
      //Loop!
      for (int i=0;i<bckT->GetEntries();i++)
	{      
	  bckT->GetEntry(i); 
	  vector<RDCluster *> clusters = evReco->GetClusters();
	  int nclusters = evReco->GetNClusters();
      
	  if (nclusters==2) //Select only events with two clusters
	    {      
	      if (clusters.at(0)->f90 > 0.2 && clusters.at(1)->f90<0.2 &&
		  clusters.at(0)->rep == 1 && clusters.at(1)->rep == 1)
		{
		  Double_t S1 = clusters.at(0)->charge;
		  Double_t S2 = clusters.at(1)->charge;
		  hbg->Fill(S1);
		  if (S1 > xlow && S1 < xhigh) //select peak
		    {
		      Double_t tdrift = 
			(clusters.at(1)->cdf_time-clusters.at(0)->cdf_time)*0.002;
		      xybck->Fill(clusters.at(1)->pos_x,clusters.at(1)->pos_y);
		      zbck->Fill(tdrift);
		      s2s1bck->Fill(S2/S1);
		    }		  
		}
	    }
	}
      cout << "Background file read " << endl;
    }
  
  //START PLOTTING

  // Energy spectra
  hsg->SetLineColor(kBlue);
  hbg->SetLineColor(kViolet);
  
  hsubt->GetXaxis()->SetTitle("Charge (PE)");
  hsubt->GetXaxis()->SetTitleOffset(1.26);
  hsubt->GetYaxis()->SetTitle("Counts");
  hsubt->GetYaxis()->SetTitleOffset(1.26);

  hsubt->Sumw2();
  hsg->Sumw2();
  hbg->Sumw2();
  
  Int_t xmin = hsg->GetXaxis()->FindBin(xhigh+100);
  Int_t xmax = hsg->GetXaxis()->FindBin(3000);
  double integral = hsg->Integral(xmin,xmax); //sorgente
  double integral1 = hbg->Integral(xmin,xmax); //bkg

  cout << "Tail: background has " << integral1 
       << " events and the source " << integral << endl;
  

  double norm = integral1/integral;
  cout << "Background normalization factor in the tail: " << 1./norm << endl;
  /*
  xmin = hsg->GetXaxis()->FindBin(xlow);
  xmax = hsg->GetXaxis()->FindBin(xhigh);
  integral = hsg->Integral(xmin,xmax); //sorgente
  integral1 = hbg->Integral(xmin,xmax); //bkg
  cout << "Peak: background has " << integral1 << " events and the source " << integral << endl;
  */
  TCanvas* c1 = new TCanvas();
  c1->cd();
  hbg->Scale(1./norm);
  hsg->Draw("hist E");
  hbg->Draw("same E");
  
  TPad *p1 = new TPad("p1","p1",.45,.3,.87,.8);
  p1->Draw();
  p1->cd();
  hsubt->GetXaxis()->SetRangeUser(0,xhigh*2);
  hsubt->SetBit(TH1::kNoTitle);
  //Scale background to source area
  hsubt->Add(hbg,-1);
  hsubt->DrawCopy("hist");

  Int_t xmin1 = hsubt->GetXaxis()->FindBin(xlow);
  Int_t xmax1 = hsubt->GetXaxis()->FindBin(xhigh);
  Double_t integral_cout = hsubt->Integral(xmin1,xmax1);   
  cout << "Events in the " << (isKr ? "Kr" : "Am") << " peak: " << integral_cout << endl;
    

  //XY distribution
  TCanvas* c2 = new TCanvas();
  c2->cd();
  cout << "Entries: background: " << xybck->GetEntries() << ", source " << xysource->GetEntries() << endl;
  if (xybck->GetEntries())    
    {
      c2->Divide(2,1);
      c2->cd(1);
    }
  else
    c2->cd();
  xysource->DrawCopy("colz text");
  if (xybck->GetEntries())
    {
      c2->cd(2);
      xybck->DrawCopy("colz text");    
      //Scale background to source area, using the same coefficient as above
      TCanvas* c3 = new TCanvas();
      c3->cd();
      xybck->Scale(1./norm);
      xysource->Add(xybck,-1);     
    }
  Double_t totalEntries = xysource->GetEntries();
  cout << "After subtraction: " << totalEntries << " entries" << endl;
  //Double_t averagePerChannel = totalEntries/1.; 
  //xysource->Scale(1./averagePerChannel);
  Double_t maxEntries = xysource->GetMaximum();
  cout << "Maximum entries in XY: " << maxEntries << endl;
  xysource->Scale(1./maxEntries);
  xysource->Draw("colz text");
  //DONE

  //NOw TDRIFT
  TCanvas* c4 = new TCanvas();
  c4->cd();
  //zsource->SetLineColor(kBlue);
  zbck->SetLineColor(kViolet);
  zbck->Scale(1./norm);
  zsource->SetMarkerStyle(21);
  zsource->SetMarkerSize(0.8);
  zsource->GetXaxis()->SetTitle("Drift time (#mu s)");
  zsource->DrawCopy("PE");
  zsource->DrawCopy("hist SAME");
  if (zbck->GetEntries())
    {
      zbck->DrawCopy("same PE");
      TCanvas* c5 = new TCanvas();
      c5->cd();
      zsource->Add(zbck,-1);
      zsource->DrawCopy("HIST EP");    
    }
  Double_t xxmax=65;
  //Double_t xxmax=130;
  Double_t xxmin = xxmax-10.;
  TF1 *fdrift = new TF1("fdrif", FuncDrift, xxmin, xxmax,3);
  fdrift->SetNpx(10000);
  double parD[3]={50.,10.,1.};
  fdrift->SetParameters(parD);
  zsource->Fit(fdrift,"QL0N","",xxmin,xxmax);
  fdrift->DrawCopy("SAME");
  TLatex l;
  l.DrawLatex(0.5*(xxmin+xxmax),200,Form("t_{drift}^{max}=%5.2lf +/- %5.2lf ",fdrift->GetParameter(0),
			  fdrift->GetParError(0)));


  /*
  TPad *p2 = new TPad("p2","p2",.45,.3,.87,.8);
  p2->Draw();
  p2->cd();
  */
 
  //And finally S2/S1
  TCanvas* c6 = new TCanvas();
  c6->cd();
  //zsource->SetLineColor(kBlue);
  s2s1bck->SetLineColor(kViolet);
  s2s1bck->Scale(1./norm);
  s2s1source->SetMarkerStyle(21);
  s2s1source->SetMarkerSize(0.8);
  s2s1source->GetXaxis()->SetTitle("S2/S1");
  s2s1source->DrawCopy("PE");
  //s2s1source->DrawCopy("hist");
  if (s2s1bck->GetEntries())
    {
      s2s1bck->DrawCopy("same PE");
      TCanvas* c7 = new TCanvas();
      c7->cd();
      s2s1source->Add(s2s1bck,-1);
      s2s1source->DrawCopy("EP");    
    }
 if (!isKr)
    s2s1source->Fit("gaus","","",6,16);
  /* 200 V/cm */
  /*
  if (!isKr)
    s2s1source->Fit("gaus","","",10,29);
  TPad *p2 = new TPad("p2","p2",.55,.55,.9,.86);
  p2->Draw();
  p2->cd();
  hsg->GetXaxis()->SetRangeUser(0,1000);
  hsg->DrawCopy("hist");
  */
  return 0;
}
