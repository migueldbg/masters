#include <iostream>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <vector>
#include <sstream>
#include <fstream>

#include <TFile.h>
#include <TF1.h>
#include <TTree.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TChain.h>
#include <TLegend.h>
#include <TCut.h>
#include <TROOT.h>

#include "TCanvas.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

void guadagno_multi(TString filename, int nbin=300, int nmin=-2000, 
		    int nmax=12000, int npeaks=5){

    
    TFile *f = new TFile(filename, "read");
    //TFile *f = new TFile(Form("rootfiles/run_%d_b01.root", run), "read");
    TTree *data = (TTree*)f->Get("reco");

// Cuts
    //TCut charge = Form("charge[%d]>-5e6 && charge[%d]<5e6",nch,nch);
    //TCut f90 = Form("f90[%d]>0 && f90[%d]<0.4",nch,nch);

//    TCut charge = "charge[0]>-5e6 && charge[0]<5e6";
//    TCut f90 = "f90[0]>0 && f90[0]<0.4";
    
    gROOT->SetStyle("Plain");
    
// Use only in case of comparison between histos
    //gStyle->SetOptStat("");

    
////////////// CHARGE ////////////////////

    TCanvas *c1 = new TCanvas();
    
    c1->Divide(5,6);
    //c1->cd();
    //c1->SetLogy();
    
    vector<TH1F*> h1;
    TCanvas *c2 = new TCanvas();
    c2->Divide(5,6);
    
    int npeak;
    TGraphErrors *gr[28];
    TSpectrum *s = new TSpectrum();

    ofstream o("ser_out.cfg");

    for (int i = 0; i < 28; i++) 
      {
	c1->cd(i + 1)->SetLogy();
	h1.push_back(new TH1F(Form("h1%d", i),Form("h1%d", i),nbin,nmin,nmax));
	data->Draw(Form("charge[%d]>>h1%d", i, i));
	
	h1[i]->SetTitle("Single Electron Response");
	h1[i]->GetXaxis()->SetTitle("Charge [ADC counts]");
        h1[i]->GetXaxis()->SetTitleOffset(1.16);
        h1[i]->GetYaxis()->SetTitle("Counts");
        h1[i]->GetYaxis()->SetTitleOffset(1.16);
        h1[i]->SetLineWidth(1);
        h1[i]->SetLineColor(kBlue);
        
        h1[i]->Rebin(2);
        
	/////////////// FIT & DRAW /////////////////////////

	//Use TSpectrum to find the peak candidates
	int nfound = s->Search(h1[i],4,"",0.0005); //(nome isto,sigma,opt,soglia)
	if (nfound<=npeaks) {
	  npeak=nfound;
	} else {
	  npeak=npeaks;
	}
    
	double *xus = s->GetPositionX(); //array with X-positions of the centroids found by TSpectrum
	double *y = s->GetPositionY();
    
	vector<double> x(xus, xus + npeak);
	sort(x.begin(), x.end());
    
	//    for (int i=0; i<x.size();i++)
	//        cout << "lista: " << x[i] << endl;
    
	float dmu = 0;
	for (int j = 0; j < npeak - 1; j++) dmu += x[j + 1] - x[j];
	dmu = dmu/npeak;
    
	//cout << "dmu: " << dmu << "\n" << endl;
    
	printf("Found %d candidate peaks to fit\n",nfound);
	printf("Found %d useful peaks to fit\n",npeak);
	//printf("Now fitting\n");
    
	//Loop on all found peaks.
    
	TF1 *g[npeak];
	for (int p=0;p<npeak;p++) {
	  g[p] = new TF1("gaus","gaus",x[p] -dmu/2, x[p] + dmu/2);
	  
	  g[p]->SetLineWidth(2);
	  g[p]->SetLineColor(kRed);
	  h1[i]->Fit(g[p],"R+");
	}
    
	TF1 *sum = new TF1("mysum",
			   "gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)",x[0] - dmu/2, x[npeak - 1] + dmu/2);
	sum->SetNpx(1000);
	for (int j=0;j<3*npeak;j++){
	  sum->SetParameter(j,g[(j-j%3)/3]->GetParameter(j%3));
	}
    
	sum->SetLineWidth(2);
	sum->SetLineColor(kBlack);
	h1[i]->Fit(sum,"R+");

	//c1->Update();
	h1[i]->DrawCopy();
      
	
	c2->cd(i + 1);
    
	gr[i] = new TGraphErrors();
	gr[i]->SetName(Form("gr%d", i));
	gr[i]->SetLineColor(kRed);
	gr[i]->SetLineWidth(2);
	gr[i]->SetMarkerSize(2);
	gr[i]->SetMarkerColor(kBlue);
	gr[i]->SetMarkerStyle(21);
	gr[i]->SetTitle("Calibration curve");
	gr[i]->GetXaxis()->SetTitle("PE");
	gr[i]->GetYaxis()->SetTitle("Charge [ADC counts]");
	gr[i]->GetYaxis()->SetTitleOffset(1.16);
	
	for (int j = 0; j < npeak; j++) {
	  gr[i]->SetPoint(j, j, sum->GetParameter(1 + j*3));
	  gr[i]->SetPointError(j, 0, sum->GetParError(2 + j*3));
	  cout <<  i << " " << j << " " << sum->GetParameter(1 + j*3) << endl;
	}
	TF1 *myline = new TF1("myline",
			      "pol1(0)");
	gr[i]->Fit(myline);
	gr[i]->Draw("AP");
        o << myline->GetParameter(1) << std::endl; 
      }

   for (int i = 0; i < 4; i++) o << "-1" << std::endl; // WARNING add -1 for the last four channels
   o.close();
}
