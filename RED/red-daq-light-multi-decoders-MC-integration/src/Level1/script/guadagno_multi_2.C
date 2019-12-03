#include <iostream>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <vector>
#include <sstream>

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

#include "red-daq/EvRec0.hh"
#include "TCanvas.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

void guadagno_multi_2(TString filename, int nbin=300, int nmin=-2000, 
		    int nmax=12000, int npeaks=5)
{   
    TFile *f = new TFile(filename, "read");
    if (!(f->IsOpen()))
    {
      cout << "could not open file: " << filename << endl;
      return;
    }
    TTree *data = (TTree*)f->Get("reco");
    EvRec0* evReco = new EvRec0();
    data->SetBranchAddress("recoevent",&evReco);
    gROOT->SetStyle("Plain");
    
    data->GetEntry(0);
    size_t nchannels = evReco->GetCharge().size();
    
    ////////////// CHARGE ////////////////////
    TCanvas *c1 = new TCanvas();
    
    c1->Divide(4,4);
    //c1->cd();
    //c1->SetLogy();
    
    vector<TH1F*> h1(nchannels,0);
    
    TCanvas *c2 = new TCanvas();
    c2->Divide(4,4);
    
    int npeak;
    TGraphErrors *gr[16];
    TSpectrum *s = new TSpectrum();

    for (size_t i = 0; i < nchannels; i++) 
      {	
	h1.at(i) = new TH1F(Form("h1%d", i),Form("h1%d", i),nbin,nmin,nmax);
	//data->Draw(Form("charge[%d]>>h1%d", i, i));
	
	h1[i]->SetTitle("Single Electron Response");
	h1[i]->GetXaxis()->SetTitle("Charge [ADC counts]");
        h1[i]->GetXaxis()->SetTitleOffset(1.16);
        h1[i]->GetYaxis()->SetTitle("Counts");
        h1[i]->GetYaxis()->SetTitleOffset(1.16);
        h1[i]->SetLineWidth(1);
        h1[i]->SetLineColor(kBlue);
      }
    
    for (int iev=0;iev<data->GetEntries();iev++)
      {
	data->GetEntry(iev);
	for (size_t ii=0;ii<evReco->GetCharge().size();ii++)
	  {
	    //Here one can do any cut by using the other parameters
	    //
	    h1.at(ii)->Fill(evReco->GetCharge().at(ii));
	  }
      }

    //
    ofstream ofile("ser.out");
    for (size_t i = 0; i < nchannels; i++) 
      {
	c1->cd(i + 1)->SetLogy();
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
    
	//printf("Found %d candidate peaks to fit\n",nfound);
	//printf("Found %d useful peaks to fit\n",npeak);
	//printf("Now fitting\n");
    
	//Loop on all found peaks.
    
	TF1 *g[npeak];
	for (int p=0;p<npeak;p++) {
	  g[p] = new TF1("gaus","gaus",x[p] -dmu/2, x[p] + dmu/2);
	  
	  g[p]->SetLineWidth(2);
	  g[p]->SetLineColor(kRed);
	  h1[i]->Fit(g[p],"R+Q");
	}
    
	TF1 *sum = new TF1("mysum",
			   "gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)",x[0] - dmu/2, x[npeak - 1] + dmu/2);
	sum->SetNpx(1000);
	for (int j=0;j<3*npeak;j++){
	  sum->SetParameter(j,g[(j-j%3)/3]->GetParameter(j%3));
	}
    
	sum->SetLineWidth(2);
	sum->SetLineColor(kBlack);
	h1[i]->Fit(sum,"R+Q");

	for (int ij=0;ij<4;ij++)
	  {
	    //Distance between one mean and the following
	    Double_t ref = sum->GetParameter(4+ij*3) - sum->GetParameter(1+ij*3);	 
	    cout << "Canale: " << i << ". Distanza fra picco " << ij+1 << " and " << ij << " --> " << 
	      ref << endl;
	  }	

	//c1->Update();
	h1[i]->DrawCopy();
      
	
	c2->cd(i + 1);
    
	gr[i] = new TGraphErrors();
	gr[i]->SetName(Form("gr%d", i));
	gr[i]->SetLineColor(kRed);
	gr[i]->SetLineWidth(1);
	//gr[i]->SetMarkerSize(2);
	gr[i]->SetMarkerColor(kBlue);
	gr[i]->SetMarkerStyle(21);
	gr[i]->SetTitle("Calibration curve");
	gr[i]->GetXaxis()->SetTitle("PE");
	gr[i]->GetYaxis()->SetTitle("Charge [ADC counts]");
	gr[i]->GetYaxis()->SetTitleOffset(1.16);
	
	for (int j = 0; j < npeak; j++) {
	  gr[i]->SetPoint(j, j, sum->GetParameter(1 + j*3));
	  gr[i]->SetPointError(j, 0, sum->GetParError(2 + j*3));
	  //cout <<  i << " " << j << " " << sum->GetParameter(1 + j*3) << endl;
	}
	TF1 *myline = new TF1("myline",
			      "pol1(0)");
	gr[i]->Fit(myline,"Q");
	gr[i]->Draw("AP");	
	ofile << i << " " << myline->GetParameter(0) << " " << 
	  myline->GetParameter(1) << endl;
	cout << i << " " << myline->GetParameter(0) << " " << 
	  myline->GetParameter(1) << endl;
      }

    

}
