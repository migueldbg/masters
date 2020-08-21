/////////////////////////////////
//Created by Simone on 2018_07_05
// Rev. by Simone on 2019_01_21
/////////////////////////////////
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

void serrator2(int run, int nbin=300, int nmin=-1500, int nmax=15000, int npeaks=5){

    TFile *f = new TFile(Form("/storage/DATA-02/darkside/red/reco/rm3reco/naples/laser/run_%d.root", run), "read"); // change this path in order to make LNS SERs
    TTree *data = (TTree*)f->Get("reco");
    
    gROOT->SetStyle("Plain");

    
////////////// Extraction of the SERs ////////////////////

    TCanvas *c1 = new TCanvas("c1","c1",10,10,1800,1200);
    
    //c1->Divide(8,4); //config. LNS
    c1->Divide(7,4); //config. Napoli Dic. 2018
    
    TCanvas *c2 = new TCanvas("c2","c2",10,10,1800,1200);
    
    //c2->Divide(8,4); //config. LNS
    c2->Divide(7,4); //config. Naples Dec. 2018
    
    TH1F *h1[32];
    TF1 *sum;
    for (int i = 0; i < 29; i++) {
        c1->cd(i + 1);/*->SetLogy();*/
        h1[i] = new TH1F(Form("h1%d", i),Form("h1%d", i),nbin,nmin,nmax);
        //data->Draw(Form("charge[%d]>>h1%d", i+16, i)); // config. LNS
        data->Draw(Form("charge[%d]>>h1%d", i, i)); // config. Naples December 2018
        //data->Draw(Form("baseline_mean[%d]-ymin[%d]>>h1%d", i, i, i)); // SER in amplitude

        h1[i]->SetTitle("Single Electron Response");
        h1[i]->GetXaxis()->SetTitle("Charge (ADC*Sample)");
        //h1[i]->GetXaxis()->SetTitle("Amplitude (ADC*Sample)");
        h1[i]->GetXaxis()->SetTitleOffset(1.16);
        h1[i]->GetYaxis()->SetTitle("Counts");
        h1[i]->GetYaxis()->SetTitleOffset(1.16);
        h1[i]->SetLineWidth(1);
        h1[i]->SetLineColor(kBlue);

    }

    ofstream ofile(Form("ser_%d.cfg", run));
    ofstream ofile2(Form("peak_dist_%d.cfg", run));
    //for (int i = 0; i < 48; i++) // config. LNS
    for (int i = 0; i < 29; i++) // config. Naples December 2018
    {

      if (h1[i]->GetEntries()==0) {
            continue;
        }
        
        c1->cd(i + 1);/*->SetLogy();*/
        h1[i]->Smooth(2);
        h1[i]->Rebin(2);

/////////////// FIT & DRAW /////////////////////////

    //Use TSpectrum to find the peak candidates
    int npeak;
    TSpectrum *s = new TSpectrum();
    int nfound = s->Search(h1[i],4,"goff",0.0005); //(histo,sigma,opt,threshold)
        if (nfound<=npeaks) {
            npeak=nfound;
        } else {
            npeak=npeaks;
        }
 
    double *xus = (double *)s->GetPositionX(); //array with X-positions of the centroids found by TSpectrum
    double *y = (double *)s->GetPositionY();
 
    vector<double> x(xus, xus + npeak);
    sort(x.begin(), x.end());

    float dmu = 0;
    for (int j = 0; j < npeak - 1; j++) dmu += x[j + 1] - x[j];
    dmu = dmu/npeak;
    
    //cout << "dmu: " << dmu << "\n" << endl;

    printf("Found %d candidate peaks to fit\n",nfound);
    printf("Found %d useful peaks to fit\n",npeak);
    printf("Now fitting\n");
    
    //Loop on all found peaks.

     TF1 *g[npeak];
     for (int p=0;p<npeak;p++) {
        g[p] = new TF1("gaus","gaus",x[p] -dmu/3, x[p] + dmu/3);
 
        g[p]->SetLineWidth(2);
        g[p]->SetLineColor(kRed);
        h1[i]->Fit(g[p],"R+Q+0");
    }

    string sgaus = "gaus(0) ";
        for (int ss = 1; ss < npeak; ss++) sgaus += Form("+ gaus(%d) ", 3*ss);
        
    sum = new TF1("mysum",sgaus.c_str(),x[0] - dmu/2, x[npeak - 1] + dmu/2);
    sum->SetNpx(1000);
    for (int k=0;k<3*npeak;k++){
       sum->SetParameter(k,g[(k-k%3)/3]->GetParameter(k%3));
       if(!(k-1)%3) sum->FixParameter(k,g[(k-k%3)/3]->GetParameter(k%3));
       if(!(k-1)%3) sum->SetParLimits(k, sum->GetParameter(k) - dmu/3,sum->GetParameter(k) + dmu/3);
      }

    h1[i]->Fit(sum,"R+Q");

    for (int k=0;k<3*npeak;k++){
       sum->SetParameter(k,g[(k-k%3)/3]->GetParameter(k%3));
       if(!(k-1)%3) sum->ReleaseParameter(k);
    }

    sum->SetLineWidth(2);
    sum->SetLineColor(kBlack);
    h1[i]->Fit(sum,"R+Q");


    for (int ij=0;ij<npeak-1;ij++)
        {
        //Distance between one mean and the following
        Double_t ref = sum->GetParameter(4+ij*3) - sum->GetParameter(1+ij*3);
        //cout << "Channel: " << i << ". Distance from peak " << ij+1 << " and " << ij << " --> " <<
        //ref << endl;
            if(ij+1<npeak){
        
              ofile2 << "" << i << " Distance from peak " << ij+1 << " and " << ij << " --> " << ref << endl;
        
            }
            
        }
        
    h1[i]->DrawCopy();
        
    c2->cd(i + 1);
         
    TGraphErrors *gr[32];
    gr[i] = new TGraphErrors();
    gr[i]->SetName(Form("gr%d", i));
    gr[i]->SetLineColor(kRed);
    gr[i]->SetLineWidth(1);
    gr[i]->SetMarkerSize(2);
    gr[i]->SetMarkerColor(kBlue);
    gr[i]->SetMarkerStyle(21);
    gr[i]->SetTitle("Calibration curve");
    gr[i]->GetXaxis()->SetTitle("PE");
    gr[i]->GetYaxis()->SetTitle("Charge (ADC*Sample)");
    //gr[i]->GetYaxis()->SetTitle("Amplitude (ADC*Sample)");
    gr[i]->GetYaxis()->SetTitleOffset(1.16);
         
    for (int j = 0; j < npeak; j++) {
    gr[i]->SetPoint(j, j, sum->GetParameter(1 + j*3));
    gr[i]->SetPointError(j, 0, sum->GetParError(2 + j*3));
    //cout <<  i << " " << j << " " << sum->GetParameter(1 + j*3) << endl;
    }
    TF1 *myline = new TF1("myline","pol1(0)");
    gr[i]->Fit(myline,"Q");
    gStyle->SetOptFit(111);
    gr[i]->Draw("AP");
        
    ofile << myline->GetParameter(1) << endl; //write only the SER values

    }
 

    ofile.close();
    ofile2.close();

c1->cd();
c1->SaveAs(Form("ser_run_%d.png", run));
    
c2->cd();
c2->SaveAs(Form("plot/naples/ser_cal_run_%d.png", run));
    
}
