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

// This macro calculates the SERs from laser runs.
// It produces plots, an external file with the SER values per channels
// and it fills with -1 those channels whose are not from TPC.
// Plots and external file are saved on the same folder of the running macro.
//
// Usage: root -l
//        root[0] .L serrator.C
//        root[1] serrator("path_to_file", run_number, 300, -1500, 25000, 5)
//
// Simone S. on 08/02/2020

void serrator(TString filename, int run, int nbin=300, int nmin=-1500, int nmax=25000, int npeaks=5){

    TFile *f = new TFile(filename, "read");
    TTree *data = (TTree*)f->Get("reco");
    
    gROOT->SetStyle("Plain");
    
// Use only in case of comparison between histos
    //gStyle->SetOptStat("");

    
////////////// CHARGE ////////////////////

    TCanvas *c1 = new TCanvas("c1","c1",10,10,1800,1200);
    
    //c1->Divide(7,4); //config. Napoli Dic. 2018
    c1->Divide(6,4); //config. LNS Feb. 2020
    
    TCanvas *c2 = new TCanvas("c2","c2",10,10,1800,1200);
    
    //c2->Divide(7,4); //config. Napoli Dic. 2018
    c2->Divide(6,4); //config. LNS Feb. 2020

    vector<double> ser;
    
    vector<TString> *chanIDs = new vector<TString>;
    //Read metadata
    TTree* metaevent = (TTree*) f->Get("metaevent");
    metaevent->SetBranchAddress("chanID",&chanIDs);
    metaevent->GetEntry(0);
/*
    vector<TH1F*> h1;
    vector<TGraphErrors*> gr;
    
    for (int i = 0; i < 48; i++) {
        
        h1.push_back(new TH1F(Form("h1%d", i),Form("h1%d", i),nbin,nmin,nmax));
        gr.push_back(new TGraphErrors());
                     
    }
 */
    TH1F *h1[48];
    TF1 *sum;
    
    int counter = 0.;

    for (int i = 0; i < 48; i++) {
        c1->cd(i+1);
        h1[i] = new TH1F(Form("h1%d", i),Form("h1%d", i),nbin,nmin,nmax);
        data->Draw(Form("charge[%d]>>h1%d", i, i)); // LNS solo TPC completa
        //data->Draw(Form("charge[%d]>>h1%d", i, i)); // configurazione Napoli Dicembre 2018
        //data->Draw(Form("baseline_mean[%d]-ymin[%d]>>h1%d", i, i, i)); // per avere le SER in ampiezza
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
    
    for (int i = 0; i < 48; i++)
    {
        //cout << "nch: " << chanIDs->at(i) << endl;

      if (chanIDs->at(i)=="NULL") {
          //cout << "Sono vuoto" <<  endl;
          ofile << "-1" << endl;
          continue;
      }else{
          counter++;
        //cout << "Sono pieno" <<  endl;
	c1->cd(counter)->SetLogy();
        h1[i]->Smooth(1);

/////////////// FIT & DRAW /////////////////////////

    //Use TSpectrum to find the peak candidates
    int npeak;
    TSpectrum *s = new TSpectrum();
    int nfound = s->Search(h1[i],4,"goff",0.0005); //(nome isto,sigma,opt,soglia)
        if (nfound<=npeaks) {
            npeak=nfound;
        } else {
            npeak=npeaks;
        }
 
    double *xus = (double *)s->GetPositionX(); //array with X-positions of the centroids found by TSpectrum
    double *y = (double *)s->GetPositionY();
 
    vector<double> x(xus, xus + npeak);
    sort(x.begin(), x.end());
/*
    for (int i=0; i<x.size();i++)
        cout << "lista: " << x[i] << endl;
*/
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
        cout << "Canale: " << i << ". Distanza fra picco " << ij+1 << " and " << ij << " --> " <<
        ref << endl;
            
        }	

    h1[i]->DrawCopy();

    c2->cd(counter);
         
    TGraphErrors *gr[48];
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
        
    ser.push_back(myline->GetParameter(1));
      
    cout << "Io sono la SER del canale " << i << ": " << ser.back() << endl;
          
    ofile << ser.back() << endl;
        
    //ofile << myline->GetParameter(1) << endl; // scrive nel file solo il valore della ser senza canale
        
    cout << "Mean_first_peak: " << sum->GetParameter(4) << " " << "Sigma_first_peak: " << sum->GetParameter(5) << endl;
         
    cout << "Res: " << sum->GetParameter(5)/sum->GetParameter(4) << endl;
         
    cout << "" << i << " " << "SNR: " << sum->GetParameter(4)/sum->GetParameter(2) << " +/- " << (sum->GetParError(4)/sum->GetParameter(4))+(sum->GetParError(2)/sum->GetParameter(2)) << endl;
          
      }
    }
/*
    for (int i=0;i<48;i++){
        
        if (i<6 || i>29){ ofile << "-1" << endl;
        }else if (i>5 && i<30) ofile << ser.at(i-6) << endl;
    
    }
 */

    ofile.close();
    
    c1->cd();
    
    c1->SaveAs(Form("ser_run_%d.png", run));
    
    c2->cd();
    
    c2->SaveAs(Form("ser_cal_run_%d.png", run));
    
}
