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
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TChain.h>
#include <TLegend.h>
#include <TCut.h>
#include <TROOT.h>
#include <TAxis.h>
void bananator(int run, int nbin=100, int xmin=0, int xmax=20000, int ymin=0, int ymax=20000){
//run numbers, bin number, xmin, xmax, ymin, ymax

    TFile *f = new TFile(Form("runs/run_%d.root", run), "read");


    TTree *data = (TTree*)f->Get("reco");

//Style
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(2210);
    gStyle->SetOptFit();


////////////// Bananator ////////////////////

    TCanvas *c1 = new TCanvas("c1","c1",10,10,900,600);

    c1->cd();
    c1->SetLogz();

    TString title = "#DeltaE / E Spectrum ";

    TH2F *h = new TH2F("h",title,nbin,xmin,xmax,nbin,ymin,ymax);

    data->Draw("baseline_mean[30] - ymin[30]:baseline_mean[31] - ymin[31]>>h");

    h->GetXaxis()->SetTitle("E [ADC counts]");
    h->GetXaxis()->SetTitleOffset(1.16);
    h->GetYaxis()->SetTitle("#Delta E [ADC counts]");
    h->GetYaxis()->SetTitleOffset(1.36);

    h->DrawCopy("colz");


}
