#include <iostream>

#include <TTree.h>
#include <TF1.h>
#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>

double ffit(double *x, double *p)
{
	double N = p[0];
        double mu = p[1];
        double sL = p[2];
        double sR = p[3];

        double sAvg = (sL + sR) / 2.; 
                                  
        if(x[0] <= mu) return 2* N/(TMath::Sqrt(2*TMath::Pi())*sAvg)*TMath::Exp(-0.5*TMath::Power((x[0] - mu)/sL, 2));
        else return 2*N/(TMath::Sqrt(2*TMath::Pi())*sAvg)*TMath::Exp(-0.5*TMath::Power((x[0] - mu)/sR, 2));
}

using namespace std;

void asyfitter(string const& fname, int nbin = 500, double xmin = 11700, double xmax = 12500) {

   TFile* fin = new TFile(fname.c_str(), "read");
   TTree* reco = (TTree*) fin->Get("reco");
   
   TCanvas *c = new TCanvas("c", "c", 1);

   TH1F *h = new TH1F("h","Elastic Peak; Amplitude [a.u.]; counts", nbin, xmin, xmax);

   reco->Project("h","-ymin[13]+baseline_mean[13]"); 

   TF1 *g = new TF1("g","gaus", xmin, xmax);
  
   h->Fit("g","0");

   TF1 *f = new TF1("f", ffit, xmin, xmax, 4);
   f->SetParNames("N", "#mu", "#sigma_{L}","#sigma_{R}");
   f->SetNpx(500);

   f->SetParameters(g->GetParameter(0), g->GetParameter(1), g->GetParameter(2), g->GetParameter(2));
   h->Fit("f");

}

