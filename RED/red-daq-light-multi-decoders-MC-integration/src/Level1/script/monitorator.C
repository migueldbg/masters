#include<iostream>

#include <iostream>

#include <TTree.h>
#include <TF1.h>
#include <TH1F.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TText.h>

double ffit(double *x, double *p)
{
        double pi = TMath::Pi();
	double N = p[0];
        double mu = p[1];
        double sL = p[2];
        double sR = p[3];
        double m = p[4];
        double q = p[5];

        double sAvg = (sL + sR) / 2.; 
                                  
        if(x[0] <= mu) return N/sqrt(2*pi)/sAvg*exp(-0.5*pow((x[0] - mu)/sL, 2)) + m*x[0] + q;
        else return N/sqrt(2*pi)/sAvg*exp(-0.5*pow((x[0] - mu)/sR, 2)) + m*x[0] + q;
}

using namespace std;

void monitorator(string const& fname, int nbin = 200, double xmin = 9800, double xmax = 10100) {

   TFile* fin = new TFile(fname.c_str(), "read");
   TTree* reco = (TTree*) fin->Get("reco");
   
   TCanvas *c = new TCanvas("c", "c", 1);

   TH1F *h = new TH1F("h","h", nbin, xmin, xmax);

   reco->Project("h","baseline_mean[31] - ymin[31]"); 

   TF1 *g = new TF1("g","gaus", xmin, xmax);
  
   h->Fit("g","0");

   TF1 *f = new TF1("f", ffit, xmin, xmax, 6);
   f->SetParNames("N", "#mu", "#sigma_{L}","#sigma_{R}");
   f->SetNpx(500);

   f->FixParameter(4, 0);
   f->FixParameter(5, 0);


   f->SetParameters(g->GetParameter(0), g->GetParameter(1), g->GetParameter(2), g->GetParameter(2));
   h->Fit("f");

   f->ReleaseParameter(4);
   f->ReleaseParameter(5);
   h->Fit("f");

   double norm = f->GetParameter(0)/h->GetBinWidth(1);
   cout << "Total C12 evt = " << norm << endl;   

   h->SetTitle(Form("Total C12 evt = %2.2f; Amplitude [ADC]; counts", norm));


}
