#include<iostream>
#include <cstring>

#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TVirtualFitter.h>
#include <TLegend.h>
#include <TText.h>


using namespace std;

void silicalibrator(string const& fname, int nbin = 500, double xmin = 1950, double xmax = 2300, double ref_peak = 2135, double sigma_peak = 4.58, int nref = 4) {

//void silicalibrator(string const& fname, int nbin = 200, double xmin = 2150, double xmax = 2550, double ref_peak = 2205, double sigma_peak = 6.4, int nref = 1) {
   TFile* f = new TFile(fname.c_str(), "read");
   TTree* reco = (TTree*) f->Get("reco");
   
   TCanvas *c = new TCanvas("c", "c", 1);
   TH1F *h = new TH1F("h","Alpha Spectrum; Amplitude [a.u.]; counts", nbin, xmin, xmax);

   reco->Project("h", "ymax[11]-baseline_mean[11]"); // E 
   //reco->Project("h","-ymin[13]+baseline_mean[13]"); // E monitor

   h->Draw();

   double peak[] = {5144.3, 5156.59, 5388.2, 5442.8, 5485.56, 5762.7, 5804.82};
   int npeak = sizeof(peak)/sizeof(peak[0]); 

   string fun = "gaus(0)";
   for (int i = 0; i < npeak; i++) fun += Form(" + gaus(%d)", 3*i);
   
   TF1 *fit = new TF1("fit", fun.c_str(), xmin, xmax);
   fit->SetNpx(500);
   for (int i = 0; i < npeak; i++) {
      fit->FixParameter(i*3 + 1, peak[i]/peak[nref]*ref_peak);
      fit->FixParameter(i*3 + 2, sigma_peak);
      fit->SetParameter(i*3 + 0, 1);
   }
  
   h->Fit("fit", "EMR");    

   for (int i = 0; i < npeak; i++) fit->ReleaseParameter(i*3 + 2);
    

   h->Fit("fit", "EMR");

   for (int i = 0; i < npeak; i++) fit->ReleaseParameter(i*3 + 1);
    
   h->Fit("fit", "EMR");    

   TCanvas *c1 = new TCanvas("c1", "c1", 1);

   TGraphErrors *g = new TGraphErrors();
   g->SetTitle("SiE Detector Calibration; Amplitude [a.u.];  E [keV];");

   int k = 0;
   for (int i = 1; i < npeak; i++) {
      g->SetPoint(k,  peak[i], fit->GetParameter(i*3 + 1));
      g->SetPointError(k, 0,  fit->GetParError(i*3 + 1));
      k++; 
   }

   g->SetPoint(k, 28000 - 3110, 9996);
   g->SetPointError(k , 0, 0.7);
   k++;

   g->SetPoint(k, 18000 - 4418, 5754);
   g->SetPointError(k , 0, 0.5);
   k++;

   g->SetPoint(k , 26000 - 3309, 9168);
   g->SetPointError(k , 0, 0.4);
   k++;

   g->SetMarkerStyle(21);
   g->SetMarkerColor(kBlue);

   g->Draw("ap");

   TF1 *l = new TF1("l","[0]*x+[1]", 0, 28000);
   l->SetParNames("m", "q");
   l->SetLineColor(kBlack);
   g->Fit("l");
 
   
   TGraphErrors *g1 = new TGraphErrors(g->GetN());
   TGraphErrors *g2 = new TGraphErrors(g->GetN());

   for (int i = 0; i < 28000; i++) {
      g1->SetPoint(i, i, 0);
      g2->SetPoint(i, i, 0);
   }

   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(g1, 0.68);
   (TVirtualFitter::GetFitter())->GetConfidenceIntervals(g2, 0.95);
   
   g1->SetLineWidth(2);
   g1->SetLineColor(kGreen);
   g1->SetMarkerColor(kGreen);
   g1->SetFillColor(kGreen);


   g2->SetLineWidth(2);
   g2->SetLineColor(kYellow);
   g2->SetMarkerColor(kYellow);
   g2->SetFillColor(kYellow);


   TGraphErrors *ga = new TGraphErrors(6);
   for (int i = 0; i < 6; i++) ga->SetPoint(i, g->GetX()[i], g->GetY()[i]);
   ga->SetMarkerStyle(21);
   ga->SetMarkerColor(kRed);

      
   g2->SetTitle("Calibration Plot; E [keV]; Amplitude [ADC]");
   g2->Draw("APL");
   g1->Draw("LPsame");

   g->Draw("Psame");
   ga->Draw("Psame");

   l->Draw("same"); 

   TLegend *leg = new TLegend(0.1, 0.9, 0.3,0.7);
   leg->AddEntry(g, "Au elastic peaks", "p");
   leg->AddEntry(ga, "Alpha calibations", "p");
   leg->AddEntry(l, "Calibration curves", "l");
   leg->AddEntry(g1 , "Confidence Band 68%", "f");
   leg->AddEntry(g2, "Confidence Band 95%", "f");

   leg->Draw("same");

   double m = l->GetParameter(0);
   double q = l->GetParameter(1);

   TText *t = new TText(.5,.5, Form("E [keV] = (Amp [ADC] - (%2.2f))/%2.2f", q, m));
   t->SetTextColor(kRed);
   t->SetTextFont(43);
   t->SetTextSize(40);
   t->Draw("same");

} 
