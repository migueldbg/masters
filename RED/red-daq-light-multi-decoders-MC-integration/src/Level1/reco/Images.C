#include <TFile.h>
#include <TLine.h>
#include <TString.h>
#include <TTree.h>


TH2F* F90vChargeHistogram( TString file_name ){

  TFile* file = new TFile( file_name );
  TTree* reco; file -> GetObject("reco", reco);

  reco -> Draw("clusters[0].f90:clusters[0].charge >> hist", "clusters[0].f90 > 0. && clusters[0].f90 < 1.0 && clusters[0].charge > 0. && clusters[0].charge < 1000.", "goff");
  TH2F* f90vCharge = (TH2F*) gDirectory -> Get("hist");

  f90vCharge -> SetDirectory(0);
  file -> Close();

  return f90vCharge;
}

void F90vChargeCanvas( int run ){

  TH2F* f90vCharge_data  = F90vChargeHistogram( Form("run_%d.root",      run) );  f90vCharge_data  -> SetTitle("f90 vs Charge (1220); Charge (PE); f90");
  TH2F* f90vCharge_mc_er = F90vChargeHistogram( Form("run_%d_MCER.root", run) );  f90vCharge_mc_er -> SetTitle("f90 vs Charge (1220, MC ER); Charge (PE); f90");
  TH2F* f90vCharge_mc_nr = F90vChargeHistogram( Form("run_%d_MCNR.root", run) );  f90vCharge_mc_nr -> SetTitle("f90 vs Charge (1220, MC NR); Charge (PE); f90");

  f90vCharge_data  -> SetMarkerStyle(8);    f90vCharge_data  -> SetMarkerSize(0.5);
  f90vCharge_mc_er -> SetMarkerStyle(8);    f90vCharge_mc_er -> SetMarkerSize(0.3);   f90vCharge_mc_er -> SetMarkerColorAlpha(30, 0.5);
  f90vCharge_mc_nr -> SetMarkerStyle(8);    f90vCharge_mc_nr -> SetMarkerSize(0.3);   f90vCharge_mc_nr -> SetMarkerColor(46);

  TLine* f90_low = new TLine(0., 0.2,  1050., 0.2);    f90_low -> SetLineStyle(10);  f90_low -> SetLineWidth(3); f90_low -> SetLineColor(13);
  TLine* f90_mid = new TLine(0., 0.4,  1050., 0.4);    f90_mid -> SetLineStyle(10);  f90_mid -> SetLineWidth(3); f90_mid -> SetLineColor(13);
  TLine* f90_upp = new TLine(0., 0.65, 1050., 0.65);   f90_upp -> SetLineStyle(10);  f90_upp -> SetLineWidth(3); f90_upp -> SetLineColor(13);

  gROOT->SetBatch(kTRUE);

  Double_t y_size = 1000.;
  Double_t x_size = 2 * 1.61803398875 * y_size;

  TCanvas* f90vC_canvas = new TCanvas("f90vC_canvas", "f90vC_canvas", x_size, y_size);
  f90vC_canvas -> Divide(2,1);

  f90vC_canvas -> cd(1);
  f90vCharge_data -> Draw(); f90_low -> Draw("same"); f90_mid -> Draw("same"); f90_upp -> Draw("same");

  f90vC_canvas -> cd(2);
  f90vCharge_mc_nr -> Draw(); f90vCharge_mc_er -> Draw("same"); f90_low -> Draw("same"); f90_mid -> Draw("same"); f90_upp -> Draw("same");

  f90vC_canvas -> SaveAs("f90 v Charge Distribution.png"); f90vC_canvas -> SaveAs("f90 v Charge Distribution.pdf");

  gROOT->SetBatch(kFALSE);
}

void Images( int run ){

  F90vChargeCanvas( run );


}
