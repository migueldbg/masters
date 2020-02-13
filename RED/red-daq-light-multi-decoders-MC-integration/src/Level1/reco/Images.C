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

TH1F* F90ClusterHist( TString file_name, int cluster ){

  TFile* file = new TFile( file_name );
  TTree* reco; file -> GetObject("reco", reco);

  reco -> Draw( Form("clusters[%d].f90 >> hist", cluster), Form("clusters[%d].f90 > 0. && clusters[%d].f90 < 1.", cluster, cluster), "goff" );
  TH1F* f90cluster_hist = (TH1F*) gDirectory -> Get("hist");

  f90cluster_hist -> SetDirectory(0);
  file -> Close();

  return f90cluster_hist;
}

TH1F* F90SIPMHist( TString file_name, int sipm_number ){

  TFile* file = new TFile( file_name );
  TTree* reco; file -> GetObject("reco", reco);

  reco -> Draw( Form("f90[%d] >> hist", sipm_number), Form("f90[%d] > 0. && f90[%d] < 1.", sipm_number, sipm_number), "goff" );
  TH1F* f90sipm_hist = (TH1F*) gDirectory -> Get("hist");

  f90sipm_hist -> SetDirectory(0);
  file -> Close();

  return f90sipm_hist;
}

void GenerateF90vChargeCanvas( int run ){

  TH2F* f90vCharge_data  = F90vChargeHistogram( Form("run_%d.root",      run) );  f90vCharge_data  -> SetTitle("f90 vs Charge (1220); Charge (PE); f90");
  TH2F* f90vCharge_mc_er = F90vChargeHistogram( Form("run_%d_MCER.root", run) );  f90vCharge_mc_er -> SetTitle("f90 vs Charge (1220, MC ER); Charge (PE); f90");
  TH2F* f90vCharge_mc_nr = F90vChargeHistogram( Form("run_%d_MCNR.root", run) );  f90vCharge_mc_nr -> SetTitle("f90 vs Charge (1220, MC NR); Charge (PE); f90");

  f90vCharge_data  -> SetMarkerStyle(8);    f90vCharge_data  -> SetMarkerSize(0.5);
  f90vCharge_mc_er -> SetMarkerStyle(8);    f90vCharge_mc_er -> SetMarkerSize(0.3);   f90vCharge_mc_er -> SetMarkerColorAlpha(30, 0.5);
  f90vCharge_mc_nr -> SetMarkerStyle(8);    f90vCharge_mc_nr -> SetMarkerSize(0.3);   f90vCharge_mc_nr -> SetMarkerColor(46);

  TLine* f90_low = new TLine(0., 0.2,  1050., 0.2);    f90_low -> SetLineStyle(10);  f90_low -> SetLineWidth(3); f90_low -> SetLineColor(13);
  TLine* f90_mid = new TLine(0., 0.4,  1050., 0.4);    f90_mid -> SetLineStyle(10);  f90_mid -> SetLineWidth(3); f90_mid -> SetLineColor(13);
  TLine* f90_upp = new TLine(0., 0.65, 1050., 0.65);   f90_upp -> SetLineStyle(10);  f90_upp -> SetLineWidth(3); f90_upp -> SetLineColor(13);

  Double_t y_size = 1000.;
  Double_t x_size = 2 * 1.61803398875 * y_size;

  TCanvas* f90vC_canvas = new TCanvas("f90vC_canvas", "f90vC_canvas", x_size, y_size);
  f90vC_canvas -> Divide(2,1);

  f90vC_canvas -> cd(1);
  f90vCharge_data -> Draw(); f90_low -> Draw("same"); f90_mid -> Draw("same"); f90_upp -> Draw("same");

  f90vC_canvas -> cd(2);
  f90vCharge_mc_nr -> Draw(); f90vCharge_mc_er -> Draw("same"); f90_low -> Draw("same"); f90_mid -> Draw("same"); f90_upp -> Draw("same");

  f90vC_canvas -> SaveAs("f90 v Charge Distribution.png"); f90vC_canvas -> SaveAs("f90 v Charge Distribution.pdf");
}

void GenerateF90comparisonCanvas( int run, int cluster, int number_of_sipm, bool isMC = false){

  TH1F* f90cluster_hist = F90ClusterHist(Form("run_%d%s.root", run, isMC ? "_MCER":""), cluster);  f90cluster_hist -> SetTitle("f90 Distribution (Cluster); f90");

  TH1F** f90sipm_hist = new TH1F*[number_of_sipm];
  for (Int_t i = 0; i < number_of_sipm; i++){
    f90sipm_hist[i] = F90SIPMHist( Form("run_%d%s.root", run, isMC ? "_MCER":""), i );
  }

  Double_t y_size = 1000.;
  Double_t x_size = 2 * 1.61803398875 * y_size;

  TCanvas* f90comparison_canvas = new TCanvas("f90comparison_canvas", "f90comparison_canvas", x_size, y_size);
  f90comparison_canvas -> Divide(2,1);

  f90comparison_canvas -> cd(1);  f90cluster_hist -> Draw("hist");

  f90comparison_canvas -> cd(2);
  for (int i = 0; i < number_of_sipm; i++){
    f90sipm_hist[0] -> SetTitle("f90 Distribution (SIPMs); f90");
    f90sipm_hist[i] -> Draw(i == 0 ? "hist":"same");
  }

  if ( isMC ){
    f90comparison_canvas -> SaveAs("f90 Distribution Comparison (1220, MC).png");
    f90comparison_canvas -> SaveAs("f90 Distribution Comparison (1220, MC).pdf");
  } else {
    f90comparison_canvas -> SaveAs("f90 Distribution Comparison (1220, Data).png");
    f90comparison_canvas -> SaveAs("f90 Distribution Comparison (1220, Data).pdf");
  }

}

void Images( int run ){

  int cluster = 0;
  int number_of_sipm = 28;
  bool isMC = true;

  GenerateF90vChargeCanvas( run );

  GenerateF90comparisonCanvas( run, cluster, number_of_sipm );

  GenerateF90comparisonCanvas( run, cluster, number_of_sipm, isMC );

}
