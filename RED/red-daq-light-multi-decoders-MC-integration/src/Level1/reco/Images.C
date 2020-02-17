#include <TCut.h>
#include <TFile.h>
#include <TLine.h>
#include <TString.h>
#include <TTree.h>

/* TH2F* F90vChargeHistogram ( TString file_name )
 *
 *  Summary of Function:
 *
 *    Generates a 2D histogram with the y-axis being the f90 parameter and the x-axis being the total charge
 *    of the event. The data is taken from a root file. The function then returns a pointer to the generated
 *    TH2F histogram.
 *
 *  Parameters   :  file_name >> the name of the file containing the data from which to construct the hist.
 *
 *  Return Value :  TH2F* f90vCharge.
 */
TH2F* F90vChargeHistogram( TString file_name, bool quality_cuts ){

  TFile* file = new TFile( file_name );
  TTree* reco; file -> GetObject("reco", reco);

  TCut cut_all;
  if (quality_cuts){
    TCut cut_f90_min         = "clusters[0].f90 > 0.2";
    TCut cut_f90_max         = "clusters[0].f90 < 0.6";
    TCut cut_charge_min      = "clusters[0].charge > 0.";
    TCut cut_charge_max      = "clusters[0].charge < 1000.";
    TCut cut_cluster_number  = "number_of_clusters == 1";
    TCut cut_rep             = "clusters[0].rep == 1";

    cut_all = cut_f90_min && cut_f90_max && cut_charge_max && cut_charge_min && cut_cluster_number && cut_rep;
  } else {
    TCut cut_f90_min         = "clusters[0].f90 > 0.";
    TCut cut_f90_max         = "clusters[0].f90 < 1.0";
    TCut cut_charge_min      = "clusters[0].charge > 0.";
    TCut cut_charge_max      = "clusters[0].charge < 1000.";

    cut_all = cut_f90_min && cut_f90_max && cut_charge_max && cut_charge_min;
  }

  reco -> Draw("clusters[0].f90:clusters[0].charge >> hist", cut_all, "goff");
  TH2F* f90vCharge = (TH2F*) gDirectory -> Get("hist");

  f90vCharge -> SetDirectory(0);
  file -> Close();

  return f90vCharge;
}

/* F90ClusterHist ( TString file_name, int cluster )
 *
 *  Summary of Function:
 *
 *    The function generates a histogram of f90. It automatically applies quality cuts (f90 > 0 and f90 < 1).
 *    The user must also define the specific cluster from witch to construct the histogram, where each one
 *    indicates a specific part of the signal (1 = S1, 2 = S2 and 3 = S3).
 *
 *  Parameters   : file_name >> the root file containing the data to construct the histogram.
 *                 cluster   >> an integer indicating which cluster to consider.
 *
 *  Return Value : TH1F* f90cluster_hist.
 */
TH1F* F90ClusterHist( TString file_name, int cluster ){

  TFile* file = new TFile( file_name );
  TTree* reco; file -> GetObject("reco", reco);

  TCut cut_f90_min = Form("clusters[%d].f90 > 0.0",  cluster);
  TCut cut_f90_max = Form("clusters[%d].f90 < 1.0",  cluster);

  TCut cut_all = cut_f90_min && cut_f90_max;

  reco -> Draw( Form("clusters[%d].f90 >> hist", cluster), cut_all, "goff" );
  TH1F* f90cluster_hist = (TH1F*) gDirectory -> Get("hist");

  f90cluster_hist -> SetDirectory(0);
  file -> Close();

  return f90cluster_hist;
}

/* TH1F* F90SIPMHist ( TString file_name, int sipm_number )
 *
 *  Summary of Function:
 *
 *    The function creates a histogram of f90, but it considers only the signal detected by a specific sipm.
 *
 *  Parameters   : file_name   >> the root file containing the data to construc the histogram.
 *                 sipm_number >> the specific sipm for witch the f90 histogram is generated.
 *
 *  Return Value : TH1F* f90sipm_hist.
 */
TH1F* F90SIPMHist( TString file_name, int sipm_number ){

  TFile* file = new TFile( file_name );
  TTree* reco; file -> GetObject("reco", reco);

  TCut cut_f90_min = Form("f90[%d] > 0.0", sipm_number);
  TCut cut_f90_max = Form("f90[%d] < 1.0", sipm_number);

  TCut cut_all = cut_f90_min && cut_f90_max;

  reco -> Draw( Form("f90[%d] >> hist", sipm_number), cut_all, "goff" );
  TH1F* f90sipm_hist = (TH1F*) gDirectory -> Get("hist");

  f90sipm_hist -> SetDirectory(0);
  file -> Close();

  return f90sipm_hist;
}

void GenerateF90vChargeCanvas( int run, bool quality_cuts ){

  TH2F* f90vCharge_data  = F90vChargeHistogram( Form("runs/run_%d.root",      run), quality_cuts );  f90vCharge_data  -> SetTitle("f90 vs Charge (1220); Charge (PE); f90");
  TH2F* f90vCharge_mc_er = F90vChargeHistogram( Form("runs/run_%d_MCER.root", run), quality_cuts );  f90vCharge_mc_er -> SetTitle("f90 vs Charge (1220, MC ER); Charge (PE); f90");
  TH2F* f90vCharge_mc_nr = F90vChargeHistogram( Form("runs/run_%d_MCNR.root", run), quality_cuts );  f90vCharge_mc_nr -> SetTitle("f90 vs Charge (1220, MC NR); Charge (PE); f90");

  f90vCharge_data  -> SetMarkerStyle(8);    f90vCharge_data  -> SetMarkerSize(0.5);
  f90vCharge_mc_er -> SetMarkerStyle(8);    f90vCharge_mc_er -> SetMarkerSize(0.3);   f90vCharge_mc_er -> SetMarkerColorAlpha(30, 0.5);
  f90vCharge_mc_nr -> SetMarkerStyle(8);    f90vCharge_mc_nr -> SetMarkerSize(0.3);   f90vCharge_mc_nr -> SetMarkerColor(46);

  // This lines may not be necessary if I instead generate another image with all the quality cuts.
  TLine* f90_low = new TLine(0., 0.2, 1050., 0.2);  f90_low -> SetLineStyle(10);  f90_low -> SetLineWidth(3);  f90_low -> SetLineColor(13);
  TLine* f90_mid = new TLine(0., 0.4, 1050., 0.4);  f90_mid -> SetLineStyle(10);  f90_mid -> SetLineWidth(3);  f90_mid -> SetLineColor(13);
  TLine* f90_upp = new TLine(0., 0.6,  900., 0.6);  f90_upp -> SetLineStyle(10);  f90_upp -> SetLineWidth(3);  f90_upp -> SetLineColor(13);

  Double_t y_size = 1000.;
  Double_t x_size = 2 * 1.61803398875 * y_size;

  TCanvas* f90vC_canvas = new TCanvas("f90vC_canvas", "f90vC_canvas", x_size, y_size);
  f90vC_canvas -> Divide(2,1);

  f90vC_canvas -> cd(1);
  f90vCharge_data -> Draw(); f90_low -> Draw("SAME"); f90_mid -> Draw("SAME"); f90_upp -> Draw("SAME");

  f90vC_canvas -> cd(2);
  f90vCharge_mc_nr -> Draw(); f90vCharge_mc_er -> Draw("SAME"); f90_low -> Draw("SAME"); f90_mid -> Draw("SAME"); f90_upp -> Draw("SAME");

  if ( quality_cuts ){
    f90vC_canvas -> SaveAs("plots/1220/Study of Monte Carlo/Quality Cuts/f90 v Charge Distribution (Quality Cuts).png");
    f90vC_canvas -> SaveAs("plots/1220/Study of Monte Carlo/Quality Cuts/f90 v Charge Distribution (Quality Cuts).pdf");
  } else {
    f90vC_canvas -> SaveAs("plots/1220/Study of Monte Carlo/Quality Cuts/f90 v Charge Distribution (No Cuts).png");
    f90vC_canvas -> SaveAs("plots/1220/Study of Monte Carlo/Quality Cuts/f90 v Charge Distribution (No Cuts).pdf");
  }
}

void GenerateF90comparisonCanvas( int run, int cluster, int number_of_sipm, bool isMC ){

  TH1F* f90cluster_hist = F90ClusterHist(Form("runs/run_%d%s.root", run, isMC ? "_MCER":""), cluster);
  f90cluster_hist -> SetTitle("f90 Distribution (Cluster); f90");

  TH1F** f90sipm_hist = new TH1F*[number_of_sipm];
  for (Int_t i = 0; i < number_of_sipm; i++){
    f90sipm_hist[i] = F90SIPMHist( Form("runs/run_%d%s.root", run, isMC ? "_MCER":""), i );
  }

  Double_t y_size = 1000.;
  Double_t x_size = 2 * 1.61803398875 * y_size;

  TCanvas* f90comparison_canvas = new TCanvas("f90comparison_canvas", "f90comparison_canvas", x_size, y_size);
  f90comparison_canvas -> Divide(2,1);

  f90comparison_canvas -> cd(1);  f90cluster_hist -> Draw("hist");

  f90comparison_canvas -> cd(2);
  for (int i = 0; i < number_of_sipm; i++){
    f90sipm_hist[0] -> SetTitle("f90 Distribution (SIPMs); f90");
    f90sipm_hist[i] -> Draw(i == 0 ? "hist":"SAME");
  }

  if ( isMC ){
    f90comparison_canvas -> SaveAs("plots/1220/Study of Monte Carlo/f90 Distribution Behaviour/f90 Distribution Comparison (1220, MC).png");
    f90comparison_canvas -> SaveAs("plots/1220/Study of Monte Carlo/f90 Distribution Behaviour/f90 Distribution Comparison (1220, MC).pdf");
  } else {
    f90comparison_canvas -> SaveAs("plots/1220/Study of Monte Carlo/f90 Distribution Behaviour/f90 Distribution Comparison (1220, Data).png");
    f90comparison_canvas -> SaveAs("plots/1220/Study of Monte Carlo/f90 Distribution Behaviour/f90 Distribution Comparison (1220, Data).pdf");
  }

}

// ------------------------------------- Images MACRO ------------------------------------- //
void Images( int run ){

  int cluster = 0;
  int number_of_sipm = 28;

  GenerateF90vChargeCanvas( run, false );
  GenerateF90vChargeCanvas( run, true );

  GenerateF90comparisonCanvas( run, cluster, number_of_sipm, false );
  GenerateF90comparisonCanvas( run, cluster, number_of_sipm, true );
}
