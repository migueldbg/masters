<<<<<<< HEAD
=======
#include "sidutility.cc"

>>>>>>> runs_check
#include <TCut.h>
#include <TFile.h>
#include <THStack.h>
#include <TH1.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

<<<<<<< HEAD
/* TFile* CheckFile( TString path_name )
 *
 * Summary of CheckFile function:
 *
 *    The CheckFile function checks if a given file already exists and returns this information to the user.
 *    If the file exists, CheckFile returns it. If not, it creates it and then returns it.
 *
 * Parameters   : path_name >> The name of the desired file.
 *
 * Return value : TFile* file.
 */
TFile* CheckFile( TString path_name ){

  if (gSystem -> AccessPathName(path_name)) {
    std::cout << "The " << path_name << " file does not exist. Creating it..." << std::endl;
  } else {
    std::cout << "Opening file: " << path_name << std::endl;
  }

  TFile *file = new TFile(path_name, "UPDATE");
  if ( !(file -> IsOpen()) ) { std::cout << "Unable to open file." << std::endl;}

  return file;
}

TCut DefineCuts(){

  TCut f90_min            = "clusters[0].f90 >= 0.";
  TCut f90_max            = "clusters[0].f90 <= 1.";
  TCut rep                = "clusters[0].rep == 1 && clusters[1].rep == 1";
  TCut number_of_clusters = "number_of_clusters == 2";

  TCut final_cut = f90_min && f90_max && rep && number_of_clusters;

  return final_cut;
}

/* TH1* NormalizeHist( TH1* hist, Double_t norm = 1. )
 *
 * Summary of NormalizeHist function:
 *
 *    The NormalizeHist function takes a TH1* histogram, makes a copy and normalizes the copy to the
 *    desired value. It then returns the normalized copy of histogram.
 *
 * Parameters   : hist >> histogram to be normalized.
 *                norm >> value to normalize the histogram to (default = 1.0).
 *
 * Return value : TH1* normalized_hist.
 */
TH1* NormalizeHist( TH1* hist, Double_t norm = 1. ){

  Double_t scale = norm / (hist -> Integral());

  TH1* normalized_hist = (TH1*)hist -> Clone();
  normalized_hist -> Scale(scale);

  return normalized_hist;
}
=======
>>>>>>> runs_check

void DrifTime(int run, const char* recoil_type ){

  string recoil_suffix;

  if ( std::strncmp( recoil_type, "er", 2 ) == 0 || std::strncmp( recoil_type, "ER", 2 ) == 0 ) {
    recoil_suffix = "ER";
  } else if ( std::strncmp( recoil_type, "nr", 2 ) == 0 || std::strncmp( recoil_type, "NR", 2 ) == 0 ) {
    recoil_suffix = "NR";
  } else{
    std::cout << "Invalid 'recoil_type' parameter." << endl;
    exit(EXIT_FAILURE);
  }

  TFile *data_file = new TFile( Form("runs/run_%d.root", run) );
  TFile *mc_file   = new TFile( Form("runs/run_%d_MC%s.root", run, recoil_suffix.c_str()) );

  TTree *data_reco = (TTree*) data_file -> Get("reco");
  TTree *mc_reco   = (TTree*) mc_file   -> Get("reco");

  const char *drift_time = "clusters[1].start_time - clusters[0].start_time >> hist";
  TCut quality_cuts = DefineCuts();

  THStack *dt_hist = new THStack("drif_time_hist", "Drif Time");

  // Drift time histogram from data:
  TH1F *da_dt_hist = new TH1F("da_dt_hist", "Drift Time", 100, 0, 40000);
  data_reco -> Draw(drift_time, quality_cuts,"goff");

  da_dt_hist = (TH1F*) gDirectory -> Get("hist");
  da_dt_hist = (TH1F*) NormalizeHist(da_dt_hist);

  da_dt_hist -> SetLineColor(kRed);

  dt_hist -> Add(da_dt_hist);
  Double_t da_nevents = da_dt_hist -> GetEntries();

  // Drift time histogram from Monte Carlo simulation:
  TH1F *mc_dt_hist = new TH1F("mc_dt_hist", "Drift Time (MC)", 100, 0, 40000);
  mc_reco -> Draw(drift_time, quality_cuts,"goff");

  mc_dt_hist = (TH1F*) gDirectory -> Get("hist");
  mc_dt_hist = (TH1F*) NormalizeHist(mc_dt_hist);

  mc_dt_hist -> SetLineColor(kBlue);

  dt_hist -> Add(mc_dt_hist);
  Double_t mc_nevents = mc_dt_hist -> GetEntries();

  TCanvas *canvas = new TCanvas("canvas", "canvas", 500*1.62, 500);

  dt_hist -> Draw("hist nostack");
  dt_hist -> GetXaxis() -> SetTitle("Drift Time (samples)");
  dt_hist -> GetYaxis() -> SetTitle("Counts (Normalized)");

  auto legend = new TLegend(0.710396,0.751579,0.888614,0.882105);
  legend -> AddEntry(da_dt_hist, "Data", "l");
  legend -> AddEntry(mc_dt_hist, "MC",   "l");
  legend -> Draw();

  canvas->Modified();

  std::cout << "Number of Events (Data): " << da_nevents << endl;
  std::cout << "Number of Events (MC): "   << mc_nevents << endl;

}
