#include <iostream>

//#include "../../Modules/RDCluster.hh"
#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>


TH1F* Createf90Dist(char *file_name, Double_t charge_max, Double_t f90_min, Double_t f90_max, Double_t division_size, int bin_number) {
  TFile *run_data = new TFile(file_name);
  TTree *reco_data; run_data -> GetObject("reco", reco_data);

  Double_t cbin_up  = bin_number * division_size;
  Double_t cbin_low = (bin_number - 1) * division_size;

  reco_data -> Draw("clusters[0].f90 >> hist",
               Form("clusters[0].f90 >= %f && clusters[0].f90 <= %f && clusters[0].charge >= %f && clusters[0].charge <= %f", f90_min, f90_max, cbin_low, cbin_up),"goff");
  TH1F* f90_dist = (TH1F*)gDirectory -> Get("hist");

  if (f90_dist -> GetSumw2N() == 0) f90_dist -> Sumw2(kTRUE);

  Double_t norm = 1.0;
  Double_t scale = norm/(f90_dist -> Integral());
  f90_dist -> Scale(scale);

  return f90_dist;
}

void CompareSimulation(int run, int number_divisions = 10, double max_charge_run = 2000){

  //double max_charge_MCER = 20000;
  //double max_charge_MCNR = 9000;


  TH1F* dist = Createf90Dist(Form("run_1220.root"), max_charge_run, 0., 1., 200., 4);
  dist -> Draw();

  /*
  TFile *run_file = new TFile(Form("run_%d.root", run));
  TFile *MCER_file = new TFile(Form("run_%dMCER.root", run));
  TFile *MCNR_file = new TFile(Form("run_%dMCNR.root", run));

  TTree *reco_run;    run_file  -> GetObject("reco", reco_run);
  TTree *reco_MCER;   MCER_file -> GetObject("reco", reco_MCER);
  TTree *reco_MCNR;   MCNR_file -> GetObject("reco", reco_MCNR);

  double division_size_run  = max_charge_run/number_divisions;
  double division_size_MCER = max_charge_MCER/number_divisions;
  double division_size_MCNR = max_charge_MCNR/number_divisions;
  */
}
