#include <iostream>
#include <array>
#include <string>

//#include "../../Modules/RDCluster.hh"
#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>


TH1F* Createf90Dist(char *file_name, Double_t charge_max, Double_t f90_min, Double_t f90_max, int number_divisions, int bin_number) {
  TFile *run_data = new TFile(file_name);
  TTree *reco_data; run_data -> GetObject("reco", reco_data);

  Double_t division_size = charge_max/number_divisions;
  Double_t cbin_up  = bin_number * division_size;
  Double_t cbin_low = (bin_number - 1) * division_size;

  reco_data -> Draw("clusters[0].f90 >> hist",
               Form("clusters[0].f90 >= %f && clusters[0].f90 <= %f && clusters[0].charge >= %f && clusters[0].charge <= %f", f90_min, f90_max, cbin_low, cbin_up),
               "goff");
  TH1F* f90_dist = (TH1F*)gDirectory -> Get("hist");

  if (f90_dist -> GetSumw2N() == 0) f90_dist -> Sumw2(kTRUE);

  Double_t norm = 1.0;
  Double_t scale = norm/(f90_dist -> Integral());
  f90_dist -> Scale(scale);

  return f90_dist;
}

void CompareSimulation(int run, int number_divisions = 10, Double_t max_charge_run = 2000., Double_t max_charge_MCER = 20000.0, Double_t max_charge_MCNR = 9000.0){

  std::array<Double_t, 3> max_charge = {max_charge_run, max_charge_MCER, max_charge_MCNR};

  TH1F** f90hist_run  = new TH1F*[number_divisions];
  TH1F** f90hist_MCER = new TH1F*[number_divisions];
  TH1F** f90hist_MCNR = new TH1F*[number_divisions];

  string file_suffix[3] = {"", "MCER", "MCNR"};

  /*for (int i = 0; i < file_suffix.size(); i++){

    for (int j = 0; j < number_divisions; j++){
      if (i == 0){
        f90hist_run[j] = Createf90Dist(Form("run_%d%s.root", run, file_suffix[i].c_str()), max_charge[i], 0., 1., number_divisions, j+1);
      } else if (i == 1) {
        f90hist_MCER[j] = Createf90Dist(Form("run_%d%s.root", run, file_suffix[i].c_str()), max_charge[i], 0., 1., number_divisions, j+1);
      } else if (i == 2) {
        f90hist_MCER[j] = Createf90Dist(Form("run_%d%s.root", run, file_suffix[i].c_str()), max_charge[i], 0., 1., number_divisions, j+1);
      }
    }

  }*/


  TH1F* dist = Createf90Dist(Form("run_%d%s.root", run, file_suffix[0].c_str()), max_charge[0], 0., 1., number_divisions, 4);
  dist -> Draw();

}
