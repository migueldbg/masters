#include <iostream>
#include <array>
#include <string>

#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>


TH1F* Createf90Dist(char *file_name, Double_t charge_max, Double_t f90_min, Double_t f90_max, int number_divisions, int bin_number) {
  TFile *file = new TFile(file_name);
  TTree *reco; file -> GetObject("reco", reco);

  Double_t division_size = charge_max/number_divisions;
  Double_t cbin_up  = bin_number * division_size;
  Double_t cbin_low = (bin_number - 1) * division_size;

  reco -> Draw("clusters[0].f90 >> hist",
               Form("clusters[0].f90 >= %f && clusters[0].f90 <= %f && clusters[0].charge >= %f && clusters[0].charge <= %f", f90_min, f90_max, cbin_low, cbin_up),
               "goff");
  TH1F* f90_dist = (TH1F*)gDirectory -> Get("hist");

  if (f90_dist -> GetSumw2N() == 0) f90_dist -> Sumw2(kTRUE);

  // Normalize the histogram for better comparision later on.
  Double_t norm = 1.0;
  Double_t scale = norm/(f90_dist -> Integral());
  f90_dist -> Scale(scale);

  // Remove the histogram from the current directory (TFile) so that I can close the file and still have access to it.
  f90_dist -> SetDirectory(0);
  file -> Close();

  return f90_dist;
}

void CompareSimulation(int run, int number_divisions = 10, Double_t max_charge_run = 2000., Double_t max_charge_MCER = 20000., Double_t max_charge_MCNR = 9000.){

  std::array<Double_t, 3>    max_charge  = {max_charge_run, max_charge_MCER, max_charge_MCNR};
  std::array<std::string, 3> file_suffix = {"", "MCER", "MCNR"};

  TH1F** f90hist_run  = new TH1F*[number_divisions];
  TH1F** f90hist_MCER = new TH1F*[number_divisions];
  TH1F** f90hist_MCNR = new TH1F*[number_divisions];

  // Cheks wether the root file already exists and tells the user. The "UPDATE" option already takes into account the possibility of the file not existing.
  if (gSystem -> AccessPathName(Form("hist_%d.root", run))){
    std::cout << "The " << Form("hist_%d.root", run) << " file does not exist. Creating it..." << std::endl;

  } else {
    std::cout << "Opening file: " << Form("hist_%d.root", run) << std::endl;
  }

  TFile *hist_file = new TFile(Form("hist_%d.root", run), "UPDATE");
  if ( !(hist_file -> IsOpen()) ) { std::cout << "Unable to open file." << std::endl;}

  for (int i = 0; i < 3; i++){
    for (int j = 0; j < number_divisions; j++){
      if (i == 0){
        f90hist_run[j] = Createf90Dist(Form("run_%d%s.root", run, file_suffix[i].c_str()), max_charge[i], 0., 1., number_divisions, j+1);
        f90hist_run[j] -> SetName(Form("f90_distribution"));
        f90hist_run[j] -> SetTitle("f90 Distrbution; f90");
        f90hist_run[j] -> Write(Form("f90_distribution"), TObject::kOverwrite);
      } else if (i == 1) {
        f90hist_MCER[j] = Createf90Dist(Form("run_%d%s.root", run, file_suffix[i].c_str()), max_charge[i], 0., 1., number_divisions, j+1);
        f90hist_MCER[j] -> SetName(Form("f90_distribution_%s", file_suffix[i].c_str()));
        f90hist_MCER[j] -> SetTitle("f90 Distrbution; f90");
        f90hist_MCER[j] -> Write(Form("f90_distribution"), TObject::kOverwrite);
      } else if (i == 2) {
        f90hist_MCNR[j] = Createf90Dist(Form("run_%d%s.root", run, file_suffix[i].c_str()), max_charge[i], 0., 1., number_divisions, j+1);
        f90hist_MCNR[j] -> SetName(Form("f90_distribution_%s", file_suffix[i].c_str()));
        f90hist_MCNR[j] -> SetTitle("f90 Distrbution; f90");
        f90hist_MCNR[j] -> Write(Form("f90_distribution"), TObject::kOverwrite);
      }
    }
  }

  f90hist_MCNR[3] -> Draw("HIST");
  hist_file -> Close();

}
