#include <iostream>
#include <array>
#include <string>

#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

TH1F* Generatef90Hist(char *file_name, Double_t charge_max, Double_t f90_min, Double_t f90_max, int number_divisions, int bin_number) {

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

// A function to create a directory that check wether the histogram already exists. If it does, it returns a pointer to it. If not, it creates one with there
// desired name and title and returns a pointer to it.
TDirectory* MakeDirectory(string dir_name, string dir_title){

  const char *dir1 = dir_name.c_str();
  const char *dir2 = dir_title.c_str();
  TDirectory *dir = gDirectory -> GetDirectory(dir1);
  if (!dir) dir = gDirectory -> mkdir(dir1, dir2);

  return dir;
}

int GenerateAllHist(int run, int number_divisions, Double_t f90_min, Double_t f90_mid, Double_t f90_max, Double_t max_charge_run,
                    Double_t max_charge_MC_ER, Double_t max_charge_MC_NR){

  std::array<Double_t, 3>    max_charge  = {max_charge_run, max_charge_MC_ER, max_charge_MC_NR};
  std::array<std::string, 3> file_suffix = {"", "MC_ER", "MC_NR"};
  std::array<std::string, 6> dir_name = {"f90_histograms", "data", "monte_carlo", "both", "ER", "NR"};

  TH1F** f90hist_run    = new TH1F*[number_divisions];
  TH1F** f90hist_run_ER = new TH1F*[number_divisions];
  TH1F** f90hist_run_NR = new TH1F*[number_divisions];
  TH1F** f90hist_MC_ER  = new TH1F*[number_divisions];
  TH1F** f90hist_MC_NR  = new TH1F*[number_divisions];

  // Cheks wether the root file already exists and tells the user. The "UPDATE" option already takes into account the possibility of the file not existing.
  if (gSystem -> AccessPathName(Form("hist_%d.root", run))){
    std::cout << "The " << Form("hist_%d.root", run) << " file does not exist. Creating it..." << std::endl;
  } else {
    std::cout << "Opening file: " << Form("hist_%d.root", run) << std::endl;
  }

  TFile *hist_file = new TFile(Form("hist_%d.root", run), "UPDATE");
  if ( !(hist_file -> IsOpen()) ) { std::cout << "Unable to open file." << std::endl;}

  // Checks wether a folder for saving the f90 histograms exists. If not, create it. Then repeat the process for all sub-directories.
  TDirectory *f90_dir = MakeDirectory(dir_name[0], dir_name[0]);

  f90_dir -> cd();
  TDirectory *data_dir    = MakeDirectory(dir_name[1], dir_name[1]);
  TDirectory *MC_dir      = MakeDirectory(dir_name[2], dir_name[2]);

  data_dir -> cd();
  TDirectory *data_dir_both = MakeDirectory(dir_name[3], dir_name[3]);
  TDirectory *data_dir_ER   = MakeDirectory(dir_name[4], dir_name[4]);
  TDirectory *data_dir_NR   = MakeDirectory(dir_name[5], dir_name[5]);

  MC_dir -> cd();
  TDirectory *MC_dir_ER   = MakeDirectory(dir_name[4], dir_name[4]);
  TDirectory *MC_dir_NR   = MakeDirectory(dir_name[5], dir_name[5]);


  // Generate all the desire f90 histograms (run, run er, run nr, er only and nr only).
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < number_divisions; j++){
      if (i == 0){
        f90hist_run[j] = Generatef90Hist(Form("run_%d%s.root", run, file_suffix[i].c_str()), max_charge[i], f90_min, f90_max, number_divisions, j+1);
        f90hist_run[j] -> SetName(Form("f90_distribution%d", j+1));
        f90hist_run[j] -> SetTitle(Form("f90 Distribution (Bin Number: %d); f90", j+1));

        data_dir_both -> WriteObject(f90hist_run[j], Form("f90_distribution%d", j+1), "OverWrite");


        f90hist_run_ER[j] = Generatef90Hist(Form("run_%d%s.root", run, file_suffix[i].c_str()), max_charge[i], f90_min, f90_mid, number_divisions, j+1);
        f90hist_run_ER[j] -> SetName(Form("f90_distribution_ER%d", j+1));
        f90hist_run_ER[j] -> SetTitle(Form("f90 Distribution (ER, Bin Number: %d); f09", j+1));

        data_dir_ER -> WriteObject(f90hist_run_ER[j], Form("f90_distribution_ER%d", j+1), "OverWrite");


        f90hist_run_NR[j] = Generatef90Hist(Form("run_%d%s.root", run, file_suffix[i].c_str()), max_charge[i], f90_mid, f90_max, number_divisions, j+1);
        f90hist_run_NR[j] -> SetName(Form("f90_distribution_NR%d", j+1));
        f90hist_run_NR[j] -> SetTitle(Form("f90 Distribution (NR, Bin Number: %d); f09", j+1));

        data_dir_NR -> WriteObject(f90hist_run_NR[j], Form("f90_distribution_NR%d", j+1), "OverWrite");

      } else if (i == 1) {
        f90hist_MC_ER[j] = Generatef90Hist(Form("run_%d%s.root", run, file_suffix[i].c_str()), max_charge[i], f90_min, f90_mid, number_divisions, j+1);
        f90hist_MC_ER[j] -> SetName(Form("f90_distribution_%s%d", file_suffix[i].c_str(), j+1));
        f90hist_MC_ER[j] -> SetTitle(Form("f90 Distrbution (MC, Bin Number: %d); f90", j+1));

        MC_dir_ER -> WriteObject(f90hist_MC_ER[j], Form("f90_distribution_ER%d", j+1), "OverWrite");

      } else if (i == 2) {
        f90hist_MC_NR[j] = Generatef90Hist(Form("run_%d%s.root", run, file_suffix[i].c_str()), max_charge[i], f90_mid, f90_max, number_divisions, j+1);
        f90hist_MC_NR[j] -> SetName(Form("f90_distribution_%s%d", file_suffix[i].c_str(), j+1));
        f90hist_MC_NR[j] -> SetTitle(Form("f90 Distrbution (MC, Bin Number: %d); f90", j+1));

        MC_dir_NR -> WriteObject(f90hist_MC_NR[j], Form("f90_distribution_NR%d", j+1), "OverWrite");
      }
    }
  }

  return 0;
}

void CompareSimulation(int run, int number_divisions, bool hist_exist = false, Double_t f90_min = 0., Double_t f90_mid = 0.4, Double_t f90_max = 1.,
                       Double_t max_charge_run = 2000., Double_t max_charge_MC_ER = 20000., Double_t max_charge_MC_NR = 9000.){

  /* If the histograms to be studied have yet to be generated (hist_exist = false), the function bellow is called, generating them. If the f90_histograms
  already exist (his_exist = true), then the function bellow is called and instead the histograms are obtained from the root file. */
  if (!hist_exist){GenerateAllHist(run, number_divisions, f90_min, f90_mid, f90_max, max_charge_run, max_charge_MC_ER, max_charge_MC_NR);}

  TH1F** f90hist_run_ER = new TH1F*[number_divisions];
  TH1F** f90hist_run_NR = new TH1F*[number_divisions];
  TH1F** f90hist_MC_ER  = new TH1F*[number_divisions];
  TH1F** f90hist_MC_NR  = new TH1F*[number_divisions];

  std::array<std::string, 6> dir_name = {"f90_histograms", "data", "monte_carlo", "both", "ER", "NR"};

  TFile *hist_file = new TFile("hist_1220.root", "UPDATE");
  TDirectory *f90_dir = MakeDirectory(dir_name[0], dir_name[0]);

  f90_dir -> cd();
  TDirectory *data_dir    = MakeDirectory(dir_name[1], dir_name[1]);
  TDirectory *MC_dir      = MakeDirectory(dir_name[2], dir_name[2]);

  data_dir -> cd();
  TDirectory *data_dir_both = MakeDirectory(dir_name[3], dir_name[3]);
  TDirectory *data_dir_ER   = MakeDirectory(dir_name[4], dir_name[4]);
  TDirectory *data_dir_NR   = MakeDirectory(dir_name[5], dir_name[5]);

  MC_dir -> cd();
  TDirectory *MC_dir_ER   = MakeDirectory(dir_name[4], dir_name[4]);
  TDirectory *MC_dir_NR   = MakeDirectory(dir_name[5], dir_name[5]);

  for (int i = 0; i < number_divisions; i++){
    f90hist_run_ER[i] = (TH1F *)data_dir_ER -> Get(Form("f90_distribution_ER%d", i+1));
  }
}
