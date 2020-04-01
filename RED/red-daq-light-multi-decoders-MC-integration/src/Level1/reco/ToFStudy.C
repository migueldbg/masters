#include <TCut.h>
#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

void study(int run, Double_t charge_min = 20, Double_t charge_max = 1000){

  TFile *file = new TFile(Form("runs/run_%d.root", run));
  TTree *reco = (TTree*) file->Get("reco");

  Double_t x_size = 500.;
  Double_t y_size = 2 * 1.61803398875 * x_size;

  TCanvas *canvas = new TCanvas("canvas", "canvas", x_size, y_size);
  canvas -> Divide(2,1);

  TH2F *hist_nocut = new TH2F("hist_nocut", "No ToF Cut", 100, charge_min, charge_max, 100, 0, 1);
  reco -> Draw("clusters[0].f90:clusters[0].charge >> hist_nocut", "clusters[0].f90 > 0 && clusters[0].f90 < 1. && clusters[0].charge > 20 && clusters[0].charge < 1000 && number_of_clusters == 2 && clusters[0].rep == 1 && clusters[1].rep == 1");
  //hist_nocut -> SetDirectory(0);

  TH2F *hist_tofcut = new TH2F("hist_tofcut", "ToF Cut", 100, charge_min, charge_max, 100, 0, 1);
  reco -> Draw("clusters[0].f90:clusters[0].charge >> hist_nocut", "clusters[0].f90 > 0 && clusters[0].f90 < 1. && clusters[0].charge > 20 && clusters[0].charge < 1000 && xmin[30] - clusters[0].min_x < 100 && xmin[30] - clusters[0].min_x > 0 && number_of_clusters == 2 && clusters[0].rep == 1 && clusters[1].rep == 1");
  //hist_tofcut -> SetDirectory(0);

  canvas -> cd(1);
  hist_nocut -> Draw();

  canvas -> cd(2);
  hist_tofcut -> Draw();

  //file -> Close();
}
