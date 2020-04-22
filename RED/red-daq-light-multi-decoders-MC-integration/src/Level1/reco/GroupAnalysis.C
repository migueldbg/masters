#include "sidutility.cc"

#include <TCut.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TLine.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <iterator>
#include <vector>

void MergeRuns(){

  std::vector<int> runs = {1501, 1502, 1503, 1504, 1506, 1508, 1509, 1510, 1511, 1512, 1513, 1514, 1515, 1516, 1519, 1520, 1521};

  TChain* reco_chain = new TChain("reco");
  TString file_name;
  for (Int_t i = 0; i < runs.size(); i++){
    file_name = Form("runs/run_%d.root", runs.at(i));
    reco_chain -> Add(file_name.Data());
  }

  TString output_file_name = Form("runs/run_%d%d.root", runs.front(), runs.back());
  reco_chain -> Merge(output_file_name.Data(), "fast");
}


TCut LowBeCut( Int_t run, Int_t num_bins, Int_t E_min, Int_t E_max, Int_t dE_min, Int_t dE_max ){

  TH2F* banana_hist = Bananator(run, num_bins, E_min, E_max, dE_min, dE_max);
  //                            15011521, 400, 1500, 4000, 400, 1200

  // Boundaries to enclose the low energy Berillium blob. This values are calibrated to the 1501-1521 runs.
  Double_t E_low  = 2250;  Double_t E_upp  = 2725;
  Double_t dE_low = 700;   Double_t dE_upp = 850;

  TLine* E_low_L  = new TLine(E_low, dE_low, E_low, dE_upp);
  TLine* E_upp_L  = new TLine(E_upp, dE_low, E_upp, dE_upp);
  TLine* dE_low_L = new TLine(E_low, dE_low, E_upp, dE_low);
  TLine* dE_upp_L = new TLine(E_low, dE_upp, E_upp, dE_upp);

  E_low_L  -> SetLineWidth(3);  E_low_L  -> SetLineColor(kRed);
  E_upp_L  -> SetLineWidth(3);  E_upp_L  -> SetLineColor(kRed);
  dE_low_L -> SetLineWidth(3);  dE_low_L -> SetLineColor(kRed);
  dE_upp_L -> SetLineWidth(3);  dE_upp_L -> SetLineColor(kRed);

  Double_t height = 600;
  Double_t width  = 1.618*height;
  TCanvas* canvas = new TCanvas("canvas", "canvas", width, height);

  gStyle -> SetPalette(kSunset);
  banana_hist -> SetOption("COLZ");
  banana_hist -> Draw();
  E_low_L -> Draw("SAME");  E_upp_L -> Draw("SAME");  dE_low_L -> Draw("SAME");  dE_upp_L -> Draw("SAME");

  TCut E_low_cut  = Form("baseline_mean[31] - ymin[31] > %f", E_low);
  TCut E_upp_cut  = Form("baseline_mean[31] - ymin[31] < %f", E_upp);
  TCut dE_low_cut = Form("baseline_mean[30] - ymin[30] > %f", dE_low);
  TCut dE_upp_cut = Form("baseline_mean[30] - ymin[30] < %f", dE_upp);

  TCut lowBe_cut = E_low_cut && E_upp_cut && dE_low_cut && dE_upp_cut;
  return lowBe_cut;
}


void TimeOfFlight(Int_t run){

  TString runFile_name = Form("runs/run_%d.root", run);
  bool normalize = false;

  Int_t cfg = 2;
  Double_t s1_min  = 50.;  Double_t s1_max  = 5000.;
  Double_t s2_min  = 50.;  Double_t s2_max  = 10000.;
  Double_t f90_min = 0.;   Double_t f90_max = 1.;
  Double_t tof_min = 0.;   Double_t tof_max = 100.;

  TCut histogram_cuts = DefineCuts(cfg, f90_min, f90_max, s1_min, s1_max, s2_min, s2_max, tof_min, tof_min);
  TCut lowBe_cut = LowBeCut( run, 400, 1500, 4000, 400, 1200 );

  TCut total_cut = histogram_cuts && lowBe_cut;

  TH1F* tof_histogram = GenerateToFHistogram(runFile_name, total_cut, 50, tof_min, tof_max, normalize);

  tof_histogram -> Draw();
}

void F90vToF(Int_t run){

  TString runFile_name = Form("runs/run_%d.root", run);
  TFile* file = CheckFile(runFile_name);
  TTree* reco;  file -> GetObject("reco", reco);

  Int_t cfg = 2;
  Double_t f90_min = 0.;   Double_t f90_max = 1.;
  Double_t tof_min = -100.;   Double_t tof_max = 100.;

  TCut histogram_cuts = DefineCuts(cfg, f90_min, f90_max, 0, 0, 0, 0, tof_min, tof_min);
  TCut lowBe_cut = LowBeCut( run, 400, 1500, 4000, 400, 1200 );

  TCut total_cut = histogram_cuts && lowBe_cut;

  reco -> Draw("clusters[0].f90:2*(xmin[30] - clusters[0].min_x) >> hist(100, -100, 100, 120, 0, 1)", total_cut, "colz");

}
