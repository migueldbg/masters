#include "sidutility.cc"

#include <TFile.h>
#include <TH2.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <vector>

/* ************************************************************************************************************************* *
 * File: TimeOfFlight.C (ROOT macro).                                                                                        *
 *                                                                                                                           *
 * Author: Miguel Del Ben Galdiano (miguel.galdiano@gmail.com).                                                              *
 * Date of creation : June 2 2020.                                                                                           *
 *                                                                                                                           *
 * Summary of File:                                                                                                          *
 *                                                                                                                           *
 *    This macro's goal is to do a study on the different time of flight parameters involved in the Litium beam experiment   *
 *    configuration. There is three such parameters: SiTel-TPC, SiTel-LSci and TPC-LSci. This macro contains a number of     *
 *    different analysis considering.                                                                                        *
 *                                                                                                                           *
 *    Am Runs: 1480.                                                                                                         *
 * ************************************************************************************************************************* */

Int_t    expCfg = 2;
Double_t f90Min = 0.0;
Double_t f90Max = 1.0;
Double_t s1Min  = 40.0;
Double_t s1Max  = 1000.0;

const double gdRatio = 1.61803398875;

void F90vS1( Int_t run ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> cd();

  TString file_name = Form("runs/run_%d.root", run);
  TFile* file = CheckFile(file_name);

  TTree* reco;  file -> GetObject("reco", reco);


  Double_t tpcToFMin = 20;    Double_t tpcToFMax = 30;
  Double_t lsciToFMin = -30;  Double_t lsciToFMax = -20;

  TCut tpcToF = Form("2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time) >= %f && 2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time) <= %f", tpcToFMin, tpcToFMax);
  TCut lsciToF = Form("2*(0.5*(start_time[30] + start_time[31] - 7.45) - start_time[0]) >= %f && 2*(0.5*(start_time[30] + start_time[31] - 7.45) - start_time[0]) <= %f", lsciToFMin, lsciToFMax);
  TCut qualityCut = DefineCuts(expCfg, f90Min, f90Max, s1Min, s1Max);
  TCutG* lowBeCut = LowBeGraphCut(run, "LowBeCut");
  TCut tpcCut = qualityCut && "LowBeCut";

  TH2F* f90vS1Hist[3];
  const char* histNames[3] = {"tpc_cut", "tpc_cut + tpcToF", "tpc_cut + tpcToF + lsciToF"};
  TCut cuts[3] = {tpcCut, tpcCut && tpcToF, tpcCut && tpcToF && lsciToF};

  Double_t height = 500; Double_t width = gdRatio * height;
  TCanvas* canvas = new TCanvas("canvas", "canvas", 3*width, height);
  canvas -> Divide(3,1);

  for (Int_t i = 0; i < 3; i++){

    f90vS1Hist[i] = new TH2F(histNames[i], "f90 v S1", 100, s1Min, s1Max, 100, f90Min, f90Max);
    reco -> Project(histNames[i], "clusters[0].f90:clusters[0].charge", cuts[i]);

    canvas -> cd(i + 1);
    if (i < 2){
      f90vS1Hist[i] -> Draw("COLZ");
    } else {
      f90vS1Hist[i] -> SetMarkerStyle(20);
      f90vS1Hist[i] -> SetMarkerSize(1);

      f90vS1Hist[i] -> Draw();
    }
  }

}


void ProjectX( TH1* hist ){



}

TStyle* SetSidStyle(){
  auto sidStyle = new TStyle("sidStyle", "Sid's Style");
  sidStyle -> SetPalette(kSunset);
  sidStyle -> SetLabelFont(102, "xyz");
  sidStyle -> SetTitleFont(102, "xyz");
  sidStyle -> SetTitleFont(102, "t");
  sidStyle -> SetCanvasBorderMode(0);
  sidStyle -> SetPadBorderMode(0);
  sidStyle -> SetTitleX(0.5);
  sidStyle -> SetTitleBorderSize(0);
  sidStyle -> SetTitleAlign(23);
  sidStyle -> SetOptStat(0);
  return sidStyle;
}

void SiTelTPCToFvS1( int run ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> cd();

  TString file_name = Form("runs/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);

  Double_t tofMin = 0;   Double_t tofMax = 50;
  Double_t tofBinSize = 0.25;
  Double_t tofBinNumber = (tofMax - tofMin)/tofBinSize;

  Double_t s1BinSize = 5;
  Double_t s1BinNumber = (s1Max - s1Min)/s1BinSize;

  TCut tpcCut     = DefineCuts(expCfg, f90Min, f90Max, s1Min, s1Max);
  TCutG* lowBeCut = LowBeGraphCut(run);
  TCut totalCut   = tpcCut && "lowBeCut";

  TH2F* ToFvS1Hist   = new TH2F("ToFvS1_hist", "SiTel-TPC ToF vs S1", s1BinNumber, s1Min, s1Max, tofBinNumber, tofMin, tofMax);
  TH2F* ToFvS1HistBe = new TH2F("ToFvS1_histBe", "SiTel-TPC ToF vs S1 (Low Be)", s1BinNumber, s1Min, s1Max, tofBinNumber, tofMin, tofMax);

  reco -> Project("ToFvS1_hist", "2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time):clusters[0].charge", tpcCut);
  reco -> Project("ToFvS1_histBe", "2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time):clusters[0].charge", totalCut);

  Double_t height = 500; Double_t width = gdRatio * height ;
  TCanvas* canvas = new TCanvas("canvas", "canvas", 2*width, height);
  canvas -> Divide(2,1);

  canvas -> cd(1);
  ToFvS1Hist -> GetXaxis() -> SetTitle("S1[PE]");
  ToFvS1Hist -> GetYaxis() -> SetTitle("ToF (SiTel - TPC)[ns]");
  ToFvS1Hist -> Draw("COLZ");

  canvas -> cd(2);
  ToFvS1HistBe -> GetXaxis() -> SetTitle("S1[PE]");
  ToFvS1HistBe ->  Draw("COLZ");

  ToFvS1HistBe -> ProjectionY("_py", -1, 0) -> Draw();

  std::cout << ToFvS1Hist -> GetNbinsX() << std::endl;
  std::cout << ToFvS1HistBe -> GetNbinsX() << std::endl;

  file -> Close();
}

void LSciTPCToF( int run ){

  gStyle -> SetLabelFont(102, "xyz");
  gStyle -> SetTitleFont(102, "xyz");
  gStyle -> SetTitleFont(102, "t");

  Double_t lsToFMin = -500;
  Double_t lsToFMax = 500;


  TString file_name = Form("runs/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);

  /*TCut histogramCuts = DefineCuts(exp_cfg, f90Min, f90Max, s1Min, s1Max);

  TCut E_low_cut  = Form("baseline_mean[31] - ymin[31] > %f", E_low);
  TCut E_upp_cut  = Form("baseline_mean[31] - ymin[31] < %f", E_upp);
  TCut dE_low_cut = Form("baseline_mean[30] - ymin[30] > %f", dE_low);
  TCut dE_upp_cut = Form("baseline_mean[30] - ymin[30] < %f", dE_upp);
  TCut lowBe_cut = E_low_cut && E_upp_cut && dE_low_cut && dE_upp_cut;*/

  //reco -> Project("hist", "2*(start_time[30] - start_time[1])", histogramCuts && lowBe_cut);

  reco -> Draw("(start_time[0] - 0.5*(start_time[31] + start_time[30] - 7.83))*2.", "charge[1]>2000 && Iteration$==0", "");
}
