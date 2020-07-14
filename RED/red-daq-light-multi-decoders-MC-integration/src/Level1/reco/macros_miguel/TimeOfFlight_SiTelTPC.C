#include "sidutility.cc"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

/* ************************************************************************************************************************* *
 * File: TimeOfFlight_SiTelTPC.C (ROOT macro).                                                                               *
 *                                                                                                                           *
 * Author: Miguel Del Ben Galdiano (miguel.galdiano@gmail.com).                                                              *
 * Date of creation : June 18 2020.                                                                                          *
 *                                                                                                                           *
 * Summary of File:                                                                                                          *
 *                                                                                                                           *
 *    This macro should contain all functions that relate to the time of flight between the SiTel and TPC. The purpouse of   *
 *    each function will be properly explained in a comment located right above the function. Some concepts are used by all  *
 *    functions, and as such will be explained here.                                                                         *
 *                                                                                                                           *
 *    > TIME OF FLIGHT PARAMETER >                                                                                           *
 *      Definition: 2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time).                                 *
 *                                                                                                                           *
 *      The parameters start_time[]/cdf_time are used instead of xmin[]/min_x because the later are quantized in 2ns steps,  *
 *      while the former are continuous. This allows for a better histogram. Also, instead of using only the timing of the   *
 *      thinner silicon detector (start_time[30]), we used a modified avarage of the thinner and larger silicon detectors.   *
 *      The value subtracted from start_time[31] is obtained from constructing a histogram of the difference between the two *
 *      start_time[] values and fitting it with a gaussian function and taking its mean. This subtraction is made to match   *
 *      timing of the larger silicon detector, that happens later, with that of the smaller silicon detector.                *
 *                                                                                                                           *
 * ************************************************************************************************************************* */

using namespace std;

void ToFHistogram( int run ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> cd();
  sidStyle -> SetOptStat(0);

  TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
  TTree* reco;
  TFile* file = CheckFile(file_name);  file -> GetObject("reco", reco);

  Double_t height = 500;    Double_t width = gdRatio * height;
  TCanvas* canvas = new TCanvas("canvas", "Time of Flight (SiTel-TPC)", width, height);
  canvas -> SetGrid();

  Int_t cfg = 2;  // Double Phase = 2; Single Phase = 1

  Double_t f90Min   = 0.1;    Double_t f90Max   = 1.0;
  Double_t f90erMin = 0.2;    Double_t f90erMax = 0.4;
  Double_t f90nrMin = 0.4;    Double_t f90nrMax = 0.6;
  Double_t s1Min    = 50;     Double_t s1Max    = 1E3;

  TCut generalCut = DefineCuts(2, f90Min,   f90Max,   s1Min, s1Max);
  TCut nrCut      = DefineCuts(2, f90nrMin, f90nrMax, s1Min, s1Max);
  TCut erCut      = DefineCuts(2, f90erMin, f90erMax, s1Min, s1Max);
  TCutG* lowBeCut = LowBeGraphCut(run, "LowBeCut");

  Double_t tofMin = 0;    Double_t tofMax = 50;
  Double_t tofBinSize   = 0.25;
  Double_t tofBinNumber = (tofMax - tofMin)/tofBinSize;

  const Int_t numOfHists = 3;

  TH1F*       tofHistograms[numOfHists];
  //TH1F*       tofLogHistograms[numOfHists];
  TCut        cuts[numOfHists]      = { erCut && "LowBeCut", nrCut && "LowBeCut", generalCut && "LowBeCut"};
  const char* histNames[numOfHists] = {"tof", "tof_nr", "tof_er"};
  EColor      colors[numOfHists] = {kBlue, kRed, kBlack};

  TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
  const char* legendLabel[numOfHists] = {"ER", "NR", "Ambos"};

  for (Int_t i = 0; i < numOfHists; i++){

    tofHistograms[i] = new TH1F(histNames[i], "", tofBinNumber, tofMin, tofMax);
    reco -> Project(histNames[i], "2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time)", cuts[i]);

    tofHistograms[i] -> SetAxisRange(0.1, 900., "Y");
    tofHistograms[i] -> SetLineWidth(2);
    tofHistograms[i] -> SetLineColor(colors[i]);

    tofHistograms[i] -> GetXaxis() -> SetTitle("ToF (SiTel-TPC) [ns]");
    tofHistograms[i] -> GetYaxis() -> SetTitle("Eventos / 0.25 ns");

    legend -> AddEntry(tofHistograms[i], legendLabel[i]);

    canvas -> cd();
    if (i == 0) {
      tofHistograms[i] -> Draw("HIST");

    } else {
      tofHistograms[i] -> Draw("SAME HIST");
    }

  }
  legend -> Draw();
}

void F90vsToF( int run ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> cd();
  sidStyle -> SetOptStat(0);

  TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);

  Int_t expCfg = 2;
  Double_t s1Min = 40;    Double_t s1Max = 1000;
  Double_t f90Min = 0.;   Double_t f90Max = 1.0;
  Double_t tofMin = 0;   Double_t tofMax = 50;
  Double_t tofBinSize = 0.125;
  Double_t tofBinNumber = (tofMax - tofMin)/tofBinSize;

  Double_t f90BinSize = 5E-3;
  Double_t f90BinNumber = (f90Max - f90Min)/f90BinSize;

  TCut   tpcCut   = DefineCuts(expCfg, f90Min, f90Max, s1Min, s1Max);
  TCutG* lowBeCut = LowBeGraphCut(run, "lowBeCut");
  TCut   totalCut = tpcCut && "lowBeCut";

  TH2F* ToFvS1Hist   = new TH2F("ToFvS1_hist", "", tofBinNumber, tofMin, tofMax, f90BinNumber, f90Min, f90Max);

  reco -> Project("ToFvS1_hist", "clusters[0].f90:2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time)", totalCut);

  Double_t height = 500; Double_t width = gdRatio * height;
  TCanvas* canvas = new TCanvas("canvas", "canvas", width, height);
  canvas -> SetGrid();

  ToFvS1Hist -> GetXaxis() -> SetTitle("ToF (SiTel-TPC) [ns]");
  ToFvS1Hist -> GetYaxis() -> SetTitle("f_{prompt}");
  ToFvS1Hist -> Draw("COLZ");

  file -> Close();
}


void ToFHistLSciCoincidence( int run ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> cd();

  TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
  TTree* reco;
  TFile* file = CheckFile(file_name);  file -> GetObject("reco", reco);

  Double_t f90Min     = 0.1;    Double_t f90Max     = 1.0;
  Double_t s1Min      = 50;     Double_t s1Max      = 1E3;
  Double_t spLsciToFMin = -40;    Double_t spLsciToFMax = -20;
  Double_t lpLsciToFMin = 550;    Double_t lpLsciToFMax = 600;

  Double_t tofMin = -100;    Double_t tofMax = 100;
  Double_t tofBinSize   = 0.25;
  Double_t tofBinNumber = (tofMax - tofMin)/tofBinSize;

  const Int_t numOfHists = 6;

  Int_t chanIDs[numOfHists] = {0, 1, 3, 4, 5, 8};

  TCut tpcCut = DefineCuts(2, f90Min, f90Max, s1Min, s1Max);
  TCutG* lowBeCut = LowBeGraphCut(run, "LowBeCut");

  TCut spLsciToF[numOfHists];
  TCut lpLsciToF[numOfHists];
  for (Int_t i = 0; i < numOfHists; i++){
    spLsciToF[i] = Form("2*(0.5*(start_time[30] + start_time[31] - 7.45) - start_time[%d]) >= %f && 2*(0.5*(start_time[30] + start_time[31] - 7.45) - start_time[%d]) <= %f", chanIDs[i], spLsciToFMin, chanIDs[i], spLsciToFMax);
    lpLsciToF[i] = Form("2*(0.5*(start_time[30] + start_time[31] - 7.45) - start_time[%d]) >= %f && 2*(0.5*(start_time[30] + start_time[31] - 7.45) - start_time[%d]) <= %f", chanIDs[i], lpLsciToFMin, chanIDs[i], lpLsciToFMax);
  }

  TH1F*       spToFHistograms[numOfHists];
  TH1F*       lpToFHistograms[numOfHists];
  const char* spHistNames[numOfHists] = {"LSci0_sp", "LSci1_sp", "LSci3&7_sp", "LSci4&6_sp", "LSci5_sp", "LSci8_sp"};
  const char* lpHistNames[numOfHists] = {"LSci0_lp", "LSci1_lp", "LSci3&7_lp", "LSci4&6_lp", "LSci5_lp", "LSci8_lp"};
    Int_t       colors[numOfHists] = {1179, 1230, 1281, 1332, 1383, 1433};

  THStack*    spStackedToF = new THStack("tof_stack_sp", "Time of Flight (SiTel-TPC)");
  THStack*    lpStackedToF = new THStack("tof_stack_lp", "Time of Flight (SiTel-TPC)");

  for (Int_t i = 0; i < numOfHists; i++){

    spToFHistograms[i] = new TH1F(spHistNames[i], "", tofBinNumber, tofMin, tofMax);
    lpToFHistograms[i] = new TH1F(lpHistNames[i], "", tofBinNumber, tofMin, tofMax);

    reco -> Project(spHistNames[i], "2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time)", spLsciToF[i] && tpcCut && "LowBeCut");
    reco -> Project(lpHistNames[i], "2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time)", lpLsciToF[i] && tpcCut && "LowBeCut");

    spToFHistograms[i] -> SetFillColor(colors[i]);    spToFHistograms[i] -> SetLineColor(colors[i]);
    lpToFHistograms[i] -> SetFillColor(colors[i]);    lpToFHistograms[i] -> SetLineColor(colors[i]);
    spStackedToF -> Add(spToFHistograms[i]);
    lpStackedToF -> Add(lpToFHistograms[i]);

  }

  Double_t height = 500;    Double_t width = gdRatio * height;
  TCanvas* canvas = new TCanvas("canvas", "Time of Flight (SiTel-TPC)", 2*width, height);
  canvas -> Divide(2,1);

  canvas -> cd(1);    spStackedToF -> Draw();
  spStackedToF -> GetXaxis() -> SetTitle("ToF (SiTel-TPC) [ns]");
  spStackedToF -> GetYaxis() -> SetTitle("Counts/0.25 ns");

  canvas -> cd(2);    lpStackedToF -> Draw();
  lpStackedToF -> GetXaxis() -> SetTitle("ToF (SiTel-TPC) [ns]");
  lpStackedToF -> GetYaxis() -> SetTitle("Counts/0.25 ns");


  TLegend*legend = new TLegend(0.85,0.74,0.95,0.95);
  legend -> SetTextFont(102);
  legend -> SetTextSize(0.03);
  for (Int_t i = 0; i < numOfHists; i++){
    legend -> AddEntry(spToFHistograms[i], spHistNames[i], "l");
  }
  legend -> Draw();

}
