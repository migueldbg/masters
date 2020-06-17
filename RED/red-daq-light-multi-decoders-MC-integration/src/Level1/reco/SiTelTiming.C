#include "sidutility.cc"

#include <TFile.h>
#include <TH2.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <vector>

Int_t    exp_cfg = 2;
Double_t f90Min = 0.0;
Double_t f90Max = 1.0;
Double_t s1Min  = 50.0;
Double_t s1Max  = 1000.0;
Double_t s2Min  = 0.0;
Double_t s2Max  = 0.0;

const double gdRatio = 1.61803398875;

void SiTelTiming (int run){

  gStyle -> SetLabelFont(102, "xyz");
  gStyle -> SetTitleFont(102, "xyz");
  gStyle -> SetTitleFont(102, "t");
  gStyle -> SetOptStat(0);


  TString file_name = Form("runs/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);

  Double_t tofMin = 0.0;  Double_t tofMax = 40.0;
  Double_t binSize = 0.25;
  Double_t binNumber = (tofMax - tofMin)/binSize;

  TCut   tpcCut   = DefineCuts(exp_cfg, f90Min, f90Max, s1Min, s1Max);
  TCutG* lowBeCut = LowBeGraphCut(run, "LowBeCut");
  TCut   totalCut = tpcCut && "LowBeCut";

  TH1F* siTeldTHist      = new TH1F("SiTeldTHist",   "SiTel Time Difference", binNumber, tofMin, tofMax);
  TH1F* siTeldTHistBe    = new TH1F("SiTeldTHistBe", "", binNumber, tofMin, tofMax);
  TH1F* siTeldTLogHist   = new TH1F("SiTeldTLogHist",   "", binNumber, tofMin, tofMax);
  TH1F* siTeldTLogHistBe = new TH1F("SiTeldTLogHistBe", "", binNumber, tofMin, tofMax);

  reco -> Project("SiTeldTHist", "2*(start_time[31] - start_time[30])", tpcCut);
  reco -> Project("SiTeldTHistBe", "2*(start_time[31] - start_time[30])", totalCut);

  siTeldTHist   -> SetLineColor(kBlue);
  siTeldTHist   -> SetLineWidth(2);
  siTeldTHistBe -> SetLineColor(kRed);
  siTeldTHistBe -> SetLineWidth(2);

  siTeldTLogHist   = (TH1F*) siTeldTHist -> Clone();
  siTeldTLogHistBe = (TH1F*) siTeldTHistBe -> Clone();

  Double_t height = 500; Double_t width = gdRatio * height ;
  TCanvas* canvas = new TCanvas("canvas", "canvas", width, height);

  siTeldTHist -> GetXaxis() -> SetTitle("dT (E - #DeltaE)[ns]");
  siTeldTHist -> GetYaxis() -> SetTitle("events/0.25 ns");
  siTeldTHist -> Draw();
  siTeldTHistBe -> Draw("SAME");

  Double_t pointX = 0.5;    Double_t pointY = 0.3;    Double_t padRatio = 0.4;

  TPad* logPad = new TPad("logPad", "", pointX, pointY, pointX + padRatio, pointY + padRatio);
  logPad -> SetLogy();   logPad -> Draw();   logPad -> cd();

  siTeldTLogHist -> SetTitle("");
  siTeldTLogHist -> Draw();
  siTeldTLogHistBe -> Draw("SAME");

  canvas -> cd();

  TLegend* legend = new TLegend(0.85, 0.74, 0.95, 0.95);
  legend -> SetTextFont(102);
  legend -> SetTextSize(0.03);
  legend -> AddEntry(siTeldTHist,   "TPC Cuts",          "l");
  legend -> AddEntry(siTeldTHistBe, "TPC Cuts + Low Be", "l");
  legend -> Draw();

  Double_t entries = siTeldTHist -> GetEntries();
  Double_t entriesBe = siTeldTHistBe -> GetEntries();

  std::cout << "Nº Events (TPC Cuts): " << entries << std::endl;
  std::cout << "Nº Events (TPC Cuts + Low Be): " << entriesBe << std::endl;

  file -> Close();
  // Distribution peak position = 14.89 us.
  // Difference between SiTel components (in samples) = 7.45
}

void TimingDiff( int run ){

  gStyle -> SetLabelFont(102, "xyz");
  gStyle -> SetTitleFont(102, "xyz");
  gStyle -> SetTitleFont(102, "t");
  gStyle -> SetOptStat(0);

  TString file_name = Form("runs/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);

  Double_t tofMin = 0.0;  Double_t tofMax = 50.0;
  Double_t binSize = 0.25;
  Double_t binNumber = (tofMax - tofMin)/binSize;

  TCut   tpcCut   = DefineCuts(exp_cfg, f90Min, f90Max, s1Min, s1Max);
  TCutG* lowBeCut = LowBeGraphCut(run);
  TCut   totalCut = tpcCut && "lowBeCut";

  TH1F* tofHistdE  = new TH1F("dEHist",    "", binNumber, tofMin, tofMax);
  TH1F* tofHistSum = new TH1F("sumHist",   "", binNumber, tofMin, tofMax);

  reco -> Project("dEHist",    "2*(start_time[30] - clusters[0].cdf_time)", totalCut);
  reco -> Project("sumHist",   "2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time)", totalCut);

  tofHistdE -> SetTitle("TPC ToF Comparison");
  tofHistdE -> GetXaxis() -> SetTitle("ToF [ns]");
  tofHistdE -> GetYaxis() -> SetTitle("Counts/0.25 ns");

  tofHistdE  -> SetLineColor(kBlue);   tofHistdE  -> SetLineWidth(2);
  tofHistSum -> SetLineColor(kRed);    tofHistSum -> SetLineWidth(2);

  tofHistdE  -> Draw("HIST");
  tofHistSum -> Draw("SAME");

  TLegend *legend = new TLegend(0.85,0.74,0.95,0.95);
  legend -> SetTextFont(102);
  legend -> SetTextSize(0.05);
  legend -> AddEntry(tofHistdE,  "#DeltaE",     "l");
  legend -> AddEntry(tofHistSum, "#DeltaE + E", "l");
  legend -> Draw();

  file -> Close();
  /* Mean Values  >> dE = 24.381 +/- 0.016
   *                 dE + E = 24.377 +/- 0.014
   * Sigma Values >> dE = 1.369 +/- 0.015 (5.6%)
   *                 dE + E = 1.158 +/- 0.014 (4.8%)
   *
   * Values were obtained selecting the low energy Be.
  */
}
