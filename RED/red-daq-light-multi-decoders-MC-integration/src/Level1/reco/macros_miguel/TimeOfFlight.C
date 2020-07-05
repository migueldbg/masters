#include "sidutility.cc"

#include <TFile.h>
#include <TH2.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <vector>

using namespace std;

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


void F90vS1( Int_t run ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> cd();

  TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(file_name);

  TTree* reco;  file -> GetObject("reco", reco);


  Double_t tpcToFMin = 20;    Double_t tpcToFMax = 30;
  Double_t lsciToFMin = -40;  Double_t lsciToFMax = -20;

  TCut tpcToF = Form("2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time) >= %f && 2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time) <= %f", tpcToFMin, tpcToFMax);
  TCut lsciToF = Form("2*(0.5*(start_time[30] + start_time[31] - 7.45) - start_time[0]) >= %f && 2*(0.5*(start_time[30] + start_time[31] - 7.45) - start_time[0]) <= %f", lsciToFMin, lsciToFMax);
  TCut qualityCut = DefineCuts(expCfg, f90Min, f90Max, s1Min, s1Max);
  TCutG* lowBeCut = LowBeGraphCut(run, "LowBeCut");
  TCut tpcCut = qualityCut && "LowBeCut";

  TH2F* f90vS1Hist[3];
  const char* histNames[3] = {"tpc_cut", "tpc_cut+tpcToF", "tpc_cut+tpcToF+lsciToF"};
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

  TCut   tpcCut   = DefineCuts(expCfg, f90Min, f90Max, s1Min, s1Max);
  TCutG* lowBeCut = LowBeGraphCut(run, "lowBeCut");
  TCut   totalCut = tpcCut && "lowBeCut";

  TH2F* ToFvS1Hist   = new TH2F("ToFvS1_hist", "SiTel-TPC ToF vs S1", s1BinNumber, s1Min, s1Max, tofBinNumber, tofMin, tofMax);
  TH2F* ToFvS1HistBe = new TH2F("ToFvS1_histBe", "SiTel-TPC ToF vs S1 (Low Be)", s1BinNumber, s1Min, s1Max, tofBinNumber, tofMin, tofMax);

  reco -> Project("ToFvS1_hist", "2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time):clusters[0].charge", tpcCut);
  reco -> Project("ToFvS1_histBe", "2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time):clusters[0].charge", totalCut);

  Double_t height = 500; Double_t width = gdRatio * height;
  TCanvas* canvas = new TCanvas("canvas", "canvas", 2*width, height);
  canvas -> Divide(2,1);

  canvas -> cd(1);
  ToFvS1Hist -> GetXaxis() -> SetTitle("S1[PE]");
  ToFvS1Hist -> GetYaxis() -> SetTitle("ToF (SiTel - TPC)[ns]");
  ToFvS1Hist -> Draw("COLZ");

  canvas -> cd(2);
  ToFvS1HistBe -> GetXaxis() -> SetTitle("S1[PE]");
  ToFvS1HistBe ->  Draw("COLZ");

  //ToFvS1HistBe -> ProjectionY("_py", 20, 22) -> Draw();

  std::cout << ToFvS1Hist -> GetNbinsX() << std::endl;
  std::cout << ToFvS1HistBe -> GetNbinsX() << std::endl;

  file -> Close();
}


void TPCLSciToF( Int_t run, Int_t chanID ){

  TStyle* sidStyle = SetSidStyle();   sidStyle -> cd();

  TString file_name = Form("runs/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);

  Int_t cfg = 2;
  Double_t s1Min = 50;    Double_t s1Max = 1000;
  Double_t f90Min = 0.1;    Double_t f90Max = 1.0;
  Double_t tofMin = 300;   Double_t tofMax = 700;

  Double_t erRange = 0.3;

  Double_t tofBinSize = 1;
  Double_t tofBinNumber = (tofMax - tofMin)/tofBinSize;

  TCut tpcCut = DefineCuts(cfg, f90Min, f90Max, s1Min, s1Max);
  TCut tpcCutER = DefineCuts(cfg, f90Min, f90Min + erRange, s1Min, s1Max);
  TCut tpcCutNR = DefineCuts(cfg, f90Min + erRange, f90Max, s1Min, s1Max);
  TCutG* lowBeCut = LowBeGraphCut(run, "LowBeCut");

  TH1F* tpcLSciTof[2];   TH1F* tpcLSciTofER[2];   TH1F* tpcLSciTofNR[2];

  tpcLSciTof[0]   = new TH1F("tpcLSciTof", "No Be Selection", tofBinNumber, tofMin, tofMax);
  tpcLSciTofER[0] = new TH1F("tpcLSciTofER", "", tofBinNumber, tofMin, tofMax);
  tpcLSciTofNR[0] = new TH1F("tpcLSciTofNR", "", tofBinNumber, tofMin, tofMax);

  tpcLSciTof[1]   = new TH1F("tpcLSciTof_LowBe", "Low Be Selection", tofBinNumber, tofMin, tofMax);
  tpcLSciTofER[1] = new TH1F("tpcLSciTofER_LowBe", "", tofBinNumber, tofMin, tofMax);
  tpcLSciTofNR[1] = new TH1F("tpcLSciTofNR_LowBe", "", tofBinNumber, tofMin, tofMax);


  const char* histName[6] = {"tpcLSciTof", "tpcLSciTofER", "tpcLSciTofNR",
                             "tpcLSciTof_LowBe", "tpcLSciTofER_LowBe", "tpcLSciTofNR_LowBe"};

  TCut cut[6] = {tpcCut, tpcCutER, tpcCutNR, tpcCut && "LowBeCut", tpcCutER && "LowBeCut", tpcCutNR && "LowBeCut"};

  for (Int_t i = 0; i < 6; i++){
    reco -> Project(histName[i], Form("(clusters[0].cdf_time - start_time[%d])*2.", chanID), cut[i]);
  }

  for (Int_t i = 0; i < 2; i++){

    tpcLSciTof[i]   -> SetLineWidth(2);
    tpcLSciTofER[i] -> SetLineWidth(2);    tpcLSciTofER[i] -> SetLineColor(kRed);
    tpcLSciTofNR[i] -> SetLineWidth(2);    tpcLSciTofNR[i] -> SetLineColor(kBlue);

    tpcLSciTof[i] -> GetXaxis() -> SetTitle("ToF (TPC-LSci) [ns]");
    tpcLSciTof[i] -> GetYaxis() -> SetTitle("Counts/1 ns");
  }

  Double_t height = 500; Double_t width = gdRatio * height;
  TCanvas* canvas = new TCanvas("canvas", "Time of Flight (TPC - LSci)", 2*width, height);
  canvas -> Divide(2,1);

  TLegend*legend = new TLegend(0.85,0.74,0.95,0.95);
  legend -> SetTextFont(102);
  legend -> SetTextSize(0.03);
  legend -> AddEntry(tpcLSciTof[1], Form("%2.1f #leq f90 #geq %2.1f", f90Min, f90Max), "l");
  legend -> AddEntry(tpcLSciTofER[1], Form("%2.1f #leq f90 #geq %2.1f", f90Min, f90Min + erRange), "l");
  legend -> AddEntry(tpcLSciTofNR[1], Form("%2.1f #leq f90 #geq %2.1f", f90Min + erRange, f90Max), "l");

  for (Int_t i = 0; i < 2; i++){
    canvas -> cd(i+1);
    tpcLSciTof[i] -> Draw("HIST");
    tpcLSciTofER[i] -> Draw("HIST SAME");
    tpcLSciTofNR[i] -> Draw("HIST SAME");
    legend -> Draw();
  }

  Double_t lowBound = 500;    Double_t uppBound = 560;
  Int_t lowBin = (Int_t)(lowBound - tofMin)/tofBinSize;
  Int_t uppBin = (Int_t)(uppBound - tofMin)/tofBinSize;

  Double_t peak[2];   Double_t erPeak[2];   Double_t nrPeak[2];
  for (Int_t i = 0; i < 2; i++){
    peak[i]   = tpcLSciTof[i] -> Integral(lowBin, uppBin);
    erPeak[i] = tpcLSciTofER[i] -> Integral(lowBin, uppBin);
    nrPeak[i] = tpcLSciTofNR[i] -> Integral(lowBin, uppBin);

    cout << "Nº Events at Peak Region (Total): " << peak[i] << endl;
    cout << "Nº Events at Peak Region (Nuclear Recoil): " << nrPeak[i] << " (" << (nrPeak[i]/peak[i])*100 << "%)" << endl;
    cout << "Nº Events at Peak Region (Electron Recoil): " << erPeak[i] << " (" << (erPeak[i]/peak[i])*100 << "%)" << endl;
  }
}

void TPCLSciToFComparison( Int_t run ){

  TStyle* sidStyle = SetSidStyle();   sidStyle -> cd();

  TString file_name = Form("runs/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);

  Int_t cfg = 2;
  Double_t s1Min = 50;    Double_t s1Max = 1000;
  Double_t f90Min = 0.1;    Double_t f90Max = 1.0;
  Double_t tofMin = 300;   Double_t tofMax = 700;

  Double_t tofBinSize = 0.5;
  Double_t tofBinNumber = (tofMax - tofMin)/tofBinSize;

  TCut tpcCut = DefineCuts(cfg, f90Min, f90Max, s1Min, s1Max);
  TCutG* lowBeCut = LowBeGraphCut(run, "LowBeCut");

  Double_t height = 500; Double_t width = gdRatio * height;
  TCanvas* canvas = new TCanvas("canvas", "Time of Flight (TPC - LSci)", width, height);

  TH1F* tpcLSciToF[6];

  Int_t chanIDs[6] = {0, 1, 3, 4, 5, 8};
  const char* histNames[6] = {"LSci_0", "LSci_1", "LSci_3&LSci_7", "LSci_4&LSci_6", "LSci_5", "LSci_8"};

  for (Int_t i = 0; i < 6; i++){
    tpcLSciToF[i] = new TH1F(histNames[i], "", tofBinNumber, tofMin, tofMax);
    reco -> Project(histNames[i], Form("(clusters[0].cdf_time - start_time[%d])*2.", chanIDs[i]), tpcCut && "LowBeCut");
    tpcLSciToF[i] -> SetLineWidth(2);
    tpcLSciToF[i] -> SetLineColor(2+i);

    if (i == 0) {
      tpcLSciToF[i] -> SetTitle("ToF TPC-LSci (for all Scintilators);ToF (TPC-LSci) [ns];Counts/1 ns");
      tpcLSciToF[i] -> Draw("HIST");
    } else {
      tpcLSciToF[i] -> Draw("HIST SAME");
    }
  }

}


void SiTelLSciToF( Int_t run, Int_t chanID ){

  TStyle* sidStyle = SetSidStyle();   sidStyle -> cd();

  TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);

  Int_t cfg = 2;
  Double_t s1Min = 50;    Double_t s1Max = 1000;
  Double_t f90Min = 0.1;    Double_t f90Max = 1.0;
  Double_t tofMin = -1000;   Double_t tofMax = 1000 ;


  Double_t tofBinSize = 0.5;
  Double_t tofBinNumber = (tofMax - tofMin)/tofBinSize;

  TCut tpcCut = DefineCuts(cfg, f90Min, f90Max, s1Min, s1Max);
  TCutG* lowBeCut = LowBeGraphCut(run, "LowBeCut");

  TH1F* siTelLSciTof[3];
  const char* histNames[3] = {"noCut", "tpcCut", "tpc+lowBeCut"};
  TCut cuts[3] = {"", tpcCut, tpcCut && "LowBeCut"};
  EColor colors[3] = {kBlack, kBlue, kRed};

  for (Int_t i = 0; i < 3; i++){

    siTelLSciTof[i] = new TH1F(histNames[i], "Time of Flight (SiTel-LSci)", tofBinNumber, tofMin, tofMax);
    reco -> Project(histNames[i], Form("2*(0.5*(start_time[30] + start_time[31] - 7.45) - start_time[%d])", chanID), cuts[i]);

    siTelLSciTof[i] -> SetLineColor(colors[i]);
    siTelLSciTof[i] -> SetLineWidth(2);

    if (i == 0) {
      siTelLSciTof[i] -> GetXaxis() -> SetTitle("ToF (SiTel-LSci) [ns]");
      siTelLSciTof[i] -> GetXaxis() -> SetTitle("Counts/0.5 ns");
      siTelLSciTof[i] -> Draw("HIST");
    } else {
      siTelLSciTof[i] -> Draw("HIST SAME");
    }

  }

}

void SiTelLSciToFNeutron( Int_t run, Int_t chanID ){

  TStyle* sidStyle = SetSidStyle();   sidStyle -> cd();

  TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);


  Double_t psdMin      = 0.0;   Double_t psdMax      = 1.0;
  Double_t lsChargeMin = 0;     Double_t lsChargeMax = 200e3;
  Double_t tofMin      = -200;  Double_t tofMax      = 800 ;

  Double_t tofBinSize = 0.5;
  Double_t tofBinNumber = (tofMax - tofMin)/tofBinSize;


  TH1F* siTelLSciTof[3];
  const char* histNames[3] = {"noCut", "neutronCut", "gammaCut"};
  TCut cuts[3] = {"", Form("lsci_psd_tot[%d] >= 0.15 && lsci_psd_tot[%d] <= 1.0 && charge[%d] > 400", chanID, chanID, chanID),
                      Form("lsci_psd_tot[%d] >= 0.0 && lsci_psd_tot[%d] <= 0.15 && charge[%d] > 400", chanID, chanID, chanID)};
  EColor colors[3] = {kBlack, kBlue, kRed};

  for (Int_t i = 0; i < 3; i++){

    siTelLSciTof[i] = new TH1F(histNames[i], "Time of Flight (SiTel-LSci)", tofBinNumber, tofMin, tofMax);
    reco -> Project(histNames[i], Form("2*(0.5*(start_time[30] + start_time[31] - 7.45) - start_time[%d])", chanID), cuts[i]);

    siTelLSciTof[i] -> SetLineColor(colors[i]);
    siTelLSciTof[i] -> SetLineWidth(2);

    if (i == 0) {
      siTelLSciTof[i] -> GetXaxis() -> SetTitle("ToF (SiTel-LSci) [ns]");
      siTelLSciTof[i] -> GetXaxis() -> SetTitle("Counts/0.5 ns");
      siTelLSciTof[i] -> Draw("HIST");
    } else {
      siTelLSciTof[i] -> Draw("HIST SAME");
    }
  }

  TLegend*legend = new TLegend(0.85,0.74,0.95,0.95);
  legend -> SetTextFont(102);
  legend -> SetTextSize(0.03);
  legend -> AddEntry(siTelLSciTof[0], "No Cuts", "l");
  legend -> AddEntry(siTelLSciTof[1], "NR Selection", "l");
  legend -> AddEntry(siTelLSciTof[2], "ER Selection", "l");
  legend -> Draw();

}

void LSciPSD( Int_t run, Int_t chanID ){

  TStyle* sidStyle = SetSidStyle();   sidStyle -> cd();

  TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);


  Double_t psdMin       = 0;   Double_t psdMax       = 1.0;
  //Double_t psdMinZ      = 0;   Double_t psdMaxZ      = 1.0;
  Double_t lsChargeMin  = 0;   Double_t lsChargeMax  = 200e3;
  Double_t lsChargeMinZ = 0;   Double_t lsChargeMaxZ = 2e3;

  Double_t psdBinSize       = 0.002;  Int_t psdBinNumber       = (psdMax - psdMin)/psdBinSize;
  Double_t lsChargeBinSize  = 200;    Int_t lsChargeBinNumber  = (lsChargeMax - lsChargeMin)/lsChargeBinSize;
  Double_t lsChargeBinSizeZ = 10;     Int_t lsChargeBinNumberZ = (lsChargeMaxZ - lsChargeMinZ)/lsChargeBinSizeZ;

  TCut lsChargeCut  = Form("charge[%d] >= %f && charge[%d] <= %f", chanID, lsChargeMin,  chanID, lsChargeMax);
  TCut lsChargeCutZ = Form("charge[%d] >= %f && charge[%d] <= %f", chanID, lsChargeMinZ, chanID, lsChargeMaxZ);
  TCut psdCut = Form("f90[%d] >= %f && f90[%d] <= %f", chanID, psdMin, chanID, psdMax);

  TH2F* PSDvCharge = new TH2F("PSDvCharge", "PSD vs Charge (LSci 0)",
                              lsChargeBinNumber, lsChargeMin, lsChargeMax, psdBinNumber, psdMin, psdMax);
  TH2F* PSDvChargeZ = new TH2F("PSDvChargeZ", "PSD vs Charge (LSci 0)",
                               lsChargeBinNumberZ, lsChargeMinZ, lsChargeMaxZ, psdBinNumber, psdMin, psdMax);

  PSDvCharge  -> GetXaxis() -> SetTitle("Charge");    PSDvCharge  -> GetYaxis() -> SetTitle("PSD");
  PSDvChargeZ -> GetXaxis() -> SetTitle("Charge");    PSDvChargeZ -> GetYaxis() -> SetTitle("PSD");

  reco -> Project("PSDvCharge",  Form("f90[%d]:charge[%d]", chanID, chanID), lsChargeCut  && psdCut);
  reco -> Project("PSDvChargeZ", Form("f90[%d]:charge[%d]", chanID, chanID), lsChargeCutZ && psdCut);


  Double_t height = 500; Double_t width = gdRatio * height;
  TCanvas* canvas = new TCanvas("canvas", "canvas", 2*width, height);
  canvas -> Divide(2,1);

  canvas -> cd(1);  PSDvCharge  -> Draw("HIST COLZ");
  canvas -> cd(2);  PSDvChargeZ -> Draw("HIST COLZ");
}
