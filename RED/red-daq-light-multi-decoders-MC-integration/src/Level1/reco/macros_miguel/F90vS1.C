#include "sidutility.cc"

#include <TFile.h>
#include <TH2.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

/* ************************************************************************************************************************* *
 * File: TimeOfFlight_SiTelTPC.C (ROOT macro).                                                                               *
 *                                                                                                                           *
 * Author: Miguel Del Ben Galdiano (miguel.galdiano@gmail.com).                                                              *
 * Date of creation : June 29 2020.                                                                                          *
 *                                                                                                                           *
 * Summary of File:                                                                                                          *
 *                                                                                                                           *
 *    This macro should contain all functions that relate to the time of flight between the SiTel and LSci. The purpouse of  *
 *    each function will be properly explained in a comment located right above the function. Some concepts are used by all  *
 *    functions, and as such will be explained here.                                                                         *
 *                                                                                                                           *
 *    > TIME OF FLIGHT PARAMETER >                                                                                           *
 *      Definition: 2*(start_time[] - 0.5*(start_time[30] + start_time[31] - 7.45)).                                         *
 *                                                                                                                           *
 *      The parameter start_time[] is used instead of xmin[] because the later is quantized in 2ns steps, while the former   *
 *      is continuous. This allows for a better histogram. Also, instead of using only the timing of the thinner silicon     *
 *      detector (start_time[30]), we used a modified avarage of the thinner and larger silicon detectors. The value         *
 *      subtracted from start_time[31] is obtained from constructing a histogram of the difference between the two           *
 *      start_time[] values and fitting it with a gaussian function and taking its mean. This subtraction is made to match   *
 *      the timing of the larger silicon detector, that happens later, with that of the smaller silicon detector.            *
 *                                                                                                                           *
 *    > LSci CHANNEL MAPPING >                                                                                               *
 *      LSci:       0    1    3 & 7    4 & 6    5    8                                                                       *
 *      Channel:    0    1      3        4      5    8                                                                       *
 * ************************************************************************************************************************* */

using namespace std;

TCut DefineSiTelTPCToFCut(Double_t tofMin, Double_t tofMax);
TCut DefineSiTelLSciToFCut(Int_t chanID, Double_t tofMin, Double_t tofMax);


Int_t expCfg = 2;

Double_t f90Min = 0.0;    Double_t f90Max = 1.0;
Double_t s1Min  = 60.0;   Double_t s1Max  = 1000.0;

Double_t f90BinSize = 0.01;   Double_t f90BinNumber = (f90Max - f90Min)/f90BinSize;
Double_t s1BinSize  = 10;     Double_t s1BinNumber  = (s1Max - s1Min)/s1BinSize;


void F90vS1( Int_t run ){

  TStyle* sidStyle = SetSidStyle();
  //sidStyle -> SetOptStat(0);
  sidStyle -> cd();

  TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);

  // Defining the selection cuts to be applied:
  Double_t tpcToFMin  = 22;   Double_t tpcToFMax  = 27;
  Double_t lsciToFMin = 20;   Double_t lsciToFMax = 40;

  const Int_t chanTotal = 6;
  Int_t chanID[chanTotal]     = { 0   , 1  ,  3  ,  4  ,  5  ,  8};
  Int_t chanColors[chanTotal] = {1179, 1230, 1281, 1332, 1383, 1433};

  TCut lsciCuts[chanTotal];
  for (Int_t i = 0; i < chanTotal; i++){
    lsciCuts[i] = DefineSiTelLSciToFCut(chanID[i], lsciToFMin, lsciToFMax);
  }

  TCut   tpcToF     = DefineSiTelTPCToFCut(tpcToFMin, tpcToFMax);
  TCut   tpcCut     = DefineCuts(expCfg, f90Min, f90Max, s1Min, s1Max);
  TCutG* lowBeCut   = LowBeGraphCut(run, "LowBeCut");
  TCut   qualityCut = tpcCut && "LowBeCut";

  TCut generalCuts[2] = {qualityCut, qualityCut && tpcToF};

  // Histograms that will contain the results after the application of the selection cuts.
  TH2F* f90vS1Hist[2];
  TH2F* lsciF90vS1Hist[chanTotal];
  TString histNames[2]             = {"tpcCut", "tpc+tofCut"};
  TString lsciHistNames[chanTotal] = {"LSci 0", "LSci 1", "LSci 3 & 7", "LSci 4 & 6", "LSci 5", "LSci 8"};
  TString histTitles[3]            = {"f90 vs S1 (TPC Cuts)",
                                      "f90 vs S1 (TPC + SiTel-TPC ToF Cuts)",
                                      "f90 vs S1 (TPC + SiTel-TPC ToF + SiTel-LSci ToF Cuts)"};
  Int_t lsciEvents = 0;

  Double_t height = 400;    Double_t width = gdRatio * height;
  TCanvas* canvas = new TCanvas("canvas", "canvas", 3*width, height);
  canvas -> Divide(3,1);

  canvas -> cd(3);
  TLegend* legend = new TLegend(0.2, 0.2);

  for (Int_t i = 0; i < 2; i++){

    f90vS1Hist[i] = new TH2F(histNames[i], histTitles[i], s1BinNumber, s1Min, s1Max, f90BinNumber, f90Min, f90Max);
    reco -> Project(histNames[i], "clusters[0].f90:clusters[0].charge", generalCuts[i]);

    canvas -> cd(i + 1);
    f90vS1Hist[i] -> GetXaxis() -> SetTitle("S1 [PE]");
    f90vS1Hist[i] -> GetYaxis() -> SetTitle("f90");
    f90vS1Hist[i] -> Draw("COLZ");

    if (i == 1){
      for (Int_t j = 0; j < chanTotal; j++){
        lsciF90vS1Hist[j] = new TH2F(lsciHistNames[j], histTitles[i+1], s1BinNumber, s1Min, s1Max, f90BinNumber, f90Min, f90Max);
        reco -> Project(lsciHistNames[j], "clusters[0].f90:clusters[0].charge", lsciCuts[j] && generalCuts[i]);

        lsciF90vS1Hist[j] -> SetMarkerStyle(20);
        lsciF90vS1Hist[j] -> SetMarkerSize(1);
        lsciF90vS1Hist[j] -> SetMarkerColor(chanColors[j]);

        canvas -> cd(3);
        lsciF90vS1Hist[j] -> GetXaxis() -> SetTitle("S1 [PE]");
        lsciF90vS1Hist[j] -> GetYaxis() -> SetTitle("f90");
        lsciF90vS1Hist[j] -> Draw("SAME");

        legend -> AddEntry(lsciF90vS1Hist[j], lsciHistNames[j], "lep");

        lsciEvents = lsciEvents + (Int_t) lsciF90vS1Hist[j] -> GetEntries();
      }
      legend -> Draw();
    }
  }

  cout << "Number of events after LSci timing cut: " << lsciEvents << endl;

}


TCut DefineSiTelTPCToFCut(Double_t tofMin, Double_t tofMax){

  TCut tofMinCut = Form("2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time) >= %f", tofMin);
  TCut tofMaxCut = Form("2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time) <= %f", tofMax);
  TCut tofCut = tofMinCut && tofMaxCut;

  return tofCut;
}

TCut DefineSiTelLSciToFCut(Int_t chanID, Double_t tofMin, Double_t tofMax){

  TCut tofMinCut = Form("2*(start_time[%d] - 0.5*(start_time[30] + start_time[31] - 7.45)) >= %f", chanID, tofMin);
  TCut tofMaxCut = Form("2*(start_time[%d] - 0.5*(start_time[30] + start_time[31] - 7.45)) <= %f", chanID, tofMax);
  TCut tofCut = tofMinCut && tofMaxCut;

  return tofCut;
}
