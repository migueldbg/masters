#include "sidutility.cc"

#include <TFile.h>
#include <TH2.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>


/* ************************************************************************************************************************* *
 * File: AmRun.C (ROOT macro).                                                                                               *
 *                                                                                                                           *
 * Author: Miguel Del Ben Galdiano (miguel.galdiano@gmail.com).                                                              *
 * Date of creation : May 8 2020.                                                                                            *
 *                                                                                                                           *
 * Summary of File:                                                                                                          *
 *                                                                                                                           *
 *    This macro uses a Am run to define the 'Electron Recoil Region' on the F90 vs S1 2D histogram. This region will then   *
 *    be used to quantify how effective a given set of cuts is at removing ER events while preserving NR events. To define   *
 *    the region, the 'CreateF90Histograms.C' macro will be used to create f90 histograms for different S1 bins. For each    *
 *    histogram a gaussian fit will be applied to obtain the necessary parameters such as to construct the diferent regions  *
 *    of interset (68%, 95%, etc). The region will be constructed as two curves, one defining the lower boundary and one the *
 *    upper boundary. The macro then saves these curves.                                                                     *                *
 *                                                                                                                           *
 *    Am Runs: 1480.                                                                                                         *
 * ************************************************************************************************************************* */


Int_t    exp_cfg = 2;
Double_t f90_min = 0.1;
Double_t f90_max = 1.0;
Double_t s1_min  = 0.0;
Double_t s1_max  = 1000.0;


void F90vS1( Int_t run ){

  TString file_name = Form("runs/run_%d.root", run);
  TFile* file = CheckFile(file_name);

  TTree* reco;  file -> GetObject("reco", reco);

  TH2F* f90vS1Hist = new TH2F("f90vS1_hist", "f90vS1_hist", 100, s1_min, s1_max, 100, f90_min, f90_max);
  TCut histogram_cuts = DefineCuts(exp_cfg, f90_min, f90_max, s1_min, s1_max);

  reco -> Project("f90vS1_hist", "clusters[0].f90:clusters[0].charge", histogram_cuts);

  f90vS1Hist -> DrawCopy("COLZ");
}

void ERRegion( int run ){

  TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);

  const Int_t nq = 2;
  Double_t xq[nq] = {0.03, 0.97};
  Double_t yq[nq];

  Int_t binsNumber = 50;
  Double_t binSize = (s1_max - s1_min)/binsNumber;
  Double_t s1LowBound;
  Double_t s1UppBound;
  Double_t s1Center;

  Double_t qMin[binsNumber];  TGraph* gr05 = new TGraph(binsNumber);
  Double_t qMax[binsNumber];  TGraph* gr95 = new TGraph(binsNumber);

  TH1F* f90Hist = new TH1F("f90Hist", "f90 Histogram", 100, f90_min, f90_max);
  TCut histogramCuts;

  for (Int_t i = 0; i < binsNumber; i++){

    s1LowBound = i*binSize;
    s1UppBound = (i + 1)*binSize;
    s1Center =  (s1LowBound + s1UppBound)/2;

    histogramCuts = DefineCuts(exp_cfg, f90_min, f90_max, s1LowBound, s1UppBound);

    reco -> Project("f90Hist", "clusters[0].f90", histogramCuts);
    f90Hist -> GetQuantiles(nq, yq, xq);

    gr05 -> SetPoint(i, s1Center, yq[0]);
    gr95 -> SetPoint(i, s1Center, yq[1]);

  }
  histogramCuts = DefineCuts(exp_cfg, 0, 1., s1_min, s1_max);
  TCanvas* c1 = new TCanvas("c1", "test", 900, 600);
  gPad -> DrawFrame(0,0, s1_max, 1.0);
  reco -> Draw("clusters[0].f90:clusters[0].charge >> hist(100,0,1000,100,0,1)", histogramCuts,  "colz");
  gr05 -> Draw("same");
  gr95 -> Draw("same");
}
