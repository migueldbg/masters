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

/* Low Energy Berilium Blob Bi-Gaussian Fit Results (fitted to merget 1501 to 1521 runs)

    > Entries:               210184
    > Mean x:        2747 +/- 1.037
    > Mean y:      661.5 +/- 0.2445
    > Std Dev x:   461.3 +/- 0.7331
    > Std Dev y:   108.8 +/- 0.1729
    > Chi2/ndf:           6336/5352
    > Constant:  7.353e+05 +/- 3168
    > MeanX:        2499 +/- 0.5052
    > SigmaX:      108.2 +/- 0.4641
    > MeanY:       772.5 +/- 0.1581
    > SigmaY:      33.09 +/- 0.1434
    > Rho:      -0.4482 +/- 0.00423
 */

const double gdRatio = 1.61803398875;

TCut LowBeCut( Int_t run, Int_t num_bins, Int_t EMin, Int_t EMax, Int_t dEMin, Int_t dEMax ){

  TH2F* bananaHist = Bananator(run, num_bins, EMin, EMax, dEMin, dEMax);
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
  bananaHist -> SetOption("COLZ");
  bananaHist -> Draw("CONT Z LIST");
  E_low_L -> Draw("SAME");  E_upp_L -> Draw("SAME");  dE_low_L -> Draw("SAME");  dE_upp_L -> Draw("SAME");

  TCut E_low_cut  = Form("baseline_mean[31] - ymin[31] > %f", E_low);
  TCut E_upp_cut  = Form("baseline_mean[31] - ymin[31] < %f", E_upp);
  TCut dE_low_cut = Form("baseline_mean[30] - ymin[30] > %f", dE_low);
  TCut dE_upp_cut = Form("baseline_mean[30] - ymin[30] < %f", dE_upp);

  TCut lowBe_cut = E_low_cut && E_upp_cut && dE_low_cut && dE_upp_cut;
  return lowBe_cut;
}

void LowBeContourCut( Int_t run, Int_t contourNumber, Int_t contourSize ){

  Int_t num_bins = 400;
  Int_t EMin  = 1500;  Int_t EMax  = 4000;
  Int_t dEMin = 400;   Int_t dEMax = 1200;

  TH2F* bananaHist = Bananator(run, num_bins, EMin, EMax, dEMin, dEMax);

  Double_t goldRatio = 1.61803398875;
  Double_t height = 500;    Double_t width = goldRatio * height;

  TCanvas* canvas1 = new TCanvas("canvas1", "Delta E/E Contours", width, height);

  Double_t contours[contourNumber];

  for (Int_t i = 0; i < contourNumber; i++){
    contours[i] = (i+1)*contourSize;
  }

  bananaHist -> SetContour(contourNumber, contours);

  bananaHist -> Draw("CONT Z LIST");
  canvas1 -> Update();

  // Get Contours
  TObjArray* contourList    = (TObjArray*)gROOT -> GetListOfSpecials() -> FindObject("contours");
  TList*     contourLevel   = NULL;
  TGraph*    curve          = NULL;
  TGraph*    curveClone     = NULL;
  TGraph*    maxCurve       = NULL;

  Int_t contourTotal = 0;

  if ( contourList == NULL ){
    std::cout << "No contours were extracted!" << std::endl;
    contourTotal = 0;
    return 0;
  } else {
    contourTotal = contourList -> GetSize();
  }

  std::cout << "Total Contours:" << contourTotal << std::endl;

  TCanvas* canvas2 = new TCanvas("canvas2", "Delta E/E Contours", width, height);
  canvas2 -> Divide(2,1);
  TH2F* hr = new TH2F("hr", "Contours", 2, EMin, EMax, 2, dEMin, dEMax);

  canvas2 -> cd(1);
  hr->Draw();
  Double_t x0, y0, z0;
  Int_t maxEntries = 0;

  for(Int_t i = 0; i < contourTotal; i++){
    contourLevel = (TList*)contourList->At(i);
    z0 = contours[i];

    // Get first graph from list on curves on this level
    curve = (TGraph*) contourLevel -> First();
    for(Int_t j = 0; j < contourLevel -> GetSize(); j++){
      curve -> GetPoint(0, x0, y0);
      curve -> SetLineColor(2 + i );

      if ( curve -> GetN() > maxEntries ) {
        maxCurve = curve;
        maxEntries = maxCurve -> GetN();
      }

      // Draw clones of the graphs to avoid deletions in case the 1st pad is redrawn.
      curveClone = (TGraph*) curve -> Clone();
      curveClone -> Draw("C");

      curve = (TGraph*) contourLevel -> After(curve); // Get Next graph
    }
  }
  canvas2 -> cd(2);
  maxCurve -> Draw();
  canvas2 -> Update();

  std::cout << "Amount of entries in curve with largest number of entries:" << maxCurve -> GetN() << std::endl;
  gStyle->SetTitleW(0.);
  gStyle->SetTitleH(0.);

  TFile* output_file = new TFile("LowBeCut.root", "UPDATE");
  output_file -> WriteObject(maxCurve, Form("lowBecut_%d", run), "OverWrite");
  output_file -> Close();
}

RooEllipse* ConfidenceEllipse( Int_t run, Double_t multiplier = 1. ){

  TF2* bigaus = new TF2("bigaus", "bigaus", 2200, 2780, 680, 880);
  bigaus -> SetParameters(7e+05, 2500, 100, 770, 30, -0.4 );

  TH2F* bananaHist = Bananator(run, 400, 1500, 4000, 400, 1200);
  bananaHist -> Fit("bigaus", "R0");

  Double_t meanX = bigaus -> GetParameter("MeanX");   Double_t sigmaX = bigaus -> GetParameter("SigmaX");
  Double_t meanY = bigaus -> GetParameter("MeanY");   Double_t sigmaY = bigaus -> GetParameter("SigmaY");
  Double_t rho   = bigaus -> GetParameter("Rho");

  auto* confidenceEllipse = new RooEllipse(Form("%2.1fSigmaEllipse", multiplier), meanX, meanY, multiplier*sigmaX, multiplier*sigmaY, rho);

  return confidenceEllipse;
}

void PlotEllipses( Int_t run ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> cd();

  TH2F* bananaHist = Bananator(run, 400, 1500, 4000, 400, 1200);

  TGraph* oneSigma = (TGraph*) ConfidenceEllipse(run, 1);    oneSigma -> SetLineWidth(3);    oneSigma -> SetLineColor(kRed-4);
  TGraph* twoSigma = (TGraph*) ConfidenceEllipse(run, 2);    twoSigma -> SetLineWidth(3);    twoSigma -> SetLineColor(kRed-4);

  Double_t height = 500; Double_t width = gdRatio * height;
  TCanvas* c1 = new TCanvas("c1", "Banana Histogram", width, height);
  c1 -> cd();

  bananaHist -> Draw("COLZ");   oneSigma -> Draw("SAME");  twoSigma -> Draw("SAME");

  oneSigma -> SetTitle("1#sigma Confidence Ellipse");
  twoSigma -> SetTitle("2#sigma Confidence Ellipse");

  TFile* output_file = new TFile("LowBeCut.root", "UPDATE");
  output_file -> WriteObject(oneSigma, Form("LowBeOneSigmaCut_%d", run), "OverWrite");
  output_file -> WriteObject(twoSigma, Form("LowBeTwoSigmaCut_%d", run), "OverWrite");
  output_file -> Close();
}

void TimeOfFlight(Int_t run){

  TStyle* sidStyle = SetSidStyle();   sidStyle -> cd();

  TString runFile_name = Form("runs/run_%d.root", run);  TFile* runFile = CheckFile(runFile_name);
  TTree* reco;    runFile -> GetObject("reco", reco);

  bool normalize = false;
  Int_t cfg = 2;
  Double_t s1Min  = 50.;  Double_t s1Max  = 5000.;
  Double_t s2Min  = 50.;  Double_t s2Max  = 10000.;
  Double_t f90Min = 0.;   Double_t f90Max = 1.;
  Double_t tofMin = 0.;   Double_t tofMax = 40.;

  Double_t tofBinSize = 0.25;
  Double_t tofBinNumber = (tofMax - tofMin)/tofBinSize;

  TCut tpcCut = DefineCuts(cfg, f90Min, f90Max, s1Min, s1Max, s2Min, s2Max, 0, 0);
  TCutG* lowBeCut = LowBeGraphCut(run, "LowBeCut");
  TCut combinedCut = tpcCut && "LowBeCut";

  TH1F* tofHistogram       = GenerateToFHistogram(runFile_name, tpcCut, tofBinNumber, tofMin, tofMax, normalize);
  TH1F* tofBeHistogram     = GenerateToFHistogram(runFile_name, combinedCut, tofBinNumber, tofMin, tofMax, normalize);
  TH1F* tofLogHistogram    = new TH1F("tofLogHist", "", tofBinNumber, tofMin, tofMax);
  TH1F* tofLogBeHistogram  = new TH1F("tofLogBeHist", "", tofBinNumber, tofMin, tofMax);

  tofHistogram   -> SetLineColor(kBlue);    tofHistogram   -> SetLineWidth(2);
  tofBeHistogram -> SetLineColor(kRed);     tofBeHistogram -> SetLineWidth(2);

  tofLogHistogram   = (TH1F*) tofHistogram -> Clone();
  tofLogBeHistogram = (TH1F*) tofBeHistogram -> Clone();

  tofHistogram   -> SetName("tofHist");
  tofBeHistogram -> SetName("tofBeHist");
  tofHistogram -> GetYaxis() -> SetTitle(Form("Counts/%3.2f ns", tofBinSize));
  tofHistogram -> GetXaxis() -> SetTitle("ns");

  tofHistogram   -> Draw();
  tofBeHistogram -> Draw("SAME HIST");

  TLegend *legend = new TLegend(0.85,0.74,0.95,0.95);
  legend -> SetTextFont(102);
  legend -> SetTextSize(0.03);
  legend -> AddEntry(tofHistogram, "No Be Selection", "l");
  legend -> AddEntry(tofBeHistogram, "Be Selection", "l");
  legend -> Draw();

  Double_t pointX = 0.15;    Double_t pointY = 0.3;    Double_t padRatio = 0.4;

  TPad* logPad = new TPad("logPad", "", pointX, pointY, pointX + padRatio, pointY + padRatio);
  logPad -> SetLogy();   logPad -> Draw();   logPad -> cd();

  tofLogHistogram -> SetTitle("");
  tofLogHistogram -> GetXaxis() -> SetTitle("");
  tofLogHistogram -> Draw();

  tofLogBeHistogram -> SetTitle("");
  tofLogBeHistogram -> GetXaxis() -> SetTitle("");
  tofLogBeHistogram -> Draw("SAME HIST");
}

void F90vToF(Int_t run){

  TString runFile_name = Form("runs/run_%d.root", run);
  TFile* file = CheckFile(runFile_name);
  TTree* reco;  file -> GetObject("reco", reco);

  Int_t cfg = 2;
  Double_t f90Min = 0.;   Double_t f90Max = 1.;
  Double_t tofMin = -100.;   Double_t tofMax = 100.;

  Double_t tofBinSize = 0.25;
  Double_t tofBinNumber = (tofMax - tofMin)/tofBinSize;
  Double_t f90BinSize = 0.005;
  Double_t f90BinNumber = (f90Max - f90Min)/f90BinSize;

  TCut tpcCut = DefineCuts(cfg, f90Min, f90Max, 0, 0, 0, 0, tofMin, tofMax);
  TCutG* lowBeCut = LowBeGraphCut(run);
  TCut combinedCut = tpcCut && "lowBeCut";

  TH2F* f90ToF_hist = new TH2F("f90ToF_hist", "F90 vs Time of Flight (TPC and SiTEL); ToF (ns); f90", tofBinNumber, tofMin, tofMax, f90BinNumber, f90Min, f90Max);
  reco -> Project( "f90ToF_hist", "clusters[0].f90:2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time)", combinedCut );

  f90ToF_hist -> Draw("COLZ");
}
