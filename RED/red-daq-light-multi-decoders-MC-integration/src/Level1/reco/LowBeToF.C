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


TCut LowBeCut( Int_t run, Int_t num_bins, Int_t E_min, Int_t E_max, Int_t dE_min, Int_t dE_max ){

  TH2F* bananaHist = Bananator(run, num_bins, E_min, E_max, dE_min, dE_max);
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
  Int_t E_min  = 1500;  Int_t E_max  = 4000;
  Int_t dE_min = 400;   Int_t dE_max = 1200;

  TH2F* bananaHist = Bananator(run, num_bins, E_min, E_max, dE_min, dE_max);

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
  TH2F* hr = new TH2F("hr", "Contours", 2, E_min, E_max, 2, dE_min, dE_max);

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


TCutG* LowBeGraphCut( int run ){

  TFile*  bCutFile    = new TFile("LowBeCut.root", "UPDATE");
  TCutG*  bGraphCut   = new TCutG("lowBe_cut");
  TString bCutName    = Form("lowBecut_%d", run);
  TGraph* sourceGraph = new TGraph();

  if ( bCutFile -> IsOpen() ){
    sourceGraph = (TGraph*)bCutFile -> Get(bCutName);
    bCutFile -> Close();
  }

  double* x = sourceGraph -> GetX();
  double* y = sourceGraph -> GetY();

  for (Int_t i = 0; i < sourceGraph -> GetN(); i++){
    bGraphCut -> SetPoint(i, x[i], y[i]);
  }

  bGraphCut -> SetVarX("baseline_mean[31] - ymin[31]");
  bGraphCut -> SetVarY("baseline_mean[30] - ymin[30]");

  /*
  TFile* data_file = new TFile(Form("runs/run_%d.root", run), "UPDATE");
  TTree* reco;    data_file -> GetObject("reco", reco);


  reco -> Draw("baseline_mean[30] - ymin[30]:baseline_mean[31] - ymin[31] >> hist", "gcut", "colz");
  */

  return bGraphCut;
}

void TimeOfFlight(Int_t run){

  TString runFile_name = Form("runs/run_%d.root", run);  TFile* runFile = CheckFile(runFile_name);
  TTree* reco;    runFile -> GetObject("reco", reco);

  bool normalize = false;
  Int_t cfg = 2;
  Double_t s1_min  = 50.;  Double_t s1_max  = 5000.;
  Double_t s2_min  = 50.;  Double_t s2_max  = 10000.;
  Double_t f90_min = 0.;   Double_t f90_max = 1.;
  Double_t tof_min = 0.;   Double_t tof_max = 40.;

  Double_t binSize = 0.25;
  Double_t binNumber = (tof_max - tof_min)/binSize;

  TCut histogram_cuts = DefineCuts(cfg, f90_min, f90_max, s1_min, s1_max, s2_min, s2_max, 0, 0);
  TCutG* lowBe_cut = LowBeGraphCut(run);
  TCut combinedCut = histogram_cuts && "lowBe_cut";

  TH1F* tofHistogram   = GenerateToFHistogram(runFile_name, histogram_cuts, binNumber, tof_min, tof_max, normalize);
  TH1F* tofBeHistogram = GenerateToFHistogram(runFile_name, combinedCut, binNumber, tof_min, tof_max, normalize);

  gStyle -> SetLabelFont(102, "xyz");
  gStyle -> SetTitleFont(102, "xyz");
  gStyle -> SetTitleFont(102, "t");

  tofHistogram -> SetName("tofHist");
  tofHistogram -> GetYaxis() -> SetTitle(Form("Counts/%3.2f ns", binSize));
  tofHistogram -> GetXaxis() -> SetTitle("ns");
  tofHistogram -> SetLineColor(kBlue);
  tofHistogram -> Draw();

  tofBeHistogram -> SetName("tofBeHist");
  tofBeHistogram -> SetLineColor(kRed);
  tofBeHistogram -> Draw("SAME HIST");

  TLegend *legend = new TLegend(0.85,0.74,0.95,0.95);
  legend -> SetTextFont(102);
  legend -> SetTextSize(0.03);
  legend -> AddEntry(tofHistogram, "No Be Selection", "l");
  legend -> AddEntry(tofBeHistogram, "Be Selection", "l");
  legend -> Draw();
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

  TH2F* f90ToF_hist = new TH2F("f90ToF_hist", "F90 vs Time of Flight (TPC and SiTEL); ToF (ns); f90", 400, tof_min, tof_max, 120, f90_min, f90_max);
  reco -> Project( "f90ToF_hist", "clusters[0].f90:2*(start_time[30] - clusters[0].cdf_time)", total_cut );

  f90ToF_hist -> Draw("COLZ");
}
