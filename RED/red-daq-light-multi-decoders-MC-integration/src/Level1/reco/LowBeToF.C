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
  banana_hist -> Draw("CONT Z LIST");
  E_low_L -> Draw("SAME");  E_upp_L -> Draw("SAME");  dE_low_L -> Draw("SAME");  dE_upp_L -> Draw("SAME");

  TCut E_low_cut  = Form("baseline_mean[31] - ymin[31] > %f", E_low);
  TCut E_upp_cut  = Form("baseline_mean[31] - ymin[31] < %f", E_upp);
  TCut dE_low_cut = Form("baseline_mean[30] - ymin[30] > %f", dE_low);
  TCut dE_upp_cut = Form("baseline_mean[30] - ymin[30] < %f", dE_upp);

  TCut lowBe_cut = E_low_cut && E_upp_cut && dE_low_cut && dE_upp_cut;
  return lowBe_cut;
}

void LowBeContourCut( Int_t run, Int_t num_bins, Int_t E_min, Int_t E_max, Int_t dE_min, Int_t dE_max ){

  TH2F* banana_hist = Bananator(run, num_bins, E_min, E_max, dE_min, dE_max);

  Double_t gold_ratio = 1.61803398875;
  Double_t height = 500;    Double_t width = gold_ratio * height;

  TCanvas* canvas1 = new TCanvas("canvas1", "Delta E/E Contours", width, height);

  Int_t number_of_contours = 1;   Double_t size_of_contours = 10/number_of_contours;
  Double_t contours[number_of_contours];

  for (Int_t i = 0; i < number_of_contours; i++){
    contours[i] = (i+1)*size_of_contours;
  }

  banana_hist -> SetContour(number_of_contours, contours);

  banana_hist -> Draw("CONT Z LIST");
  canvas1 -> Update();

  // Get Contours
  TObjArray *conts     = (TObjArray*)gROOT -> GetListOfSpecials() -> FindObject("contours");
  TList* contour_level = NULL;
  TGraph* curve        = NULL;
  TGraph* curve_clone  = NULL;
  TGraph* max_curve    = NULL;

  Int_t number_of_graphs = 0;
  Int_t total_contours   = 0;

  if ( conts == NULL ){
    std::cout << "No contours were extracted!" << std::endl;
    total_contours = 0;
    return 0;
  } else {
    total_contours = conts -> GetSize();
  }

  std::cout << "Total Contours:" << total_contours << std::endl;

  for (Int_t i = 0; i < total_contours; i++){
    contour_level = (TList*) conts -> At(i);
    std::cout << "Contour " << i << " has " << contour_level -> GetSize() << " graphs." << std::endl;
    number_of_graphs += contour_level -> GetSize();
  }

  number_of_graphs = 0;

  TCanvas* canvas2 = new TCanvas("canvas2", "Delta E/E Contours", width, height);
  canvas2 -> Divide(2,1);
  TH2F* hr = new TH2F("hr", "Contours", 2, E_min, E_max, 2, dE_min, dE_max);

  canvas2 -> cd(1);
  hr->Draw();
  Double_t x0, y0, z0;
  Int_t max_entries = 0;

  for(Int_t i = 0; i < total_contours; i++){
    contour_level = (TList*)conts->At(i);
    z0 = contours[i];
    std::cout << "Z-Level Passed in as:  Z = " << z0 << std::endl;

    // Get first graph from list on curves on this level
    curve = (TGraph*)contour_level->First();
    for(Int_t j = 0; j < contour_level -> GetSize(); j++){
      curve -> GetPoint(0, x0, y0);
      curve -> SetLineColor(2 + i );
      number_of_graphs ++;
      std::cout << "Graph: " << number_of_graphs << " -- " << curve -> GetN() << " Elements" << std::endl;

      if ( curve -> GetN() > max_entries ) {
        max_curve = curve;
        max_entries = max_curve -> GetN();
      }

      // Draw clones of the graphs to avoid deletions in case the 1st pad is redrawn.
      curve_clone = (TGraph*)curve->Clone();
      curve_clone->Draw("C");

      curve = (TGraph*)contour_level->After(curve); // Get Next graph
    }
  }
  canvas2 -> cd(2);
  max_curve -> Draw();

  canvas2->Update();
  std::cout << "Amount of entries in curve with largest number of entries:" << max_curve -> GetN();
  printf("\n\n\tExtracted %d Contours and %d Graphs \n", total_contours, number_of_graphs );
  gStyle->SetTitleW(0.);
  gStyle->SetTitleH(0.);

  TFile* output_file = new TFile("LowBeCut.root", "UPDATE");
  output_file -> WriteObject(max_curve, Form("lowBe_cut_%d", run), "OverWrite");
}

void TimeOfFlight(Int_t run){

  TString runFile_name = Form("runs/run_%d.root", run);  TFile* runFile = CheckFile(runFile_name);
  TTree* reco;    runFile -> GetObject("reco", reco);
  bool normalize = false;

  Int_t cfg = 2;
  Double_t s1_min  = 50.;  Double_t s1_max  = 5000.;
  Double_t s2_min  = 50.;  Double_t s2_max  = 10000.;
  Double_t f90_min = 0.;   Double_t f90_max = 1.;
  Double_t tof_min = 0.;   Double_t tof_max = 80.;

  TCut histogram_cuts = DefineCuts(cfg, f90_min, f90_max, s1_min, s1_max, s2_min, s2_max, tof_min, tof_min);
  TCut lowBe_cut = LowBeCut( run, 400, 1500, 4000, 400, 1200 );
  TCut total_cut = histogram_cuts && lowBe_cut;

  reco -> Draw("2*(start_time[30] - clusters[0].cdf_time) >> hist(160, 0, 80)", total_cut);
  /*TH1F* tof_histogram = GenerateToFHistogram(runFile_name, total_cut, 160, tof_min, tof_max, normalize);

  tof_histogram -> GetYaxis() -> SetTitle("Counts/2 ns");
  tof_histogram -> GetXaxis() -> SetTitle("ns");
  tof_histogram -> Draw();*/
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

  TH2F* f90ToF_hist = new TH2F("f90ToF_hist", "F90 vs Time of Flight (TPC and SiTEL); ToF (ns); f90", 100, tof_min, tof_max, 120, f90_min, f90_max);
  reco -> Project( "f90ToF_hist", "clusters[0].f90:2*(xmin[30] - clusters[0].min_x)", total_cut );

  f90ToF_hist -> Draw("COLZ");
}
