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

Int_t CountOf2DHistograms(TDirectory* directory);
TCut SiTelTPCToFCut(Double_t tofMin, Double_t tofMax);
TCut SiTelLSciToFCut(Int_t chanID, Double_t tofMin, Double_t tofMax);
TCut SiTelLSciToFCutAllChan(Int_t *chanID, Int_t chanCount, Double_t tofMin, Double_t tofMax);


Int_t expCfg = 2;

Double_t f90Min = 0.0;    Double_t f90Max = 1.0;
Double_t s1Min  = 60.0;   Double_t s1Max  = 1000.0;

Double_t f90BinSize = 5E-3;   Double_t f90BinCount = (f90Max - f90Min)/f90BinSize;
Double_t s1BinSize  = 2.5;     Double_t s1BinCount  = (s1Max - s1Min)/s1BinSize;


void F90vS1( Int_t run ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> SetOptStat(0);
  sidStyle -> cd();

  TString runFileName = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(runFileName);
  TTree* reco;  file -> GetObject("reco", reco);

  // Defining the selection cuts to be applied:
  Double_t tpcToFMin  = 22;   Double_t tpcToFMax  = 27;
  Double_t lsciToFMin = 20;   Double_t lsciToFMax = 40;

  const Int_t chanTotal = 6;
  Int_t chanID[chanTotal]     = { 0   , 1  ,  3  ,  4  ,  5  ,  8};
  Int_t chanColors[chanTotal] = {1179, 1230, 1281, 1332, 1383, 1433};

  TCut lsciCuts[chanTotal];
  for (Int_t i = 0; i < chanTotal; i++){
    lsciCuts[i] = SiTelLSciToFCut(chanID[i], lsciToFMin, lsciToFMax);
  }

  TCut   tpcToF     = SiTelTPCToFCut(tpcToFMin, tpcToFMax);
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

    f90vS1Hist[i] = new TH2F(histNames[i], histTitles[i], s1BinCount, s1Min, s1Max, f90BinCount, f90Min, f90Max);
    reco -> Project(histNames[i], "clusters[0].f90:clusters[0].charge", generalCuts[i]);

    canvas -> cd(i + 1);
    f90vS1Hist[i] -> GetXaxis() -> SetTitle("S1 [PE]");
    f90vS1Hist[i] -> GetYaxis() -> SetTitle("f90");
    f90vS1Hist[i] -> Draw("COLZ");

    if (i == 1){
      for (Int_t j = 0; j < chanTotal; j++){
        lsciF90vS1Hist[j] = new TH2F(lsciHistNames[j], histTitles[i+1], s1BinCount, s1Min, s1Max, f90BinCount, f90Min, f90Max);
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

  cout << "Count of events after LSci timing cut: " << lsciEvents << endl;

}

/* void F90vS1ProgressiveCuts( Int_t run, bool write = false, bool draw = true )
 *
 *  Summary of Function:
 *
 *    The function construncs 4 2D histograms with x-axis equal to S1 charge and y-axis equal to f90. Each
 *    histogram, while taken from the same raw dataset define by the run Count, has different selection cuts
 *    applied to it. Each new cut also applies the previous ones to it. The are: no cuts, TPC signal quality
 *    cuts (Count of clusters, complete signal reconstruction), SiTel-TPC time of flight cuts and SiTel-LSci
 *    time of flight cuts. As such, each histogram has progressively stricter cut selections, which should
 *    result in a greater percentage of nuclear recoil events with respect to background, as is desired.
 *
 *    The code, by default, then draws the generated histograms side by side. The user may opt to not draw the
 *    histograms if so desired. The user may also opt to save the generated histograms in a root file named
 *    'analysis_<run>.root', were '<run>' is the run Count.
 *
 *  Parameters   : run   >> the run containing the reconstructed data.
 *                 write >> tells the code to save the generated histograms if set to 'true'.
 *                 draw  >> tells the code to draw the generated histograms if set to 'true'.
 *
 *  Return Value :
 */
void F90vS1ProgressiveCuts( Int_t run, bool write = false, bool draw = true ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> cd();
  sidStyle -> SetOptStat(0);


  TString runFileName = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* runFile = CheckFile(runFileName);
  TTree* reco;  runFile -> GetObject("reco", reco);


  // Defining the selection cuts to be applied:
  const Int_t chanCount = 6;
  Int_t chanID[chanCount] = {0, 1, 3, 4, 5, 8};

  Double_t tpcToFMin  = 22;   Double_t tpcToFMax  = 27;
  Double_t lsciToFMin = 20;   Double_t lsciToFMax = 40;

  TCut tpcCut     = DefineCuts(expCfg, f90Min, f90Max, s1Min, s1Max);
  TCut tpcToFCut  = SiTelTPCToFCut(tpcToFMin, tpcToFMax);
  TCut lsciToFCut = SiTelLSciToFCutAllChan(chanID, chanCount, lsciToFMin, lsciToFMax);


  const Int_t histCount = 4;
  TH2F*   f90vS1Hist[histCount];
  TString histName[histCount]  = {"f90vs1_histogram_nocut"    , "f90vs1_histogram_tpccut",
                                   "f90vs1_histogram_tpctofcut", "f90vs1_histogram_lscitofcut"};
  TString histTitle[histCount] = {"F90 vs S1 (No Cuts)",
                                   "F90 vs S1 (TPC Cuts)",
                                   "F90 vs S1 (TPC + SiTel-TPC ToF Cuts)",
                                   "F90 vs S1 (TPC + SiTel-TPC ToF + SiTel-LSci ToF Cuts)"};

  TCut cut[histCount] = {"", tpcCut, tpcCut && tpcToFCut, tpcCut && tpcToFCut && lsciToFCut};

  for (Int_t i = 0; i < histCount; i++){
    f90vS1Hist[i] = new TH2F(histName[i], histTitle[i], s1BinCount, s1Min, s1Max, f90BinCount, f90Min, f90Max);
    reco -> Project(histName[i], "clusters[0].f90:clusters[0].charge", cut[i]);
  }

  // Saves the generated histograms to root file for later use.
  if (write){

    TString outputFileName = analysisDirectoryPath + Form("/analysis_%d.root", run);
    TFile* outputFile = CheckFile(outputFileName);

    TDirectory* histDir = MakeDirectory("histograms", "histograms");
    histDir -> cd();
    TDirectory* f90vs1Dir = MakeDirectory("f90vS1", "f90vS1");

    for (Int_t i = 0; i < histCount; i++){
      f90vs1Dir -> WriteObject(f90vS1Hist[i], histName[i]);
    }
    outputFile -> Close();
  }

  // Draws the generated histograms for the user.
  if (draw){

    Double_t height = 300;    Double_t width = gdRatio*height;
    TCanvas* canvas = new TCanvas("canvas", "F90 vs S1", 2*width, 2*height);
    canvas -> Divide(2,2);

    for (Int_t i = 0; i < histCount; i++){

      canvas -> cd(i + 1);

      if (i == 3 ){
        f90vS1Hist[i-1] -> SetMarkerColorAlpha(kGray, 0.01);
        f90vS1Hist[i-1] -> Draw("");
        f90vS1Hist[i] -> SetMarkerStyle(20);
        f90vS1Hist[i] -> SetMarkerSize(1);
        f90vS1Hist[i] -> Draw("SAME");
      } else {
        f90vS1Hist[i] -> Draw("COLZ");
      }
    }

  }

}

void CutEffectAnalysis( Int_t run, Double_t s1BinSize = 20 ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> cd();
  sidStyle -> SetOptStat(0);


  TString analysisFileName = analysisDirectoryPath + Form("/analysis_%d.root", run);
  TFile* analysisFile = CheckFile(analysisFileName);

  TDirectory* histDir = MakeDirectory("histograms", "histograms");
  histDir -> cd();
  TDirectory* f90vs1Dir = MakeDirectory("f90vS1", "f90vS1");


  const Int_t histCount = 4;
  TString histName[histCount] = {"f90vs1_histogram_nocut"    , "f90vs1_histogram_tpccut",
                                  "f90vs1_histogram_tpctofcut", "f90vs1_histogram_lscitofcut"};

  TH2F* f90vs1Hist[histCount];
  for (Int_t i = 0; i < histCount; i++){
    f90vs1Hist[i] =  (TH2F*) f90vs1Dir -> Get(histName[i]);
  }


  Double_t s1Min = f90vs1Hist[0] -> GetXaxis() -> GetXmin();
  Double_t s1Max = f90vs1Hist[0] -> GetXaxis() -> GetXmax();

  const Int_t s1BinCount = (Int_t) (s1Max - s1Min)/s1BinSize;

  // Setting up the values that will be used to impose the region of interest when counting the number of events.
  Double_t nrEventCount[histCount][s1BinCount];
  Double_t erEventCount[histCount][s1BinCount];
  Double_t s1[s1BinCount];

  Double_t nrF90Min   = 0.4;    Double_t nrF90Max    = 0.6;
  Double_t erF90Min   = 0.2;    Double_t erF90Max    = 0.4;
  Double_t s1LowBound = 0.;     Double_t s1UppBound0 = 0.;

  Int_t nrF90MinBin;      Int_t nrF90MaxBin;
  Int_t erF90MinBin;      Int_t erF90MaxBin;
  Int_t s1LowBoundBin;    Int_t s1UppBoundBin;

  for (Int_t i = 0; i < histCount; i++){

    nrF90MinBin = f90vs1Hist[i] -> GetYaxis() -> FindBin(nrF90Min);
    nrF90MaxBin = f90vs1Hist[i] -> GetYaxis() -> FindBin(nrF90Max);

    erF90MinBin = f90vs1Hist[i] -> GetYaxis() -> FindBin(erF90Min);
    erF90MaxBin = f90vs1Hist[i] -> GetYaxis() -> FindBin(erF90Max);

    cout << "NR Bin Range: " << nrF90MinBin << " to " << nrF90MaxBin << endl;
    cout << "ER Bin Range: " << erF90MinBin << " to " << erF90MaxBin << endl;

    Int_t j = 0;
    for (j = 0; j < s1BinCount; j++){

      cout << j << endl;

      //nrEventCount[i][j] = Integral();
    }
  }

}
/*
void CutEffectAnalysis( Int_t run ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> cd();
  sidStyle -> SetOptStat(0);

  TString runFileName = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(runFileName);
  TTree* reco;  file -> GetObject("reco", reco);


  // Defining the selection cuts to be applied:
  const Int_t chanCount = 6;
  Int_t chanID[chanCount] = {0, 1, 3, 4, 5, 8};

  Double_t tpcToFMin  = 22;   Double_t tpcToFMax  = 27;
  Double_t lsciToFMin = 20;   Double_t lsciToFMax = 40;

  TCut tpcCut    = DefineCuts(expCfg, f90Min, f90Max, s1Min, s1Max);
  TCut tpcToFCut = SiTelTPCToFCut(tpcToFMin, tpcToFMax);
  TCut lsciToFCut = SiTelLSciToFCutAllChan(chanID);
  //

  TString histName[2] = {"cut1_hist", "cut2_hist"};
  TString histTitle[2] = {"F90 vs S1 (TPC Cuts)", "F90 vs S1 (TPC + SiTel-TPC ToF Cuts)"};

  TH2F* firstCutHist  = new TH2F("cut1_hist", "F90 vs S1 (TPC Cuts)", s1BinCount, s1Min, s1Max, f90BinCount, f90Min, f90Max);
  TH2F* secondCutHist = new TH2F("cut2_hist", "F90 vs S1 (TPC Cuts)", s1BinCount, s1Min, s1Max, f90BinCount, f90Min, f90Max);

  Double_t nrF90Min = 0.4;   Double_t nrF90Max = 0.6;
  Double_t erF90Min = 0.2;   Double_t erF90Max = 0.4;



  // Cuts to select the desired histogram slice.
  TCut nrHistSlice;   TCut erHistSlice;

  Double_t s1Size = 94;   const Int_t pointsCount = (s1Max - s1Min)/s1Size;
  Double_t s1LowBound;    Double_t s1UppBound;
  Double_t s1[pointsCount];


  Double_t nrCount[pointsCount];
  Double_t erCount[pointsCount];
  Double_t total[pointsCount];

  for (Int_t i = 0; i < pointsCount; i++){

    s1LowBound = s1Min + i*s1Size;
    s1UppBound = s1Min + (i + 1)*s1Size;
    s1[i] = (s1LowBound + s1UppBound) / 2;

    nrHistSlice = DefineF90Range(nrF90Min, nrF90Max) && DefineS1Range(s1LowBound, s1UppBound);
    erHistSlice = DefineF90Range(erF90Min, erF90Max) && DefineS1Range(s1LowBound, s1UppBound);

    reco -> Project("nr_hist", "clusters[0].f90", nrHistSlice && qualityCut && tpcToFCut);
    reco -> Project("er_hist", "clusters[0].f90", erHistSlice && qualityCut && tpcToFCut);

    total[i] = nrHist -> GetEntries() + erHist -> GetEntries();
    nrCount[i] = (nrHist -> GetEntries())/total[i];
    erCount[i] = (erHist -> GetEntries())/total[i];
  }


  TGraph* nrGraph = new TGraph(pointsCount, s1, nrCount);
  TGraph* erGraph = new TGraph(pointsCount, s1, erCount);

  erGraph -> SetLineColor(kRed);

  TCanvas* c1 = new TCanvas("c1", "c1", 500*gdRatio, 500);
  erGraph -> Draw();

  TCanvas* c2 = new TCanvas("c2", "c2", 500*gdRatio, 500);
  nrGraph -> Draw();
}
*/

// =================================== AUXILIARY FUNCTIONS USED IN THE MACRO =================================== //

/* Int_t CountOf2DHistograms ( TDirectory* directory )
 *
 * Summary of function:
 *
 *    The CountofHistograms function counts the Count of two dimenstional histograms in a given directory
 *    by interating over all objects located in it and checking which inherit from TH2. It then returns the
 *    Count of 2D histograms in the directory.
 *
 * Parameters   : directory >> a pointer to the the directory where the histograms are located.
 *
 * Return value : Int_t histogramCount.
 */
Int_t CountOf2DHistograms( TDirectory* directory ){

  Int_t histogramCount = 0;
  TKey* key;

  TIter next( (TList *)directory -> GetListOfKeys() );

  while ( (key = (TKey *)next()) ) {
    TClass *object_class = gROOT -> GetClass( key->GetClassName() );
    if ( object_class -> InheritsFrom("TH2") ) {
      histogramCount++;
    }
  }

  return histogramCount;
}

TCut SiTelTPCToFCut(Double_t tofMin, Double_t tofMax){

  TCut tofMinCut = Form("2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time) >= %f", tofMin);
  TCut tofMaxCut = Form("2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time) <= %f", tofMax);
  TCut tofCut = tofMinCut && tofMaxCut;

  return tofCut;
}

TCut SiTelLSciToFCut(Int_t chanID, Double_t tofMin, Double_t tofMax){

  TCut tofMinCut = Form("2*(start_time[%d] - 0.5*(start_time[30] + start_time[31] - 7.45)) >= %f", chanID, tofMin);
  TCut tofMaxCut = Form("2*(start_time[%d] - 0.5*(start_time[30] + start_time[31] - 7.45)) <= %f", chanID, tofMax);
  TCut tofCut = tofMinCut && tofMaxCut;

  return tofCut;
}

TCut SiTelLSciToFCutAllChan(Int_t *chanID, Int_t chanCount,Double_t tofMin, Double_t tofMax){

  TCut tofLSciCut[chanCount];
  for (Int_t i = 0; i < chanCount; i++) tofLSciCut[i] = SiTelLSciToFCut(chanID[i], tofMin, tofMax);

  TCut tofCut = tofLSciCut[0] || tofLSciCut[1] || tofLSciCut[2] || tofLSciCut[3] || tofLSciCut[4] || tofLSciCut[5];

  return tofCut;
}
