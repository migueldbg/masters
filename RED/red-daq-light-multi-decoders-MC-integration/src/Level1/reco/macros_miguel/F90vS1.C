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


// <code-fold> AUXILIARY FUNCTIONS DECLARATION.
Int_t CountOf2DHistograms(TDirectory* directory);
TGraph* DivideGraphs(TGraph* numGraph, TGraph* denGraph);
TCut SiTelTPCToFCut(Double_t tofMin, Double_t tofMax);
TCut SiTelLSciToFCut(Int_t chanID, Double_t tofMin, Double_t tofMax);
TCut SiTelLSciToFCutAllChan(Int_t *chanID, Int_t chanCount, Double_t tofMin, Double_t tofMax);
// </code-fold>


// <code-fold> DECLARATION OF VARIABLES TO BE USED IN ALL FUNCTIONS.
Int_t expCfg = 2;

Double_t f90Min = 0.0;    Double_t f90Max = 1.0;
Double_t s1Min  = 50.0;   Double_t s1Max  = 1000.0;

Double_t f90BinSize = 5E-3;   Double_t f90BinCount = (f90Max - f90Min)/f90BinSize;
Double_t s1BinSize  = 2.5;     Double_t s1BinCount  = (s1Max - s1Min)/s1BinSize;
// </code-fold>

//TODO: write the intendend workflow when using this macro.

void tempF90vS1( Int_t run ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> SetOptStat(111);
  sidStyle -> cd();

  TString runFileName = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(runFileName);
  TTree* reco;  file -> GetObject("reco", reco);

  TCut tpcCut = DefineCuts(1, f90Min, f90Max, s1Min, s1Max);
  TCutG* lowBeCut = LowBeGraphCut(run, "LowBeCut");

  // Histograms that will contain the results after the application of the selection cuts.
  TH2F* f90vS1Hist = new TH2F("f_{prompt} vs S1", "; S1 [PE]; #it{f}_{prompt}", s1BinCount, s1Min, s1Max, f90BinCount, f90Min, f90Max);

  Double_t height = 500;    Double_t width = gdRatio * height;
  TCanvas* canvas = new TCanvas("canvas", "canvas", width, height);



  reco -> Project("f_{prompt} vs S1", "clusters[0].f90:clusters[0].charge", "LowBeCut");

  f90vS1Hist -> Draw("COLZ");
}

void F90vS1LSciCuts( Int_t run ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> SetOptStat(0);
  sidStyle -> cd();

  TString runFileName = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(runFileName);
  TTree* reco;  file -> GetObject("reco", reco);

  //Defining the selection cuts to be applied:
  Double_t tpcToFMin  = 22;   Double_t tpcToFMax  = 27;
  Double_t lsciToFMin = 20;   Double_t lsciToFMax = 40;

  const Int_t chanTotal = 6;
  Int_t chanID[chanTotal]     = { 0   , 1  ,  3  ,  4  ,  5  ,  8};
  Int_t chanColors[chanTotal] = {1179, 1230, 1281, 1332, 1383, 1403};

  TCut lsciCuts[chanTotal];
  for (Int_t i = 0; i < chanTotal; i++){
    lsciCuts[i] = SiTelLSciToFCut(chanID[i], lsciToFMin, lsciToFMax);
  }

  TCut tpcToF = SiTelTPCToFCut(tpcToFMin, tpcToFMax);
  TCut tpcCut = DefineCuts(expCfg, f90Min, f90Max, s1Min, s1Max);
  TCutG* lowBeCut = LowBeGraphCut(run, "LowBeCut");
  TCut qualityCut = tpcCut && "LowBeCut" && tpcToF;

  // Histograms that will contain the results after the application of the selection cuts.
  TH2F* lsciF90vS1Hist[chanTotal];
  TString lsciHistNames[chanTotal] = {"LSci 0", "LSci 1", "LSci 3 & 7", "LSci 4 & 6", "LSci 5", "LSci 8"};

  Int_t lsciEvent[chanTotal];
  Int_t lsciEventTotal = 0;

  Double_t height = 500;    Double_t width = gdRatio * height;
  TCanvas* canvas = new TCanvas("canvas", "canvas", width, height);
  canvas -> SetGrid();

  TLegend* legend = new TLegend(0.2, 0.2);


  for (Int_t j = 0; j < chanTotal; j++){
    lsciF90vS1Hist[j] = new TH2F(lsciHistNames[j], "", s1BinCount, s1Min, s1Max, f90BinCount, f90Min, f90Max);
    reco -> Project(lsciHistNames[j], "clusters[0].f90:clusters[0].charge", lsciCuts[j] && qualityCut);

    lsciF90vS1Hist[j] -> SetMarkerStyle(20);
    lsciF90vS1Hist[j] -> SetMarkerSize(2);
    lsciF90vS1Hist[j] -> SetMarkerColor(chanColors[j]);

    lsciF90vS1Hist[j] -> GetXaxis() -> SetTitle("S1 [PE]");
    lsciF90vS1Hist[j] -> GetYaxis() -> SetTitle("f_{prompt}");
    lsciF90vS1Hist[j] -> Draw("SAME");

    lsciEvent[j] = (Int_t) lsciF90vS1Hist[j] -> GetEntries();
    lsciEventTotal = lsciEventTotal + (Int_t) lsciF90vS1Hist[j] -> GetEntries();

    legend -> AddEntry(lsciF90vS1Hist[j], lsciHistNames[j] + Form("(%d events)", lsciEvent[j]), "lep");
  }
  legend -> Draw();

  cout << "Count of events after LSci timing cut: " << lsciEventTotal << endl;

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
void ProgressiveCutsF90vS1( Int_t run, bool write = false, bool draw = true ){

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

  TCutG* lowBeCut = LowBeGraphCut(run, "LowBeCut");
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

  TCut cut[histCount] = {"", tpcCut && "LowBeCut", tpcCut && tpcToFCut && "LowBeCut", tpcCut && tpcToFCut && lsciToFCut && "LowBeCut"};

  for (Int_t i = 0; i < histCount; i++){
    f90vS1Hist[i] = new TH2F(histName[i], "", s1BinCount, s1Min, s1Max, f90BinCount, f90Min, f90Max);
    reco -> Project(histName[i], "clusters[0].f90:clusters[0].charge", cut[i]);

    f90vS1Hist[i] -> GetXaxis() -> SetTitle("S1 [PE]");
    f90vS1Hist[i] -> GetYaxis() -> SetTitle("f_{prompt}");
  }

  // Saves the generated histograms to root file for later use.
  if (write){

    TString outputFileName = analysisDirectoryPath + Form("/analysis_%d.root", run);
    TFile* outputFile = CheckFile(outputFileName);

    TDirectory* histDir = MakeDirectory("histograms", "histograms");
    histDir -> cd();
    TDirectory* f90vs1Dir = MakeDirectory("f90vS1", "f90vS1");

    for (Int_t i = 0; i < histCount; i++){
      f90vs1Dir -> WriteObject(f90vS1Hist[i], histName[i], "OverWrite");
    }
    outputFile -> Close();
  }

  // Draws the generated histograms for the user.
  if (draw){

    Double_t height = 300;    Double_t width = gdRatio*height;
    TCanvas* canvas[histCount];

    for (Int_t i = 0; i < histCount; i++){

      canvas[i] = new TCanvas(Form("canvas_%d", i), "F90 vs S1", 2*width, 2*height);

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

/* void ProgressiveCutsGraphs( Int_t run, Double_t s1BinSize = 20 )
 *
 *  Summary of Functions:
 *
 *    This functions takes the 2D f90 vs S1 histograms generated by 'ProgressiveCutsF90vS1()' and generates two
 *    TGraph objects from each. It does by slicing up the histogram in slices of equal S1 size. For each slice
 *    the program counts how many events fall within the NR region (0.4 <= f90 <= 0.6)and how many fall within
 *    the ER region (0.2 <= f90 <=). It then divides the counts obtained for each region by the total number of
 *    events in the S1 slice (0.0 <= f90 <= 1.0), resulting in a percentage NR and ER events. The code then
 *    creats 'percentage graphs' withi x-axis equal to the s1 value and y-axis equal to the event percentage.
 *    The graphs are then saved in a foldere within the 'analysis' root file corresponding to the run set used.
 *
 *    Note: the f90 values used to determine the NR and ER regions are constant for all values of S1 and were
 *    chosen based on a visual analysis of the distributions of the two populations in the f90 vs S1 histogram.
 *    A more robust definition of such region, through the use of NR and ER exclusive sources, would be preferable
 *    and will be sought after.
 *
 *  Parameters   : run >> the run containing the reconstructed data.
 *                 s1BinSize >> the size of the s1 slices used.
 */
void ProgressiveCutsGraphs( Int_t run, Double_t s1BinSize = 20 ){

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
  for (Int_t i = 0; i < histCount; i++) f90vs1Hist[i] =  (TH2F*) f90vs1Dir -> Get(histName[i]);


  // <code-fold> SETTING UP ARRAYS AND BOUNDARY VALUES.

  Double_t s1Min = f90vs1Hist[0] -> GetXaxis() -> GetXmin();
  Double_t s1Max = f90vs1Hist[0] -> GetXaxis() -> GetXmax();

  const Int_t s1BinCount = (Int_t) (s1Max - s1Min)/s1BinSize;

  // These arrays will hould the number of events within each regioSetting up the values that will be used to impose the region of interest when counting the number of events.n defined
  Double_t nrEventCount[histCount][s1BinCount];   Double_t nrPercentage[histCount][s1BinCount];
  Double_t erEventCount[histCount][s1BinCount];   Double_t erPercentage[histCount][s1BinCount];
  Double_t s1[s1BinCount];
  Double_t totalEventCount[histCount][s1BinCount];


  Double_t nrF90Min   = 0.4;    Double_t nrF90Max    = 0.6;
  Double_t erF90Min   = 0.2;    Double_t erF90Max    = 0.4;
  Double_t s1LowBound = 0.;     Double_t s1UppBound  = 0.;

  // The f90 range bin numbers will be equal for all histograms, since they all have the same number of bins.
  Int_t nrF90MinBin = f90vs1Hist[0] -> GetYaxis() -> FindBin(nrF90Min);
  Int_t nrF90MaxBin = f90vs1Hist[0] -> GetYaxis() -> FindBin(nrF90Max);

  Int_t erF90MinBin = f90vs1Hist[0] -> GetYaxis() -> FindBin(erF90Min);
  Int_t erF90MaxBin = f90vs1Hist[0] -> GetYaxis() -> FindBin(erF90Max);

  Int_t f90MinBin = f90vs1Hist[0] -> GetYaxis() -> FindBin(f90Min);
  Int_t f90MaxBin = f90vs1Hist[0] -> GetYaxis() -> FindBin(f90Max);

  // The same applies for the s1 boundaries, though they will change as they run through the values of s1.
  Int_t s1LowBoundBin;    Int_t s1UppBoundBin;
  // </code-fold>

  // Getting the count of NR and ER events per S1 bin.
  for (Int_t i = 0; i < histCount; i++){

    Int_t j = 0;
    for (j = 0; j < s1BinCount; j++){

      s1LowBound = j * s1BinSize + s1Min;
      s1UppBound = (j + 1)* s1BinSize + s1Min;

      s1[j] = (s1LowBound + s1UppBound)/2;

      s1LowBoundBin = f90vs1Hist[i] -> GetXaxis() -> FindBin(s1LowBound);
      s1UppBoundBin = f90vs1Hist[i] -> GetXaxis() -> FindBin(s1UppBound);

      if (s1UppBound == 1000) s1UppBoundBin = 376;

      nrEventCount[i][j] = f90vs1Hist[i] -> Integral(s1LowBoundBin, s1UppBoundBin, nrF90MinBin, nrF90MaxBin);
      erEventCount[i][j] = f90vs1Hist[i] -> Integral(s1LowBoundBin, s1UppBoundBin, erF90MinBin, erF90MaxBin);

      totalEventCount[i][j] = f90vs1Hist[i] -> Integral(s1LowBoundBin, s1UppBoundBin, f90MinBin, f90MaxBin);
      nrPercentage[i][j] = nrEventCount[i][j]/totalEventCount[i][j];
      erPercentage[i][j] = erEventCount[i][j]/totalEventCount[i][j];
    }
  }


  // Directory were the TGraph objects will be saved.
  analysisFile -> cd();
  TDirectory* graphDir = MakeDirectory("graphs", "graphs");
  graphDir -> cd();
  TDirectory* cutEffDir = MakeDirectory("cut effects", "cut effects");
  cutEffDir -> cd();

  TGraph* percGraph[6];
  // <code-fold> SETTING UP GRAPH PARAMETERS.
  TString graphName[6]  = {"nr_graph_nocut", "nr_graph_tpccut", "nr_graph_tpctofcut",
                           "er_graph_nocut", "er_graph_tpccut", "er_graph_tpctofcut"};
  TString graphTitle[6] = {"Nenhum Corte (NR)", "TPC + Be (NR)", "TPC + ToF + Be (NR)",
                           "Nenhum Corte (ER)", "TPC + Be (ER)", "TPC + ToF + Be (ER)"};

  Int_t graphColor[3] = {1179, 1230, 1281};

  Int_t nrMarker = 23;
  Int_t erMarker = 22;
  // </code-fold>

  // Define and write graphs in 'cutEffDir'.
  for (Int_t i = 0; i < 6; i++){

    if (i < 3) {
      percGraph[i] = new TGraph(s1BinCount, s1, nrPercentage[i]);
      percGraph[i] -> SetMarkerStyle(nrMarker);
    } else {
      percGraph[i] = new TGraph(s1BinCount, s1, erPercentage[i%3]);
      percGraph[i] -> SetMarkerStyle(erMarker);
    }

    percGraph[i] -> SetName(graphName[i]);
    percGraph[i] -> SetTitle(graphTitle[i]);
    percGraph[i] -> SetLineColor(graphColor[i % 3]);
    percGraph[i] -> SetMarkerColor(graphColor[i % 3]);

    cutEffDir -> WriteObject(percGraph[i], graphName[i], "OverWrite");
  }
  analysisFile -> Close();
}

/* void ProgressiveCutsEffectAnalysis( Int_t run, bool singleEventAllCut = true, bool bothEventOneCut = false )
 *
 *  Summary of Function:
 *
 *    The function takes the graphs constructed by 'ProgresiveCutsGraphs()' and plots them in up to two ways,
 *    depending on the parameters passed. The SiTel-LSci ToF isn't considered due to the small number of data
 *    points that survive the cut.
 *
 *    The first constitutes in creating two separate graphs, one for NR and one for ER. Each plot constitutes
 *    of three event percentage graphs, each corresponding to different set of cuts. The ratio of the graphs
 *    with respect to the graph with no cuts is also plotted, such as to more explicitely show how the cuts
 *    affect the percentage of events.
 *
 *    The second constitutes of three plots, one for each cut applied. Each plot has the NR and ER event
 *    percentage graphs plotted together. The ration of the two graphs is also plotted.
 *
 *  Parameters   : run >> the run containing the reconstructed data.
 *                 singleEventAllCut >> tells the code to construct the first set of plots.
 *                 bothEventOneCut >> tells the code to construct the second set of plots.
 *
 *  Return Value : void.
 */
void ProgressiveCutsEffectAnalysis( Int_t run, bool singleEventAllCut = true, bool bothEventOneCut = false ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> cd();
  sidStyle -> SetOptStat(0);

  TString analysisFileName = analysisDirectoryPath + Form("/analysis_%d.root", run);
  TFile* analysisFile = CheckFile(analysisFileName);

  TDirectory* graphDir = MakeDirectory("graphs", "graphs");
  graphDir -> cd();
  TDirectory* cutEffDir = MakeDirectory("cut effects", "cut effects");

  // <code-fold> GETTING THE TGRAPHS OBJECTS GENERATED BY 'ProgressiveCutsGraphs()'.
  const Int_t graphCount = 3;
  TGraph* nrGraph[graphCount];
  TGraph* erGraph[graphCount];

  TString nrGraphName[graphCount] = {"nr_graph_nocut", "nr_graph_tpccut", "nr_graph_tpctofcut"};
  TString erGraphName[graphCount] = {"er_graph_nocut", "er_graph_tpccut", "er_graph_tpctofcut"};

  for (Int_t i = 0; i < graphCount; i++){
    nrGraph[i] = (TGraph*) cutEffDir -> Get(nrGraphName[i]);
    erGraph[i] = (TGraph*) cutEffDir -> Get(erGraphName[i]);

    nrGraph[i] -> SetLineWidth(3);    nrGraph[i] -> SetMarkerSize(1);
    erGraph[i] -> SetLineWidth(3);    erGraph[i] -> SetMarkerSize(1);   erGraph[i] -> SetLineStyle(2);

    //TODO: Settle on the esthetic parameters of the graphs.
  }

  analysisFile -> Close();
  // </code-fold>

  Double_t height = 500;    Double_t width = gdRatio*height;

  // Takes all NR percentage graphs and plot them together to show the effect of the cuts. The same is done for ER.
  if (singleEventAllCut) {

    TMultiGraph* nrAllGraph = new TMultiGraph("nr_all_cuts", "; ; Fração de Eventos");
    TMultiGraph* erAllGraph = new TMultiGraph("er_all_cuts", "; ; Fração de Eventos");

    TGraph* nrRatio[graphCount];
    TGraph* erRatio[graphCount];

    TMultiGraph* nrAllRatio = new TMultiGraph("nr_ratio", "; S1[PE]; Ratio");
    TMultiGraph* erAllRatio = new TMultiGraph("er_ratio", "; S1[PE]; Ratio");

    // Getting the ratio graphs and adding all to the TMultiGraphs.
    for (Int_t i = 0; i < graphCount; i++){
      nrAllGraph -> Add(nrGraph[i]);
      erAllGraph -> Add(erGraph[i]);

      if (i > 0){
        nrRatio[i] = DivideGraphs(nrGraph[i], nrGraph[0]);
        erRatio[i] = DivideGraphs(erGraph[i], erGraph[0]);

        nrRatio[i] -> SetLineWidth(2);
        erRatio[i] -> SetLineWidth(2);

        nrAllRatio -> Add(nrRatio[i]);
        erAllRatio -> Add(erRatio[i]);
      }
    }


    TMultiGraph* allGraph[2] = {nrAllGraph, erAllGraph};
    TMultiGraph* allRatio[2] = {nrAllRatio, erAllRatio};

    TCanvas* allCanvas[2];
    TPad* mainPad[2];
    TPad* ratioPad[2];

    // Set up the canvas and plot the graphs.
    for (Int_t i = 0; i < 2; i++){

      allCanvas[i] = new TCanvas(Form("all_cut_canvas%d", i+1), "Cut Effect Analysis", width, height);

      mainPad[i]  = new TPad("mainPad" , "", 0.0, 0.2, 1.0, 1.0);
      ratioPad[i] = new TPad("ratioPad", "", 0.0, 0.0, 1.0, 0.25);

      mainPad[i] -> Draw();   ratioPad[i] -> Draw();

      mainPad[i] -> cd();
      mainPad[i] -> SetMargin(0.1, 0.05, 0.08, 0.02);
      mainPad[i] -> SetGrid();
      allGraph[i] -> Draw("alp");
      allGraph[i] -> GetXaxis() -> SetLabelOffset(999);
      mainPad[i] -> BuildLegend();

      ratioPad[i] -> cd();
      ratioPad[i] -> SetMargin(0.1, 0.05, 0.25, 0.05);
      ratioPad[i] -> SetGrid();
      allRatio[i] -> Draw("alp");
      allRatio[i] -> GetXaxis() -> SetLabelSize(0.11);
      allRatio[i] -> GetXaxis() -> SetTitleSize(0.1);
      allRatio[i] -> GetYaxis() -> SetLabelSize(0.1);
      allRatio[i] -> GetYaxis() -> SetTitleSize(0.12);
      allRatio[i] -> GetYaxis() -> SetTitleOffset(0.2);
    }
  }

  // For each cut, creats a plot of with both the NR and ER percentage. Also plots a ratio of the two.
  if (bothEventOneCut) {

    TMultiGraph* cutGraph[graphCount];
    TGraph*      ratioGraph[graphCount];

    TCanvas* canvas[graphCount];
    TPad* mainPad[graphCount];
    TPad* ratioPad[graphCount];

    for (Int_t i = 0; i < graphCount; i++){

      cutGraph[i] = new TMultiGraph(Form("cut%d_graph", i), "; ; Fração de Eventos");
      cutGraph[i] -> Add(nrGraph[i]);
      cutGraph[i] -> Add(erGraph[i]);

      ratioGraph[i] = DivideGraphs(nrGraph[i], erGraph[i]);
      ratioGraph[i] -> SetTitle(";S1 [PE]; Razão");
      ratioGraph[i] -> SetLineWidth(2);

      // <code-fold> SETTING UP THE CANVAS CONFIGURATIONS FOR PLOTTING THE GRAPHS.
      canvas[i] = new TCanvas(Form("nrer_cut%d_graph", i), "", width, height);

      mainPad[i]  = new TPad("mainPad" , "", 0.0, 0.2, 1.0, 1.0);
      ratioPad[i] = new TPad("ratioPad", "", 0.0, 0.0, 1.0, 0.25);

      mainPad[i] -> Draw();   ratioPad[i] -> Draw();

      mainPad[i] -> cd();
      mainPad[i] -> SetMargin(0.1, 0.05, 0.08, 0.02);
      mainPad[i] -> SetGrid();
      cutGraph[i] -> Draw("alp");
      cutGraph[i] -> GetXaxis() -> SetLabelOffset(999);
      mainPad[i] -> BuildLegend();

      ratioPad[i] -> cd();
      ratioPad[i] -> SetMargin(0.1, 0.05, 0.25, 0.05);
      ratioPad[i] -> SetGrid();
      ratioGraph[i] -> Draw("alp");
      ratioGraph[i] -> GetXaxis() -> SetLabelSize(0.11);
      ratioGraph[i] -> GetXaxis() -> SetTitleSize(0.1);
      ratioGraph[i] -> GetYaxis() -> SetLabelSize(0.1);
      ratioGraph[i] -> GetYaxis() -> SetTitleSize(0.12);
      ratioGraph[i] -> GetYaxis() -> SetTitleOffset(0.3);
      // </code-fold>
    }
  }
}


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
Int_t CountOf2DHistograms(TDirectory* directory){

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

TGraph* DivideGraphs(TGraph* numGraph, TGraph* denGraph){

  Double_t* numY = numGraph -> GetY();
  Double_t* denY = denGraph -> GetY();
  Double_t* numX = numGraph -> GetX();

  Int_t numPoints = numGraph -> GetN();
  Int_t denPoints = denGraph -> GetN();

  Double_t divY[numPoints];

  if (numPoints != denPoints){
    cout << "Invalid operation. TGraphs posses different number of points";
    exit(EXIT_FAILURE);
  } else if (numPoints == denPoints){
    for (Int_t i = 0; i < numPoints; i++)
      divY[i] = numY[i]/denY[i];
  }

  TGraph* divGraph = new TGraph(numPoints, numX, divY);
  divGraph -> SetLineColor(numGraph -> GetLineColor());
  divGraph -> SetMarkerColor(numGraph -> GetMarkerColor());

  return divGraph;
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
