#include "sidutility.h"

#include <TChain.h>
#include <TCut.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>
#include <vector>


TStyle* SetSidStyle(){
  auto sidStyle = new TStyle("sidStyle", "Sid's Style");

  sidStyle -> SetCanvasBorderMode(0);
  sidStyle -> SetCanvasColor(0);

  sidStyle -> SetLabelFont(102, "xyz");

  sidStyle -> SetPadBorderMode(0);
  sidStyle -> SetPadColor(0);

  sidStyle -> SetPalette(kSunset);

  sidStyle -> SetStatColor(0);

  sidStyle -> SetTitleAlign(23);
  sidStyle -> SetTitleBorderSize(0);
  sidStyle -> SetTitleColor(0, "t");
  sidStyle -> SetTitleFont(102, "xyz");
  sidStyle -> SetTitleFont(102, "t");
  sidStyle -> SetTitleX(0.5);

  //sidStyle -> SetOptStat(0);

  return sidStyle;
}

/* TFile* CheckFile()
 *
 * Summary of CheckFile function:
 *
 *    The CheckFile function checks if a given file already exists and returns this information to the user.
 *    If the file exists, CheckFile returns it. If not, it creates it and then returns it.
 *
 * Parameters   : path_name >> The name of the desired file.
 *
 * Return value : TFile* file.
 */
TFile* CheckFile( TString path_name ){

  if (gSystem -> AccessPathName(path_name)) {
    std::cout << "The " << path_name << " file does not exist. Creating it..." << std::endl;
  } else {
    std::cout << "Opening file: " << path_name << std::endl;
  }

  TFile *file = new TFile(path_name, "UPDATE");
  if ( !(file -> IsOpen()) ) { std::cout << "Unable to open file." << std::endl;}

  return file;
}


/* TDirectory* MakeDirectory()
  *
  * Summary of MakeDirectory function:
  *
  *    The MakeDirectory function checks wether a given directory currently exists. If it does, it returns a
  *    pointer to it. If it doesn't, the function creates the directory with the given name in the current
  *    path and returns a pointer to it.
  *
  * Parameters   : dir_name  >> name of the desired directory.
  *                dir_title >> title of the directory, for bookeeping purpouses.
  *
  * Return value : TDirectory* directory.
  */
TDirectory* MakeDirectory( const char* dir_name, const char* dir_title ){

  TDirectory *directory = gDirectory -> GetDirectory(dir_name);
  if (!directory) directory = gDirectory -> mkdir(dir_name, dir_title);

  return directory;
}


/* TFile* MergeRuns( std::vector<int> runs )
 *
 *  Summary of Function:
 *
 *    The function takes a set of root files and merges the 'reco' tree of each one into a single
 *    TTree object. The resulting merged tree is then saved into a root file.
 *
 *  Parameters   : runs >> a vector of integers containing the number of each run to be merged.
 *
 *  Return Value : TFile* file.
 */
TFile* MergeRuns( std::vector<int> runs ){

  TChain* reco_chain = new TChain("reco");
  TString file_name;
  for (Int_t i = 0; i < runs.size(); i++){
    file_name = Form("runs/run_%d.root", runs.at(i));
    reco_chain -> Add(file_name.Data());
  }

  TString output_file_name = Form("runs/run_%d%d.root", runs.front(), runs.back());
  reco_chain -> Merge(output_file_name.Data(), "fast");

  TFile* file = new TFile(output_file_name, "update");
  return file;
}


TCut DefineF90Range( Double_t f90Low, Double_t f90Upp ){

  TCut f90Range;

  if ( f90Low > f90Upp ) {
    std::cout << "Invalid f90 range: lower boundary is greater than upper boundary." << std::endl;
    exit(EXIT_FAILURE);
  } else if ( f90Low == f90Upp ) {
    f90Range = "";
  } else {
    f90Range = Form( "clusters[0].f90 >= %f && clusters[0].f90 <= %f", f90Low, f90Upp );
  }

  return f90Range;
}


TCut DefineS1Range( Double_t s1Low, Double_t s1Upp ){

  TCut s1Range;

  if ( s1Low > s1Upp ) {
    std::cout << "Invalid S1 charge range: lower boundary is greater than upper boundary." << std::endl;
    exit(EXIT_FAILURE);
  } else if ( s1Low == s1Upp ) {
    s1Range = "";
  } else {
    s1Range = Form( "clusters[0].charge >= %f && clusters[0].charge <= %f", s1Low, s1Upp );
  }

  return s1Range;
}


TCut DefineS2Range( Double_t s2Low, Double_t s2Upp ){

  TCut s2Range;

  if ( s2Low > s2Upp ) {
    std::cout << "Invalid S2 charge range: lower boundary is greater than upper boundary." << std::endl;
    exit(EXIT_FAILURE);
  } else if ( s2Low == s2Upp ) {
    s2Range = "";
  } else {
    s2Range = Form( "clusters[1].charge >= %f && clusters[1].charge <= %f", s2Low, s2Upp );
  }

  return s2Range;
}


TCut DefineSiTelTPCToFRange( Double_t tofLow, Double_t tofUpp ){

  TCut tofRange;

  if ( tofLow > tofUpp ) {
    std::cout << "Invalid time of flight range: lower boundary is greater than upper boundary." << std::endl;
    exit(EXIT_FAILURE);
  } else if ( tofLow == tofUpp ) {
    tofRange = "";
  } else {
    tofRange = Form( "2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time) >= %f && 2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time) <= %f", tofLow, tofUpp );
  }

  return tofRange;
}


TCut DefineQualityCuts( Int_t experiment_cfg ){

  TCut number_of_clusters = "";
  TCut rep = "";

  if ( experiment_cfg == 1 ){
    number_of_clusters = "number_of_clusters >= 1";
    rep = "clusters[0].rep == 1";
  } else if ( experiment_cfg == 2 ){
    number_of_clusters = "number_of_clusters >= 2";
    rep = "clusters[0].rep == 1 && clusters[1].rep == 1";
  }

  TCut quality_cut = number_of_clusters && rep;

  return quality_cut;
}


/* DefineCuts()
 *
 *  Summary of Function:
 *
 *    This function expects that the following functions are defined: DefineF90Range(), DefineS1ChargeRange(), DefineToFRange() and
 *    DefineQualityCuts(). It then call each one with the appropriate parameters and constructs a TCut object that is the sum of
 *    all cuts. The function then return the TCut.
 *
 *  Parameters   : cfg     >> indicates if the data set is from a single phase (1) or dual phase (2) run.
 *                 f90Low >> lower boundary for the f90 parameter.
 *                 f90Upp  >> upper boundary fo the f90 parameter.
 *                 s1Low  >> lower boundary for the S1 charge.
 *                 s1Upp   >> upper boundary for the S1 charge.
 *                 s2Low  >> lower boundary for the S2 charge.
 *                 s2Upp   >> upper boundary for the S2 charge.
 *                 tofLow >> lower boundary for the time of flight parameter.
 *                 tofUpp  >> upper boundary for the time of flight parameter.
 *
 *  Return Value : TCut final_cut.
 */
TCut DefineCuts( Int_t cfg, Double_t f90Low = 0, Double_t f90Upp = 0, Double_t s1Low = 0, Double_t s1Upp = 0, Double_t s2Low = 0, Double_t s2Upp = 0, Double_t tofLow = 0., Double_t tofUpp = 0. ){

  TCut qualityCut = DefineQualityCuts( cfg );
  TCut s1Range     = DefineS1Range( s1Low, s1Upp );
  TCut s2Range     = DefineS2Range( s2Low, s2Upp );
  TCut f90Range    = DefineF90Range( f90Low, f90Upp );
  TCut tofRange    = DefineSiTelTPCToFRange( tofLow, tofUpp );

  TCut cut = qualityCut && s1Range && s2Range && f90Range && tofRange;

  return cut;
}


/* TCutG* LowBeGraphCut ( Int_t run, const char* name, Int_t sigma )
 *
 * Summary of LowBeGraphCut function:
 *
 *    This function takes the TGraph object saved in a root file (LowBeCut.root) and converts it into
 *    a TCutG* object. This object is then returned. The graph saved in LowBeCut.root is generated by
 *    the LowBeStudy.C macro. There is two graphs generated for each set of runs considered: one for
 *    the one-sigma confidence ellipse and another for the two-sigma confidence ellipse.
 *
 * Parameters   : run >> indicates with respect to witch run the TGraph object was created.
 *                name >> the name of the TCutG* object returned.
 *                sigma >> indicates witch graph to return (one-sigma or two-sigma)
 *
 * Return value : TCutG* bGraphCut
 */
TCutG* LowBeGraphCut( Int_t run, const char* name = "LowBeCut", Int_t sigma = 2 ){

  TFile*  bCutFile    = new TFile("LowBeCut.root", "UPDATE");
  TCutG*  bGraphCut   = new TCutG(name);
  TGraph* sourceGraph = new TGraph();

  TString bCutName;
  if ( sigma == 1 ){
    bCutName = Form("LowBeOneSigmaCut_%d", run);
  } else if ( sigma == 2 ){
    bCutName = Form("LowBeTwoSigmaCut_%d", run);
  }

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

  return bGraphCut;
}


TH2F* Bananator(int run, int nbin=100, int xmin=0, int xmax=20000, int ymin=0, int ymax=20000){

  TFile *file = new TFile( runsDirectoryPath + Form("/run_%d.root", run), "read" );

  TTree *reco = (TTree*)file -> Get("reco");

  // Style
  gStyle->SetOptFit();


  // ************** Bananator ************** //

  TString hist_title = "#DeltaE / E Spectrum ";

  TH2F *hist = new TH2F("hist", hist_title, nbin, xmin, xmax, nbin, ymin, ymax);

  reco -> Draw("baseline_mean[30] - ymin[30]:baseline_mean[31] - ymin[31] >> hist", "", "goff");

  hist -> GetXaxis() -> SetTitle("E [ADC counts]");
  hist -> GetXaxis() -> SetTitleOffset(1.16);
  hist -> GetYaxis() -> SetTitle("#Delta E [ADC counts]");
  hist -> GetYaxis() -> SetTitleOffset(1.36);

  //hist -> DrawCopy("colz");

  hist -> SetDirectory(0);
  return hist;
}



/* TH1* NormalizeHist()
 *
 * Summary of NormalizeHist function:
 *
 *    The NormalizeHist function takes a TH1* histogram, makes a copy and normalizes the copy to the
 *    desired value. It then returns the normalized copy of histogram.
 *
 * Parameters   : hist >> histogram to be normalized.
 *                norm >> value to normalize the histogram to (default = 1.0).
 *
 * Return value : TH1* normalized_hist.
 */
TH1* NormalizeHistogram( TH1* hist, Double_t norm = 1. ){

  Double_t scale = norm / (hist -> Integral());

  TH1* normalized_hist = (TH1*)hist -> Clone();
  normalized_hist -> Scale(scale);

  return normalized_hist;
}
