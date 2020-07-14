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


TCut DefineF90Range( Double_t f90_low, Double_t f90_up ){

  TCut f90_range;

  if ( f90_low > f90_up ) {
    std::cout << "Invalid f90 range: lower boundary is greater than upper boundary." << std::endl;
    exit(EXIT_FAILURE);
  } else if ( f90_low == f90_up ) {
    f90_range = "";
  } else {
    f90_range = Form( "clusters[0].f90 >= %f && clusters[0].f90 <= %f", f90_low, f90_up );
  }

  return f90_range;
}


TCut DefineS1Range( Double_t s1_low, Double_t s1_up ){

  TCut s1_range;

  if ( s1_low > s1_up ) {
    std::cout << "Invalid S1 charge range: lower boundary is greater than upper boundary." << std::endl;
    exit(EXIT_FAILURE);
  } else if ( s1_low == s1_up ) {
    s1_range = "";
  } else {
    s1_range = Form( "clusters[0].charge >= %f && clusters[0].charge <= %f", s1_low, s1_up );
  }

  return s1_range;
}


TCut DefineS2Range( Double_t s2_low, Double_t s2_up ){

  TCut s2_range;

  if ( s2_low > s2_up ) {
    std::cout << "Invalid S2 charge range: lower boundary is greater than upper boundary." << std::endl;
    exit(EXIT_FAILURE);
  } else if ( s2_low == s2_up ) {
    s2_range = "";
  } else {
    s2_range = Form( "clusters[1].charge >= %f && clusters[1].charge <= %f", s2_low, s2_up );
  }

  return s2_range;
}


TCut DefineSiTelTPCToFRange( Double_t tof_low, Double_t tof_up ){

  TCut tof_range;

  if ( tof_low > tof_up ) {
    std::cout << "Invalid time of flight range: lower boundary is greater than upper boundary." << std::endl;
    exit(EXIT_FAILURE);
  } else if ( tof_low == tof_up ) {
    tof_range = "";
  } else {
    tof_range = Form( "2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time) >= %f && 2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time) <= %f", tof_low, tof_up );
  }

  return tof_range;
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
 *                 f90_low >> lower boundary for the f90 parameter.
 *                 f90_up  >> upper boundary fo the f90 parameter.
 *                 s1_low  >> lower boundary for the S1 charge.
 *                 s1_up   >> upper boundary for the S1 charge.
 *                 s2_low  >> lower boundary for the S2 charge.
 *                 s2_up   >> upper boundary for the S2 charge.
 *                 tof_low >> lower boundary for the time of flight parameter.
 *                 tof_up  >> upper boundary for the time of flight parameter.
 *
 *  Return Value : TCut final_cut.
 */
TCut DefineCuts( Int_t cfg, Double_t f90_low, Double_t f90_up, Double_t s1_low, Double_t s1_up, Double_t s2_low = 0, Double_t s2_up = 0, Double_t tof_low = 0., Double_t tof_up = 0. ){

  TCut quality_cuts = DefineQualityCuts( cfg );
  TCut s1_range     = DefineS1Range( s1_low, s1_up );
  TCut s2_range     = DefineS2Range( s2_low, s2_up );
  TCut f90_range    = DefineF90Range( f90_low, f90_up );
  TCut tof_range    = DefineSiTelTPCToFRange( tof_low, tof_up );

  TCut total_cut = quality_cuts && s1_range && s2_range && f90_range && tof_range;

  return total_cut;
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


/* TH1F* GenerateF90Histogram()
 *
 * Summary of Function:
 *
 *    The Generatef90Histogram creates a histogram of f90, with the events taken from a TTree in a file decided by the
 *    user. The user also provides the parameters of the histogram, such as number of bins, minimum and maximum
 *    f90 value. The histogram is constructed using the cuts provided and the function then returns it.
 *
 * Parameters   : file_name      >> the name of the file containing the events from which to construct the histogram.
 *                histogram_cuts >> contain the cuts to applied when generating the histogram.
 *                num_bins       >> number of bins of the histogram;
 *                f90_min        >> minimum f90 value.
 *                f90_max        >> maximum f90 value.
 *                normalize      >> indicates if the output histogram is to be normalized.
 *
 * Return value : TH1F* f90_hist
 */
TH1F* GenerateF90Histogram( TString file_name, TCut hist_cuts, Int_t num_bins, Double_t f90_min, Double_t f90_max, bool normalize ){

    TFile* file = new TFile(file_name);
    TTree* reco; file -> GetObject("reco", reco);

    TH1F* f90_hist = new TH1F("f90_hist", "f90; f90", num_bins, f90_min, f90_max);
    reco -> Project( "f90_hist", "clusters[0].f90", hist_cuts);

    if (f90_hist -> GetSumw2N() == 0) f90_hist -> Sumw2(kTRUE);

    if ( normalize ) f90_hist = (TH1F*) NormalizeHistogram(f90_hist, 1.0);

    f90_hist -> SetOption("HIST");
    f90_hist -> SetDirectory(0);
    file -> Close();

    return f90_hist;
}


/* void WriteF90Histogram ()
 *
 *  Summary of Function:
 *
 *    The function receives a run number, from witch a f90 histogram is generated using the GenerateF()Histogram() function.
 *    The resulting histogram then receives a name, a title and is then saved into the specified directory.
 *
 *  Parameters   : save_dir       >> directory were the generated histogram will be written to.
 *                 file_name      >> name of the file were the data used to generate the histogram is located.
 *                 numb_bins      >> number of bins of the histogram.
 *                 f90_min        >> minimum f90 value.
 *                 f90_max        >> maximum f90 value.
 *                 normalize      >> indicates if the output histogram is to be normalized.
 *                 histogram_cuts >> cuts to be used when constructing the histogram.
 *                 hist_name      >> name of the histogram to be generated.
 *                 hist_title     >> title of the histogram to be generated.
 *
 *  Return Value : void.
 *
 */
void WriteF90Histogram( TDirectory* save_dir, TString file_name, TCut histogram_cuts, Int_t num_bins, Double_t f90_min,
                       Double_t f90_max, bool normalize, TString hist_name, TString hist_title ){

  TH1F* hist = GenerateF90Histogram( file_name, histogram_cuts, num_bins, f90_min, f90_max, normalize );

  hist -> SetName(hist_name);
  hist -> SetTitle(hist_title);

  save_dir -> WriteObject( hist, hist_name, "OverWrite" );
}


/* TH1F* GenerateS1Histogram()
 *
 * Summary of Function:
 *
 *    The GenerateS1Histogram creates a histogram of S1, with the events taken from a TTree in a file decided by the
 *    user. The user also provides the parameters of the histogram, such as number of bins, minimum charge and
 *    maximum charge. The histogram is constructed using the cuts provided and the function then returns it.
 *
 * Parameters   : file_name      >> the name of the file containing the events from which to construct the histogram.
 *                histogram_cuts >> contain the cuts to applied when generating the histogram.
 *                num_bins       >> number of bins of the histogram;
 *                s1_min         >> minimum S1 value.
 *                s1_max         >> maximum S1 value.
 *                normalize      >> indicates if the output histogram is to be normalized.
 *
 * Return value : TH1F* s1_hist
 */
TH1F* GenerateS1Histogram( TString file_name, TCut hist_cuts, Int_t num_bins, Double_t s1_min, Double_t s1_max, bool normalize ){

    TFile* file = new TFile(file_name);
    TTree* reco; file -> GetObject("reco", reco);

    TH1F* s1_hist = new TH1F("s1_hist", "S1 Charge; S1 (PE)", num_bins, s1_min, s1_max);
    reco -> Project( "s1_hist", "clusters[0].charge", hist_cuts);

    if (s1_hist -> GetSumw2N() == 0) s1_hist -> Sumw2(kTRUE);

    if ( normalize ) s1_hist = (TH1F*) NormalizeHistogram(s1_hist, 1.0);

    s1_hist -> SetOption("HIST");
    s1_hist -> SetDirectory(0);
    file -> Close();

    return s1_hist;
}


/* void WriteS1Histogram ()
 *
 *  Summary of Function:
 *
 *    The function receives a run number, from witch a s1 histogram is generated using the GenerateS1Histogram() function.
 *    The resulting histogram then receives a name, a title and is then saved into the specified directory.
 *
 *  Parameters   : save_dir       >> directory were the generated histogram will be written to.
 *                 file_name      >> name of the file were the data used to generate the histogram is located.
 *                 numb_bins      >> number of bins of the histogram.
 *                 s1_min         >> minimum S1 value.
 *                 s1_max         >> maximum S1 value.
 *                 normalize      >> indicates if the output histogram is to be normalized.
 *                 histogram_cuts >> cuts to be used when constructing the histogram.
 *                 hist_name      >> name of the histogram to be generated.
 *                 hist_title     >> title of the histogram to be generated.
 *
 *  Return Value : void.
 *
 */
void WriteS1Histogram( TDirectory* save_dir, TString file_name, TCut histogram_cuts, Int_t num_bins, Double_t s1_min,
                       Double_t s1_max, bool normalize, TString hist_name, TString hist_title ){

  TH1F* hist = GenerateS1Histogram( file_name, histogram_cuts, num_bins, s1_min, s1_max, normalize );

  hist -> SetName(hist_name);
  hist -> SetTitle(hist_title);

  save_dir -> WriteObject( hist, hist_name, "OverWrite" );
}


/* TH1F* GenerateS2Histogram()
 *
 * Summary of Function:
 *
 *    The GenerateS1Histogram creates a histogram of S2, with the events taken from a TTree in a file decided by the
 *    user. The user also provides the parameters of the histogram, such as number of bins, minimum charge and
 *    maximum charge. The histogram is constructed using the cuts provided and the function then returns it.
 *
 * Parameters   : file_name      >> the name of the file containing the events from which to construct the histogram.
 *                histogram_cuts >> contain the cuts to applied when generating the histogram.
 *                num_bins       >> number of bins of the histogram;
 *                s2_min         >> minimum S1 value.
 *                s2_max         >> maximum S1 value.
 *                normalize      >> indicates if the output histogram is to be normalized.
 *
 * Return value : TH1F* s2_hist
 */
TH1F* GenerateS2Histogram( TString file_name, TCut hist_cuts, Int_t num_bins, Double_t s2_min, Double_t s2_max, bool normalize ){

    TFile* file = new TFile(file_name);
    TTree* reco; file -> GetObject("reco", reco);

    TH1F* s2_hist = new TH1F("s2_hist", "S2 Charge; S2 (PE)", num_bins, s2_min, s2_max);
    reco -> Project( "s2_hist", "clusters[1].charge", hist_cuts);

    if (s2_hist -> GetSumw2N() == 0) s2_hist -> Sumw2(kTRUE);

    if ( normalize ) s2_hist = (TH1F*) NormalizeHistogram(s2_hist, 1.0);

    s2_hist -> SetOption("HIST");
    s2_hist -> SetDirectory(0);
    file -> Close();

    return s2_hist;
}


/* void WriteS2Histogram ()
 *
 *  Summary of Function:
 *
 *    The function receives a run number, from witch a S2 histogram is generated using the GenerateS2Histogram() function.
 *    The resulting histogram then receives a name, a title and is then saved into the specified directory.
 *
 *  Parameters   : save_dir       >> directory were the generated histogram will be written to.
 *                 file_name      >> name of the file were the data used to generate the histogram is located.
 *                 numb_bins      >> number of bins of the histogram.
 *                 s2_min         >> minimum S2 value.
 *                 s2_max         >> maximum S2 value.
 *                 normalize      >> indicates if the output histogram is to be normalized.
 *                 histogram_cuts >> cuts to be used when constructing the histogram.
 *                 hist_name      >> name of the histogram to be generated.
 *                 hist_title     >> title of the histogram to be generated.
 *
 *  Return Value : void.
 *
 */
void WriteS2Histogram( TDirectory* save_dir, TString file_name, TCut histogram_cuts, Int_t num_bins, Double_t s2_min,
                       Double_t s2_max, bool normalize, TString hist_name, TString hist_title ){

  TH1F* hist = GenerateS2Histogram( file_name, histogram_cuts, num_bins, s2_min, s2_max, normalize );

  hist -> SetName(hist_name);
  hist -> SetTitle(hist_title);

  save_dir -> WriteObject( hist, hist_name, "OverWrite" );
}


/* TH1F* GenerateDriftTimeHistogram()
 *
 * Summary of Function:
 *
 *    The GenerateS1Histogram creates a histogram of drift time, with the events taken from a TTree in a file
 *    decided by the user. The user also provides the parameters of the histogram, such as number of bins,
 *    minimum and maximum drift time. The histogram is constructed using the cuts provided and the function then
 *    returns it.
 *
 * Parameters   : file_name      >> the name of the file containing the events from which to construct the histogram.
 *                histogram_cuts >> contain the cuts to applied when generating the histogram.
 *                num_bins       >> number of bins of the histogram;
 *                dt_min         >> minimum drift time value (in microseconds).
 *                dt_max         >> maximum drift time value (in microseconds).
 *                normalize      >> indicates if the output histogram is to be normalized.
 *
 * Return value : TH1F* dt_hist
 */
TH1F* GenerateDriftTimeHistogram( TString file_name, TCut hist_cuts, Int_t num_bins, Double_t dt_min, Double_t dt_max, bool normalize ){

    TFile* file = new TFile(file_name);
    TTree* reco; file -> GetObject("reco", reco);

    TH1F* dt_hist = new TH1F("dt_hist", "Drift Time; Time (#mus)", num_bins, dt_min, dt_max);
    reco -> Project( "dt_hist", "(clusters[1].start_time - clusters[0].start_time)*(2/1000)", hist_cuts);

    if (dt_hist -> GetSumw2N() == 0) dt_hist -> Sumw2(kTRUE);

    if ( normalize ) dt_hist = (TH1F*) NormalizeHistogram(dt_hist, 1.0);

    dt_hist -> SetOption("HIST");
    dt_hist -> SetDirectory(0);
    file -> Close();

    return dt_hist;
}


/* void WriteDriftTimeHistogram ()
 *
 *  Summary of Function:
 *
 *    The function receives a run number, from witch a drift time histogram is generated using the GenerateDriftTimeHistogram()
 *    function. The resulting histogram then receives a name, a title and is then saved into the specified directory.
 *
 *  Parameters   : save_dir       >> directory were the generated histogram will be written to.
 *                 file_name      >> name of the file were the data used to generate the histogram is located.
 *                 numb_bins      >> number of bins of the histogram.
 *                 dt_min         >> minimum drift time value.
 *                 dt_max         >> maximum drift time value.
 *                 normalize      >> indicates if the output histogram is to be normalized.
 *                 histogram_cuts >> cuts to be used when constructing the histogram.
 *                 hist_name      >> name of the histogram to be generated.
 *                 hist_title     >> title of the histogram to be generated.
 *
 *  Return Value : void.
 *
 */
void WriteDriftTimeHistogram( TDirectory* save_dir, TString file_name, TCut histogram_cuts, Int_t num_bins, Double_t dt_min,
                       Double_t dt_max, bool normalize, TString hist_name, TString hist_title ){

  TH1F* hist = GenerateDriftTimeHistogram( file_name, histogram_cuts, num_bins, dt_min, dt_max, normalize );

  hist -> SetName(hist_name);
  hist -> SetTitle(hist_title);

  save_dir -> WriteObject( hist, hist_name, "OverWrite" );
}


/* TH1F* GenerateToFHistogram()
 *
 * Summary of Function:
 *
 *    The GenerateToFHistogram creates a histogram of time of flighy, with the events taken from a TTree in a file
 *    decided by the user. The user also provides the parameters of the histogram, such as number of bins,
 *    minimum and maximum time of flight. The histogram is constructed using the cuts provided and the function then
 *    returns it.
 *
 * Parameters   : file_name      >> the name of the file containing the events from which to construct the histogram.
 *                histogram_cuts >> contain the cuts to applied when generating the histogram.
 *                num_bins       >> number of bins of the histogram;
 *                tof_min        >> minimum drift time value (in microseconds).
 *                tof_max        >> maximum drift time value (in microseconds).
 *                normalize      >> indicates if the output histogram is to be normalized.
 *
 * Return value : TH1F* dt_hist
 */
 TH1F* GenerateToFHistogram( TString file_name, TCut hist_cuts, Int_t num_bins, Double_t tof_min, Double_t tof_max, bool normalize ){

  TFile* file = new TFile(file_name);
  TTree* reco; file -> GetObject("reco", reco);

  TH1F* tof_hist = new TH1F("tof_hist", "Time of Flight (TPC and SiTel); Time (ns)", num_bins, tof_min, tof_max);
  reco -> Project( "tof_hist", "2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time)", hist_cuts);

  if (tof_hist -> GetSumw2N() == 0) tof_hist -> Sumw2(kTRUE);

  if ( normalize ) tof_hist = (TH1F*) NormalizeHistogram(tof_hist, 1.0);

  tof_hist -> SetOption("HIST");
  tof_hist -> SetDirectory(0);
  file -> Close();

  return tof_hist;
}


/* void WriteToFHistogram ()
 *
 *  Summary of Function:
 *
 *    The function receives a run number, from witch a time of flight histogram is generated using the GenerateToFHistogram()
 *    function. The resulting histogram then receives a name, a title and is then saved into the specified directory.
 *
 *  Parameters   : save_dir       >> directory were the generated histogram will be written to.
 *                 file_name      >> name of the file were the data used to generate the histogram is located.
 *                 numb_bins      >> number of bins of the histogram.
 *                 tof_min        >> minimum time of flight value.
 *                 tog_max        >> maximum time of flight value.
 *                 normalize      >> indicates if the output histogram is to be normalized.
 *                 histogram_cuts >> cuts to be used when constructing the histogram.
 *                 hist_name      >> name of the histogram to be generated.
 *                 hist_title     >> title of the histogram to be generated.
 *
 *  Return Value : void.
 *
 */
void WriteToFHistogram( TDirectory* save_dir, TString file_name, TCut histogram_cuts, Int_t num_bins, Double_t tof_min,
                       Double_t tof_max, bool normalize, TString hist_name, TString hist_title ){

  TH1F* hist = GenerateToFHistogram( file_name, histogram_cuts, num_bins, tof_min, tof_max, normalize );

  hist -> SetName(hist_name);
  hist -> SetTitle(hist_title);

  save_dir -> WriteObject( hist, hist_name, "OverWrite" );
}
