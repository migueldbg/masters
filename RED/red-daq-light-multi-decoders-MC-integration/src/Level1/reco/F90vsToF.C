#include <TCut.h>
#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

/* TFile* CheckFile( TString path_name )
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

/* TDirectory* MakeDirectory( const char* dir_name, const char* dir_title )
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



TCut DefineF90Range( Double_t f90_low, Double_t f90_up ){

  TCut f90_range;

  if ( f90_low > f90_up ) {
    std::cout << "Invalid f90 range: lower boundary is greater than upper boundary.";
    exit(EXIT_FAILURE);
  } else if ( f90_low == f90_up ) {
    f90_range = "";
  } else {
    f90_range = Form( "clusters[0].f90 >= %f && clusters[0].f90 <= %f", f90_low, f90_up );
  }

  return f90_range;
}

TCut DefineS1ChargeRange( Double_t charge_low, Double_t charge_up ){

  TCut charge_range;

  if ( charge_low > charge_up ) {
    std::cout << "Invalid S1 charge range: lower boundary is greater than upper boundary.";
    exit(EXIT_FAILURE);
  } else if ( charge_low == charge_up ) {
    charge_range = "";
  } else {
    charge_range = Form( "clusters[0].charge >= %f && clusters[0].charge <= %f", charge_low, charge_up );
  }

  return charge_range;
}

TCut DefineToFRange( Double_t tof_low, Double_t tof_up ){

  TCut tof_range;

  if ( tof_low > tof_up ) {
    std::cout << "Invalid time of flight range: lower boundary is greater than upper boundary.";
    exit(EXIT_FAILURE);
  } else if ( tof_low == tof_up ) {
    tof_range = "";
  } else {
    tof_range = Form( "2*(xmin[30] - clusters[0].min_x) >= %f && 2*(xmin[30] - clusters[0].min_x) <= %f", tof_low, tof_up );
  }

  return tof_range;
}

TCut DefineQualityCuts( Int_t experiment_cfg ){

  TCut number_of_clusters = "";
  TCut rep = "";

  if ( experiment_cfg == 1 ){
    number_of_clusters = "number_of_clusters == 1";
    rep = "clusters[0].rep == 1";
  } else if ( experiment_cfg == 2 ){
    number_of_clusters = "number_of_clusters == 2";
    rep = "clusters[0].rep == 1 && clusters[1].rep == 1";
  }

  TCut quality_cut = number_of_clusters && rep;

  return quality_cut;
}

/* DefineCuts( Int_t cfg, Double_t f90_low, Double_t f90_up, Double_t charge_low, Double_t charge_up, Double_t tof_low, Double_t tof_up )
 *
 *  Summary of Function:
 *
 *    This function expects that the following functions are defined: DefineF90Range(), DefineS1ChargeRange(), DefineToFRange() and
 *    DefineQualityCuts(). It then call each one with the appropriate parameters and constructs a TCut object that is the sum of
 *    all cuts. The function then return the TCut.
 *
 *  Parameters   : cfg        >> indicates if the data set is from a single phase (1) or dual phase (2) run.
 *                 f90_low    >> lower boundary for the f90 parameter.
 *                 f90_up     >> upper boundary fo the f90 parameter.
 *                 charge_low >> lower boundary for the S1 charge.
 *                 charge_up  >> upper boundary for the S1 charge.
 *                 tof_low    >> lower boundary for the time of flight parameter.
 *                 tof_up     >> upper boundary for the time of flight parameter.
 *
 *  Return Value : TCut final_cut.
 */
TCut DefineCuts( Int_t cfg, Double_t f90_low, Double_t f90_up, Double_t charge_low, Double_t charge_up, Double_t tof_low = 0., Double_t tof_up = 0. ){

  TCut quality_cuts = DefineQualityCuts( cfg );
  TCut charge_range = DefineS1ChargeRange( charge_low, charge_up );
  TCut f90_range    = DefineF90Range( f90_low, f90_up );
  TCut tof_range    = DefineToFRange( tof_low, tof_up );

  TCut total_cut = quality_cuts && charge_range && f90_range && tof_range;

  return total_cut;
}



/* TH1F* GenerateF90Hist( TString file_name, TCut histogram_cuts )
 *
 * Summary of GenerateF90Hist function:
 *
 *    The GenerateF90Hist creates a histogram of f90, with the events taken from a TTree in a file decided by the
 *    user. The user also determines the boundaries of the f90 range to be considered, as well as the boundaries of
 *    the charge values to be considered. The function then returns a normalized histogram.
 *
 * Parameters   : file_name      >> the name of the file containing the events from which to construct the histogram.
 *                histogram_cuts >> contain the cuts to applied when generating the histogram.
 *
 * Return value : TH1F* f90_hist
 */
TH1F* GenerateF90Hist( TString file_name, TCut histogram_cuts ) {

  TFile* file = new TFile(file_name);
  TTree* reco; file -> GetObject("reco", reco);

  reco -> Draw("clusters[0].f90 >> hist", histogram_cuts, "goff");
  TH1F* f90_hist = (TH1F*)gDirectory -> Get("hist");

  if (f90_hist -> GetSumw2N() == 0) f90_hist -> Sumw2(kTRUE);

  //f90_hist = (TH1F*)NormalizeHist(f90_hist, 1.0);

  f90_hist -> SetOption("HIST");

  // Remove the histogram from the current directory (TFile) so that I can close the file and still have access to it.
  f90_hist -> SetDirectory(0);
  file -> Close();

  return f90_hist;
}

/* void WriteF90Hist ( TDirectory* save_dir, TString file_name, TCut histogram_cuts, TString hist_name, TString hist_title )
 *
 *  Summary of Function:
 *
 *    The function receives a run number, from witch a f90 histogram is generated using the GenerateF90Hist() function.
 *    The resulting histogram then receives a name, a title and is then saved into the specified directory.
 *
 *  Parameters   : save_dir   >> directory were the generated histogram will be written to.
 *                 file_name  >> name of the file were the data used to generate the histogram is located.
 *                 histogram_cuts >> cuts to be used when constructing the histogram.
 *                 hist_name  >> name of the histogram to be generated.
 *                 hist_title >> title of the histogram to be generated.
 *
 *  Return Value : void.
 *
 */
void WriteF90Hist( TDirectory* save_dir, TString file_name, TCut histogram_cuts, TString hist_name, TString hist_title ){

  TH1F* hist = GenerateF90Hist( file_name, histogram_cuts );

  hist -> SetName(hist_name);
  hist -> SetTitle(hist_title);

  save_dir -> WriteObject( hist, hist_name, "OverWrite" );
}



void F90vsToF( int run, Int_t number_of_bins, Double_t tof_min = 20, Double_t tof_max = 60 ){

  Int_t    cfg         = 2;
  Double_t f90_min     = 0.;
  Double_t f90_max     = 1.;
  Double_t charge_min  = 50.;
  Double_t charge_max  = 1000.;

  TFile* output_file = CheckFile( Form("analysis_%d.root", run) );

  output_file -> cd();
  TDirectory* histograms_dir = MakeDirectory("histograms", "histograms");
  histograms_dir -> cd();
  TDirectory* tof_dir = MakeDirectory("time_of_flight", "time_of_flight");
  tof_dir -> cd();
  TDirectory* f90_dir = MakeDirectory("f90", "f90");
  f90_dir -> Delete("*;*");

  TString file_name = Form("runs/run_%d.root", run );

  TString hist_name; TString hist_title; TCut histogram_cuts;

  Double_t bin_size = (tof_max - tof_min) / number_of_bins;
  Double_t tof_low; Double_t tof_up;

  for ( Int_t i = 0; i < number_of_bins; i++ ){

    tof_low = ( i      * bin_size) + tof_min;
    tof_up  = ((i + 1) * bin_size) + tof_min;

    histogram_cuts = DefineCuts(cfg, f90_min, f90_max, charge_min, charge_max, tof_low, tof_up);

    hist_name  = Form("f90_histogram_tof%d", i+1);
    hist_title = Form("f90 Distribution (ToF Interval: %4.2f - %4.2f ns); f90", tof_low, tof_up);

    WriteF90Hist(f90_dir, file_name, histogram_cuts, hist_name, hist_title);
  }
}
