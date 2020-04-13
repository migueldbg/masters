/* File: CreateF90Histograms.C (ROOT macro).
 *
 * Author: Miguel Del Ben Galdiano (miguel.galdiano@gmail.com).
 * Date of creation : February 4 2020.
 *
 * Summary of File:
 *
 *    The goal of this macro is to generate a number of f90 histograms for a given run. We utilize three different data sets:
 *    actual data, a electron recoil (ER) only Monte Carlo (MC) simulation and a neutron recoil (NR) only MC simulation. The
 *    user defines the max total charge value of the events to be considered and also determines the bin size (in PE) to be
 *    utilized when constructing the histograms. The program then generates a number of histograms equal to the max charge
 *    divided by the size of the bin, each one moving along the total charge (starts at 0, goes up to max total charge in
 *    incriments defined by the bin size). The number of histograms must be a whole integer, so this must be taken into account
 *    when deciding the values of max charge and bin size.The generated histograms are then saved in a .root file, separated in
 *    different folders to facilitate later access, if necessary.
 */

#include <array>
#include <cstring>

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
    tof_range = Form( "xmin[30] - clusters[0].min_x >= %f && xmin[30] - clusters[0].min_x <= %f", tof_low, tof_up );
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

/* TH1* NormalizeHist( TH1* hist, Double_t norm = 1. )
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
TH1* NormalizeHist( TH1* hist, Double_t norm = 1. ){

  Double_t scale = norm / (hist -> Integral());

  TH1* normalized_hist = (TH1*)hist -> Clone();
  normalized_hist -> Scale(scale);

  return normalized_hist;
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

/* TDirectory* MakeDirStruct( TFile* file, bool isMC, Int_t ERorNR )
 *
 *  Summary of Function:
 *
 *    The function creates the desire directory structure to save the histograms to be generated. It considers
 *    if the source data is from a MC simulation or not, and also creates different folder names depending if
 *    the selected events are electron or nuclear recoils.
 *
 *  Parameters   : file >> file were the structure will be created in.
 *                 isMC >> indicates if the source data comes from a MC simulation or not.
 *                 ERorNR >> indicates the recoil type (0 for ER, 1 for NR).
 *
 *  Return value : TDirectory* directory
 */
TDirectory* MakeDirStruct( TFile* file, bool isMC, Int_t ERorNR ){

  file -> cd();
  TDirectory* histograms_dir  = MakeDirectory("histograms", "histograms");

  histograms_dir -> cd();
  TDirectory* f90_histograms_dir = MakeDirectory("f90", "f90");

  f90_histograms_dir -> cd();
  TDirectory* type_dir;

  if ( !isMC ) {
    type_dir = MakeDirectory( "data", "data" );
  } else if ( isMC ){
    type_dir = MakeDirectory( "monte_carlo", "monte_carlo" );
  }

  type_dir -> cd();
  TDirectory* directory;
  if ( ERorNR == 0 ) { directory = MakeDirectory("ER","ER"); }
  else if ( ERorNR == 1 ) { directory = MakeDirectory("NR","NR"); }

  directory -> Delete("*;*"); // Clears the directory in case there are objects already there.

  return directory;
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



/* void CreateF90ERHistograms ( int run, Int_t exp_cfg, bool isMC = false, bool tof_cut = false, Int_t number_of_bins = 50, Double_t min_charge = 0., Double_t max_charge = 1000.)
 *
 *  Summary of Function:
 *
 *    This function constructs a number of electron recoil f90 histograms within a total S1 charge range defined by the user. The
 *    determination of wether an event is classified as an ER can be done by f90 cuts (tof_cut = false) or ToF cuts (tof_cut = true).
 *    Changing both of these values must be done within the code. The function can also take into consideration wether the run is
 *    single phase or dual phase and call for the appropriate quality cuts. The resulting histograms are then saved in a root file
 *    named analysis_#run number#.root.
 *
 *  Parameters   : run >> run number.
 *                 data_type >> in dicates if the source is taken from data or MC simulation.
 *                 exp_cfg >> indicates if the run was single phase or dual phase.
 *                 tof_cut >> indicates wheter ToF cuts are to be applied for ER event selection.
 *                 number_of_bins >> total number of histograms that will be generated, each for a different charge range.
 *                 min_charge >> the lower boundary of the S1 total charge range to be considered.
 *                 max_charge >> the upper boundary of the S1 total charge range to be considered.
 *
 *  Return Value : void.
 */
void CreateF90ERHistograms ( int run, Int_t exp_cfg, bool isMC = false, bool tof_cut = false, Int_t number_of_bins = 50,
                             Double_t min_charge = 0., Double_t max_charge = 1000.){

  TFile* output_file = CheckFile( Form("analysis_%d.root", run) );

  TDirectory* dir = MakeDirStruct( output_file, isMC, 0 );

  TString run_file_name;
  if ( !isMC ) { run_file_name = Form("runs/run_%d.root", run); }
  else if ( isMC ) { run_file_name = Form("runs/run_%d_MCER.root", run); }

  // ELECTRON RECOIL EVENT SELECTION PARAMETERS //
  Double_t f90_min; Double_t f90_max;
  Double_t tof_min; Double_t tof_max; // If the ToF border values are equal, the code doesn't apply a ToF cut.
  // ------------------------------------------ //

  if ( tof_cut ){
    f90_min = 0.0;  f90_max = 1.0;
    tof_min = 0.0;  tof_max = 50;
  } else {
    f90_min = 0.2;  f90_max = 0.4;
    tof_min = 0.0;  tof_max = 0.0;
  }

  Double_t bin_size = ( max_charge - min_charge ) / number_of_bins;
  Double_t charge_low;   Double_t charge_up;
  TString hist_name;     TString hist_title;
  TCut histogram_cuts;

  for (int i = 0; i < number_of_bins; i++){

    // Calculate the boundaries of the current bin.
    charge_low = ( i       * bin_size ) + min_charge;
    charge_up  = ( (i + 1) * bin_size ) + min_charge;

    histogram_cuts = DefineCuts(exp_cfg, f90_min, f90_max, charge_low, charge_up, tof_min, tof_max);

    hist_name  = Form( "f90_histogram_%ser_%d", (isMC)?"mc":"", i+1 );
    hist_title = Form( "f90 Distribution (%s, Charge Interval: %d - %d PE); f90", (isMC)?"MC ER":"ER", (int) charge_low, (int) charge_up );

    WriteF90Hist( dir, run_file_name, histogram_cuts, hist_name, hist_title );
  }

  output_file -> Close();
}

/* void CreateF90NRHistograms ( int run, Int_t exp_cfg, bool isMC = false, bool tof_cut = false, Int_t number_of_bins = 50, Double_t min_charge = 0., Double_t max_charge = 1000.)
 *
 *  Summary of Function:
 *
 *    This function constructs a number of nuclear recoil f90 histograms within a total S1 charge range defined by the user. The
 *    determination of whether an event is classified as an NR can be done by f90 cuts (tof_cut = false) or ToF cuts (tof_cut = true).
 *    Changing both of these values must be done within the code. The function can also take into consideration wether the run is
 *    single phase or dual phase and call for the appropriate quality cuts. The resulting histograms are then saved in a root file
 *    named analysis_#run number#.root.
 *
 *  Parameters   : run >> run number.
 *                 data_type >> in dicates if the source is taken from data or MC simulation.
 *                 exp_cfg >> indicates if the run was single phase or dual phase.
 *                 tof_cut >> indicates wheter ToF cuts are to be applied for ER event selection.
 *                 number_of_bins >> total number of histograms that will be generated, each for a different charge range.
 *                 min_charge >> the lower boundary of the S1 total charge range to be considered.
 *                 max_charge >> the upper boundary of the S1 total charge range to be considered.
 *
 *  Return Value : void.
 */
void CreateF90NRHistograms ( int run, Int_t exp_cfg, bool isMC = false, bool tof_cut = false, Int_t number_of_bins = 50,
                             Double_t min_charge = 0., Double_t max_charge = 1000.){

  TFile* output_file = CheckFile( Form("analysis_%d.root", run) );

  TDirectory* dir = MakeDirStruct( output_file, isMC, 1 );

  TString run_file_name;
  if ( !isMC ) { run_file_name = Form("runs/run_%d.root", run); }
  else if ( isMC ) { run_file_name = Form("runs/run_%d_MCNR.root", run); }

  // NUCLEAR RECOIL EVENT SELECTION PARAMETERS //
  Double_t f90_min; Double_t f90_max;
  Double_t tof_min; Double_t tof_max; // If the ToF border values are equal, the code doesn't apply a ToF cut.
  // ----------------------------------------- //

  if ( tof_cut ){
    f90_min = 0.0;  f90_max = 1.0;
    tof_min = 0.0;  tof_max = 50;
  } else if ( !tof_cut ) {
    f90_min = 0.4;  f90_max = 0.7;
    tof_min = 0.0;  tof_max = 0.0;
  }

  Double_t bin_size = ( max_charge - min_charge ) / number_of_bins;
  Double_t charge_low;   Double_t charge_up;
  TString hist_name;     TString hist_title;
  TCut histogram_cuts;

  for (int i = 0; i < number_of_bins; i++){

    // Calculate the boundaries of the current bin.
    charge_low = ( i       * bin_size ) + min_charge;
    charge_up  = ( (i + 1) * bin_size ) + min_charge;

    histogram_cuts = DefineCuts(exp_cfg, f90_min, f90_max, charge_low, charge_up, tof_min, tof_max);

    hist_name  = Form( "f90_histogram_%snr_%d", (isMC)?"mc":"", i+1 );
    hist_title = Form( "f90 Distribution (%s, Charge Interval: %d - %d PE); f90", (isMC)?"MC NR":"NR", (int) charge_low, (int) charge_up );

    WriteF90Hist( dir, run_file_name, histogram_cuts, hist_name, hist_title );
  }

  output_file -> Close();
}

void CreateF90Histograms (int run, Double_t bin_size = 20., Double_t max_charge = 1000., Double_t min_charge = 20.){

  TString file_name = Form("analysis_%d.root", run);
  TFile* hist_file = CheckFile(file_name);

  // --------------------- CREATING NECESSARY DIRECTORIES ---------------------- //
  hist_file -> cd();
  TDirectory* histograms_dir     = MakeDirectory("histograms", "histograms");

  histograms_dir -> cd();
  TDirectory* f90_histograms_dir = MakeDirectory("f90", "f90");

  f90_histograms_dir -> cd();
  TDirectory* data_dir           = MakeDirectory("data", "data");
  TDirectory* monte_carlo_dir    = MakeDirectory("monte_carlo", "monte_carlo");

  data_dir -> cd();
  TDirectory* da_both_dir         = MakeDirectory("both","both");
  TDirectory* da_er_dir           = MakeDirectory("ER","ER");
  TDirectory* da_nr_dir           = MakeDirectory("NR","NR");

  monte_carlo_dir -> cd();
  TDirectory* mc_er_dir           = MakeDirectory("ER","ER");
  TDirectory* mc_nr_dir           = MakeDirectory("NR","NR");
  // --------------------------------------------------------------------------- //

  // These values define the ranges used to obtain the total f90 distribution, the ER only peak and the NR only peak.
  // Total -> f90_min - f90_max; ER only -> f90_min - f90_mid; NR only -> f90_mid - f90_max.
  Double_t f90_min = 0.2;
  Double_t f90_mid = 0.4;
  Double_t f90_max = 0.6;

  // These variables are used to determine the maximum number of bins, as well as saving the charge boundaries of a given bin.
  int number_of_divisions = (int) (max_charge - min_charge) / bin_size;
  Double_t charge_low = 0;
  Double_t charge_up  = 0;

  for (int i = 0; i < number_of_divisions; i++){

    // Calculate the boundaries of the current bin.
    charge_low = ( i       * bin_size ) + min_charge;
    charge_up  = ( (i + 1) * bin_size ) + min_charge;
    /*
    WriteF90Hist(run, da_both_dir, charge_low, charge_up, f90_min, f90_max, Form("f90_histogram_%d", i+1),      Form("f90 Distribution (Charge Interval: %d - %d PE); f90", (int) charge_low, (int) charge_up) );

    WriteF90Hist(run, da_er_dir,   charge_low, charge_up, f90_min, f90_max, Form("f90_histogram_er_%d", i+1),   Form("f90 Distribution (ER, Charge Interval: %d - %d PE); f90", (int) charge_low, (int) charge_up) );

    WriteF90Hist(run, da_nr_dir,   charge_low, charge_up, f90_min, f90_max, Form("f90_histogram_nr_%d", i+1),   Form("f90 Distribution (NR, Charge Interval: %d - %d PE); f90", (int) charge_low, (int) charge_up));

    WriteF90Hist(run, mc_er_dir,   charge_low, charge_up, f90_min, f90_max, Form("f90_histogram_mcer_%d", i+1), Form("f90 Distribution (MC ER, Charge Interval: %d - %d PE); f90", (int) charge_low, (int) charge_up))

    WriteF90Hist(run, mc_nr_dir,   charge_low, charge_up, f90_min, f90_max, Form("f90_histogram_mcnr_%d", i+1), Form("f90 Distribution (MC NR, Charge Interval: %d - %d PE); f90", (int) charge_low, (int) charge_up))
    */
  }

  hist_file -> Close();
}
