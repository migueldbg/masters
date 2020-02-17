/* File: CreateF90Histograms.C (ROOT macro).
 *
 * Author: Miguel Del Ben Galdiano (miguel.galdiano@gmail.com).
 * Date of creation : February 4 2020.
 * Last update      : February 5 2020.
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

/* TH1F* Generatef90Hist( TString file_name, Double_t charge_low, Double_t charge_up, Double_t f90_low, Double_t f90_up )
 *
 * Summary of Generatef90Hist function:
 *
 *    The Generatef90Hist creates a histogram of f90, with the events taken from a TTree in a file decided by the
 *    user. The user also determines the boundaries of the f90 range to be considered, as well as the boundaries of
 *    the charge values to be considered. The function then returns a normalized histogram.
 *
 * Parameters   : file_name  >> the name of the file containing the events from which to construct the histogram.
 *                charge_low >> the lower boundary of the charge value (in PE).
 *                charge_up >> the upper boundary of the charge value (in PE).
 *                f90_low    >> the lower boundary of the f90 value.
 *                f90_up    >> the upper boundary of the f90 value.
 *
 * Return value : TH1F* f90_hist
 */
TH1F* Generatef90Hist( TString file_name, Double_t charge_low, Double_t charge_up, Double_t f90_low, Double_t f90_up ) {

  TFile *file = new TFile(file_name);
  TTree *reco; file -> GetObject("reco", reco);

  TCut cut_f90_min         = Form("clusters[0].f90 >= %f", f90_low);
  TCut cut_f90_max         = Form("clusters[0].f90 <= %f", f90_up);
  TCut cut_charge_min      = Form("clusters[0].charge >= %f", charge_low);
  TCut cut_charge_max      = Form("clusters[0].charge <= %f", charge_up);
  TCut cut_cluster_number  = "number_of_clusters == 1";
  TCut cut_rep             = "clusters[0].rep == 1";

  TCut cut_all = cut_f90_min && cut_f90_max && cut_charge_max && cut_charge_min && cut_cluster_number && cut_rep;

  reco -> Draw("clusters[0].f90 >> hist", cut_all, "goff");
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

// WriteHistogram() function?

void CreateF90Histograms (int run, Double_t bin_size = 20., Double_t max_charge = 1000., Double_t min_charge = 20.){

  TString file_name = Form("hist_%d.root", run);
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

  TH1F** f90hist_data    = new TH1F*[number_of_divisions];
  TH1F** f90hist_data_er = new TH1F*[number_of_divisions];    TH1F** f90hist_data_nr = new TH1F*[number_of_divisions];
  TH1F** f90hist_mc_er   = new TH1F*[number_of_divisions];    TH1F** f90hist_mc_nr   = new TH1F*[number_of_divisions];

  for (int i = 0; i < number_of_divisions; i++){
    // Calculate the boundaries of the current bin.
    charge_low = ( i       * bin_size ) + min_charge;
    charge_up  = ( (i + 1) * bin_size ) + min_charge;

    //--------------------------------------------------------------------------------------------------------------------------------------//

    f90hist_data[i] = Generatef90Hist( Form("runs/run_%d.root", run), charge_low, charge_up, f90_min, f90_max );
    f90hist_data[i] -> SetName ( Form("f90_histogram_%d", i+1) );
    f90hist_data[i] -> SetTitle( Form("f90 Distribution (Charge Interval: %d - %d PE); f90", (int) charge_low, (int) charge_up) );

    da_both_dir -> WriteObject( f90hist_data[i], Form("f90_histogram_%d", i+1), "OverWrite" );

    //--------------------------------------------------------------------------------------------------------------------------------------//

    f90hist_data_er[i] = Generatef90Hist( Form("runs/run_%d.root", run), charge_low, charge_up, f90_min, f90_mid );
    f90hist_data_er[i] -> SetName ( Form("f90_histogram_er_%d", i+1) );
    f90hist_data_er[i] -> SetTitle( Form("f90 Distribution (ER, Charge Interval: %d - %d PE); f90", (int) charge_low, (int) charge_up) );

    da_er_dir -> WriteObject( f90hist_data_er[i], Form("f90_histogram_er_%d", i+1), "OverWrite" );

    //--------------------------------------------------------------------------------------------------------------------------------------//

    f90hist_data_nr[i] = Generatef90Hist( Form("runs/run_%d.root", run), charge_low, charge_up, f90_mid, f90_max );
    f90hist_data_nr[i] -> SetName ( Form("f90_histogram_nr_%d", i+1) );
    f90hist_data_nr[i] -> SetTitle( Form("f90 Distribution (NR, Charge Interval: %d - %d PE); f90", (int) charge_low, (int) charge_up) );

    da_nr_dir -> WriteObject( f90hist_data_nr[i], Form("f90_histogram_nr_%d", i+1), "OverWrite" );

    //--------------------------------------------------------------------------------------------------------------------------------------//

    f90hist_mc_er[i] = Generatef90Hist( Form("runs/run_%d_MCER.root", run), charge_low, charge_up, f90_min, f90_mid );
    f90hist_mc_er[i] -> SetName ( Form("f90_histogram_mcer_%d", i+1) );
    f90hist_mc_er[i] -> SetTitle( Form("f90 Distribution (MC ER, Charge Interval: %d - %d PE); f90", (int) charge_low, (int) charge_up) );

    mc_er_dir -> WriteObject( f90hist_mc_er[i], Form("f90_histogram_mcer_%d", i+1), "OverWrite" );

    //--------------------------------------------------------------------------------------------------------------------------------------//

    f90hist_mc_nr[i] = Generatef90Hist( Form("runs/run_%d_MCNR.root", run), charge_low, charge_up, f90_mid, f90_max );
    f90hist_mc_nr[i] -> SetName ( Form("f90_histogram_mcnr_%d", i+1) );
    f90hist_mc_nr[i] -> SetTitle( Form("f90 Distribution (MC NR, Charge Interval: %d - %d PE); f90", (int) charge_low, (int) charge_up) );

    mc_nr_dir -> WriteObject( f90hist_mc_nr[i], Form("f90_histogram_mcnr_%d", i+1), "OverWrite" );

  }

  hist_file -> Close();
}
