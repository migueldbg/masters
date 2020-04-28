#include "sidutility.cc"

#include <TCut.h>
#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <array>
#include <cstring>

/* ************************************************************************************************************************* *
 * File: CreateF90Histograms.C (ROOT macro).                                                                                 *
 *                                                                                                                           *
 * Author: Miguel Del Ben Galdiano (miguel.galdiano@gmail.com).                                                              *
 * Date of creation : February 4 2020.                                                                                       *
 *                                                                                                                           *
 * Summary of File:                                                                                                          *
 *                                                                                                                           *
 *    The goal of this macro is to generate a number of f90 histograms for a given run. We utilize three different data      *
 *    sets: actual data, an electron recoil (ER) only Monte Carlo (MC) simulation and a neutron recoil (NR) only MC          *
 *    simulation. The user defines the max total charge value of the events to be considered and also determines the number  *
 *    of bins utilized when constructing the histograms. The program then generates a number of histograms equal to the      *
 *    number of bins provided by the size of the bin, each one moving along the total charge: starts at 0, goes up to max    *
 *    total charge in incriments defined by the bin size. The bin size is defined as the total charge range divide by the    *
 *    number of bins. The generated histograms are then saved in a root file, separated in different folders to facilitate   *
 *    later access, if necessary.                                                                                            *
 *                                                                                                                           *
 * ************************************************************************************************************************* */



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
  if      ( !isMC ) { run_file_name = Form("runs/run_%d.root", run); }
  else if ( isMC )  { run_file_name = Form("runs/run_%d_MCER.root", run); }

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
  Double_t s1_charge_low;   Double_t s1_charge_up;
  TString hist_name;     TString hist_title;

  TCut histogram_cuts;

  for (int i = 0; i < number_of_bins; i++){

    // Calculate the boundaries of the current bin.
    s1_charge_low = ( i      * bin_size ) + min_charge;
    s1_charge_up  = ((i + 1) * bin_size ) + min_charge;

    histogram_cuts = DefineCuts(exp_cfg, f90_min, f90_max, s1_charge_low, s1_charge_up, tof_min, tof_max);

    hist_name  = Form( "f90_histogram_%ser_%d", (isMC)?"mc":"", i+1 );
    hist_title = Form( "f90 Distribution (%s, Charge Interval: %d - %d PE); f90", (isMC)?"MC ER":"ER", (int) s1_charge_low, (int) s1_charge_up );

    WriteF90Histogram( dir, run_file_name, histogram_cuts, 100, f90_min, f90_max, false, hist_name, hist_title );
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
  Double_t f90_min;  Double_t f90_max;
  Double_t tof_min;  Double_t tof_max; // If the ToF border values are equal, the code doesn't apply a ToF cut.

  if ( tof_cut ){
    f90_min = 0.0;  f90_max = 1.0;
    tof_min = 30;  tof_max = 42;  // based on runs 1501 to 1521.
  } else if ( !tof_cut ) {
    f90_min = 0.4;  f90_max = 0.7;
    tof_min = 0.0;  tof_max = 0.0;
  }
  // ----------------------------------------- //

  // LOW ENERGY BERILIUM BLOB SELECTION (basend on runs 1501 to 1521) //
  Double_t E_min  = 2250.; Double_t E_max  = 2725.;
  Double_t dE_min = 700.;  Double_t dE_max = 850.;

  TCut E_min_cut  = Form("baseline_mean[31] - ymin[31] > %f", E_min);
  TCut E_max_cut  = Form("baseline_mean[31] - ymin[31] < %f", E_max);
  TCut dE_min_cut = Form("baseline_mean[30] - ymin[30] > %f", dE_min);
  TCut dE_max_cut = Form("baseline_mean[30] - ymin[30] < %f", dE_max);

  TCut lowBe_cut = E_min_cut && E_max_cut && dE_min_cut && dE_max_cut;
  // ---------------------------------------------------------------- //

  Double_t bin_size = ( max_charge - min_charge ) / number_of_bins;
  Double_t s1_charge_low;  Double_t s1_charge_up;

  TString  hist_name;      TString  hist_title;
  TCut histogram_cuts;  TCut total_cuts;

  for (Int_t i = 0; i < number_of_bins; i++){

    // Calculate the boundaries of the current bin.
    s1_charge_low = ( i       * bin_size ) + min_charge;
    s1_charge_up  = ( (i + 1) * bin_size ) + min_charge;

    histogram_cuts = DefineCuts(exp_cfg, f90_min, f90_max, s1_charge_low, s1_charge_up, 0, 50000, tof_min, tof_max);
    total_cuts = histogram_cuts && lowBe_cut;

    hist_name  = Form( "f90_histogram_%snr_%d", (isMC)?"mc":"", i+1 );
    hist_title = Form( "f90 Distribution (%s, Charge Interval: %d - %d PE); f90", (isMC)?"MC NR":"NR", (int) s1_charge_low, (int) s1_charge_up );

    WriteF90Histogram( dir, run_file_name, total_cuts, 100, f90_min, f90_max, false, hist_name, hist_title);
  }

  output_file -> Close();
}



void CreateF90Histograms ( int run, Int_t exp_cfg, bool isMC = false, bool tof_cut = false, Int_t ERorNR = 10, Int_t number_of_bins = 50, Double_t min_charge = 20., Double_t max_charge = 1000. ){

  switch(ERorNR){
    case 0:
        CreateF90ERHistograms( run, exp_cfg, isMC, tof_cut, number_of_bins, min_charge, max_charge );
        break;
    case 1:
        CreateF90NRHistograms( run, exp_cfg, isMC, tof_cut, number_of_bins, min_charge, max_charge );
        break;
    case 10:
        CreateF90ERHistograms( run, exp_cfg, isMC, tof_cut, number_of_bins, min_charge, max_charge );
        CreateF90NRHistograms( run, exp_cfg, isMC, tof_cut, number_of_bins, min_charge, max_charge );
        break;
    default:
        std::cout << "Invalid 'ERorNR' parameter value. Plese enter one of the following options: 0 (ER), 1 (NR) or 10 (ER & NR)." << std::endl;
  }
}
