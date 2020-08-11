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


TDirectory* MakeDirStruct(TFile* file, bool isMC, Int_t ERorNR);



/* void CreateF90ERHistograms ( int run, Int_t expCfg, bool isMC = false, bool cutTof = false, Double_t binSize = 20, Double_t s1Min = 0., Double_t s1Max = 1000.)
 *
 *  Summary of Function:
 *
 *    This function constructs a number of electron recoil f90 histograms within a total S1 charge range defined by the user. The
 *    determination of wether an event is classified as an ER can be done by f90 cuts (cutTof = false) or ToF cuts (cutTof = true).
 *    Changing both of these values must be done within the code. The function can also take into consideration wether the run is
 *    single phase or dual phase and call for the appropriate quality cuts. The resulting histograms are then saved in a root file
 *    named analysis_#run number#.root.
 *
 *  Parameters   : run       >> run number.
 *                 data_type >> in dicates if the source is taken from data or MC simulation.
 *                 expCfg    >> indicates if the run was single phase or dual phase.
 *                 cutTof    >> indicates wheter ToF cuts are to be applied for ER event selection.
 *                 binSize   >> size of the S1 slices for which a f90 histogram will be generated.
 *                 s1Min     >> the maximum value of S1 total charge to be considered.
 *                 s1Max     >> the minimum value of S1 total charge to be considered.
 *
 *  Return Value : void.
 */
void CreateF90ERHistograms ( int run, Int_t expCfg, bool isMC = false, bool cutTof = false, Double_t binSize = 20, Double_t s1Min = 0.,
                             Double_t s1Max = 1000.){

  TString analysisFileName = analysisDirectoryPath + Form("/analysis_%d.root", run);
  TFile*  analysisFile = CheckFile(analysisFileName);

  TDirectory* dir = MakeDirStruct(analysisFile, isMC, 0); // 0 creates the 'ER' folder.

  TString runFileName;
  if ( !isMC ) { runFileName = runsDirectoryPath + Form("/run_%d.root", run); }
  else if ( isMC ) { runFileName = runsDirectoryPath + Form("/run_%d_MCER.root", run); }

  TFile* runFile = CheckFile(runFileName);
  TTree* reco;  runFile -> GetObject("reco", reco);


  // ELECTRON RECOIL EVENT SELECTION PARAMETERS //
  Double_t f90Min;  Double_t f90Max;
  Double_t tofMin;  Double_t tofMax; // If the ToF border values are equal, the code doesn't apply a ToF cut.

  if ( cutTof ){
    f90Min = 0.0;  f90Max = 1.0;
    tofMin = 30;   tofMax = 42;  // based on runs 1501 to 1521.
  } else if ( !cutTof ) {
    f90Min = 0.1;  f90Max = 0.4; //TODO: change f90 limits.
    tofMin = 0.0;  tofMax = 0.0;
  }


  Int_t binNumber = ( s1Max - s1Min ) / binSize;
  Double_t s1LowBound;
  Double_t s1UppBound;

  TString histName;     TString histTitle;

  TCut tpcCut;

  for (int i = 0; i < binNumber; i++){

    // Calculate the boundaries of the current bin.
    s1LowBound = ( i      * binSize ) + s1Min;
    s1UppBound = ((i + 1) * binSize ) + s1Min;

    tpcCut = DefineCuts(expCfg, f90Min, f90Max, s1LowBound, s1UppBound, 0, 0, tofMin, tofMax);

    histName  = Form( "f90_histogram_%ser_%d", (isMC)?"mc":"", i+1 );
    histTitle = Form( "f90 Distribution (%s, Charge Interval: %d - %d PE); f90", (isMC)?"MC ER":"ER", (int) s1LowBound, (int) s1UppBound );

    TH1F* f90Hist = new TH1F(histName, histTitle, 100, f90Min, f90Max);
    reco -> Project(histName, "clusters[0].f90", tpcCut);

    dir -> WriteObject(f90Hist, histName, "OverWrite");
  }

  runFile -> Close();
  analysisFile -> Close();
}

/* void CreateF90NRHistograms ( int run, Int_t expCfg, bool isMC = false, bool cutTof = false, Double_t binSize = 20, Double_t s1Min = 0., Double_t s1Max = 1000.)
 *
 *  Summary of Function:
 *
 *    This function constructs a number of nuclear recoil f90 histograms within a total S1 charge range defined by the user. The
 *    determination of whether an event is classified as an NR can be done by f90 cuts (cutTof = false) or ToF cuts (cutTof = true).
 *    Changing both of these values must be done within the code. The function can also take into consideration wether the run is
 *    single phase or dual phase and call for the appropriate quality cuts. The resulting histograms are then saved in a root file
 *    named analysis_#run number#.root.
 *
 *  Parameters   : run       >> run number.
 *                 data_type >> in dicates if the source is taken from data or MC simulation.
 *                 expCfg    >> indicates if the run was single phase or dual phase.
 *                 cutTof    >> indicates wheter ToF cuts are to be applied for ER event selection.
 *                 binSize   >> size of the S1 slices for which a f90 histogram will be generated.
 *                 s1Min     >> the maximum value of S1 total charge to be considered.
 *                 s1Max     >> the minimum value of S1 total charge to be considered.
 *
 *  Return Value : void.
 */
void CreateF90NRHistograms ( int run, Int_t expCfg, bool isMC = false, bool cutTof = false, Double_t binSize = 20, Double_t s1Min = 0.,
                             Double_t s1Max = 1000.){

  TString analysisFileName = analysisDirectoryPath + Form("/analysis_%d.root", run);
  TFile*  analysisFile = CheckFile(analysisFileName);

  TDirectory* dir = MakeDirStruct(analysisFile, isMC, 1);

  TString runFileName;
  if ( !isMC ) { runFileName = runsDirectoryPath + Form("/run_%d.root", run); }
  else if ( isMC ) { runFileName = runsDirectoryPath + Form("/run_%d_MCNR.root", run); }

  TFile* runFile = CheckFile(runFileName);
  TTree* reco;  runFile -> GetObject("reco", reco);


  // NUCLEAR RECOIL EVENT SELECTION PARAMETERS //
  Double_t f90Min;  Double_t f90Max;
  Double_t tofMin;  Double_t tofMax; // If the ToF border values are equal, the code doesn't apply a ToF cut.

  if ( cutTof ){
    f90Min = 0.0;  f90Max = 1.0;
    tofMin = 30;   tofMax = 42;
  } else if ( !cutTof ) {
    f90Min = 0.4;  f90Max = 0.7;
    tofMin = 0.0;  tofMax = 0.0;
  }


  Int_t binNumber = ( s1Max - s1Min ) / binSize;
  Double_t s1LowBound;
  Double_t s1UppBound;

  TString  histName;      TString  histTitle;
  TCut tpcCut;

  for (Int_t i = 0; i < binNumber; i++){

    // Calculate the boundaries of the current bin.
    s1LowBound = ( i       * binSize ) + s1Min;
    s1UppBound = ( (i + 1) * binSize ) + s1Min;

    tpcCut = DefineCuts(expCfg, f90Min, f90Max, s1LowBound, s1UppBound, 0, 0, tofMin, tofMax);

    histName  = Form( "f90_histogram_%snr_%d", (isMC)?"mc":"", i+1 );
    histTitle = Form( "f90 Distribution (%s, Charge Interval: %d - %d PE); f90", (isMC)?"MC NR":"NR", (int) s1LowBound, (int) s1UppBound );

    TH1F* f90Hist = new TH1F(histName, histTitle, 100, f90Min, f90Max);
    reco -> Project(histName, "clusters[0].f90", tpcCut);

    dir -> WriteObject(f90Hist, histName, "OverWrite");
  }

  analysisFile -> Close();
}


void CreateF90Histograms ( int run, Int_t expCfg, bool isMC = false, bool cutTof = false, Int_t ERorNR = 10, Double_t binSize = 20, Double_t s1Min = 20., Double_t s1Max = 1000. ){

  switch(ERorNR){
    case 0:
        CreateF90ERHistograms(run, expCfg, isMC, cutTof, binSize, s1Min, s1Max);
        break;
    case 1:
        CreateF90NRHistograms(run, expCfg, isMC, cutTof, binSize, s1Min, s1Max);
        break;
    case 10:
        CreateF90ERHistograms(run, expCfg, isMC, cutTof, binSize, s1Min, s1Max);
        CreateF90NRHistograms(run, expCfg, isMC, cutTof, binSize, s1Min, s1Max);
        break;
    default:
        std::cout << "Invalid 'ERorNR' parameter value. Plese enter one of the following options: 0 (ER), 1 (NR) or 10 (ER & NR)." << std::endl;
  }
}


// =================================== AUXILIARY FUNCTIONS USED IN THE MACRO =================================== //

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
TDirectory* MakeDirStruct(TFile* file, bool isMC, Int_t ERorNR){

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
