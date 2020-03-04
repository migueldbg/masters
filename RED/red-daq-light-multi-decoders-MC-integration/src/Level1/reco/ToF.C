/* File: ToF.C
 *
 * Author: Miguel Del Ben Galdiano]
 * Date of Creation: March 3 2020.
 *
 * Summary of File:
 *
 *    The goal of this macro is to create a 2D histogram with y-axis equal to f90 and x-axis equal to a parameter to be
 *    defined in the code, the so called "time of flight" (TOF). This parameter will be used to apply a cut, such as to
 *    extract nuclear recoils from background events.
 */

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

TH2* GenerateTOFHist( TString file_name ){

  TFile* file = new TFile(file_name);
  TTree* reco; file -> GetObject("reco", reco);

  TCut cut_f90_min = "clusters[0].f90 >= 0.0";
  TCut cut_f90_max = "clusters[0].f90 <= 1.0";
  TCut cut_rep     = "clusters[0].rep == 1";

  TCut all_cuts = cut_f90_min + cut_f90_max + cut_rep;

  const char* expression = "clusters[0].f90:clusters[0].start_time - start_time[30] >> hist(120, -2000, 40000, 120, 0, 1)";

  reco -> Draw(expression, all_cuts, "goff");
  TH2* TOFvf90_hist = (TH2*) gDirectory -> Get("hist");

  TOFvf90_hist -> SetOption("colz");
  TOFvf90_hist -> SetDirectory(0);

  file -> Close();

  return TOFvf90_hist;
}

// ---------------------------------------------------- MACRO::ToF ---------------------------------------------------- //

void ToF ( int run ){

  TString file_name = Form("analysis_%d.root", run);
  TFile* analysis_file = CheckFile(file_name);

  // ------------------------ CREATING NECESSARY DIRECTORIES ------------------------ //
  analysis_file -> cd();
  TDirectory* histograms_dir = MakeDirectory("histograms", "histograms");

  histograms_dir -> cd();
  TDirectory* tof_dir = MakeDirectory("time_of_flight", "time_of_flight");
  // -------------------------------------------------------------------------------- //

  TH2F* hist = (TH2F*) GenerateTOFHist( Form("runs/run_%d.root", run) );

  hist -> Draw();

}
