#include <TCut.h>
#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

/* ************************************************************************************************************************
 * File: ToF.C                                                                                                            *
 *                                                                                                                        *
 * Author: Miguel Del Ben Galdiano                                                                                        *
 * Date of Creation: March 3 2020.                                                                                        *
 *                                                                                                                        *
 * Summary of File:                                                                                                       *
 *                                                                                                                        *
 *    The goal of this macro is to create a 2D histogram with y-axis equal to f90 and x-axis equal to a parameter to be   *
 *    defined in the code, the so called "time of flight" (TOF). This parameter will be used to apply a cut, such as to   *
 *    extract nuclear recoils from background events.                                                                     *
 * ************************************************************************************************************************ /

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


TH1* GenerateTOFHist( TString file_name ){

  TFile* file = new TFile(file_name);
  TTree* reco; file -> GetObject("reco", reco);

  TCut f90_cut     = DefineF90Range(0.2, 1.);
  TCut quality_cut = DefineQualityCuts(2);
  TCut final_cut = f90_cut + quality_cut;

  std::string ToF = "2*(xmin[30] - clusters[0].min_x)";
  std::string histogram = " >> hist(100, -100, 100)";
  std::string expression = ToF + histogram;

  reco -> Draw(expression.c_str(), final_cut, "goff");
  TH1* ToF_hist = (TH1*) gDirectory -> Get("hist");

  //ToF_hist -> SetOption("colz");
  ToF_hist -> SetDirectory(0);

  file -> Close();

  return ToF_hist;
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

  TH1F* hist = (TH1F*) GenerateTOFHist( Form("runs/run_%d.root", run) );

  hist -> Draw();

}
