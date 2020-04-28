#include "sidutility.cc"

#include <TCut.h>
#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

/* ************************************************************************************************************************* *
 * File: ToF.C                                                                                                               *
 *                                                                                                                           *
 * Author: Miguel Del Ben Galdiano                                                                                           *
 * Date of Creation: March 3 2020.                                                                                           *
 *                                                                                                                           *
 * Summary of File:                                                                                                          *
 *                                                                                                                           *
 *    The goal of this macro is to create a 2D histogram with y-axis equal to f90 and x-axis equal to a parameter to be      *
 *    defined in the code, the so called "time of flight" (TOF). This parameter will be used to apply a cut, such as to      *
 *    extract nuclear recoils from background events.                                                                        *
 * ************************************************************************************************************************* */

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

  hist -> SetTitle(Form("Time of Flight Distribution; ToF (ns)"));
  hist -> Draw();

}
