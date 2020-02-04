/*
 * File: PlotRMS.C (ROOT macro)
 *
 * Author: Miguel Del Ben Galdiano (miguel.galdiano@gmail.com)
 * Date: February 4 2020
 *
 * Summary of File:
 *
 *  The goal of this macro is to compare the root mean square (RMS) of the f90 distribution histograms obtained from actual data and simulations
 *  We shall consider three different data sets, saved in three different .root files: actual data, an electron recoil (ER) only simulation and
 *  a neutron recoil (NR) only simulation. We consider signal events with charge up to 1000 photoelectrons (PE). We divide this range in bins of
 *  20 PE, and for each bin we construct a f90 histogram from each data set. We then obtain the RMS of each distribution and construct a plot of
 *  the RMS as a function of charge. This final plot is the result we wish to obtain.
 */

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>

/* TFile* CheckFile( TString path_name )
 *
 * Summary of CheckFile function:
 *
 *    The CheckFile function checks if a given file already exists and returns this information to the user.
 *    If the file exists, CheckFile returns it. If not, it creates it and then returns it.
 *
 * Parameters   : TString > containing the name of the desired file.
 *
 * Return Value : TFile* file.
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

void PlotRMS (int run, Double_t bin_size = 20., Double_t max_charge = 1000){

  TFile* hist_file;
  TString file_name = Form("hist_%d.root", run);

  hist_file = CheckFile(file_name);
}
