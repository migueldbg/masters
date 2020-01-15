#include <iostream>

#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <ROOT/TDataFrame.hxx>

using namespace ROOT::Experimental;

void CompareSimulation(int run, int number_divisions){

  TFile *run_file = new TFile(Form("run_%d.root", run));
  TFile *MCER_file = new TFile(Form("run_%dMCER.root", run));
  TFile *MCNR_file = new TFile(Form("run_%dMCNR.root", run));

  /*TTree *tree_reco_run;
  TTree *tree_reco_MCER;
  TTree *tree_reco_MCNR;
  */

  TDataFrame d_run("reco", run_file);
  TDataFrame d_MCER("reco", MCER_file);
  TDataFrame d_MCNR("reco", MCNR_file);

}
