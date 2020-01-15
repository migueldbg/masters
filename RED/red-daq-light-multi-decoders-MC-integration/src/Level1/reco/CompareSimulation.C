#include <iostream>
//#include "../../Modules/RDCluster.hh"
#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>

void CompareSimulation(int run, int number_divisions = 10, double max_charge = 2000){

  TFile *run_file = new TFile(Form("run_%d.root", run));
  //TFile *MCER_file = new TFile(Form("run_%dMCER.root", run));
  //TFile *MCNR_file = new TFile(Form("run_%dMCNR.root", run));

  TTree *tree_reco_run;
  TTree *tree_reco_MCER;
  TTree *tree_reco_MCNR;

  run_file  -> GetObject("reco", reco_run);
  MCER_file -> GetObject("reco", reco_MCER);
  MCNR_file -> GetObject("reco", reco_MCNR);

  

}
