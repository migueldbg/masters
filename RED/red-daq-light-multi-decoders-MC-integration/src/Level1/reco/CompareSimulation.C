#include <iostream>
//#include "../../Modules/RDCluster.hh"
#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <ROOT/RDataFrame.hxx>


TH1F* MyMacro(int count) {
  auto tmp = new TH1F("tmp_hist", "Test Histogram", count, 0, count);

  for(int i = 0; i < count; i++) {
    tmp->Fill(gRandom->Rndm(1)*count);
  }

  return(tmp);
}

void CompareSimulation(int run, int number_divisions = 10, double max_charge_run = 2000){

  //double max_charge_MCER = 20000;
  //double max_charge_MCNR = 9000;
  TH1F* hist = MyMacro(10000);
  hist -> Draw();


  /*
  TFile *run_file = new TFile(Form("run_%d.root", run));
  TFile *MCER_file = new TFile(Form("run_%dMCER.root", run));
  TFile *MCNR_file = new TFile(Form("run_%dMCNR.root", run));

  TTree *reco_run;    run_file  -> GetObject("reco", reco_run);
  TTree *reco_MCER;   MCER_file -> GetObject("reco", reco_MCER);
  TTree *reco_MCNR;   MCNR_file -> GetObject("reco", reco_MCNR);

  double division_size_run  = max_charge_run/number_divisions;
  double division_size_MCER = max_charge_MCER/number_divisions;
  double division_size_MCNR = max_charge_MCNR/number_divisions;
  */
}
