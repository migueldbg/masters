#include <iostream>

#include "TCut.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"

void vis_1520(){


  TFile *run_file = new TFile("runs/run_1520.root");

  TTree *reco; reco = (TTree*)run_file -> Get("reco");

  gStyle->SetOptStat(2210);
  gStyle->SetOptFit();
  //TCut *tof_cut = "xmin[] - clusters.min_x[0]";

  reco -> Draw("clusters.f90[0]:clusters[0].charge >> h", "clusters.f90[0] > 0 && clusters.f90[0] < 1. && xmin[30] - clusters.min_x[0] > 0. && xmin[30] - clusters.min_x[0] < 50. && clusters.charge[0] < 2000");
  
  run_file->Close();

}
