#include <iostream>

#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>


void GetChargeDistribution(int run, bool isMC = false, int ERorNR = 0)
{
  // This line is simply to make it so that the canvas generated isn't show. The last line reverts this operation.
  gROOT->SetBatch(kTRUE);

    TFile *run_file = new TFile(Form("run_%d%s%s%s.root",run, (isMC)?"MC":"", (ERorNR == 1)?"ER":"", (ERorNR == 2)?"NR":""));
    TTree *tree_reco;

    run_file  -> GetObject("reco", tree_reco);
    tree_reco -> SetMarkerStyle(20);
    tree_reco -> SetMarkerSize(0.5);
    tree_reco -> Draw("clusters[0].charge >> charge_dist", "clusters[0].f90 > 0 && clusters[0].f90 < 1");

    TH1F *charge_dist = (TH1F*)gDirectory -> Get("charge_dist");
    charge_dist -> SetName(Form("charge_distribution%s%s%s", (isMC)?"MC":"", (ERorNR == 1)?"ER":"", (ERorNR == 2)?"NR":""));
    charge_dist -> SetTitle("Charge Distribution; Charge (in PE)");

  // Remove the histogram from the current directory so that there isn't an extra copy in it.
  charge_dist -> SetDirectory(0);

    // Cheks wether the root file already exists.
    if (gSystem -> AccessPathName(Form("hist_%d.root", run))){
      std::cout << "The " << Form("hist_%d.root", run) << " file does not exist. Creating it..." << std::endl;
      TFile *hist_file = new TFile(Form("hist_%d.root", run), "CREATE");

    } else {
        std::cout << "Opening the " << Form("hist_%d.root", run) << " file." << std::endl;
        TFile *hist_file = new TFile(Form("hist_%d.root", run), "UPDATE");

        if ( !(hist_file -> IsOpen()) ) { std::cout << "Unable to open file." << std::endl;}
    }

    charge_dist -> Write(Form("charge_distribution%s%s%s", (isMC)?"MC":"", (ERorNR == 1)?"ER":"", (ERorNR == 2)?"NR":""));*/
  delete charge_dist;
  gROOT->SetBatch(kFALSE);
}
