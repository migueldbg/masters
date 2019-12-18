#include <iostream>

#include <TFile.h>
#include <TH1.h>
#include <TSystem.h>

TTree *tree_reco;

void GetChargeDistribution(int run, bool isMC = true)
{
    // Cheks wether the root file already exists.
    /*if (gSystem -> AccessPathName(Form("hist_%d.root", run))){
      std::cout << "The " << Form("hist_%d.root", run) << " file does not exist. Creating it:" << std::endl;
      TFile *hist_file = new TFile (Form("hist_%d.root", run), "CREATE");

    } else {
      std::cout << "Opening the " << Form("hist_%d.root", run) << " file." << std::endl;
      TFile *hist_file = new TFile (Form("hist_%d.root", run), "UPDATE");
    }*/

    TFile *run_file = new TFile(Form("run_%d%s.root",run, (isMC)?"MC":""));
    run_file -> GetObject("reco", tree_reco);
    tree_reco ->SetMarkerStyle(20);
    tree_reco ->SetMarkerSize(0.5);

    tree_reco -> Draw("clusters[0].charge >> charge_distribution", "clusters[0].f90 > 0 && clusters[0].f90 < 1");

}
