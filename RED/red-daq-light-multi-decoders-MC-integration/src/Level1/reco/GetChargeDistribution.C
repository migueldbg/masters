#include <iostream>

#include <TFile.h>
#include <TH1.h>
#include <TSystem.h>

TTree *tree_reco;

void GetChargeDistribution(int run, bool isMC = false, int ERorNR = 0)
{
    // Cheks wether the root file already exists.
    if (gSystem -> AccessPathName(Form("hist_%d.root", run))){
      std::cout << "The " << Form("hist_%d.root", run) << " file does not exist. Creating it..." << std::endl;
      TFile *hist_file = new TFile (Form("hist_%d.root", run), "CREATE");

    } else {
      std::cout << "Opening the " << Form("hist_%d.root", run) << " file." << std::endl;
      TFile *hist_file = new TFile (Form("hist_%d.root", run), "UPDATE");
    }

    // Might need to change the placemente of the code abover. Or just find a way to change the currenct directory, so that the histogram is saved in the
    // correct file ("hist_%d.root").

    TFile *run_file = new TFile(Form("run_%d%s.root",run, (isMC)?"MC":""));
    run_file -> GetObject("reco", tree_reco);
    tree_reco ->SetMarkerStyle(20);
    tree_reco ->SetMarkerSize(0.5);

    tree_reco -> Draw("clusters[0].charge >> charge_distribution", "clusters[0].f90 > 0 && clusters[0].f90 < 1");
    TH1F *charge_dist = (TH1F*)gDirectory -> Get("charge_distribution");
    charge_dist -> SetName(Form("charge_distribution%s%s%s", (isMC)?"MC":"", (ERorNR == 1)?"ER":"", (ERorNR == 2)?"NR":""))

    // Might add some way to rename the created histogram depending on the original root file. Maybe adding the possibility of adding MC, ER and NR to the end.
    // Maybe just add the possibility of the user himself changing the name? It opens the possibility of less cosistency though.
    // Now must save the created histogram
}
