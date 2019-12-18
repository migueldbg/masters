#include <iostream>

#include <TFile.h>
#include <TH1.h>
#include <TSystem.h>

void GetChargeDistribution(int run, bool isMC)
{
    // Cheks wether the root file already exists.
    if (gSystem -> AccessPathName(Form("hist_%d.root", run))){
      std::cout << "The " << Form("hist_%d.root", run) << " file does not exist. Creating it:" << std::endl;
      TFile *hist_file = new TFile (Form("hist_%d.root", run), "CREATE");

    } else {
      std::cout << "Opening the " << Form("hist_%d.root", run) << " file." << std::endl;
      TFile *hist_file = new TFile (Form("hist_%d.root", run), "UPDATE");
    }



}
