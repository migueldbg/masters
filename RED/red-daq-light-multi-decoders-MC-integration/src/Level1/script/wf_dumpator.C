#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "red-daq/EvRec0.hh"
#include "red-daq/EvRaw0.hh"
#include <TVirtualFFT.h>
#include <TGraph.h>
#include <TMath.h>
#include <TH1D.h>
#include <TF1.h>
#include <TString.h>
#include <TH1.h>

using namespace std;

void fourierator(int tot_ev = 100) {
  TFile* f = new TFile("/storage/local/home/darkside/nrossi/Level1/run_359.root");

  TTree* reco = (TTree*) f->Get("reco");
  TTree* raw = (TTree*) f->Get("raw");

  EvRec0* evReco = new EvRec0();
  EvRaw0* evRaw = new EvRaw0();

  reco->SetBranchAddress("recoevent",&evReco);
  raw->SetBranchAddress("rawevent",&evRaw);
  
  TFile *fout = new TFile("fourier.root", "recreate");

  vector<double> v_top, v_bot;
  TTree* t = new TTree();
  t->SetName("wf");
  t->Branch("v_top", &v_top);
  t->Branch("v_bot", &v_bot);

  int ev = 0;
  for (int i = 0; i < reco->GetEntries(); i++) { 
     raw->GetEntry(i);
     reco->GetEntry(i);
         
     if (evReco->GetChargeTot() < 20) {
        vector<double>* wf_top = evRaw->GetWF(24);
        vector<double>* wf_bot = evRaw->GetWF(14);

        for (size_t j = 0; j< wf_top->size(); j++) {
           v_top.push_back(wf_top->at(j)); 
           v_bot.push_back(wf_bot->at(j)); 
        }   
        cout << "event " << ev << " found at " << i  << endl; 
        t->Fill();
        ++ev;
     }
     if (ev == tot_ev) break;
   }     
 
   fout->cd();
   t->Write();
   fout->Close();

}
