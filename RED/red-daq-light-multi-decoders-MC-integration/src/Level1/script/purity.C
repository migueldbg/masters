#include <iostream>
#include <cstring>

#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TProfile.h>
#include <TCut.h>
#include <TTree.h>
using namespace std;

void purity(string const& fname) {
   
   string path = "/storage/DATA-02/darkside/red/reco/rm3reco/";
   path += fname;

   TFile* f = new TFile(path.c_str(), "read");
   TTree* reco = (TTree* )f->Get("reco");

   // data selection
   TCut raw = "number_of_clusters==2";
   TCut quality = "";
    
    
   TH2F* h2 = new TH2F("h2","Purity; drfit time [us]; S2/S1", 100, 0, 40000/500, 50, -15, 15);

   reco->Project("h2","clusters[1].charge/clusters[0].charge:(clusters[1].start_time-clusters[0].start_time)/500", raw);

   TCanvas* c0 = new TCanvas("c0", "c0", 1);
   h2->Draw();

   TProfile *p = h2->ProfileX();
   p->Draw("same");

}
