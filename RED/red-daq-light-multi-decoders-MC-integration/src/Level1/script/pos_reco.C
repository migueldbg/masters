#include <iostream>
#include <cstring>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>

#include "red-daq/EvRec0.hh"
#include "red-daq/RDCluster.hh"
#include "red-daq/RDconfig.h"

using namespace std;

void pos_reco(TString filename, int ev = 0)
{   

    gStyle->SetOptStat(0);
    gStyle->SetPalette(55);


    double A[6][4] = {{0, 1, 2, 3},
                     {5, 6, 4, 8},
                     {10, 11, 7, 9},
                     {2+12, 3+12, 1+12, 0+12},
                     {7+12, 8+12, 6+12, 5+12},
                     {9+12, 10+12, 4+12, 11+12}};

    TFile *f = new TFile(filename, "read");
    if (!(f->IsOpen()))
    {
      cout << "could not open file: " << filename << endl;
      return;
    }

    TTree *data = (TTree*)f->Get("reco");
    EvRec0* evReco = new EvRec0();
    data->SetBranchAddress("recoevent",&evReco);
    
    data->GetEntry(ev);
   
    TCanvas* c1 = new TCanvas("c1", "c1", 1000,1000); 
    c1->Divide(2,2);

    vector<TH2F*> h2; 
    vector<TGraph*> g; 
    vector<TGraph*> gm; 
    vector<TGraph*> gs;
        
    size_t nchannels = evReco->GetCharge().size();
    data->GetEntry(ev);
    vector<RDCluster*> clusters = evReco->GetClusters();

    size_t nc = clusters.size();	
    cout << "Found n. " << nc << " clusters" << endl;

                
    for (size_t ch = 0; ch < nc && ch < 4; ch++) {

       c1->cd(ch + 1);
       h2.push_back(new TH2F(Form("h2_%d", (int) ch), Form("Cluster[%d] X-Y", (int) ch), 4, 0, 5, 6, 0, 5));
       //h2[ch]->SetContour(30);
       g.push_back(new TGraph());
       g[ch]->SetMarkerStyle(29);
       g[ch]->SetMarkerSize(10);
       g[ch]->SetMarkerColor(kMagenta);

       gm.push_back(new TGraph());
       gm[ch]->SetMarkerStyle(28);
       gm[ch]->SetMarkerSize(10);
       gm[ch]->SetMarkerColor(kMagenta);

       gs.push_back(new TGraph());
       gs[ch]->SetMarkerStyle(34);
       gs[ch]->SetMarkerSize(10);
       gs[ch]->SetMarkerColor(kMagenta);
            
       vector<double> ch1 = clusters.at(ch)->charge_top;
       double ch_top = clusters.at(ch)->tot_charge_top;
              
       double pos_x(0), pos_y(0);
       double pos_x2(0), pos_y2(0);

       for (int iy = 0; iy < 6; iy++) {
          for (int ix = 0; ix < 4; ix++) {  
             double ch_norm = ch1.at(A[iy][ix])/ch_top;
             h2[ch]->SetBinContent(ix + 1, iy + 1, ch_norm);
             pos_x += ((0.5 + ix)*5/4)*ch_norm; // in cm
             pos_y += ((0.5 + iy)*5/6)*ch_norm;

          }
       }

       gm[ch]->SetPoint(0, h2[ch]->GetMean(1), h2[ch]->GetMean(2));             
       g[ch]->SetPoint(0, pos_x, pos_y); 

       
       int x, y, z;     
       h2[ch]->Smooth(); 
       h2[ch]->GetMaximumBin(x, y, z);
       gs[ch]->SetPoint(0, h2[ch]->GetXaxis()->GetBinCenter(x), h2[ch]->GetYaxis()->GetBinCenter(y));     

       h2[ch]->Draw("colz");
       g[ch]->Draw("Psame");
       gm[ch]->Draw("Psame");
       gs[ch]->Draw("Psame");

    }
    c1->Modified();
    c1->Update(); 
} 



