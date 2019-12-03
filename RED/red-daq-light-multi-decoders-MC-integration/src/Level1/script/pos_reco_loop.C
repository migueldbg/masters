#include <iostream>
#include <cstring>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TF1.h>

#include "red-daq/EvRec0.hh"
#include "red-daq/RDCluster.hh"
#include "red-daq/RDconfig.h"

using namespace std;

void pos_reco_loop(TString filename)
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
 
    TH2F* bar = new TH2F("bar", "Barycenter; x[cm]; y[cm]", 100, 1, 4, 100, 1, 4); 

    TH1F* charge1 = new TH1F("charge1", "Corner; charge_top [PE]; Counts", 400, 0, 4000);
    TH1F* charge2 = new TH1F("charge2", "Center; charge_top [PE]; Counts", 400, 0, 4000);
    TH2F* charge_x = new TH2F("charge_x", "Charge_radius; charge_top [PE]; x [cm]", 200, 0, 1000, 100, 0, 5);

    TF1* r = new TF1("f","-x+4.2",0, 5);

    for (int ev = 0; ev < data->GetEntries(); ev++) {                
       data->GetEntry(ev);
       vector<RDCluster*> clusters = evReco->GetClusters();

       size_t nc = clusters.size();	
        
       if (nc==2 && clusters.at(0)->rep == 1) {
     
          vector<double> ch1 = clusters.at(1)->charge_top;
          double ch_top = clusters.at(1)->tot_charge_top;
              
          double pos_x(0), pos_y(0);

          for (int iy = 0; iy < 6; iy++) {
             for (int ix = 0; ix < 4; ix++) {  
                double ch_norm = ch1.at(A[iy][ix])/ch_top;
                pos_x += ((0.5 + ix)*5/4)*ch_norm; // in cm
                pos_y += ((0.5 + iy)*5/6)*ch_norm;
             }
          }

          bar->Fill(pos_x, pos_y);
          charge_x->Fill(clusters.at(0)->tot_charge_top, pos_x);
        
          //double dr1 = sqrt(pow(pos_x - 1.9, 2) + pow(pos_x - 1.9, 2));
          //double dr2 = sqrt(pow(pos_x - 2.5, 2) + pow(pos_x - 2.4, 2));
          
          if(pos_y < r->Eval(pos_x)) charge1->Fill(clusters.at(0)->tot_charge_top);
          else charge2->Fill(clusters.at(0)->tot_charge_top);

          //if(dr1 < 0.2) charge1->Fill(clusters.at(0)->tot_charge_top);
          //if(dr2 < 0.25) charge2->Fill(clusters.at(0)->tot_charge_top);

  
       } 

       if(ev%1000==0) cout << "Event " << ev << " processed" << endl; 
    }
    

/////////////////////////////////////////////////////

    TH2F* bar1 = new TH2F("bar1", "Barycenter Left Peak; x[cm]; y[cm]", 100, 1, 4, 100, 1, 4); 




    for (int ev = 0; ev < data->GetEntries(); ev++) {                
       data->GetEntry(ev);
       vector<RDCluster*> clusters = evReco->GetClusters();

       size_t nc = clusters.size();	
        
       if (nc==2 && clusters.at(0)->rep == 1 && clusters.at(0)->tot_charge_top > 220 && clusters.at(0)->tot_charge_top < 260) {
     
          vector<double> ch1 = clusters.at(1)->charge_top;
          double ch_top = clusters.at(1)->tot_charge_top;
              
          double pos_x(0), pos_y(0);

          for (int iy = 0; iy < 6; iy++) {
             for (int ix = 0; ix < 4; ix++) {  
                double ch_norm = ch1.at(A[iy][ix])/ch_top;
                pos_x += ((0.5 + ix)*5/4)*ch_norm; // in cm
                pos_y += ((0.5 + iy)*5/6)*ch_norm;
             }
          }

          bar1->Fill(pos_x, pos_y);

       } 

       if(ev%1000==0) cout << "Event " << ev << " processed" << endl; 
    }
    


/////////////////////////////////////////////////////
    TH2F* bar2 = new TH2F("bar2", "Barycenter Right Peak; x[cm]; y[cm]", 100, 1, 4, 100, 1, 4); 


    for (int ev = 0; ev < data->GetEntries(); ev++) {                
       data->GetEntry(ev);
       vector<RDCluster*> clusters = evReco->GetClusters();

       size_t nc = clusters.size();	
        
       if (nc==2 && clusters.at(0)->rep == 1 && clusters.at(0)->tot_charge_top > 320 && clusters.at(0)->tot_charge_top < 380) {
     
          vector<double> ch1 = clusters.at(1)->charge_top;
          double ch_top = clusters.at(1)->tot_charge_top;
              
          double pos_x(0), pos_y(0);

          for (int iy = 0; iy < 6; iy++) {
             for (int ix = 0; ix < 4; ix++) {  
                double ch_norm = ch1.at(A[iy][ix])/ch_top;
                pos_x += ((0.5 + ix)*5/4)*ch_norm; // in cm
                pos_y += ((0.5 + iy)*5/6)*ch_norm;
             }
          }

          bar2->Fill(pos_x, pos_y);
        
       } 

       if(ev%1000==0) cout << "Event " << ev << " processed" << endl; 
    }
    
    TCanvas* c1 = new TCanvas("c1", "c1", 1000, 1000);
    bar->Draw("colz");

    r->SetLineColor(kRed);
    r->SetLineWidth(3);
    r->Draw("same");


    TCanvas* c2 = new TCanvas("c2", "c2", 1000, 1000);
    
    bar1->SetMarkerStyle(7); bar1->SetMarkerColor(kRed);
    bar2->SetMarkerStyle(7); bar2->SetMarkerColor(kBlue);

    bar1->Draw();
    bar2->Draw("same");

    TCanvas* c3 = new TCanvas("c3", "c3", 1000, 1000);
    charge1->SetLineColor(kRed);
    charge2->SetLineColor(kBlue);
    charge2->Draw();
    charge1->Draw("same");

    TCanvas* c4 = new TCanvas("c4", "c4", 1000, 1000);
    charge_x->Draw("colz");
} 



