#include "TROOT.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include <TGraph.h>
#include <iostream>
#include <algorithm>
#include "red-daq/EvRec0.hh"
#include "red-daq/EvRaw0.hh"

#include <numeric>
#include <vector>

using namespace std;

void filter(TString filename, int run)
{
  TFile* f = new TFile(filename);
  if (!(f->IsOpen()))
    {
      cout << "could not open file: " << filename << endl;
      return ;
    }
  TTree* reco = (TTree*) f->Get("reco");
  TTree* metaevent = (TTree*) f->Get("metaevent");
  TTree* raw = (TTree*) f->Get("raw");
  if (!raw)
    {
      cout << "TTree of rawdata not found! " << endl;
      cout << "Please re-process with RedLevel1 -w " << endl;
      return ;
    }
  EvRec0* evReco = new EvRec0();
  EvRaw0* evRaw = new EvRaw0();
  int chanmap[100];
  reco->SetBranchAddress("recoevent",&evReco);
  raw->SetBranchAddress("rawevent",&evRaw);
  metaevent->SetBranchAddress("chanmap",&chanmap);

  raw->GetEntry(0);
  metaevent->GetEntry(0);
  int imax = evRaw->GetWFs().rbegin()->first;
  cout << "Max channel is: " << imax << endl;
  vector<TH1D*> avg(imax+1,nullptr);
  vector<TH1D*> ker(imax+1,nullptr);
  //vector<double>* wf0 = evRaw->GetWF(16); //bottom channel
  //vector<double>* wf0 = evRaw->GetWF(21); //top channel

  vector<int> nevents(imax+1,0);
  //int index=0;
  // ***** AVERAGE WAVEFORM *****

  for (int ev = 0; ev < reco->GetEntries(); ev++) {
     reco->GetEntry(ev);
     raw->GetEntry(ev);
     //index = 0;
     map<int,vector<double>* > wfs = evRaw->GetWFs();
     for (map<int,vector<double>* >::iterator it=wfs.begin(); it!=wfs.end(); ++it) {
       int i = it->first; 
          
       double aCharge = evReco->GetCharge().at(i);
       double aTime = evReco->GetStartTime().at(i);
       double aAmpli = evReco->GetBaseMean().at(i)-evReco->GetYmin().at(i);
       double  aBase = evReco->GetBaseMean().at(i);
       vector<double>* wf = evRaw->GetWF(i);
       //cout << i << "Vector SIZE: " << wf->size() << endl;

       if (!avg.at(i))
       {
           TString name;
           avg.at(i) = new TH1D(Form("h%d",i),Form("Channel_%d",i),
                                wf->size(),0,wf->size()-1);
           ker.at(i) = new TH1D(Form("k%d",i),Form("Kernel_ch%d",i),
                                1300,0,1299); 
      }
      
      bool isBottom = (chanmap[i]==2) ? true : false;
      //cout << i << " " << isBottom << endl; 
      double amplow = isBottom ? 22 : 22;
      double amphigh = isBottom ? 34 : 30;   
      double chargelow = isBottom ? 2500 : 1700;
      double chargehigh = isBottom ? 4300 : 3120;

      if(aAmpli<amphigh && aAmpli>amplow && aCharge>chargelow && aCharge<chargehigh && aTime<3000 && aTime>2950) {
            for (size_t j = 0; j < wf->size(); j++) {
                Double_t val = avg.at(i)->GetBinContent(j+1);
                val += (wf->at(j) - aBase)/aCharge;
                avg.at(i)->SetBinContent(j+1,val);
            }
        nevents[i]++;
         
      }
     } //channel loop

    if (ev%500 == 0) cout << "Event " << ev << " processed" <<endl;
  } //event loop

  TH1F *h_avg = new TH1F("h_avg","",1300,0,1299);
 //cout << "qui" << endl;
  
  for (size_t i=0; i<ker.size();i++)
  {
    if (ker.at(i) && nevents[i]) //at least one event averaged 
     {
     for (size_t j = 1; j <= ker.at(i)->GetNbinsX(); j++) {
          ker.at(i)->SetBinContent(j,avg.at(i)->GetBinContent(2900+j)/nevents[i]);  
       if (i==16)
         h_avg->SetBinContent(j,ker.at(i)->GetBinContent(j));  
     }
     double norm = 0;
     for (size_t j = 1; j <= ker.at(i)->GetNbinsX(); j++) 
	norm -=  ker.at(i)->GetBinContent(j);
     ker.at(i)->Scale(1./norm);
  }
 }
 gROOT->SetStyle("Plain");

  TCanvas *c0 = new TCanvas("c0", "c0", 10, 10, 900, 600);
  h_avg->SetTitle("Kernel");
  h_avg->GetXaxis()->SetTitle("Sample (arb.)");
  h_avg->GetYaxis()->SetTitle("Counts (arb.)");
  h_avg->GetXaxis()->SetTitleOffset(1.10);
  h_avg->GetYaxis()->SetTitleOffset(1.40);
  //h_avg->DrawCopy();
  //return;
  //h_avg->SaveAs(Form("rootfiles/kernel_run_%d.root",run));
  //h_avg->SaveAs(Form("plots/kernel_run_%d.C",run));
  //c0->SaveAs(Form("rootfiles/c0_run_%d.root",run));
  //c0->SaveAs(Form("rootfiles/c0_run_%d.C",run));
  //c0->SaveAs(Form("plots/kernel_run_%d.png",run));
    
// ***** CROSS CORRELATION *****
    /*
     vector<TCanvas*> *c = new vector<TCanvas*>;
     for(int i=0; i<8; i++){
     c->push_back(new TCanvas(Form("c_%d", i), Form("c_%d", i), 10, 10, 900, 600));
     }
     for(size_t j=0; j<c->size();j++){
     c->at(j)->Divide(2,2);
     }
     */
    //TH1D *h_cc = nullptr;
    
    
  vector<TH1D*> h_charge(imax+1,nullptr);
  vector<TH1D*> h_cc(imax+1,nullptr);
    
    
  for (int ev = 0; ev < reco->GetEntries(); ev++) {
      reco->GetEntry(ev);
      raw->GetEntry(ev);
        
      map<int,vector<double>* > wfs = evRaw->GetWFs();
      for (map<int,vector<double>* >::iterator it=wfs.begin(); it!=wfs.end(); ++it) {
          int i = it->first;
            
          double aCharge = evReco->GetCharge().at(i);
          double aTime = evReco->GetStartTime().at(i);
          double aAmpli = evReco->GetBaseMean().at(i)-evReco->GetYmin().at(i);
          double  aBase = evReco->GetBaseMean().at(i);
          vector<double>* wf = evRaw->GetWF(i);
            
          //if(!h_cc) h_cc = new TH1D(Form("h%d",i),Form("Channel_%d",i),wf->size(),0,wf->size()-1);
            
            
          if (ker.at(i)==nullptr)
          {
             ker.at(i) = (TH1D*) f1->Get(Form("k%d",i));
          }
            
          if (h_charge.at(i)==nullptr)
          {
              h_charge.at(i) = new TH1D(Form("SER_%d",i),Form("charge_%d",i), 2000, -10, 200);
          }
            
          if (h_cc.at(i)==nullptr)
          {
              h_cc.at(i) = new TH1D(Form("Wff%d",i),Form("Wff_%d",i), wf->size(), 0, wf->size()-1);
          }
            
            
          //cout << i << "sono qui" << endl;
          //if (i != ch) continue;
            
            
          double charge = -100;
          double base = 0;
          
          for (size_t l = 0; l < wf->size() - ker.at(i)->GetNbinsX(); l++) {
              double sum = 0;
              for (int k = 1; k <= ker.at(i)->GetNbinsX(); k++) {
                  sum += ker.at(i)->GetBinContent(k)*(wf->at(l + k-1)-aBase);
                  //cout << i << " " << k << " " <<  ker.at(i)->GetBinContent(k) << " " << rwf[l + k] <<  " " << aBase << " " << wf->at(l) <<  " " << sum << endl;
              }
              //Find maximum
              if (l > 100 && l < 8000)
                  if (sum > charge) charge = sum;
              //if (l < 1000) base += sum;
                
              //if(ev == 0)
              h_cc->SetBinContent(l , sum);
          }
          h_charge.at(i)->Fill(charge);
          h_cc.at(i)->Fill(sum);
          //ker.at(i)->DrawCopy();
          //cc->cd(2);
          //h_cc->Draw();
            
      }
        
      if (ev%100 == 0) cout << "Event " << ev << " serred" <<endl;
  }
  cout << "Loop closed, now plotting" << endl;
  //int contatore=0;
  for(size_t k=0; k<h_charge.size(); k++){
        
      if (h_charge[k] == nullptr)
          continue;
      //cout << k << " " << contatore << " " << h_charge.size() << " " <<
      // c->size() << endl;
      //c->at(contatore/4)->cd(contatore%4+1);
      //contatore++;
      h_charge[k]->SetTitle("Single Electron Response");
      h_charge[k]->GetXaxis()->SetTitle("Amplitude (arb.)");
      h_charge[k]->GetYaxis()->SetTitle("Counts (arb.)");
      h_charge[k]->GetXaxis()->SetRangeUser(0,40);
      h_charge[k]->GetXaxis()->SetTitleOffset(1.10);
      h_charge[k]->GetYaxis()->SetTitleOffset(1.30);
        
      //h_charge.at(k)->Draw();
  }
    
  for(size_t l=0; l<h_cc.size(); l++){
        
      if (h_cc[l] == nullptr)
          continue;
      h_cc[l]->SetTitle("Filtered waveform");
      h_cc[l]->GetXaxis()->SetTitle("Sample (arb.)");
      h_cc[l]->GetYaxis()->SetTitle("Counts (arb.)");
      h_cc[l]->GetXaxis()->SetTitleOffset(1.10);
      h_cc[l]->GetYaxis()->SetTitleOffset(1.30);
        
      //h_cc.at(l)->Draw();
  }
    
  //h_charge->SaveAs(Form("rootfiles/ampli_filtered_run_%d.root",run));
  //h_charge->SaveAs(Form("rootfiles/ampli_filtered_run_%d.C",run));
  //c2->SaveAs(Form("rootfiles/c2_run_%d.root",run));
  //c2->SaveAs(Form("rootfiles/c2_run_%d.C",run));
  //c2->SaveAs(Form("plots/ampli_filtered_run_%d.png",run));
    
  TFile* ofile = new TFile(Form("rootfiles/kernel_run_%d.root",run),"RECREATE");
  ofile->cd();
  for (size_t i=0;i<ker.size();i++)
      if (ker.at(i))
          ker.at(i)->Write(ker.at(i)->GetName());
  ofile->Write();
  ofile->Close();
    
  TFile* ofile1 = new TFile(Form("rootfiles/ampli_filtered_run_%d.root",run),"RECREATE");
  ofile1->cd();
  for (size_t i=0;i<h_charge.size();i++)
      if (h_charge.at(i))
          h_charge.at(i)->Write(h_charge.at(i)->GetName());
  ofile1->Write();
  ofile1->Close();
    
  TFile* ofile2 = new TFile(Form("rootfiles/wf_filtered_run_%d.root",run),"RECREATE");
  ofile2->cd();
  for (size_t i=0;i<h_cc.size();i++)
      if (h_cc.at(i))
          h_cc.at(i)->Write(h_cc.at(i)->GetName());
  ofile2->Write();
  ofile2->Close();

  return;
}

