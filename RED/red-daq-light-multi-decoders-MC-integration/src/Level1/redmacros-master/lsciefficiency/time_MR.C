#include "TTree.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TVirtualPad.h"
#include "TChain.h"
#include "TMath.h"
#include <iostream>

Double_t DT(Double_t *x, Double_t *par);


void time_MR()
{
 
  /*
  TFile *_file0 = TFile::Open("run_74_30_int.root");
  //TFile *_file0 = TFile::Open("run_74_int.root");
  TTree* reco = (TTree*) _file0->Get("reco");

  reco->Project("h1bis","(start_time[17]-start_time[0])*2.+16*(evheader.boardtimes[1]-evheader.boardtimes[0])/2.","charge[0]/7500*60>200. && charge[17]/7500.*60.>100. && f90[0]<0.13 && f90[17]<0.13 ");
  */
  TFile *_file0 = TFile::Open("run_74_int.root.light");
  TTree* reco = (TTree*) _file0->Get("t1");

  TH1D* h1bis = new TH1D("h1bis","h1",400,-20,20);
  TH1D* h0bis = new TH1D("h0bis","h0",400,-20,20);
  reco->Project("h1bis","tof","nchan==17 && chargeFar>25000 && chargeNear>12500 && f90Far<0.13 && f90Near<0.13");


 TF1 * d = new TF1("DT",DT,-20.,20.,8);
 d->SetParameters(100.,1.2,1000.,0.5,1,0.4,10,0.5);
 d->SetParName(0,"Baseline");
 d->SetParName(1,"Offset");
 d->SetParName(2,"#Cf Photon");
 d->SetParName(3,"sigma1");
 d->SetParName(4,"sigma2");
 d->SetParName(5,"fraction");
 d->SetParName(6,"#Bk Photon");
 d->SetParName(7,"#Random Bk");
 
 h1bis->Fit(d,"LQ","",-20.,7.5);
 cout << "Ch" << 1 << "--> " << 
   "Baseline: (" << d->GetParameter(0) << " +/- " << d->GetParError(0) << 
   ") Offset: (" << d->GetParameter(1) << " +/- " << d->GetParError(1) << ")" << endl; 
 
 h1bis->DrawCopy();
		     
 //c1->Print("Fit-PMT2-PMT1.gif");
  reco->Project("h0bis","tof","nchan==16 && chargeFar>25000 && chargeNear>12500 && f90Far<0.13 && f90Near<0.13");
  //reco->Project("h0bis","(start_time[16]-start_time[0])*2.+16*(evheader.boardtimes[1]-evheader.boardtimes[0])/2.","charge[0]/7500*60>200. && charge[16]/7500.*60.>100. && f90[0]<0.13 && f90[16]<0.13 ");
 h0bis->Fit(d,"LQ","",-20.,7.5);
 cout << "Ch" << 0 << "--> " << 
   "Baseline: (" << d->GetParameter(0) << " +/- " << d->GetParError(0) << 
   ") Offset: (" << d->GetParameter(1) << " +/- " << d->GetParError(1) << ")" << endl; 
 h0bis->DrawCopy();
 //c1->Print("Fit-PMT0-PMT1.gif");
}
Double_t DT(Double_t *x, Double_t *par)
{

  Double_t DTCf = par[0]/29.98+par[1];
  Double_t NphCf = par[2];
  Double_t sigmaCf = par[3];
  Double_t sigmaCf2 = par[4];
  Double_t f1 = par[5];


  //  Double_t phCf = NphCf/TMath::Sqrt(2*3.1415)/sigmaCf*TMath::Exp(-0.5*(TMath::Power((x[0]-DTCf)/sigmaCf,2)));
  Double_t phCf = NphCf*(f1*1./TMath::Sqrt(2*3.1415)/sigmaCf*TMath::Exp(-0.5*(TMath::Power((x[0]-DTCf)/sigmaCf,2))) + 
			 (1-f1)*1./TMath::Sqrt(2*3.1415)/sigmaCf2*TMath::Exp(-0.5*(TMath::Power((x[0]-DTCf)/sigmaCf2,2))));


  Double_t NphBk = par[6];
  Double_t sigmaBk = par[7];
  Double_t DTBk = par[1]-par[0]/29.98;

  Double_t phBk = NphBk*(f1*1./TMath::Sqrt(2*3.1415)/sigmaCf*TMath::Exp(-0.5*(TMath::Power((x[0]-DTBk)/sigmaCf,2))) + 
			 (1-f1)*1./TMath::Sqrt(2*3.1415)/sigmaCf2*TMath::Exp(-0.5*(TMath::Power((x[0]-DTBk)/sigmaCf2,2))));

  //  Double_t phBk = NphBk/TMath::Sqrt(2*3.1415)/sigmaBk*TMath::Exp(-0.5*(TMath::Power((x[0]-DTBk)/sigmaBk,2)));

  //  Double_t NphCfBs = par[6];
  //  Double_t sigmaCfBs = par[7];
  //  Double_t DTCfBs = par[8];
  //
  //Double_t phCfBs = NphCfBs/TMath::Sqrt(2*3.1415)/sigmaCfBs*TMath::Exp(-0.5*(TMath::Power((x[0]-DTCfBs)/sigmaCfBs,2)));

  Double_t costBk=par[7];

 
  return phCf+phBk+costBk;


}
