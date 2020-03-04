#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"

int run_715()
{
  TFile* f1 = new TFile("run_715.root");
  TFile* f2 = new TFile("run_715_rev.root");
  
  TTree* t1 = (TTree*) f1->Get("reco");
  TTree* t2 = (TTree*) f2->Get("reco");

  TH1D* h1 = new TH1D("h1","Ref reconstruction",600,2700,3300);
  TH1D* h2 = new TH1D("h2","Rev reconstruction",600,2700,3300);
  h2->SetLineColor(kRed);

  t1->Project("h1","xmin","ymin<14000 && Iteration$<9");
  t2->Project("h2","xmin","ymin<14000 && Iteration$<9");

 

  TF1* gaus = new TF1("gaus","gaus(0)",2700,3300);
  TF1* pol0 = new TF1("pol0","pol0(0)",2700,3300);
  TF1* fun1 = new TF1("fun","gaus(0)+pol0(3)",2700,3000);
  TF1* fun2 = new TF1("fun","gaus(0)+pol0(3)",2700,3000);

  //h1
  Double_t par[4];
  h1->Fit(gaus,"","",2880,2920);
  gaus->GetParameters(&par[0]);
  h1->Fit(pol0);
  pol0->GetParameters(&par[3]);
  fun1->SetParameters(par);
  h1->Fit(fun1,"0"); //option: function not drawn

  //h2
  h2->Fit(gaus,"","",2880,2920);
  gaus->GetParameters(&par[0]);
  h2->Fit(pol0);
  pol0->GetParameters(&par[3]);
  fun2->SetParameters(par);
  h2->Fit(fun2,"0");

  //Now subtract continuum
  for (Int_t i=1;i<=h1->GetNbinsX();i++)
    {
      Double_t val = h1->GetBinContent(i);
      val -= fun1->GetParameter(3);
      h1->SetBinContent(i,val);
    }
  
  for (Int_t i=1;i<=h2->GetNbinsX();i++)
    {
      Double_t val = h2->GetBinContent(i);
      val -= fun2->GetParameter(3);
      h2->SetBinContent(i,val);
    }
  cout << h1->GetTitle() << "--> " << 
    fun1->GetParameter(0)*TMath::Sqrt(2.*TMath::Pi())*
    fun1->GetParameter(2)/h1->GetBinWidth(1) << " counts" << endl;
  cout << h2->GetTitle() << "--> " << 
    fun2->GetParameter(0)*TMath::Sqrt(2.*TMath::Pi())*
    fun2->GetParameter(2)/h2->GetBinWidth(1) << " counts" << endl;
  

  h1->Draw();
  h2->Draw("same");
  return 0;

}
