#include <iostream>
#include <cstring>

#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>

using namespace std;

void extreme_fit(string fname = "fit.root", double a = -2000, double b = 4000,  int ngaus = 3) {

   TFile* f = new TFile(fname.c_str(), "read");
   TH1F* hh = (TH1F*)f->Get("hh");
   hh->GetXaxis()->SetRangeUser(a,b);
   
   TH1F* htot = (TH1F*)hh->Clone();
   
   TCanvas *c = new TCanvas("c","c", 1);
   c->Divide(2,1);
    
   TH1F* h[ngaus];
   TF1* g[ngaus];
   for (int i = 0; i < ngaus; i++) {
      h[i] = (TH1F*)hh->Clone();
      g[i] = new TF1(Form("g%d", i), "gaus", a, b); 
   }
   
   string stot = "gaus(0)";
   for (int i = 1; i < ngaus; i++) stot += Form(" + gaus(%d)", 3*i);
   
   TF1* gtot = new TF1("gtot", stot.c_str(), a, b);  
   gtot->SetLineWidth(3);
   gtot->SetLineColor(kBlue);
   

   c->cd(1);
   
   h[0]->Draw();
   h[0]->Fit("g0", "LR", "", a, h[0]->GetBinCenter(h[0]->GetMaximumBin())); 
      
   for (int i = 1; i < ngaus; i++) {
      for (int j = 0; j < h[i - 1]->GetXaxis()->GetNbins(); j++) 
         h[i]->SetBinContent(j + 1, h[i - 1]->GetBinContent(j + 1) - g[i - 1]->Eval(h[i - 1]->GetBinCenter(j + 1)));
  
      h[i]->Draw("same");
      h[i]->Fit(Form("g%d", i), "LR", "", h[i - 1]->GetBinCenter(h[i - 1]->GetMaximumBin()), h[i]->GetBinCenter(h[i]->GetMaximumBin()));
   }  
 
   c->cd(2);
   htot->Draw();
   
   for (int i = 0; i < 3*ngaus; i++) gtot->SetParameter(i, g[i/3]->GetParameter(i%3));
   for (int i = 0; i < ngaus; i++) g[i]->Draw("same");

   htot->Fit("gtot", "LR", "", a, b);  

}
