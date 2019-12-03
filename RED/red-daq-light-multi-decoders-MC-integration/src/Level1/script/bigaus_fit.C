#include <iostream>
#include <cstdlib>
#include <cstring>
#include <algorithm> 
#include <vector> 
#include <sstream> 
#include <fstream>

#include <TFile.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TChain.h>
#include <TLegend.h>
#include <TCut.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TMath.h>

using namespace std;
/*
Double_t mybigaus(Double_t *x, Double_t *par) {
    
    Double_t mu1,mu2,sigma1,sigma2,rho;
    Double_t f,z;
    
    mu1=par[1];
    mu2=par[3];
    sigma1=par[2];
    sigma2=par[4];
    rho=par[5];
    
    z = (x[0]-mu1)*(x[0]-mu1)/sigma1/sigma1 -2.*rho*(x[0]-mu1)*(x[1]-mu2)/sigma1/sigma2
    + (x[1]-mu2)*(x[1]-mu2)/sigma2/sigma2;
    f = par[0]/(2.*TMath::Pi()*sigma1*sigma2*TMath::Sqrt(1-rho*rho));
    f *= TMath::Exp(-z/2./(1-rho*rho));
    
    return f;
    
}
*/

void bigaus_fit(TString filename, int run)
{
    TFile *file = new TFile(filename, "read");
    
    TTree *reco = (TTree*)gDirectory->Get("reco");

    //ofstream ofile("par_fit.out");

//Style
    gROOT->SetStyle("Plain");
    gStyle->SetOptFit();
    
//Canvas
    TCanvas *c1 = new TCanvas("c1","c1",10,10,900,600);

/*
//Cuts
    gROOT->ProcessLine(".x ../../rootfiles/hibe502.C");
    gROOT->ProcessLine(".x ../../rootfiles/lowbe502.C");
    gROOT->ProcessLine(".x ../../rootfiles/lip1_502.C");
    gROOT->ProcessLine(".x ../../rootfiles/lip2_502.C");
    gROOT->ProcessLine(".x ../../rootfiles/liC502.C");
*/
////////////// Beryllium loci ////////////////////
    
    TH2F *h = new TH2F("h","h",300,1000,14000,300,1000,14000);
    reco->Draw("ymax[9]-1975:ymax[10]-baseline_mean[10]>>h");
    
    h->GetXaxis()->SetTitle("E [arb.]");
    h->GetYaxis()->SetTitle("#Delta E [arb.]");
    h->GetYaxis()->SetTitleOffset(1.32);
    h->GetXaxis()->SetTitleOffset(1.12);
    h->SetMaximum(25);
    h->SetTitle("Ion energy loss spectrum");

    TF2 *f = new TF2("f","bigaus",6320,7253.33,6460,8000);
    f->SetLineColor(kRed);
    
    h->Fit("f","R+");
    h->Draw("colz");
    
    TF2 *f1 = new TF2("f1","bigaus",9000,10300,5500,6600);
    f1->SetLineColor(kGreen+2);
    
    h->Fit("f1","R+");
    h->Draw("samecolz");


////////////// Lithium loci ////////////////////
    
    TF2 *f2 = new TF2("f2","bigaus",6600,7100,4000,5000); // Li+p (1)
    f2->SetParameters(5.01e6,6.82e3,7.44e2,4.45e3,3.19e2,-8.76e-1);
    f2->SetLineColor(kYellow+1);
    
    h->Fit("f2","R+");
    h->Draw("samecolz");

    TF2 *f3 = new TF2("f3","bigaus",12000,12800,2200,3500); // Li+p (2)
    f3->SetLineColor(kViolet);
    
    h->Fit("f3","R+");
    h->Draw("samecolz");

    TF2 *f4 = new TF2("f4","bigaus",13000,13600,2200,3400); // Li+C
    f4->SetLineColor(kBlue);
    
    h->Fit("f4","R+");
    h->Draw("samecolz");
    
    
    TLegend *legend = new TLegend(0.23,0.68,0.37,0.88);
    legend->SetBorderSize(0);
    legend->SetTextFont(72);
    legend->SetTextSize(0.03);
    TLegendEntry *be1 = legend->AddEntry(f,"p(7Li,7Be)n","l");
    TLegendEntry *be2 = legend->AddEntry(f1,"p(7Li,7Be)n","l");
    TLegendEntry *li1 = legend->AddEntry(f2,"p(7Li,7Li)p'","l");
    TLegendEntry *li2 = legend->AddEntry(f3,"p(7Li,7Li)p'","l");
    TLegendEntry *li3 = legend->AddEntry(f4,"7Li(12C,12C)7Li","l");
    legend->Draw();
/*
    ofile << f->GetParameter(1) << " " << f->GetParError(1) << " " << f->GetParameter(3) << " " << f->GetParError(3) << " " << f->GetParameter(5) << endl
    
    << f1->GetParameter(1) << " " << f1->GetParError(1) << " " << f1->GetParameter(3) << " " << f1->GetParError(3) << " " << f1->GetParameter(5) << endl
    
    << f2->GetParameter(1) << " " << f2->GetParError(1) << " " << f2->GetParameter(3) << " " << f2->GetParError(3) << " " << f2->GetParameter(5) << endl
    
    << f3->GetParameter(1) << " " << f3->GetParError(1) << " " << f3->GetParameter(3) << " " << f3->GetParError(3) << " " << f3->GetParameter(5) << endl
    
    << f4->GetParameter(1) << " " << f4->GetParError(1) << " " << f4->GetParameter(3) << " " << f4->GetParError(3) << " " << f4->GetParameter(5) << endl;
*/

//    c1->SaveAs(Form("SiCal_%d.png",run));
    


}
