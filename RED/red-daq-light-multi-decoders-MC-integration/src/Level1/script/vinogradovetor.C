#include <iostream>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <vector>
#include <sstream>
#include <fstream>

#include <TFile.h>
#include <TF1.h>
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

#include "TCanvas.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

int f(int n) { return TMath::Factorial(n); }

int B(int i, int k) {
    if(i == 0 && k == 0) return 1;
    else if(i == 0 && k > 0) return 0;
    else return f(k)*f(k - 1)/f(i)/f(i - 1)/f(k - i);
}

class LL {
private:
    vector<double> n_;
    double N_;
public:
    LL (vector<double> n, double N) : n_(n), N_(N) {}
    double Evaluate(double *x, double *p) {
        double lambda = x[0];
        double P = x[1];
        
        double ll = 0;
        for (size_t k = 0; k < n_.size(); k++) {
            double nk = n_[k]/N_;
            
            double Sum = 0;
            for(int i = 0; i <= k; i++) Sum += B(i, k)*pow(lambda*(1 - P), i)*pow(P, k - i)/f(k);
            double fk = exp(-lambda)*Sum;
            
            if (fk > 0) ll += 2*(fk - nk + nk*log(nk/fk));
        }
        
        
        return ll;
    } 
    
    
};

void vinogradovetor(int run, int nbin=300, int nmin=-1500, int nmax=20000, int npeaks=5){


    TFile *f = new TFile(Form("../../rootfiles/run_%d.root", run), "read");
    TTree *data = (TTree*)f->Get("reco");
    
    gROOT->SetStyle("Plain");
    
////////////// CHARGE ////////////////////

    TCanvas *c1 = new TCanvas("c1","c1",10,10,1800,1200);
    
    //c1->Divide(8,4); //config. LNS
    c1->Divide(7,4); //config. Napoli Dic. 2018
    
    TCanvas *c2 = new TCanvas("c2","c2",10,10,1800,1200);
    
    //c2->Divide(8,4); //config. LNS
    c2->Divide(7,4); //config. Napoli Dic. 2018
    
    TCanvas *c3 = new TCanvas("c3","c3",10,10,1800,1200);
    
    //c3->Divide(8,4); //config. LNS
    c3->Divide(7,4); //config. Napoli Dic. 2018
    
    TGraphErrors *gk = new TGraphErrors();
    
    TH1F *h1[32];
    TF1 *sum;
    for (int i = 0; i < 32; i++) {
        c1->cd(i + 1);/*->SetLogy();*/
        h1[i] = new TH1F(Form("h1%d", i),Form("h1%d", i),nbin,nmin,nmax);
        //data->Draw(Form("charge[%d]>>h1%d", i+16, i));
        data->Draw(Form("charge[%d]>>h1%d", i, i)); // configurazione Napoli Dicembre 2018
        //data->Draw(Form("baseline_mean[%d]-ymin[%d]>>h1%d", i, i, i)); // per avere le SER in ampiezza

        h1[i]->SetTitle("Single Electron Response");
        h1[i]->GetXaxis()->SetTitle("Charge [ADC counts]");
        h1[i]->GetXaxis()->SetTitleOffset(1.16);
        h1[i]->GetYaxis()->SetTitle("Counts");
        h1[i]->GetYaxis()->SetTitleOffset(1.16);
        h1[i]->SetLineWidth(1);
        h1[i]->SetLineColor(kBlue);
        
    }

    
    //ofstream ofile(Form("cfg/naples/ser_%d.cfg", run)); //config. Napoli Dic. 2018
    //ofstream ofile(Form("cfg/naples/ser_%d_kdup.cfg", run)); //config. Napoli Dic. 2018
    //ofstream ofile2(Form("cfg/naples/kdup_%d_top.cfg", run)); //config. Napoli Dic. 2018
    //ofstream ofile2(Form("cfg/naples/kdup_%d_bottom.cfg", run)); //config. Napoli Dic. 2018
    
    double kdup_mean_bot = 0;
    double kdup_mean_top = 0;
    int nkdup_top = 0;
    int nkdup_bot = 0;
    
    double kdup_mean_bot_ll = 0;
    double kdup_mean_top_ll = 0;

    
    for (int i = 0; i < 32; i++)
    {

      if (h1[i]->GetEntries()==0) {
            continue;
        }
        
        c1->cd(i + 1)->SetLogy();
        h1[i]->Rebin(2);

/////////////// FIT & DRAW /////////////////////////

    //Use TSpectrum to find the peak candidates
    int npeak;
    TSpectrum *s = new TSpectrum();
    int nfound = s->Search(h1[i],4,"goff",0.0005); //(nome isto,sigma,opt,soglia)
        if (nfound<=npeaks) {
            npeak=nfound;
        } else {
            npeak=npeaks;
        }
 
    double *xus = (double *)s->GetPositionX(); //array with X-positions of the centroids found by TSpectrum
    double *y = (double *)s->GetPositionY();
 
    vector<double> x(xus, xus + npeak);
    sort(x.begin(), x.end());
    float dmu = 0;
    for (int j = 0; j < npeak - 1; j++) dmu += x[j + 1] - x[j];
    dmu = dmu/npeak;

    printf("Found %d candidate peaks to fit\n",nfound);
    printf("Found %d useful peaks to fit\n",npeak);
    printf("Now fitting\n");
    
    //Loop on all found peaks.

     TF1 *g[npeak];
     for (int p=0;p<npeak;p++) {
        g[p] = new TF1("gaus","gaus",x[p] -dmu/3, x[p] + dmu/3);
 
        g[p]->SetLineWidth(2);
        g[p]->SetLineColor(kRed);
        h1[i]->Fit(g[p],"R+Q+0");
    }

    string sgaus = "gaus(0) ";
        for (int ss = 1; ss < npeak; ss++) sgaus += Form("+ gaus(%d) ", 3*ss);
        
    sum = new TF1("mysum",sgaus.c_str(),x[0] - dmu/2, x[npeak - 1] + dmu/2);
    sum->SetNpx(1000);
    for (int k=0;k<3*npeak;k++){
       sum->SetParameter(k,g[(k-k%3)/3]->GetParameter(k%3));
       if(!(k-1)%3) sum->FixParameter(k,g[(k-k%3)/3]->GetParameter(k%3));
       if(!(k-1)%3) sum->SetParLimits(k, sum->GetParameter(k) - dmu/3,sum->GetParameter(k) + dmu/3);
      }

    h1[i]->Fit(sum,"R+Q");

    for (int k=0;k<3*npeak;k++){
       sum->SetParameter(k,g[(k-k%3)/3]->GetParameter(k%3));
       if(!(k-1)%3) sum->ReleaseParameter(k);
    }

    sum->SetLineWidth(2);
    sum->SetLineColor(kBlack);
    h1[i]->Fit(sum,"R+Q");


    h1[i]->DrawCopy();
        
    c2->cd(i + 1);
         
    TGraphErrors *gr[32];
    gr[i] = new TGraphErrors();
    gr[i]->SetName(Form("gr%d", i));
    gr[i]->SetLineColor(kRed);
    gr[i]->SetLineWidth(1);
    gr[i]->SetMarkerSize(2);
    gr[i]->SetMarkerColor(kBlue);
    gr[i]->SetMarkerStyle(21);
    gr[i]->SetTitle("Calibration curve");
    gr[i]->GetXaxis()->SetTitle("PE");
    gr[i]->GetYaxis()->SetTitle("Charge [ADC counts]");
    //gr[i]->GetYaxis()->SetTitle("Amplitude [ADC counts]");
    gr[i]->GetYaxis()->SetTitleOffset(1.16);
         
    for (int j = 0; j < npeak; j++) {
    gr[i]->SetPoint(j, j, sum->GetParameter(1 + j*3));
    gr[i]->SetPointError(j, 0, sum->GetParError(2 + j*3));
    }
    TF1 *myline = new TF1("myline","pol1(0)");
    gr[i]->Fit(myline,"Q");
    gStyle->SetOptFit(111);
    gr[i]->Draw("AP");
        
    cout << "Mean_first_peak: " << sum->GetParameter(4) << " " << "Sigma_first_peak: " << sum->GetParameter(5) << endl;
         
    cout << "Res: " << sum->GetParameter(5)/sum->GetParameter(4) << endl;
         
    cout << "" << i << " " << "SNR: " << sum->GetParameter(4)/sum->GetParameter(2) << " +/- " << (sum->GetParError(4)/sum->GetParameter(4))+(sum->GetParError(3)/sum->GetParameter(2)) << endl;
        
    vector<double>N;
    double N_tot = h1[i]->GetEntries();
    for (int l=0; l<npeak; l++){
            
        N.push_back(sum->GetParameter(3*l)*sqrt(2*TMath::Pi())*sum->GetParameter(2+3*l)/h1[i]->GetBinWidth(1));
        
        cout << N.back() << endl;
        
    }
        cout << N_tot << endl;
        
    double m(0), m2(0);
        
    for (int l=0; l<npeak; l++){
        
        m += l*N[l]/N_tot;
            
        m2 += l*l*N[l]/N_tot;
        
    }
        
    double sig2 = m2-m*m;
    double p = (sig2-m)/(sig2+m);
    double k_dup = p/(1-p);
    
    LL* oll = new LL(N, N_tot);
        
    TF2* f2 = new TF2("f2", oll, &LL::Evaluate, 0.01, 10, 0.01, 0.99, 0, "LL", "Evaluate");
        
    double lambda, P;
    f2->GetMinimumXY(lambda, P);
        
    double k_dup_ll = P/(1-P);
    
    if (i<4) {

        kdup_mean_bot += k_dup;
        nkdup_bot++;
        
        kdup_mean_bot_ll += k_dup_ll;
        
    }else{
    
        kdup_mean_top += k_dup;
        nkdup_top++;
        
        kdup_mean_top_ll += k_dup_ll;
    }
    
    gk->SetPoint(i,i,k_dup_ll);
        

    cout << "" << i << " " << "N_tot: " << N_tot << " " << "m: " << m << " " << "sig: " << sig2 << endl;
    cout << "" << i << " " << "p: " << p << " " << "k_dup: " << k_dup << " k_dup_ll: " << k_dup_ll << endl;
        
    //ofile << myline->GetParameter(1)*(1+k_dup) << endl; //scrive il valore di SER/1+k_dup
    //ofile2 << i << " " << k_dup << endl; //scrive il valore del canale e di k_dup
        
    delete oll, f2;
         
    }
    
    
    
    c3->cd();
    gk->Draw("ALP");
    
    cout << "kdup_mean_bot: " << kdup_mean_bot/nkdup_bot << endl;
    cout << "kdup_mean_top: " << kdup_mean_top/nkdup_top << endl;
    cout << "kdup_mean_bot_ll: " << kdup_mean_bot_ll/nkdup_bot << endl;
    cout << "kdup_mean_top_ll: " << kdup_mean_top_ll/nkdup_top << endl;
    cout << "ch. analized bot: " << nkdup_bot << endl;
    cout << "ch. analized top: " << nkdup_top << endl;

 
    
    //ofile.close();
    //ofile2.close();

//c1->cd();
    

//c1->SaveAs(Form("plot/naples/ser_run_%d.png", run)); //config. Napoli Dic. 2018
//c1->SaveAs(Form("plot/naples/ser_run_%d_kdup.png", run)); //config. Napoli Dic. 2018

    
//c2->cd();

//c2->SaveAs(Form("plot/naples/ser_cal_run_%d.png", run)); //config. Napoli Dic. 2018
//c2->SaveAs(Form("plot/naples/ser_cal_run_%d_kdup.png", run)); //config. Napoli Dic. 2018


    

}
