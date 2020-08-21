#include <iostream>
#include <ctime>
#include <vector>

#include <TVirtualFFT.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TGraph.h>
#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <TTree.h>
#include <TFile.h>


using namespace std;

void fft(int N = 1000) { 
   TFile* fin = new TFile("fourier.root", "read");
   TTree* wf = (TTree*) fin->Get("wf");

   vector<double>* wf_top = 0; 
   vector<double>* wf_bot = 0;
   wf->SetBranchAddress("v_top", &wf_top);
   wf->SetBranchAddress("v_bot", &wf_bot);

  
// *************************** TOP ***************************************
   TCanvas *c_top = new TCanvas("c_top", "c_top", 1);
   c_top->Divide(1, 2);

   wf->GetEntry(0);

   TGraph* g_top = new TGraph();
   g_top->SetLineColor(kBlue); g_top->SetLineWidth(2); 
   g_top->SetTitle("Top channel; Samples; ADC");


   int n = wf_top->size();
   double range = n*2e-9;
   double *in = new double[n];
   double ps[(int) n/2];

   for (size_t i = 0; i < n; i++) {
      g_top->SetPoint(i, i, wf_top->at(i));
      in[i] = wf_top->at(i);
   }
   c_top->cd(1);
   g_top->Draw("ALP");

   TVirtualFFT *fft = TVirtualFFT::FFT(1, &n, "R2C");
   fft->SetPoints(in);
   fft->Transform();

   double Re, Im;
   for (int i = 0; i < (int) n/2; i++) {
       fft->GetPointComplex(i, Re, Im);
       ps[i] = Re*Re + Im*Im;
   } 

   TGraph* ps_top = new TGraph();
   for (size_t i = 1; i < (int) n/2; i++) ps_top->SetPoint(i - 1, (i - 1)/range, ps[i]);

   ps_top->SetTitle("PowerSpectrum; frequency [Hz]; ADC");
   
   //c_top->cd(2)->SetLogx();
   c_top->cd(2)->SetLogy();
   ps_top->SetLineColor(kGreen+3);  ps_top->SetLineWidth(2);
   //ps_top->SetFillColor(kGreen+3);

   ps_top->Draw("ALF");   
 

 // *************************** AVERAGE LOOP ***************************************     


 
  vector<double> x_bot, x_top, y_bot, y_top; 
  //int N = wf->GetEntries(); 
 

  for (size_t i = 1; i < n/2; i++) {
     x_bot.push_back(i/range); 
     y_bot.push_back(0);
     x_top.push_back(i/range); 
     y_top.push_back(0);
  }

  for (int i = 0; i < N; i++) {
     wf->GetEntry(i);

     double *in_top = new double[n];
     double ps_top[(int) n/2];
     double *in_bot = new double[n];
     double ps_bot[(int) n/2];


     for (size_t i = 0; i < n; i++) {
        in_top[i] = wf_top->at(i);
        in_bot[i] = wf_bot->at(i);
     }


     TVirtualFFT *fft_top = TVirtualFFT::FFT(1, &n, "R2C");
     fft_top->SetPoints(in_top);
     fft_top->Transform();

     double Re_top, Im_top;
     for (int i = 1; i < (int) n/2; i++) {
        fft_top->GetPointComplex(i, Re_top, Im_top);
        y_top[i] += (Re_top*Re_top + Im_top*Im_top)/N;
     }

     TVirtualFFT *fft_bot = TVirtualFFT::FFT(1, &n, "R2C");
     fft_bot->SetPoints(in_bot);
     fft_bot->Transform();

     double Re_bot, Im_bot;
     for (int i = 1; i < (int) n/2; i++) {
        fft_bot->GetPointComplex(i, Re_bot, Im_bot);
        y_bot[i] += (Re_bot*Re_bot + Im_bot*Im_bot)/N;
     }

     if(!i%100) cout << "event " << i << " processed" << endl; 
  }


  TGraph* avg_top = new TGraph();
  TGraph* avg_bot = new TGraph();
    
  for (int i = 0; i < x_top.size(); i++) {
     avg_top->SetPoint(i, x_top[i], y_top[i]);
     avg_bot->SetPoint(i, x_bot[i], y_bot[i]);
  }


  TCanvas* c_avg = new TCanvas("c_avg", "c_avg", 1); 
  c_avg->Divide(1, 2);

  avg_top->SetTitle("PowerSpectrum - TOP; frequency [Hz]; ADC");
  avg_top->SetLineColor(kGreen+3);  avg_top->SetLineWidth(2); 

  avg_bot->SetTitle("PowerSpectrum - BOTTOM; frequency [Hz]; ADC");
  avg_bot->SetLineColor(kGreen+3);  avg_bot->SetLineWidth(2); 

  c_avg->cd(1)->SetLogy();
  avg_top->Draw("AL");
  c_avg->cd(2)->SetLogy();
  avg_bot->Draw("AL");

}


