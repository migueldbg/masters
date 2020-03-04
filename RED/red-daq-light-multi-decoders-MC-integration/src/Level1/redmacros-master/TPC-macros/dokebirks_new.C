#include "TGraph.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TPad.h"
#include "TH2F.h"
#include "TMath.h"
#include "TFile.h"
#include "TLegend.h"
#include "TF1.h"
#include "TStyle.h"
#include <vector>
#include <iostream>


using namespace std ;
using namespace TMath ;

TGraph *gloss;
TGraphErrors *grA, *grC ;
vector<double> ptx, pty, pte, ptf; 

double fdoke(double *x, double*p) {
  double E = x[0];
  double A = p[0];
  double C = p[1];
  double B = A/(1-C);
  double dEdx = gloss->Eval(x[0]);
  double reco = A*dEdx/ (1 + B*dEdx) + C;
  double alpha = 0.21; 
  double Nq    = E/19.5e-3 ; 
  double Ni    = Nq/(1.+alpha) ; 
  double Nx    = Nq-Ni ; 
  double Ng    = Nx + Ni*reco  ; 
  return Ng / Nq  ; 

}
double fdoke_field(double *x, double*p) {
  double E = x[0];
  double A = p[0];
  double C = p[1];
  double D = p[2];
  double field = p[3];
  C *= Exp(-D*field);
  double B = A/(1-C);
  double dEdx = gloss->Eval(x[0]); 
  double reco = A*dEdx/ (1 + B*dEdx) + C;
  double alpha = 0.21; 
  double Nq    = E/19.5e-3 ; 
  double Ni    = Nq/(1.+alpha) ; 
  double Nx    = Nq-Ni ; 
  double Ng    = Nx + Ni*reco  ; 
  return Ng / Nq  ; 

  

}
//////////////////////////////////////////////////////////////////////
//   G4DS recobination
//////////////////////////////////////////////////////////////////////
double GetRecoProb(double InitialKinEne) {
  double myRecoProb = 0 ;
    // new model obtained with 37Ar and 83mKr constraint - oct '15 6th
    double p0=   2.96766e-01 ;//  7.50908e-03   0.00000e+00  -2.61950e-01
    double p1=   3.95496e+00 ;// 2.77124e+00  -0.00000e+00  -1.52624e+01
    double p2=  -5.17812e-01 ;// 2.69001e-01   0.00000e+00  -1.88046e+02
    double p3=  -1.38485e-02 ;// 1.13375e-03  -0.00000e+00   4.47780e+03
    double p4=   9.12436e-01 ;// 1.34073e-02  -0.00000e+00  -9.04345e+02
    double p5=   6.61046e-01 ;// 1.59214e-04  -0.00000e+00   2.03157e+01

    myRecoProb = p0 * (1 - p1*exp( p2 * InitialKinEne )) *
              exp( (p3) * pow( InitialKinEne, (p4) ) ) +
              p5 ;
  

  return myRecoProb;
}
//////////////////////////////////////////////////////////////////////
//   LY ratio from G4DS
//////////////////////////////////////////////////////////////////////
double LYat200V (double *x , double *p) {

  double alpha = p[0] ; 
  double E     = x[0] ;   
  double Nq    = E/19.5e-3 ; 
  double Ni    = Nq/(1.+alpha) ; 
  double Nx    = Nq-Ni ; 
  double Ng    = Nx + Ni*GetRecoProb(E )  ; 
  return Ng / Nq  ; 

}
//-------------------------------------------------------------------------
//    Chi2
//-------------------------------------------------------------------------
double Chi2(double *p) {
  
  double A = p[0] ; 
  double C = p[1] ; 
  
  double chi2 = 0 ;
  cout << "aaa" << endl ;
  for(int i=0;i<int(ptx.size());++i) {

    double x0   = ptx[i];
    double err  = pte[i] ; 
    double data = pty[i] ; 
    
    cout << x0 << " " << err << " " << data << endl ;
    double model = fdoke(&x0, p) ; 
    chi2 += pow(data - model,2)/err/err ;      
  
  }  
  cout << chi2 << endl ;
  return chi2 ;

}
//-------------------------------------------------------------------------
//    Chi2
//-------------------------------------------------------------------------
double Chi2field(double *p) {
  
  double A     = p[0] ; 
  double C     = p[1] ; 
  double D     = p[2] ; 
  
  double chi2 = 0 ;
  for(int i=0;i<int(ptx.size());++i) {

    double x0    = ptx[i];
    double err   = pte[i] ; 
    double data  = pty[i] ; 
    double field = ptf[i] ; 
    p[3]         = field ;
    
    double model = fdoke_field(&x0, p) ; 
    chi2 += pow(data - model,2)/err/err ;      
  
  }  
  return chi2 ;

}

//-------------------------------------------------------------------------
//    FCN
//-------------------------------------------------------------------------
void minFunc(int& nDim, double* gout, double& result, double par[], int flg) {
  cout << "aaa" <<endl ;
  result = Chi2(par);	
}

//-------------------------------------------------------------------------
//    FCN
//-------------------------------------------------------------------------
void minFunc2(int& nDim, double* gout, double& result, double par[], int flg) {
  result = Chi2field(par);	
}




void exe(int field = 200) {
  TFile *fin = TFile::Open("eloss.root");  
  gloss = (TGraph*) fin->Get("gr_eloss");

  TFile *f2 = TFile::Open("dataER.root");  

  char cgr[50], cgram[50], cgrba[50];
  
  sprintf(cgr,"g%i",field);
  sprintf(cgram,"gam%i",field);
  sprintf(cgrba,"gba%i",field);

  TGraphErrors *gr   = (TGraphErrors*) f2->Get(cgr);
  TGraphErrors *grAm = (TGraphErrors*) f2->Get(cgram);
  TGraphErrors *grBa = (TGraphErrors*) f2->Get(cgrba);
   
  gr->Draw("ap");
  grAm->Draw("p");
  grBa->Draw("p");
  TF1 *fun = new TF1("fun",fdoke,0,500,2);
  //fun->SetParameters(1,3);
  fun->SetLineColor(2);
  //gr->Fit("fun");
  //gr->Fit("fun");
  
  ptx.clear();
  pty.clear();
  pte.clear();
  
  
  int N = gr->GetN();
  
  for(int i=0;i<N;++i) {
    ptx.push_back(gr->GetX()[i]);
    pty.push_back(gr->GetY()[i]);
    pte.push_back(gr->GetEY()[i]);
  }
  ptx.push_back(grAm->GetX()[0]);
  pty.push_back(grAm->GetY()[0]);
  pte.push_back(grAm->GetEY()[0]);
  ptx.push_back(grBa->GetX()[0]);
  pty.push_back(grBa->GetY()[0]);
  pte.push_back(grBa->GetEY()[0]);
  //////////////////////////////////////////////////////////////////////////////////////
  // Minimization
  double arglist[2];
  int ierflg=0;
  TMinuit *minLogL = new TMinuit(2);
  minLogL->SetFCN(minFunc); 
  
  //minLogL->mnexcm("SET PRINT -1", arglist, 1, ierflg);
  minLogL->Command("SET STRATEGY 1");
  
  minLogL->mnparm(0, "A",   1.0,     0.1, 0, 100,  ierflg);  
  minLogL->mnparm(1, "C",   0.5,     0.1, 0,   2,  ierflg);  
   

  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  
  double A, C, eA, eC; 
  minLogL->GetParameter(0,A,eA);
  minLogL->GetParameter(1,C,eC);
  fun->SetParameters(A,C);
  fun->Draw("same");
  int nn = grA->GetN();
  grA->SetPoint(nn, field, A);
  grA->SetPointError(nn, 0, eA);
  grC->SetPoint(nn, field, C);
  grC->SetPointError(nn, 0, eC);
}






void exe2() {
  gStyle->SetOptStat(0);
  TCanvas *cn = new TCanvas("cn","cn");
  
  TFile *fin = TFile::Open("eloss.root");  
  gloss = (TGraph*) fin->Get("gr_eloss");

  TFile *f2 = TFile::Open("dataER.root");  
  ptx.clear();
  pty.clear();
  pte.clear();
  ptf.clear();
  char cgr[50], cgram[50], cgrba[50];
  int field ;
  
  TH2F *hh = new TH2F("hh","",100,0,300,100,0.2,1.1);
  hh->GetXaxis()->SetTitle("Energy [keV_{ee}]");
  hh->GetYaxis()->SetTitle("S1 / S1_{0}");
  
  
  hh->Draw();
  TGraph *m1;
  TGraph *m2;
  TGraph *m3;
  
  for(int i =0 ;i<4;++i) {
    int col ;
    if(i == 0 ) {field = 50; col = kGreen-2; }
    if(i == 1 ) {field = 100;col = kOrange -3; }
    if(i == 2 ) {field = 200;col = kRed-3; }
    if(i == 3 ) {field = 500;col = kAzure-1; }
    
    sprintf(cgr,"g%i",field);
    sprintf(cgram,"gam%i",field);
    sprintf(cgrba,"gba%i",field);

    TGraphErrors *gr   = (TGraphErrors*) f2->Get(cgr);
    TGraphErrors *grAm = (TGraphErrors*) f2->Get(cgram);
    TGraphErrors *grBa = (TGraphErrors*) f2->Get(cgrba);
    
    gr->SetMarkerColor(col);
    gr->SetLineColor(col);
    gr->SetMarkerStyle(8);
    gr->SetMarkerSize(1.);
    
    
    grAm->SetMarkerColor(col);
    grAm->SetLineColor(col);
    grAm->SetMarkerStyle(22);
    
    grBa->SetMarkerColor(col);
    grBa->SetLineColor(col);
    grBa->SetMarkerStyle(23);
    
    
    gr->Draw("p");
    grAm->Draw("p");
    grBa->Draw("p");
    m1 = new TGraph ; m1->SetMarkerStyle(gr->GetMarkerStyle());
    m2 = new TGraph ; m2->SetMarkerStyle(grAm->GetMarkerStyle());m2->SetMarkerSize(grAm->GetMarkerSize());
    m3 = new TGraph ; m3->SetMarkerStyle(grBa->GetMarkerStyle());m3->SetMarkerSize(grBa->GetMarkerSize());
    
    int N = gr->GetN();

    for(int i=0;i<N;++i) {
      ptx.push_back(gr->GetX()[i]);
      pty.push_back(gr->GetY()[i]);
      pte.push_back(gr->GetEY()[i]);
      ptf.push_back(field);
    }
    ptx.push_back(grAm->GetX()[0]);
    pty.push_back(grAm->GetY()[0]);
    pte.push_back(grAm->GetEY()[0]);
    ptf.push_back(field);

    ptx.push_back(grBa->GetX()[0]);
    pty.push_back(grBa->GetY()[0]);
    pte.push_back(grBa->GetEY()[0]);
    ptf.push_back(field);

  }
  //////////////////////////////////////////////////////////////////////////////////////
  // Minimization
  double arglist[3];
  int ierflg=0;
  TMinuit *minLogL = new TMinuit(2);
  minLogL->SetFCN(minFunc2); 
  
  //minLogL->mnexcm("SET PRINT -1", arglist, 1, ierflg);
  minLogL->Command("SET STRATEGY 1");
  
  minLogL->mnparm(0, "A",   1.0,     0.1, 0, 100,  ierflg);  
  minLogL->mnparm(1, "C",   0.5,     0.1, 0,  10,  ierflg);  
  minLogL->mnparm(2, "D",   0.004,     0.001, 0,  20,  ierflg);  
   

  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  TF1 *fun1 = new TF1("fun",fdoke_field,0,500,4);
  TF1 *fun2 = new TF1("fun",fdoke_field,0,500,4);
  TF1 *fun3 = new TF1("fun",fdoke_field,0,500,4);
  TF1 *fun4 = new TF1("fun",fdoke_field,0,500,4);
  fun1->SetLineColor(1);
  fun2->SetLineColor(2);
  fun3->SetLineColor(3);
  fun4->SetLineColor(4);
  double A, C,D,  eA, eC, eD; 
  minLogL->GetParameter(0,A,eA);
  minLogL->GetParameter(1,C,eC);
  minLogL->GetParameter(2,D,eD);
  fun1->SetParameters(A,C,D,50);
  fun2->SetParameters(A,C,D,100);
  fun3->SetParameters(A,C,D,200);
  fun4->SetParameters(A,C,D,500);

  fun1->SetLineColor(kGreen-2);
  fun2->SetLineColor(kOrange -3);
  fun3->SetLineColor(kRed-3);
  fun4->SetLineColor(kAzure-1);


  fun1->Draw("same");
  fun2->Draw("same");
  fun3->Draw("same");
  fun4->Draw("same");
  
  TF1 *fP = new TF1("fP",LYat200V,0,350,1);
  fP->SetParameter(0,0.21);
  //fP->Draw("same");
  TLegend * leg = new TLegend (.2, .15, .7,.4) ; 
  leg->SetNColumns(2) ; 
  leg->SetBorderSize(0);
//  leg->SetTextFont(102) ; 
//  leg->SetBorderSize(1) ; 
//  leg->SetFillStyle(0) ; 
  leg->AddEntry(fun1, "50 V/cm", "l") ; 
  leg->AddEntry(m1, "Compton e^{-}", "p") ; 
  leg->AddEntry(fun2, "100 V/cm", "l") ; 
  leg->AddEntry(m2, "^{241}Am", "p") ; 
  leg->AddEntry(fun3, "200 V/cm", "l") ; 
  leg->AddEntry(m3, "^{133}Ba", "p") ; 
  leg->AddEntry(fun4, "500 V/cm", "l") ; 

  //leg->AddEntry(fre, "model 200 V/cm", "l") ; 
  leg->Draw("same") ; 
  cn->SaveAs("../figures/dokebirks.pdf") ; 
  
  

}
void exe3() {
  gStyle->SetOptStat(0);
  TCanvas *cn = new TCanvas("cn","cn");
  
  TFile *fin = TFile::Open("eloss.root");  
  gloss = (TGraph*) fin->Get("gr_eloss");

  TFile *f2 = TFile::Open("dataER.root");  
  ptx.clear();
  pty.clear();
  pte.clear();
  ptf.clear();
  char cgr[50], cgram[50], cgrba[50];
  int field ;
  
  TH2F *hh = new TH2F("hh","",100,0,300,100,0.3,1.1);
  hh->GetXaxis()->SetTitle("Energy [keV_{ee}]");
  hh->GetYaxis()->SetTitle("S1 / S1_{0}");
  
  
  hh->Draw();
  TGraph *m1;
  TGraph *m2;
  TGraph *m3;
  
  for(int i =0 ;i<4;++i) {
    int col ;
    if(i == 0 ) {field = 50; col = kGreen-2; }
    if(i == 1 ) {field = 100;col = kOrange -3; }
    if(i == 2 ) {field = 200;col = kRed-3; }
    if(i == 3 ) {field = 500;col = kAzure-1; }
    
    sprintf(cgr,"g%i",field);
    sprintf(cgram,"gam%i",field);
    sprintf(cgrba,"gba%i",field);

    TGraphErrors *gr   = (TGraphErrors*) f2->Get(cgr);
    TGraphErrors *grAm = (TGraphErrors*) f2->Get(cgram);
    TGraphErrors *grBa = (TGraphErrors*) f2->Get(cgrba);
    
    gr->SetMarkerColor(col);
    gr->SetLineColor(col);
    gr->SetMarkerStyle(8);
    gr->SetMarkerSize(1.);
    
    
    grAm->SetMarkerColor(col);
    grAm->SetLineColor(col);
    grAm->SetMarkerStyle(22);
    
    grBa->SetMarkerColor(col);
    grBa->SetLineColor(col);
    grBa->SetMarkerStyle(23);
    
    if(field == 200) {
    gr->Draw("p");
    grAm->Draw("p");
    grBa->Draw("p");
    m1 = new TGraph ; m1->SetMarkerStyle(gr->GetMarkerStyle());
    m2 = new TGraph ; m2->SetMarkerStyle(grAm->GetMarkerStyle());m2->SetMarkerSize(grAm->GetMarkerSize());
    m3 = new TGraph ; m3->SetMarkerStyle(grBa->GetMarkerStyle());m3->SetMarkerSize(grBa->GetMarkerSize());
    
    m1 = gr;
    m2 = grAm;
    m3 = grBa;
    
    //m1->SetMarkerColor(gr->GetMarkerColor());
    //m2->SetMarkerColor(gr->GetMarkerColor());
    //m3->SetMarkerColor(gr->GetMarkerColor());
    

    }
    
    int N = gr->GetN();

    for(int i=0;i<N;++i) {
      ptx.push_back(gr->GetX()[i]);
      pty.push_back(gr->GetY()[i]);
      pte.push_back(gr->GetEY()[i]);
      ptf.push_back(field);
    }
    ptx.push_back(grAm->GetX()[0]);
    pty.push_back(grAm->GetY()[0]);
    pte.push_back(grAm->GetEY()[0]);
    ptf.push_back(field);

    ptx.push_back(grBa->GetX()[0]);
    pty.push_back(grBa->GetY()[0]);
    pte.push_back(grBa->GetEY()[0]);
    ptf.push_back(field);

  }
  //////////////////////////////////////////////////////////////////////////////////////
  // Minimization
  double arglist[3];
  int ierflg=0;
  TMinuit *minLogL = new TMinuit(2);
  minLogL->SetFCN(minFunc2); 
  
  //minLogL->mnexcm("SET PRINT -1", arglist, 1, ierflg);
  minLogL->Command("SET STRATEGY 1");
  
  minLogL->mnparm(0, "A",   1.0,     0.1, 0, 100,  ierflg);  
  minLogL->mnparm(1, "C",   0.5,     0.1, 0,  10,  ierflg);  
  minLogL->mnparm(2, "D",   0.004,     0.001, 0,  20,  ierflg);  
   

  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  TF1 *fun1 = new TF1("fun",fdoke_field,0,500,4);
  TF1 *fun2 = new TF1("fun",fdoke_field,0,500,4);
  TF1 *fun3 = new TF1("fun",fdoke_field,0,500,4);
  TF1 *fun4 = new TF1("fun",fdoke_field,0,500,4);
  fun1->SetLineColor(1);
  fun2->SetLineColor(2);
  fun3->SetLineColor(3);
  fun4->SetLineColor(4);
  double A, C,D,  eA, eC, eD; 
  minLogL->GetParameter(0,A,eA);
  minLogL->GetParameter(1,C,eC);
  minLogL->GetParameter(2,D,eD);
  fun1->SetParameters(A,C,D,50);
  fun2->SetParameters(A,C,D,100);
  fun3->SetParameters(A,C,D,200);
  fun4->SetParameters(A,C,D,500);

  fun1->SetLineColor(kGreen-2);
  fun2->SetLineColor(kOrange -3);
  fun3->SetLineColor(kRed-3);
  fun4->SetLineColor(kAzure-1);
  
  fun3->SetLineColor(kOrange-3);
  fun3->SetLineStyle(1);
  

  //fun1->Draw("same");
  //fun2->Draw("same");
  fun3->Draw("same");
  //fun4->Draw("same");
  
  
  
  TF1 *fP = new TF1("fP",LYat200V,0,350,1);
  fP->SetParameter(0,0.21);
  fP->SetLineColor(kAzure -1);
  fP->Draw("same");
  TLegend * leg = new TLegend (.15, .2, .45,.45) ; 
  leg->SetNColumns(1) ; 
  leg->SetBorderSize(0);
//  leg->SetTextFont(102) ; 
//  leg->SetBorderSize(1) ; 
  leg->SetFillStyle(0) ; 
  leg->AddEntry(m1, "Compton e^{-}", "p") ; 
  leg->AddEntry(m2, "^{241}Am", "p") ; 
  leg->AddEntry(m3, "^{133}Ba", "p") ; 

  //leg->AddEntry(fre, "model 200 V/cm", "l") ; 
  leg->Draw("same") ; 

  TLegend * leg2 = new TLegend (.55, .65, .88,.88) ; 
  leg2->SetNColumns(1) ; 
  leg2->SetBorderSize(0);
  leg2->AddEntry(fP, "PARIS", "l") ; 
  leg2->AddEntry(fun3, "Doke-Birks", "l") ; 
  leg2->Draw("same") ; 

  TPad *pad1 = new TPad("pad1","",0.5,0.15,0.88,0.5);
  pad1->SetFillColor(0);
  pad1->SetLogx(1);
  pad1->Draw();
  pad1->cd();
  TH1 * h2 = pad1->DrawFrame(0.4, 0.3, 320, 1.1) ; 
  //TH2F *h2 = (TH2F*) hh->Clone("h2");
  //h2->GetXaxis()->SetRangeUser(0.1,320);
  h2->GetXaxis()->SetLabelSize(0.06);
  h2->GetYaxis()->SetLabelSize(0.06);
  h2->GetXaxis()->SetTitle("");
  h2->GetYaxis()->SetTitle("");
  
  //h2->Draw("same");
  m1->Draw("p");
  m2->Draw("p");
  m3->Draw("p");
  fun3->Draw("same");
  fP->Draw("same");
  
  
  cn->SaveAs("../figures/paris.pdf") ; 
  
  

}


double fdoke_vs_field(double *x, double*p) {
  double field = x[0];
  double A = p[0];
  double C = p[1];
  double D = p[2];
  double E = p[3];
  C *= Exp(-D*field);
  double B = A/(1-C);
  double dEdx = gloss->Eval(E); 
  double reco = A*dEdx/ (1 + B*dEdx) + C;
  double alpha = 0.21; 
  double Nq    = E/19.5e-3 ; 
  double Ni    = Nq/(1.+alpha) ; 
  double Nx    = Nq-Ni ; 
  double Ng    = Nx + Ni*reco  ; 
  return Ng / Nq  ; 
  

}

void dokebirks_vs_E() {

  TFile *fin = TFile::Open("eloss.root");  
  gloss = (TGraph*) fin->Get("gr_eloss");

  cout <<gloss->Eval(59.5)<<endl;

  //  TGraphErrors * myQ = new TGraphErrors("ReDLYratio.txt");

  TGraphErrors * myQ = new TGraphErrors("REDLY_SP.txt");

  TGraphErrors * myQn = new TGraphErrors("ReDLYratio_nom.txt");
  //  TGraphErrors * myQ = new TGraphErrors("ReDLYratio_Ecorr.txt");
  //  TGraphErrors * myQn = new TGraphErrors("ReDLYratio_nom_Ecorr.txt");
  TGraphErrors * myQr = new TGraphErrors("ReDLYratio_Ring.txt");
  // ARIS 
  TFile *faris = TFile::Open("dataER.root");  
  TGraphErrors *grAm;
  TGraphErrors *grAmP = new TGraphErrors(4);
  int E[4]={50,100,200,500};
  Double_t x;
  Double_t y;
  Double_t ex,ey;
  for (int i=0;i<4;i++) {

     if (i==0) grAm = (TGraphErrors*) faris->Get("gam50");
     if (i==1) grAm = (TGraphErrors*) faris->Get("gam100");
     if (i==2) grAm = (TGraphErrors*) faris->Get("gam200");
     if (i==3) grAm = (TGraphErrors*) faris->Get("gam500");

     grAm->GetPoint(0,x,y);ey = grAm->GetErrorY(0);
     printf ("%d %f %f %f\n",E[i],x,y, ey);
     grAmP->SetPoint(i,(Double_t) E[i],y);  grAmP->SetPointError(i,0.1,ey);

  }
  
  TF1 *fun1 = new TF1("funAris",fdoke_vs_field,0,900,4);
  TF1 *fun12 = new TF1("funReD",fdoke_vs_field,0,900,4);
  TF1 *fun2 = new TF1("fun",fdoke_field,0,500,4);

  fun2->SetLineColor(2);

  double A, C,D,  eA, eC, eD; 
  A=2.49221;
  C=7.72418e-01;
  D=3.50517e-03;
  double ene=59.6;//keV
  fun1->SetParameters(A,C,D,ene);
  fun2->SetParameters(A,C,D,200.);


  //set of "fit" with guessed fields
  //  A=2.7;
  //C=8.02418e-01;
  //D=3.50517e-03;
  //set of "fit" to ARIS+RED
  A=2.7;
  C=7.82418e-01;
  D=3.40517e-03;
  //set of "fit" with corrected fields
  //A=2.7;
  //C=8.06418e-01;
  //D=3.10517e-03;
  fun12->SetParameters(A,C,D,ene);

  TCanvas * c1 = new TCanvas();
  fun1->SetLineColor(4);
  fun1->Draw();
  fun1->GetXaxis()->SetTitle("E Field [V/cm]");
  fun1->GetYaxis()->SetTitle("S1/S1^{0}");
  fun12->SetLineColor(4);  fun12->SetLineStyle(2);fun12->Draw("same");
  //ReD
  myQ->SetMarkerStyle(20);
  myQ->SetName("myQ");
  myQn->SetMarkerStyle(20);
  myQn->SetMarkerColor(2);
  myQn->SetName("myQn");
  myQr->SetMarkerStyle(25);
  myQr->SetMarkerColor(2);
  myQr->SetMarkerSize(1.5);
  myQr->SetName("myQr");
  myQ->Draw("pesame");  myQn->Draw("pesame"); myQr->Draw("pesame");
  // Aris
  grAmP->SetMarkerStyle(20);
  grAmP->SetMarkerColor(3);
  grAmP->SetName("Aris");
  grAmP->Draw("pesame"); 

  TLegend * l = new TLegend(0.4,0.5,0.9,0.8);
  l->AddEntry("myQn","ReD nominal Field configuration","P");
  l->AddEntry("myQ","ReD with Anode=0 Ring=0 K=-815,1000,2000 V","P");
  l->AddEntry("myQr","ReD with Anode=0 R=85 K=-815 V","P");
  l->AddEntry("Aris","Aris ","P");
  l->AddEntry("funAris","Aris Fit to Doke-Birks  A=2.5;B=7.72;D=3.51","L");
  l->AddEntry("funReD","ReD+Aris Fit to Doke-Birks A=2.7;B=7.82;D=3.41","L");
  l->SetBorderSize(0);
  l->Draw();
  TCanvas * c2 = new TCanvas();
  fun2->Draw();


}

void myexe() {

  gStyle->SetOptStat(0);
  TCanvas *cn = new TCanvas("cn","cn");
  
  TFile *fin = TFile::Open("eloss.root");  
  gloss = (TGraph*) fin->Get("gr_eloss");

  TGraphErrors * grAm = new TGraphErrors("REDLY_SP.txt");

  ptx.clear();
  pty.clear();
  pte.clear();
  ptf.clear();
  char cgr[50], cgram[50], cgrba[50];
  int field ;

  grAm->SetMarkerColor(4);
  grAm->SetLineColor(4);
  grAm->SetMarkerStyle(22);

  int N = grAm->GetN();

  double redfield[9]={0,50.,100.,200.,300.,400.,500.,700.,1000.};

  for(int i=1;i<N;++i) {
    cout <<"59.5 "<< grAm->GetY()[i]<<
    " "<<grAm->GetEY()[i]<<" "<<redfield[i]<<endl;

    ptx.push_back(59.5);
    pty.push_back(grAm->GetY()[i]);
    pte.push_back(grAm->GetEY()[i]);
    ptf.push_back(redfield[i]);
  }

  double arglist[3];
  int ierflg=0;
  TMinuit *minLogL = new TMinuit(2);
  minLogL->SetFCN(minFunc2); 
  
  //minLogL->mnexcm("SET PRINT -1", arglist, 1, ierflg);
  minLogL->Command("SET STRATEGY 1");
  
  minLogL->mnparm(0, "A",   2.5,     0.1, 1, 10,  ierflg);  
  minLogL->mnparm(1, "C",   0.77,     0.1, 0,  10,  ierflg);  
  minLogL->mnparm(2, "D",   0.004,     0.001, 0,  20,  ierflg);  
   

  minLogL->mnfixp(0,ierflg);
  minLogL->mnfixp(0,ierflg);
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  minLogL->mnfree(1);
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migra
  minLogL->mnfree(0);
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  minLogL->mnexcm("MIGRAD",arglist,1,ierflg);//minimization with Migrad
  TF1 *fun1 = new TF1("funAris",fdoke_vs_field,0,1100,4);
  fun1->SetLineColor(6);
  double A, C,D,  eA, eC, eD; 
  minLogL->GetParameter(0,A,eA);
  minLogL->GetParameter(1,C,eC);
  minLogL->GetParameter(2,D,eD);
  fun1->SetParameters(A,C,D,59.5);
  
  grAm->Draw("AP");
  fun1->Draw("same");

}


void ReDQuenching_vs_E() {

  TFile *fin = TFile::Open("eloss.root");  
  gloss = (TGraph*) fin->Get("gr_eloss");

  cout <<gloss->Eval(59.5)<<endl;

  //  TGraphErrors * myQ = new TGraphErrors("ReDLYratio.txt");

  TGraphErrors * myQ = new TGraphErrors("REDLY_SP.txt");
  TGraphErrors * myQn = new TGraphErrors("REDLY_DP.txt");
  //  TGraphErrors * myQr = new TGraphErrors("ReDLYratio_Ring.txt");
  // ARIS 
  TFile *faris = TFile::Open("dataER.root");  
  TGraphErrors *grAm;
  TGraphErrors *grAmP = new TGraphErrors(4);
  int E[4]={50,100,200,500};
  Double_t x;
  Double_t y;
  Double_t ex,ey;
  for (int i=0;i<4;i++) {

     if (i==0) grAm = (TGraphErrors*) faris->Get("gam50");
     if (i==1) grAm = (TGraphErrors*) faris->Get("gam100");
     if (i==2) grAm = (TGraphErrors*) faris->Get("gam200");
     if (i==3) grAm = (TGraphErrors*) faris->Get("gam500");

     grAm->GetPoint(0,x,y);ey = grAm->GetErrorY(0);
     printf ("%d %f %f %f\n",E[i],x,y, ey);
     grAmP->SetPoint(i,(Double_t) E[i],y);  grAmP->SetPointError(i,0.1,ey);

  }
  
  TF1 *fun1 = new TF1("funAris",fdoke_vs_field,0,1100,4);
  TF1 *fun12 = new TF1("funReD",fdoke_vs_field,0,1100,4);
  TF1 *fun2 = new TF1("fun",fdoke_field,0,500,4);

  fun2->SetLineColor(2);

  double A, C,D,  eA, eC, eD; 
  A=2.49221;
  C=7.72418e-01;
  D=3.50517e-03;
  double ene=59.6;//keV
  fun1->SetParameters(A,C,D,ene);
  fun2->SetParameters(A,C,D,200.);


  //set of "fit" from ReD SP data
  //   1  A            1.51970e+00   8.26936e-03  -3.88991e-07  -5.50997e-04
  //   2  C            7.79286e-01   8.99723e-04  -1.93777e-08  -8.77658e-03
  //   3  D            2.03965e-03   1.19728e-05  -6.39571e-09  -1.49955e-01

  A=1.52;
  C=7.79e-01;
  D=2.04e-03;

  fun12->SetParameters(A,C,D,ene);

  TCanvas * c1 = new TCanvas();
  fun1->SetLineColor(4);
  fun1->Draw();
  fun1->GetXaxis()->SetTitle("E Field [V/cm]");
  fun1->GetYaxis()->SetTitle("S1/S1^{0}");
  fun1->GetYaxis()->SetRangeUser(0.6,0.97);
  fun12->SetLineColor(4);  fun12->SetLineStyle(2);fun12->Draw("same");
  //ReD
  myQ->SetMarkerStyle(20);
  myQ->SetName("myQ");
  myQn->SetMarkerStyle(24);
  myQn->SetLineWidth(2);
  myQn->SetMarkerColor(2);
  myQn->SetName("myQn");
  //  myQr->SetMarkerStyle(25);
  //  myQr->SetMarkerColor(2);
  //  myQr->SetMarkerSize(1.5);
  //  myQr->SetName("myQr");
  myQ->Draw("pesame");  myQn->Draw("pesame"); 
  //myQr->Draw("pesame");
  // Aris
  grAmP->SetMarkerStyle(23);
  grAmP->SetMarkerColor(3);
  grAmP->SetName("Aris");
  grAmP->Draw("pesame"); 

  TLegend * l = new TLegend(0.4,0.5,0.9,0.8);
  l->AddEntry("myQn","ReD Double Phase nominal Field config. A=5211V","P");
  l->AddEntry("myQ","ReD Single Phase config.","P");
  //  l->AddEntry("myQr","ReD with Anode=0 R=85 K=-815 V","P");
  l->AddEntry("Aris","Aris ","P");
  l->AddEntry("funAris","Aris Fit to Doke-Birks  A=2.5;B=7.72;D=3.51","L");
  l->AddEntry("funReD","ReD SP data Fit to Doke-Birks A=1.52;B=7.79;D=2.04","L");
  l->SetBorderSize(0);
  l->Draw();


  //  TCanvas * c2 = new TCanvas();
  //  fun2->Draw();


}


void dokebirks() {

  grA = new TGraphErrors ;
  grC = new TGraphErrors ;
  TCanvas *c1 = new TCanvas("c1","A");

  exe(50);
  exe(100);
  exe(200);
  exe(500);
  
  grA->Draw("ap");
  
  
  TCanvas *c2 = new TCanvas("c2","C");
  grC->Draw("ap");

}

Double_t Qy_Aris(double *x, double *p)
{
  double field = x[0];
  double A = p[0];
  double C = p[1];
  double D = p[2];
  double E = p[3];
  C *= Exp(-D*field);
  double B = A/(1-C);
  double dEdx = gloss->Eval(E); 
  double reco = A*dEdx/ (1 + B*dEdx) + C;
  //  cout <<"dEdx "<<dEdx<<" reco "<< reco<<endl;
  double alpha = 0.21; 
  double Nq    = E/19.5e-3 ; 
  double Ni    = Nq/(1.+alpha) ; 
  double Nx    = Nq-Ni ; 
  double Ng    = Nx + Ni*reco  ; 

  //Qy

  double Qy= 1 * Ni * (1-reco) * 1/E; 
  return Qy  ; 
  //  return Ng / Nq  ; 

}

Double_t e_Qy_Aris(double *x, double *p)
{
  double field = x[0];
  double A = p[0];
  double C = p[1];
  double D = p[2];
  double E = p[3];
  C *= Exp(-D*field);
  double B = A/(1-C);
  double dEdx = gloss->Eval(E); 
  double reco = A*dEdx/ (1 + B*dEdx) + C;
  //  cout <<"dEdx "<<dEdx<<" reco "<< reco<<endl;
  double alpha = 0.21; 
  double Nq    = E/19.5e-3 ; 
  double Ni    = Nq/(1.+alpha) ; 
  double Nx    = Nq-Ni ; 
  double Ng    = Nx + Ni*reco  ; 

  double DeltaA = p[4];
  double DeltaC = p[5];
  double DeltaD = p[6];

  //  Delta
  DeltaReco = pow( (dEdx/ (1 + B*dEdx) - A*dEdx/pow((1 + B*dEdx),2)*1/(1-C)) * DeltaA,2.);
  DeltaReco += pow( ( -1.*A*dEdx/pow((1 + B*dEdx),2)*A/pow((1-C),2.)+ 1)*Exp(-D*field) * DeltaC,2.);
  DeltaReco +=  pow( ( +1.*A*dEdx/pow((1 + B*dEdx),2)*A/pow((1-C),2.) - 1.) * p[1]*field*Exp(-D*field) * DeltaD,2);
  DeltaReco = sqrt(DeltaReco);

  double eQy = 1 * Ni * DeltaReco * 1/E; 
  return eQy  ; 

  //  return Ng / Nq  ; 

}

Double_t Qy_Aris2(double *x, double *p)
{
  double field = x[0];
  double A = p[0];
  double C = p[1];
  double D = p[2];
  double E1 = p[3];
  double E2 = p[4];
  C *= Exp(-D*field);
  double B = A/(1-C);
  double dEdx1 = gloss->Eval(E1); 
  double dEdx2 = gloss->Eval(E2); 
  double reco1 = A*dEdx1/ (1 + B*dEdx1) + C;
  double reco2 = A*dEdx2/ (1 + B*dEdx2) + C;
  //  cout <<"dEdx1 "<<dEdx1<<" reco1 "<< reco1<<"dEdx2 "<<dEdx2<<" reco2 "<< reco2<<endl;
  double alpha = 0.21; 
  double Nq1    = E1/19.5e-3 ; 
  double Nq2    = E2/19.5e-3 ; 
  double Ni1    = Nq1/(1.+alpha) ; 
  double Ni2    = Nq2/(1.+alpha) ; 
  double Nx    = Nq1-Ni1 ; 
  double Ng    = Nx + Ni1*reco1  ; 

  //Qy

  double Qy = 1*( Ni1 * (1-reco1) + Ni2 * (1-reco2)) * 1/(E1+E2); 
  return Qy; 
  //  return Ng / Nq  ; 

}

Double_t TI(Double_t *x, Double_t *par)
{
   Float_t xx =x[0];

   Double_t Ni=1760; // e-
   Double_t g2=1.0;  //to express Qy in e-/keV
   Double_t C=par[0]; //Scene 0.18+-0.03 V/cm^B/e-
   Double_t B=par[1]; //Scene 0.61+-0.03
   Double_t E=par[2]; //Scene Kr 41.5 keV, ReD Am 59.5 keV

   Double_t csi = Ni*C/TMath::Power(xx,B);
   Double_t f = g2*Ni/E/csi * TMath::Log(1+csi);

   return f;
}

void Plot (Double_t gasbelow=0.,  Bool_t shrinkage=kFALSE)
//(Double_t gasbelow=0., Double_t leakage=0., Bool_t shrinkage=kFALSE, Bool_t S2diffusion=kFALSE)
{
  Double_t EfieldScene[] = {49.5,96.4,193,293,486};
  Double_t e_EfieldScene[] = {0.1,0.1,0.1,0.1,0.1};
  Double_t QyScene[] = {5.0,6.5,8.6,10.2,12.6}; //Scene fig. 15
  Double_t e_QyScene[] = {0.5,0.5,0.5,0.5,0.5}; //Scene fig. 15

  Int_t nFieldScene = sizeof(EfieldScene)/sizeof(Double_t);

  //ReD

  Double_t EfieldReD[] = {100.,200.,400.,700.,1000.};
  Double_t e_EfieldReD[] = {0.1,0.1,0.1,0.1,0.1};
  Double_t QyReDold[] = {5.5,10.5,19.0,27.2,33.5}; //Giuseppe PE/KeV
  Double_t e_QyReDold[] = {0.5,0.5,0.5,0.5,0.5}; //Giuseppe
  Double_t S2ReD[] = {0,6.39,10.0,16.8,23.8,29.39}; //Simone PE/KeV
  Double_t e_S2ReD[] = {0.04,0.04,0.04,0.04,0.04,0.04}; //Simone
  Double_t QyReD[] = {6.39,10.0,16.8,23.8,29.39}; //Simone PE/KeV
  Double_t e_QyReD[] = {0.04,0.04,0.04,0.04,0.04}; //Simone
  Double_t g2=1.79;

  Double_t S1LY[] = {8.691639801, 8.1447211, 7.685204656, 6.864880158,6.053781076,5.453038048};
  Double_t e_S1LY[] = {0.0349797, 0.0357045, 0.0360498, 0.035404, 0.0359978, 0.0367451};
  Double_t g2Simo=1.79; // starting value

  Double_t KReD[] = {310.,815.,1820.,3330,4870};

  Int_t nFieldReD = sizeof(EfieldReD)/sizeof(Double_t);

  //Determine g2Simo
  TCanvas * c2 = new TCanvas();
  TGraphErrors * S1vsS2 = new TGraphErrors(nFieldReD+1,S2ReD,S1LY,e_S2ReD,e_S1LY);
  S1vsS2->SetMarkerStyle(20); S1vsS2->SetMarkerSize(0.6);
  S1vsS2->Draw("APE");
  S1vsS2->GetHistogram()->GetXaxis()->SetRangeUser(-1,35);
  S1vsS2->GetHistogram()->SetTitle("Drift Field Scan 0,100,200,400,700,1000 V");
  S1vsS2->GetHistogram()->GetXaxis()->SetTitle("S2_{corr} MPV [PE/keV]");
  S1vsS2->GetHistogram()->GetYaxis()->SetTitle("S1_{corr} MPV [PE/keV]");
  
  TF1 * lin = new TF1("lin","[0]/0.0195 - [0]/[1]*x",0,50);
  lin->SetParameters(0.2,g2Simo);
  lin->SetParNames("g1  [PE/ph]","g2 [PE/e^-]");
  S1vsS2->Fit("lin","","",4.,35.);
  TF1 * lin2 = new TF1(*lin);
  lin2->SetLineStyle(2);lin2->Draw("same");
  //  g2Simo=abs(lin->GetParameter(0)*0.0236/lin->GetParameter(1));//0.0236 = 23.6 eV Wph
  //  g2Simo=abs(lin->GetParameter(0)*0.0195/lin->GetParameter(1));//0.0236 = 23.6 eV Wph

  g2Simo=lin->GetParameter(1);//0.0236 = 23.6 eV Wph

  cout <<"g2 from linear fit "<<g2Simo<<" +/- "<< lin->GetParError(1)<<endl; 

  TCanvas * c1 = new TCanvas();

  Double_t NominalDrift=5.0; //cm
  Double_t PTFEcte = 86E-6;
  Double_t Drift=NominalDrift;
  if (shrinkage) {
    Drift *= (1+PTFEcte * (88.2-273.15-25.));
    cout <<"Nominal Drift: "<<NominalDrift<<" with shrinkage: "<<Drift<<endl;
  }

  for (int i=0;i<nFieldReD;i++) {
    QyReDold[i]=QyReDold[i]/g2;
    e_QyReDold[i]=e_QyReDold[i]/g2;
    QyReD[i]=QyReD[i]/g2Simo;
    e_QyReD[i]=e_QyReD[i]/g2Simo;
    // QyReD[i]=QyReD[i]/2.2;
    // e_QyReD[i]=e_QyReD[i]/2.2;
    KReD[i]/=(1.505*gasbelow/10.+(Drift-gasbelow/10.));
  }


  TGraphErrors * QY_S = new TGraphErrors(nFieldScene,EfieldScene,QyScene,e_EfieldScene,e_QyScene);
  QY_S->SetLineColor(4);QY_S->SetLineWidth(2); QY_S->SetMarkerStyle(20);

  TGraphErrors * QY_R = new TGraphErrors(nFieldReD,EfieldReD,QyReDold,e_EfieldReD,e_QyReDold);
  QY_R->SetLineColor(2);QY_R->SetLineWidth(2); QY_R->SetMarkerStyle(22);QY_R->SetMarkerColor(2);

  TGraphErrors * QY_RSimo = new TGraphErrors(nFieldReD,EfieldReD,QyReD,e_EfieldReD,e_QyReD);
  QY_RSimo->SetLineColor(2);QY_RSimo->SetLineWidth(2); QY_RSimo->SetMarkerStyle(22);QY_RSimo->SetMarkerColor(2);QY_RSimo->SetMarkerSize(1.5);

  TGraphErrors * QY_RSimo_m = new TGraphErrors(nFieldReD,KReD,QyReD,e_EfieldReD,e_QyReD);
  QY_RSimo_m->SetLineColor(6);QY_RSimo_m->SetLineWidth(2); QY_RSimo_m->SetMarkerStyle(23);


  TF1 * tiS = new TF1 ("tiS",TI,0,1200,3);
  TF1 * tiR = new TF1 ("tiR",TI,0,1200,3);
  tiS->SetParameters(0.18,0.61,41.5);
  tiR->SetParameters(0.18,0.61,59.5);
 
  tiS->SetLineColor(4);
  tiR->SetLineColor(2);

  tiS->GetHistogram()->GetXaxis()->SetTitle("E_{drift} [V/cm]");
  tiS->GetHistogram()->GetYaxis()->SetTitle("Qy [e-/keV]");
  tiS->GetHistogram()->SetTitle("Charge Yield vs Field");
  tiS->GetHistogram()->SetMinimum(0);
  tiS->GetHistogram()->SetMaximum(22);

  
  tiS->Draw("");
  //  tiR->Draw("same");
  
  //include theory from Doke-Birks and Aris fit
  TFile *fin = TFile::Open("eloss.root");  
  gloss = (TGraph*) fin->Get("gr_eloss");

  //  cout <<gloss->Eval(59.5)<<endl;

  TF1 *fun1 = new TF1("funAris",fdoke_vs_field,0,1200,4);
  TF1 *QyAris = new TF1("QyAris",Qy_Aris,0,1200,4);
  TF1 *e_QyAris = new TF1("e_QyAris",e_Qy_Aris,0,1200,7);
  TF1 *QyArisKr = new TF1("QyArisKr",Qy_Aris2,0,1200,5);
  //  TF1 *fun12 = new TF1("funReD",fdoke_vs_field,0,900,4);
  //  TF1 *fun2 = new TF1("fun",fdoke_field,0,500,4);

  //  fun2->SetLineColor(2);

  double A, C,D,  eA, eC, eD; 
  A=2.49221;   //0.2
  C=7.72418e-01; // 0.01
  D=3.50517e-03; //0.3
  //set of "fit" with guessed fields from ReD in September
  //A=2.7;
  //C=7.82418e-01;
  //D=3.40517e-03;

  double ene=59.6;//keV
  //fun1->SetParameters(A,C,D,ene);
  QyAris->SetParameters(A,C,D,ene);
  QyAris->SetParError(0,0.2);
  QyAris->SetParError(1,0.01);
  QyAris->SetParError(2,0.3E-3);
  QyAris->SetParError(3,1.);
  
  
  QyAris->SetLineColor(2);QyAris->SetLineStyle(2);QyAris->Draw("same");

  //Set Error band
  e_QyAris->SetParameters(A,C,D,ene,0.2,0.01,0.3E-3);
  // QyAris->SetParError(0,0.2);
  // QyAris->SetParError(1,0.01);
  // QyAris->SetParError(2,0.3E-3);
  // QyAris->SetParError(3,1.);

  double x[100];
  double y[100];
  double ex[100];
  double ey[100];

  for (int i=0;i<100;i++) {
    x[i]=i*15;
    ex[i]=0.1;
    y[i]=QyAris->Eval(x[i]);
    ey[i]=e_QyAris->Eval(x[i]);
  }
  TGraphErrors * ArisBand = new TGraphErrors(100,x,y,ex,ey);
  ArisBand->SetFillColor(46);   ArisBand->SetLineColor(6); ArisBand->SetLineStyle(2); 
  ArisBand->Draw("3Lsame");QyAris->Draw("same");
  

  QyArisKr->SetParameters(A,C,D,32.5,9.0);
  QyArisKr->SetLineColor(4);QyArisKr->SetLineStyle(2);QyArisKr->Draw("same");


  //Data Points

  QY_S->Draw("Psame");
  //  QY_R->Draw("Psame");
  QY_RSimo->Draw("Psame");
  //QY_RSimo_m->Draw("Psame");

  TLegend * l = new TLegend(0.42,0.12,0.85,0.42);
  l->AddEntry(QY_S,"Qy from SCENE (kr83)","PE");
  l->AddEntry(tiS,"Mod. Thomas Imel fit from SCENE on Kr83","L");
  //  l->AddEntry(QY_R,"Qy from Giuseppe","PE");
  //  l->AddEntry(QY_RSimo,"Qy from Simone","PE");
  l->AddEntry(QY_RSimo,"Qy from ReD","PE");
  //  l->AddEntry(QY_RSimo_m,"Qy from Simonee (Mauro's field)","PE");
  // l->AddEntry(tiR,"Mod. Thomas Imel fit from SCENE scaled to Am241","L");
  l->AddEntry(ArisBand,"ARIS Doke-Birks parametrization","LF");
    //  l->AddEntry(QyAris,"ARIS Doke-Birks parametrization","L");
  l->AddEntry(QyArisKr,"ARIS Doke-Birks parametrization (Kr83)","L");
  l->SetBorderSize(0);
  l->Draw();
   

}
