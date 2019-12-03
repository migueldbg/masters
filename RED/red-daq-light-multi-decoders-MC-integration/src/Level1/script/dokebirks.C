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

  TGraphErrors * myQ = new TGraphErrors("ReDLYratio.txt");
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
