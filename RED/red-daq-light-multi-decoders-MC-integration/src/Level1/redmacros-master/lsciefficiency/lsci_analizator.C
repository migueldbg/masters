#include <iostream>
#include <cstdlib>
#include <cstring>
#include <algorithm> 
#include <vector> 
#include <sstream> 

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
#include <TProfile.h>
#include <TLatex.h>

#include "EvRec0.hh"
#include "EvRaw0.hh"

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "TGaxis.h"
#include "RooClassFactory.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooCompton.h"

//
//+++COMANDI+++
//
// .L RooCompton.cxx+ per caricare la libreria "RooCompton.h"
// .L lsci_analizator.C
// lsci_analizator(run,run1,run2,run3,ch,pmt)

using namespace std;

Double_t Q(Double_t *x, Double_t *par)
{
    
    Double_t kB=par[0];
    Double_t GrauMalonda[21][7]={
        0.000,1.0,0.0,0.0,0.0,0.0,0.0,
        0.001,0.84859,0.38579,0.15528,-5.91812*1E-3,0.38996,0.15143,
        0.002,0.73916,0.33815,0.14114,-4.12402*1E-3,0.34586,0.13451,
        0.003,0.65585,0.3043,0.13198,-4.75264*1E-3,0.31382,0.12344,
        0.004,0.58917,0.26473,0.12173,-2.54649*1E-3,0.27388,0.11198,
        0.005,0.53551,0.2397,0.11401,-1.95557*1E-3,0.24956,0.10317,
        0.006,0.49041,0.216,0.10725,-1.46804*1E-3,0.22468,0.09544,
        0.007,0.45355,0.20231,0.10137,-0.792*1E-3,0.21445,0.0888,
        0.008,0.42149,0.18731,0.09624,-0.52346*1E-3,0.19977,0.08311,
        0.009,0.39439,0.17372,0.09129,-0.57124*1E-3,0.18524,0.07784,
        0.01,0.36851,0.16398,0.08816,-0.45562*1E-3,0.17582,0.07431,
        0.011,0.34793,0.15193,0.08427,-0.05526*1E-3,0.16443,0.07002,
        0.012,0.32903,0.14404,0.08059,-0.04536*1E-3,0.15623,0.06611,
        0.013,0.31231,0.13664,0.07744,-0.04306*1E-3,0.14961,0.06284,
        0.014,0.29668,0.12872,0.07477,0.15992*1E-3,0.14091,0.05978,
        0.015,0.28281,0.12259,0.07235,0.30371*1E-3,0.13539,0.05725,
        0.016,0.2702,0.11663,0.06985,0.34938*1E-3,0.12933,0.05453,
        0.017,0.25863,0.11133,0.06764,0.41188*1E-3,0.12414,0.05218,
        0.018,0.24808,0.10646,0.06576,0.46844*1E-3,0.11938,0.04998,
        0.019,0.23832,0.10192,0.06363,0.51606*1E-3,0.11493,0.04793,
        0.02,0.22933,0.09784,0.06184,0.56027*1E-3,0.11096,0.04604
    };
    
    Double_t dkb=1000.,diff;
    static int ikb=-1;
    //  if (ikb<0) {
    for (int k=0;k<21;k++) {
        diff=fabs(kB-GrauMalonda[k][0]);
        if (diff<dkb) {dkb=diff;ikb=k;}
        //      printf ("%f %d\n ",GrauMalonda[k][0],ikb);
    }
    //  }
    Double_t logE=log(x[0]);
    
    Double_t num = GrauMalonda[ikb][1] + GrauMalonda[ikb][2]*logE +
    GrauMalonda[ikb][3]*pow(logE,2) + GrauMalonda[ikb][4]*pow(logE,3);
    
    Double_t den = 1 + GrauMalonda[ikb][5]*logE +
    GrauMalonda[ikb][6]*pow(logE,2) + GrauMalonda[ikb][4]*pow(logE,3);
    
    if (den!=0) return num/den;
    else return 0;
}

Double_t CQ(Double_t *x, Double_t *par)
{
    Double_t xx[1];
    xx[0]=par[1]+par[2]*x[0];
    Double_t ret;
    ret = par[1]+(par[2]*x[0])/Q(xx,par);
    //  ret = Q(xx,par);
    return ret;
}

void lsci_analizator(int run, int run1, int run2, int run3, int nch, int pmt, int amin=5200, int amax=10000, int csmin=45000, int csmax=75000, int namin1=35000, int namax1=52000, int namin2=138000, int namax2=160000){
    
    double mu_am = 0.;
    double mu_am_err = 0.;
    double A_am = 0.;
    double A_cs = 0.;
    double A_cs_err = 0.;
    double A_na1 = 0.;
    double A_na2 = 0.;
    double A_na_err1 = 0.;
    double A_na_err2 = 0.;
    double sig_am = 0.;
    double sig_cs = 0.;
    double sig_na1 = 0.;
    double sig_na2 = 0.;
    double sig_am_err = 0.;
    double sig_cs_err = 0.;
    double sig_na_err1 = 0.;
    double sig_na_err2 = 0.;
    
    TFile *f = new TFile(Form("../../rootfiles/lsci/run_%d.root", run), "read"); //bkg
    TFile *g = new TFile(Form("../../rootfiles/lsci/run_%d.root", run1), "read"); //sorgente1
    TFile *l = new TFile(Form("../../rootfiles/lsci/run_%d.root", run2), "read"); //sorgente2
    TFile *m = new TFile(Form("../../rootfiles/lsci/run_%d.root", run3), "read"); //sorgente3

    TTree *reco = (TTree*)f->Get("reco");
    TTree *reco_1 = (TTree*)g->Get("reco");
    TTree *reco_2 = (TTree*)l->Get("reco");
    TTree *reco_3 = (TTree*)m->Get("reco");

//Aliases
//    reco->SetAlias("s1", "clusters[0].charge");
    
//Cuts
    TCut time = Form("(start_time[%d]-500)<50",nch); //da usare nel caso delle singole!
    //TCut f90 = Form("f90[%d]<0.15",nch);
    
//Styles
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(2210);
    gStyle->SetOptFit(111);
    
//Titles
    TString title = Form("Background Spectrum on PMT %d",pmt);
    TString title_1 = Form("241Am Spectrum on PMT %d",pmt);
    TString title_2 = Form("137Cs Spectrum on PMT %d",pmt);
    TString title_3 = Form("22Na Spectrum on PMT %d",pmt);
    TString title_4 = Form("Energy Calibration Curve for PMT %d",pmt);
    TString title_5 = Form("Resolution Curve for PMT %d",pmt);
    TString title_6 = Form("Residuals for PMT %d",pmt);
    
//Canvases
    TCanvas *c1 = new TCanvas("c1","c1",10,10,900,600);
    TCanvas *c2 = new TCanvas("c2","c2",10,10,900,600);
    TCanvas *c3 = new TCanvas("c3","c3",10,10,900,600);
    TCanvas *c4 = new TCanvas("c4","c4",10,10,900,600);
    TCanvas *c5 = new TCanvas("c5","c5",10,10,900,600);
    TCanvas *c6 = new TCanvas("c6","c6",10,10,900,600);
    TCanvas *c7 = new TCanvas("c7","c7",10,10,900,600);
    //TCanvas *c8 = new TCanvas("c8","c8",10,10,900,600);
    
/*
//Legends
    TLegend *l1 = new TLegend(0.6648107,0.6869565,0.8864143,0.8434783,NULL,"brNDC");
    l1->SetBorderSize(0);

    TLegend *l2 = new TLegend(0.6648107,0.6869565,0.8864143,0.8434783,NULL,"brNDC");
    l2->SetBorderSize(0);
*/

////////////// BKG Spectrum ////////////////////
    c1->cd()->SetLogy();
    //gStyle->SetOptStat(0);
    
    TH1F *h = new TH1F("h",title,300,0,300000);
    
    reco->Draw(Form("charge[%d]>>h",nch),time);
    
    h->GetXaxis()->SetTitle("Energy (ADC*Sample)");
    h->GetXaxis()->SetTitleOffset(1.16);
    h->GetYaxis()->SetTitle("Events");
    h->GetYaxis()->SetTitleOffset(1.16);
    h->SetLineWidth(1);
    h->SetLineColor(kBlue);
    
    h->Draw();

////////////// 241Am Spectrum ////////////////////
    c2->cd();//->SetLogy();
    gStyle->SetOptFit(111);
    
    TH1F *h1 = new TH1F("h1",title_1,300,200,50000);
    
    reco_1->Draw(Form("charge[%d]>>h1",nch),time);

    h1->GetXaxis()->SetTitle("Energy (ADC*Sample)");
    h1->GetXaxis()->SetTitleOffset(1.16);
    h1->GetYaxis()->SetTitle("Counts (arb.)");
    h1->GetYaxis()->SetTitleOffset(1.16);
    h1->SetLineWidth(1);
    h1->SetLineColor(kBlack);
    
    //h1->Add(h,-1);

    //h1->Draw();

////////////// Fitting 241Am Spectrum using RooFit ////////////////////
    RooRealVar x("x","x", 200.,50000.); //range dell'histo cambiare in caso di modello con bkg
    
    RooDataHist am("am","am",x,RooFit::Import(*h1));
//    RooDataHist am_bkg("am_bkg","am_bkg",x,RooFit::Import(*h)); //background from histo
    
    RooRealVar mam("mam","mam",7000.,200,50000);
    RooRealVar sam("sam","sam",1300., 0., 100000.);
    RooGaussian gam("gam","gam",x,mam,sam);
/*
    RooHistPdf bkg_am("bkg_am","pdf from hist",x,am_bkg,0); //background from histo
    RooRealVar bkg_frac_am("bkg_frac_am","fraction of bkg",0.1,0.,1.);
    RooAddPdf model_am("model_am","gam+bkg",RooArgList(gam,bkg_am),bkg_frac_am);
*/
    x.setRange("R1",amin,amax);
//    x.setRange("R2",amax,180000);
    
    RooPlot* fr = x.frame(RooFit::Title(title_1));
    am.plotOn(fr);
    
    gam.fitTo(am,RooFit::Range("R1"));
    gam.plotOn(fr,RooFit::Components(gam),RooFit::LineColor(kRed));
/*
    model_am.fitTo(am,RooFit::Range("R2"));
    model_am.plotOn(frame,RooFit::LineColor(kRed)); //somma Gauss + bkg
    model_am.plotOn(fr,RooFit::Components(bkg_am),RooFit::LineColor(kViolet));
*/
    //mam.Print();
    //sam.Print();
    
    mu_am = mam.getValV();
    mu_am_err = mam.getError();
    sig_am = sam.getValV();
    sig_am_err = sam.getError();
    
    fr->SetXTitle ("Energy (ADC*Sample)");
    fr->SetYTitle ("Events");
    fr->GetYaxis()->SetTitleOffset(1.2);
    fr->Draw();

////////////// 137Cs Spectrum ////////////////////
    c3->cd();//->SetLogy();
    gStyle->SetOptFit(111);
    
    TH1F *h3 = new TH1F("h3",title_2,300,0,300000); //ADC*Sample
    //TH1F *h3 = new TH1F("h3",title_2,1000,0,2000); //in keV calibrata sul 22Na
    
    reco_2->Draw(Form("charge[%d]>>h3",nch),time); //ADC*Sample
    //reco_2->Draw(Form("charge[%d]*0.00693231>>h3",nch),time); //in keV calibrata sul 22Na
    
    h3->GetXaxis()->SetTitle("Energy (ADC*Sample)");
    //h3->GetXaxis()->SetTitle("Energy (keV)");
    h3->GetXaxis()->SetTitleOffset(1.16);
    h3->GetYaxis()->SetTitle("Counts (arb.)");
    h3->GetYaxis()->SetTitleOffset(1.16);
    h3->SetLineWidth(1);
    h3->SetLineColor(kBlack);
    
    //h3->Add(h,-1);
    
    //h3->Draw();

////////////// Fitting 137Cs Spectrum using RooFit ////////////////////
    RooRealVar t("t","t", 500.,200000.); //range dell'histo
    RooRealVar hv("hv","hv", 677000); //energia picco fotoelettrico in eV
    
    RooDataHist cs("cs","cs",t,RooFit::Import(*h3));
    RooDataHist cs_bkg("cs_bkg","cs_bkg",t,RooFit::Import(*h)); //background from histo
    
    //Setup component pdfs (1 peak)
    RooRealVar A("A", "A",6.7,1.,100.);
    
    RooCompton compton("compton", "compton",t,A,hv);
    
    RooRealVar mg("mg","mg",0.);
    RooRealVar sg("sg","sg",4500., 0., 100000.);
    RooGaussian gauss("gauss","gauss",t,mg,sg);
    
    t.setBins(15000,"cache");
    
    RooFFTConvPdf cxg("cxg","compton (X) gauss",t,compton,gauss);
    
    RooHistPdf bkg_cs("bkg_cs","pdf from hist",t,cs_bkg,0); //background from histo
    RooRealVar bkg_frac_cs("bkg_frac_cs","fraction of bkg",0.1,0.,1.);
    RooAddPdf model_cs("model_cs","cxg+bkg",RooArgList(cxg,bkg_cs),bkg_frac_cs);

    t.setRange("R1",csmin,csmax);
    t.setRange("R2",csmin,180000);
    
    RooPlot* frame = t.frame(RooFit::Title(title_2));
    cs.plotOn(frame);

    cxg.fitTo(cs,RooFit::Range("R1"));
    cxg.plotOn(frame,RooFit::Components(cxg),RooFit::LineColor(kGreen));
    
    model_cs.fitTo(cs,RooFit::Range("R2"));
    //model_cs.plotOn(frame,RooFit::LineColor(kRed)); //somma Compton + bkg
    model_cs.plotOn(frame,RooFit::Components(bkg_cs),RooFit::LineColor(kViolet));

    A_cs = A.getValV();
    A_cs_err = A.getError();
    sig_cs = sg.getValV();
    sig_cs_err = sg.getError();
    
    frame->Draw();
    frame->SetXTitle ("Energy (ADC*Sample)");
    frame->SetYTitle ("Events");
    
////////////// 22Na Spectrum ////////////////////
    c4->cd();//->SetLogy();
    gStyle->SetOptFit(111);
    
    TH1F *h4 = new TH1F("h4",title_3,300,0,300000); //ADC*Sample
    //TH1F *h4 = new TH1F("h4",title_3,1000,0,2000); //in keV calibrata su 241Am - 137Cs
    
    reco_3->Draw(Form("charge[%d]>>h4",nch),time); //ADC*Sample
    //reco_3->Draw(Form("charge[%d]*0.00754546>>h4",nch),time); //in keV calibrata su 241Am - 137Cs
    
    h4->GetXaxis()->SetTitle("Energy (ADC*Sample)");
    //h4->GetXaxis()->SetTitle("Energy (keV)"); //nel caso si utilizzi l'energia calibrata
    h4->GetXaxis()->SetTitleOffset(1.16);
    h4->GetYaxis()->SetTitle("Counts (arb.)");
    h4->GetYaxis()->SetTitleOffset(1.16);
    h4->SetLineWidth(1);
    h4->SetLineColor(kBlack);
    
    //h4->Add(h,-1);
    
    //h4->Draw();

////////////// Fitting 22Na Spectrum using RooFit ////////////////////
    RooRealVar t1("t1","t1", 1000.,200000.); //range dell'histo
    RooRealVar hv1("hv1","hv1", 511000); //energia picco1 fotoelettrico in eV
    RooRealVar hv2("hv2","hv2", 1274540); //energia picco2 fotoelettrico in eV
    
    RooDataHist na("na","na",t1,RooFit::Import(*h4));
    RooDataHist na_bkg("na_bkg","na_bkg",t1,RooFit::Import(*h)); //background from histo
    
    //Setup component pdfs (1 peak)
    RooRealVar A1("A1", "A1",6.7,2.,20.);
    RooRealVar A2("A2", "A2",6.7,2.,20.);
    
    RooCompton compton1("compton1", "compton1",t1,A1,hv1);
    RooCompton compton2("compton2", "compton2",t1,A2,hv2);
    
    RooRealVar mg1("mg1","mg1",0.);
    RooRealVar sg1("sg1","sg1",4500., 1000., 10000.);
    RooGaussian gauss1("gauss1","gauss1",t1,mg1,sg1);
    
    RooRealVar mg2("mg2","mg2",0.);
    RooRealVar sg2("sg2","sg2",2000, 1000., 10000.);
    RooGaussian gauss2("gauss2","gauss2",t1,mg2,sg2);
    
    t1.setBins(15000,"cache");
    
    RooFFTConvPdf cxg1("cxg1","compton1 (X) gauss1",t1,compton1,gauss1);
    RooFFTConvPdf cxg2("cxg2","compton2 (X) gauss2",t1,compton2,gauss2);
    
    RooRealVar sigfrac("sigfrac","fraction of signal",0.5,0.,1.);
    RooAddPdf model_na1("model_na1","cxg1+cxg2",RooArgList(cxg1,cxg2),sigfrac);
    
    RooHistPdf bkg_na("bkg_na","pdf from hist",t1,na_bkg,0); //background from histo
    RooRealVar bkg_frac_na("bkg_frac_na","fraction of bkg",0.1,0.,1.);
    RooAddPdf model_na2("model_na2","cxg1+cxg2+bkg",RooArgList(model_na1,bkg_na),bkg_frac_na);
    
    t1.setRange("R1",namin1,namax1);
    t1.setRange("R2",namin2,namax2);
    t1.setRange("R3",30000,190000);
    
    
    RooPlot* frame1 = t1.frame(RooFit::Title(title_3));
    na.plotOn(frame1);
    
//Model to fit single peaks

    cxg1.fitTo(na,RooFit::Range("R1"));
    cxg1.plotOn(frame1,RooFit::Components(cxg1),RooFit::LineColor(kGreen));
    
    cxg2.fitTo(na,RooFit::Range("R2"));
    cxg2.plotOn(frame1,RooFit::Components(cxg2),RooFit::LineColor(kBlue));

//Model to fit peak1+peak2
    model_na1.fitTo(na,RooFit::Range("R3"));
//    model_na1.plotOn(frame1, RooFit::LineColor(kRed));
//    model_na1.plotOn(frame1,RooFit::Components(cxg1),RooFit::LineColor(kViolet));
//    model_na1.plotOn(frame1,RooFit::Components(cxg2),RooFit::LineColor(kPink));

//Model to fit peak1+peak2+bkg

    model_na2.fitTo(na,RooFit::Range("R3"));
//    model_na2.plotOn(frame1,RooFit::LineColor(kRed));
//    model_na2.plotOn(frame1,RooFit::Components(model_na1),RooFit::LineColor(kYellow)); //somma delle due componenti Compton
    model_na2.plotOn(frame1,RooFit::Components(bkg_na),RooFit::LineColor(kViolet)); //bkg

    A_na1 = A1.getValV();
    A_na2 = A2.getValV();
    A_na_err1 = A1.getError();
    A_na_err2 = A2.getError();
    sig_na1 = sg1.getValV();
    sig_na2 = sg2.getValV();
    sig_na_err1 = sg1.getError();
    sig_na_err2 = sg2.getError();
    
    frame1->Draw();
    frame1->SetXTitle ("Energy (ADC*Sample)");
    frame1->SetYTitle ("Events");
    
////////////// Calibration Curve ////////////////////
    
    A_am = 59.54/mu_am;
    double energy[4] {59.54,340.67,477.34,1061.71}; //in keV Am: 59.54 22Na-secondo: 1061.71
    double k[4] {A_am,A_na1/1E3,A_cs/1E3,A_na2/1E3}; //valori di A calcolati dai fit precedenti ed espressi in keV/ADC A_am,A_na2/1E3
    double k_err[4] {mu_am_err,A_na_err1/1E3,A_cs_err/1E3,A_na_err2/1E3}; //mu_am_err, A_na_err2/1E3
    double sis[4] {mu_am_err,1849.008,2672.83,2672.83}; //sistematiche calcolate dai test di Marco
    
    double ene_am[1] {59.54};
    double k_am[1] {A_am};
    double ene_err[1] {1.};
    double sis_am[1] {mu_am_err};
    
    c5->cd();
    //TH2F *grframe = new TH2F("grframe","",300,0,1100,300,0,200000);
    TH2F *grframe = new TH2F("grframe","",300,0,200000,300,0,1200);
    gStyle->SetOptStat(0);
    TGraphErrors *gr = new TGraphErrors();
        
        for (int i=1;i<4;i++){
    
            //gr->SetPoint(i,energy[i],energy[i]/k[i]); //con sistematiche
            //gr->SetPointError(i,1.,sis[i]); //con sistematiche
            gr->SetPoint(i,energy[i]/k[i],energy[i]);
            gr->SetPointError(i,sis[i],1.);
    
        }
    
    grframe->SetTitle(title_4);
    grframe->GetXaxis()->SetTitle("Ch. (ADC*Sample)");
    grframe->GetYaxis()->SetTitle("Energy (keV)");
    grframe->GetYaxis()->SetTitleOffset(1.2);
    
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.7);
    //gr->Draw("PE");

//fit ad intercetta nulla
    TF1 *mypol = new TF1("mypol","[0]+[1]*x",0,170000); //fit ad intercetta nulla
    mypol->FixParameter(0,0.);
    mypol->SetNpx(10000);
    mypol->SetLineWidth(2);
    mypol->SetLineColor(kRed);
/*
    TGraphErrors *gram = new TGraphErrors();
    
    for (int i=0;i<1;i++){
        
        gram->SetPoint(i,energy[i]/k[i],energy[i]);
        gram->SetPointError(i,sis[i],1.);
    }
    
    gram->SetMarkerStyle(20);
    gram->SetMarkerSize(0.7);
*/
    
    gr->Fit("mypol","EMR");
    grframe->Draw();
    gr->Draw("PEsame");
    mypol->Draw("same");
//    gram->Draw("PEsame");

//fit a parametri liberi
    c6->cd();
    TH2F *grframe1 = new TH2F("grframe1","",300,0,200000,300,0,1200);
    gStyle->SetOptStat(0);
    TGraphErrors *gr1 = new TGraphErrors();
    
    for (int i=1;i<4;i++){
        
        //gr1->SetPoint(i,energy[i],energy[i]/k[i]); //con sistematiche
        //gr1->SetPointError(i,1.,sis[i]); //con sistematiche
        gr1->SetPoint(i,energy[i]/k[i],energy[i]);
        gr1->SetPointError(i,sis[i],1.);
        
    }
    
    grframe1->SetTitle(title_4);
    grframe1->GetXaxis()->SetTitle("Ch. (ADC*Sample)");
    grframe1->GetYaxis()->SetTitle("Energy (keV)");
    grframe1->GetYaxis()->SetTitleOffset(1.2);
    
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(0.7);
    //gr->Draw("PE");
    
    TF1 *mypol1 = new TF1("mypol1","[0]+[1]*x",0,170000); //fit a parametri liberi
    mypol1->SetNpx(10000);
    mypol1->SetLineWidth(2);
    mypol1->SetLineColor(kBlue);
/*
    TGraphErrors *gram1 = new TGraphErrors();
    
    for (int i=0;i<1;i++){
        
        gram1->SetPoint(i,energy[i]/k[i],energy[i]);
        gram1->SetPointError(i,sis[i],1.);
    }
    
    gram1->SetMarkerStyle(20);
    gram1->SetMarkerSize(0.7);
*/
    gr1->Fit(mypol1,"EMR");
    grframe1->Draw();
    gr1->Draw("PEsame");
    mypol1->Draw("same");
//    gram1->Draw("PEsame");

//fit Birks
    c7->cd();
    TH2F *grframe2 = new TH2F("grframe2","",300,0,200000,300,0,1200);
    gStyle->SetOptStat(0);
    TGraphErrors *gr2 = new TGraphErrors();
    
    for (int i=0;i<4;i++){
        
        //gr2->SetPoint(i,energy[i],energy[i]/k[i]); //con sistematiche
        //gr2->SetPointError(i,1.,sis[i]); //con sistematiche
        gr2->SetPoint(i,energy[i]/k[i],energy[i]);
        gr2->SetPointError(i,sis[i],1.);
        
    }
    
    grframe2->SetTitle(title_4);
    grframe2->GetXaxis()->SetTitle("Ch. (ADC*Sample)");
    grframe2->GetYaxis()->SetTitle("Energy (keV)");
    grframe2->GetYaxis()->SetTitleOffset(1.2);
    
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(0.7);
    //gr2->Draw("PE");
    
    TF1 * cqe = new TF1("cqe",CQ,0.,170000.,3); //fit Birks
    cqe->SetParameter(0,0.012); //0.012; 0.02 per il pmt8
    cqe->FixParameter(1,0.);
    cqe->SetParameter(2,mypol->GetParameter(1));
    cqe->SetLineWidth(2);
    cqe->SetLineColor(kViolet);
    cqe->SetParNames("kB","c0","c1");
    gr2->Fit("cqe","","",0.,170000.);
/*
    TGraphErrors *gram2 = new TGraphErrors();
    
    for (int i=0;i<1;i++){
        
        gram2->SetPoint(i,energy[i]/k[i],energy[i]);
        gram2->SetPointError(i,sis[i],1.);
    }
    
    gram2->SetMarkerStyle(20);
    gram2->SetMarkerSize(0.7);
*/
    
    grframe2->Draw();
    gr2->Draw("PEsame");
    //printf ("%f\n",cqe->GetChisquare());
    cqe->Draw("same");
//    gram2->Draw("PEsame");

//Output
    double cal = mypol1->GetParameter(1);
    double cal_err = mypol1->GetParError(1);
    double q = mypol1->GetParameter(0);
    double q_err = mypol1->GetParError(0);

    cout << "Calibration Parameter: " << mypol1->GetParameter(1) << " +/- " << mypol1->GetParError(1) << endl;
    
    cout << "Am: " << mu_am << " +/- " << mu_am_err << "\n"
         << "Cs: " << 477.34/A_cs*1E3 << " +/- " << A_cs_err*1E3 << "\n"
         << "Na1: " << 340.67/A_na1*1E3 << " +/- " << A_na_err1*1E3 << "\n"
         << "Na2: " << 1061.71/A_na2*1E3 << " +/- " << A_na_err2*1E3 << endl;
    
    cout << "\n" << "Calibrated Energies: " << endl;

    //in caso di retta "normale"
    cout << "Am: " << (mu_am*cal)+q << " +/- " << mu_am_err << "\n"
         << "Cs: " << (477.34/A_cs*1E3*cal)+q << " +/- " << A_cs_err*1E3 << "\n"
         << "Na1: " << (340.67/A_na1*1E3*cal)+q << " +/- " << A_na_err1*1E3 << "\n"
         << "Na2: " << (1061.71/A_na2*1E3*cal)+q << " +/- " << A_na_err2*1E3 << endl;
/*
    //in caso di retta "inversa"
    cout << "Am: " << (mu_am-q)/cal << " +/- " << mu_am_err << "\n"
         << "Cs: " << (477.34/A_cs*1E3-q)/cal << " +/- " << A_cs_err*1E3 << "\n"
         << "Na1: " << (340.67/A_na1*1E3-q)/cal << " +/- " << A_na_err1*1E3 << "\n"
         << "Na2: " << (1061.71/A_na2*1E3-q)/cal << " +/- " << A_na_err2*1E3 << endl;
*/

/*
    ofile << mu_am << " " << mu_am_err << "\n"
          << 477.34/A_cs*1E3 << " " << A_cs_err*1E3 << "\n"
          << 340.67/A_na1*1E3 << " " << A_na_err1*1E3 << "\n"
          << 1061.71/A_na2*1E3 << " " << A_na_err2*1E3 << "\n"
          << mypol1->GetParameter(1) << " " << mypol1->GetParError(1) << endl;


    ofile2 << "Am: " << mu_am*cal << " +/- " << mu_am_err << "\n"
           << "Cs: " << A_cs*1E3*cal << " +/- " << A_cs_err*1E3 << "\n"
           << "Na1: " << A_na1*1E3*cal << " +/- " << A_na_err1*1E3 << "\n"
           << "Na2: " << A_na2*1E3*cal << " +/- " << A_na_err2*1E3 << endl;
    
*/
 /*
////////////// Diff. vs True Energy Curve (Residual) ////////////////////
    c6->cd();
    
    double e_th[4]{59.54,340.67,477.34,1061.71};
    double e_cal[4]{(mu_am*cal)+q,(340.67/A_na1*1E3*cal)+q,(477.34/A_cs*1E3*cal)+q,(1061.71/A_na2*1E3*cal)+q};
    double e_cal_err[4]{mu_am_err,A_na_err1*1E3,A_cs_err*1E3,A_na_err2*1E3};
 
    TGraphErrors *grk = new TGraphErrors();
    
    for (int i=0;i<4;i++){
 
        grk->SetPoint(i,e_th[i],e_th[i]-e_cal[i]);
        grk->SetPointError(i,1.,e_cal_err[i]);
 
    }
 
    grk->SetTitle(title_6);
    grk->GetYaxis()->SetTitle("Res.");
    grk->GetXaxis()->SetTitle("Energy (keV)");
    grk->GetYaxis()->SetTitleOffset(1.2);
    grk->SetMarkerStyle(20);
    grk->SetMarkerSize(0.2);
    grk->Draw("APE");
 
*/
/*
////////////// Resolution Curve ////////////////////
    c8->cd();
    
    double ene[4] {59.54,340.67,477.34,1061.71}; //in keV Am: 59.54 22Na-secondo: 1061.71
    double res[4] {sig_am/59540*100,sig_na1/340670*100,sig_cs/477340*100,sig_na2/1061710*100};
    double res_err[4] {sig_am_err/sig_am,sig_na_err1/sig_na1,sig_cs_err/sig_cs,sig_na_err2/sig_na2};
    
    TGraphErrors *gres = new TGraphErrors();
    for (int i=0;i<4;i++){
        
        gres->SetPoint(i,ene[i],res[i]);
        gres->SetPointError(i,1.,res_err[i]);
        
    }
    
    gres->SetTitle(title_5);
    gres->GetXaxis()->SetTitle("Energy (keV)");
    gres->GetYaxis()->SetTitle("Resolution (%)");
    gres->GetYaxis()->SetTitleOffset(1.1);
    gres->SetMarkerStyle(20);
    gres->SetMarkerSize(0.7);
    //gres->Draw("APE");
    
    TF1 *mypol2 = new TF1("mypol2","pol2",0,1100);
    mypol2->SetLineWidth(2);
    mypol2->SetLineColor(kRed);
    gres->Fit(mypol2,"","",0,1100);
    gres->Draw("APE");
    mypol2->Draw("same");
*/
    
    
/*
////////////// A vs Energy Curve ////////////////////
    TCanvas *c7 = new TCanvas("c7","c7",10,10,900,600);
    
    TGraphErrors *grk = new TGraphErrors();
    for (int i=0;i<4;i++){
        
        grk->SetPoint(i,energy[i],k[i]);
        //grk->SetPointError(i,0.,k_err[i]);
        
    }
    
    grk->SetTitle(title_6);
    grk->GetYaxis()->SetTitle("A (keV/ADC)");
    grk->GetXaxis()->SetTitle("Energy (keV)");
    grk->GetYaxis()->SetTitleOffset(1.2);
    grk->SetMarkerStyle(21);
    grk->SetMarkerSize(1.5);
    grk->Draw("APE");
*/
    
}
