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
#include <TH2F.h>
#include <TProfile.h>
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
#include <TH1I.h>

#include "red-daq/EvRec0.hh"
#include "red-daq/RDCluster.hh"

/*
  This macro is meant to extract informations about S1 and S2 yields together
  with the electron lifetime for the double-phase runs taken with Am241 source.
  
  It accepts the name of the root file as input together with some other
  parameters such as, in order: the run number (used in case of an autmatic
  save of spectra), the value of drift field, number of bins per each histogram,
  axis ranges of the s1 vs tdrift, s2/s1 vs tdrift, s1 and s2 corrected plots and finally
  fitting range for the above plots.

  S1 is corrected starting from the s1 vs tdrift spectrum. An order three polinomial
  is used as fitting function, and then an evaluation of the drift time in the middle of
  the chamber is done. Finally, S1 is normalized by the ratio bewtween the above mentioned
  drift time over the total value. This is done per each value of drift field. 

  S1 corrected is then used to plot S2/S1 ratio versus the drift time. So from this
  an estimation of the electron lifetime is extracted (tau).
  
  Once the above tau value is known, a correction in S2 is also performed as TMath::Exp(-tdrift/tau). 

  Finally, corrected S1 and S2 are fitted by using a Monte Carlo function provided by Davide Franco
  that accounts also for the gaussin smearing of the detector response.

  Error bars (when available) are purely statistical.
  
  It works under ROOT, as:
  .L ly_v2_ph2.C+
   ly_v2_ph2("root_filename",run_number,200,300,500,800,7,25,0,1400,0,25000,400,850,5000,16000,0,80,16,66,14,66)
  in case of a 200 V/cm run taken in Catania.
  
  Last edit by Simone S. on 31 Jan. 2020

 */

Double_t fun(Double_t *x, Double_t *par)
{
  
  return par[0]+par[1]*x[0]+par[2]*TMath::Power(x[0],2)+par[3]*TMath::Power(x[0],3);
  
}

Double_t fun2(Double_t *x, Double_t *par)
{
  
  return TMath::Exp(par[0]+par[1]*x[0]);
  //return TMath::Exp(-par[0]/par[1]*x[0]);
  
}

using namespace std ;
using namespace TMath ;

TH1F *ham ;
double fam241(double *x, double *p) {
  double LY  = p[0];
  double sigma0 = p[1] ;
  double fano   = p[2] ;
  double width = ham->GetBinWidth(3);
  double tot = 0 ;
  for(int i=1;i<ham->GetXaxis()->GetNbins();++i) {
    double ene = ham->GetBinCenter(i);
    if(ene > 60) continue; 
    double val   = ham->GetBinContent(i);
    double npe   = LY*ene ;
    double sigma = sqrt(pow(sigma0,2)+npe*fano);
    //cout << npe << " " << sigma << " " << fano << endl;
    tot += val*Gaus(x[0],npe,sigma,1)/width;    
  }
  
  return p[3]*tot ;
}

void ly_v2_ph2(TString filename, int run, int Ed=200, int nbin=300, int s1tmin=500, int s1tmax=800, int s2s1tmin=10, int s2s1tmax=30, int rmin=0, int rmax=1400, int r_s2_min=0, int r_s2_max=25000, int mins1=400, int maxs1=800, int mins2=9000, int maxs2=15800, int r_min=0, int r_max=80, int pmin1=15, int pmax1=58, int pmin2=15, int pmax2=58)
{

  TFile *f = new TFile(filename,"read");
    
//Titles
  
  TString title = "^{241}Am Spectrum";
  TString title_drift_tot = "S1 vs Drift Time";
  TString title2_drift_tot = "S2/S1 vs Drift Time";

//Canvases and Style
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1800,1200);
  c1->Divide(3,2);
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(111);
  gStyle->SetOptStat("");
 
  //let's start 
  TTree* reco = (TTree*) f->Get("reco");
  EvRec0* evReco = new EvRec0();
  reco->SetBranchAddress("recoevent",&evReco);
  cout << "reco entries: " << reco->GetEntries() << endl;
  
  //S1 vs tdrift
  TH2F *h_drift_tot = new TH2F("h_drift_tot",title_drift_tot,300,0,200,300,s1tmin,s1tmax); //0,5000
  
  //S1 and S1_corr
  TH1F *h = new TH1F("h",title,nbin,rmin,rmax);
  TH1F *h_corr = new TH1F("h_corr",title,nbin,rmin,rmax);//s1 totale

  //S2/S1_corr vs tdrift
  TH2F *h2_drift_tot = new TH2F("h2_drift_tot",title2_drift_tot,300,0,200,300,s2s1tmin,s2s1tmax);
  
  //S2 and S2_corr
  TH1F *h2 = new TH1F("h2",title,nbin,r_s2_min,r_s2_max);
  TH1F *h2_corr = new TH1F("h2_corr",title,nbin,r_s2_min,r_s2_max);
  
  //Corrections distributions
  TH1F *h_dist_tot = new TH1F("h_dist_tot","",300,0,2);
  TH1F *h2_dist_tot = new TH1F("h2_dist_top","",300,0,2);
  
  //Tdrift
  TH1F *h_time = new TH1F("h_time","",300,0,100);
    
   //S1 vs tdrift
  //loop sul numero di eventi
  for (Int_t i=0;i<reco->GetEntries();i++)
    {            
    
      reco->GetEntry(i); 
      
      vector<double> f90 = evReco->GetF90();
      vector<double> start_time = evReco->GetStartTime();
      vector<double> charge = evReco->GetCharge();
      vector<RDCluster *> clusters = evReco->GetClusters();
      int nclusters = evReco->GetNClusters();
      
      if (nclusters==2){
      
        double f90_1 = clusters.at(0)->f90;
        double f90_2 = clusters.at(1)->f90;
      
        double charge = clusters.at(0)->charge;
      	double charge_bottom = clusters.at(0)->tot_charge_bottom;
	double charge_top = clusters.at(0)->tot_charge_top;
      	double tdrift = (clusters.at(1)->cdf_time-clusters.at(0)->start_time)*2./1000.;
	double rep_1 = clusters.at(0)->rep;
	double rep_2 = clusters.at(1)->rep;
	
	if(Ed==100 && f90_1>0.2 && f90_2<0.2 && charge>600 && charge<800 && rep_1==1 && rep_2==1) h_drift_tot->Fill(tdrift,charge);
      
	//if(Ed==200 && f90_1>0.2 && f90_2<0.2 && charge>550 && charge<750 && rep_1==1 && rep_2==1) h_drift_tot->Fill(tdrift,charge);
	
	if(Ed==200 && f90_1>0.2 && f90_2<0.2 && charge>550 && charge<750 && rep_1==1 && rep_2==1) h_drift_tot->Fill(tdrift,charge); // configurazione LNS Dic. 2019
	
	if(Ed==400 && f90_1>0.2 && f90_2<0.2 && charge>470 && charge<680 && rep_1==1 && rep_2==1) h_drift_tot->Fill(tdrift,charge);
      
	if(Ed==700 && f90_1>0.2 && f90_2<0.2 && charge>400 && charge<610 && rep_1==1 && rep_2==1) h_drift_tot->Fill(tdrift,charge);
      
	if(Ed==1000 && f90_1>0.2 && f90_2<0.2 && charge>354 && charge<556 && rep_1==1 && rep_2==1) h_drift_tot->Fill(tdrift,charge);
	
	if(f90_1>0.2 && f90_2<0.2 && rep_1==1 && rep_2==1) h_time->Fill(tdrift);
	
	}
      
     }

  gStyle->SetOptStat("");
  gStyle->SetOptFit(111);
  c1->cd(1);
  gStyle->SetOptFit(111);
  gStyle->SetOptStat("");
  TF1 *fit_tot = new TF1("fit_tot",fun,0,80,4);    
     
     //fit
  TProfile *p = h_drift_tot->ProfileX();
  p->SetMarkerStyle(21);
  p->SetMarkerSize(0.5);
  fit_tot->SetLineColor(kRed);
  p->Fit("fit_tot","EMQ","",pmin1,pmax1); //7,26
  gStyle->SetOptFit(111);
  gStyle->SetOptStat("");
    
  h_drift_tot->GetXaxis()->SetTitle("t_{drift} (#mus)");
  h_drift_tot->GetXaxis()->SetTitleOffset(1.16);
  h_drift_tot->GetYaxis()->SetTitle("Charge (PE)");
  h_drift_tot->GetYaxis()->SetTitleOffset(1.16);
  h_drift_tot->GetXaxis()->SetRangeUser(r_min,r_max); //(12,26)
  h_drift_tot->GetYaxis()->SetRangeUser(s1tmin,s1tmax);
    
  h_drift_tot->Draw("colz");
  p->Draw("same");
   
  //calcolo le correzioni S1
  double corr_s1=0;
  double thalf_tot=0;
  double t=0;
  
  for (int i=0;i<reco->GetEntries();i++){
      	
	reco->GetEntry(i); 
      
      	vector<double> f90 = evReco->GetF90();
      	vector<double> start_time = evReco->GetStartTime();
      	vector<double> charge = evReco->GetCharge();
      	vector<RDCluster *> clusters = evReco->GetClusters();
      	int nclusters = evReco->GetNClusters();
      
      	if (nclusters==2){
      
        	double f90_1 = clusters.at(0)->f90;
        	double f90_2 = clusters.at(1)->f90;
      
      		double charge_bottom = clusters.at(0)->tot_charge_bottom;
		double charge_top = clusters.at(0)->tot_charge_top;
		double charge = clusters.at(0)->charge;
		
		double s2_bottom = clusters.at(1)->tot_charge_bottom;
		double s2_top = clusters.at(1)->tot_charge_top;
		
      		double tdrift = (clusters.at(1)->fixed_time-clusters.at(0)->start_time)*2./1000.;
		
		double rep_1 = clusters.at(0)->rep;
		double rep_2 = clusters.at(1)->rep;

		
		if (Ed==100){
		
			thalf_tot = fit_tot->Eval(40.3);
			
		}else if(Ed==200){
		
			thalf_tot = fit_tot->Eval(30.17);
			//thalf_tot = fit_tot->Eval(33.17);
			
		}else if(Ed==400){
		
			thalf_tot = fit_tot->Eval(18.67);
			
		}else if(Ed==700){
		
			thalf_tot = fit_tot->Eval(14.29);
			
		}else{
		
			thalf_tot = fit_tot->Eval(12.4);
			
		}
		

		//cout << "io sono thalf: " << thalf_tot << endl;

		t = fit_tot->Eval(tdrift);
		corr_s1 = thalf_tot/t; //correzione in S1
		
		if(f90_1>0.2 && f90_2<0.2 && rep_1==1 && rep_2==1){
		
			 h_corr->Fill(charge*corr_s1); //riempio S1 corretto
	
			 h->Fill(charge);
			
		}
		
		
		//spettro di distribuzione correzioni S1
		h_dist_tot->Fill(corr_s1);
		
				        
      }
  }
  
  //S2/S1_corr vs tdrift
  
  for (int i=0;i<reco->GetEntries();i++){
      	
	reco->GetEntry(i); 
      
      	vector<double> f90 = evReco->GetF90();
      	vector<double> start_time = evReco->GetStartTime();
      	vector<double> charge = evReco->GetCharge();
      	vector<RDCluster *> clusters = evReco->GetClusters();
      	int nclusters = evReco->GetNClusters();
      
      	if (nclusters==2){
      
        	double f90_1 = clusters.at(0)->f90;
        	double f90_2 = clusters.at(1)->f90;
      
      		double charge_bottom = clusters.at(0)->tot_charge_bottom;
		double charge_top = clusters.at(0)->tot_charge_top;
		double charge = clusters.at(0)->charge;
		
		double s2_bottom = clusters.at(1)->tot_charge_bottom;
		double s2_top = clusters.at(1)->tot_charge_top;
		double s2 = clusters.at(1)->charge;
		
      		double tdrift = (clusters.at(1)->fixed_time-clusters.at(0)->start_time)*2./1000.;
		
		double rep_1 = clusters.at(0)->rep;
		double rep_2 = clusters.at(1)->rep;
		
		if(Ed==100 && f90_1>0.2 && f90_2<0.2 && s2>5560 && s2<9700 && rep_1==1 && rep_2==1) h2_drift_tot->Fill(tdrift,s2/(charge*corr_s1));
		
		//if(Ed==200 && f90_1>0.2 && f90_2<0.2 && s2>10000 && s2<15500 && rep_1==1 && rep_2==1) h2_drift_tot->Fill(tdrift,s2/(charge*corr_s1));
		
		if(Ed==200 && f90_1>0.2 && f90_2<0.2 && s2>5000 && s2<11000 && charge>550 && charge<750 && rep_1==1 && rep_2==1) h2_drift_tot->Fill(tdrift,s2/(charge*corr_s1)); //configurazione LNS Dic.19
		
		//if(Ed==200 && f90_1>0.2 && f90_2<0.2 && s2>0 && s2<2000 && rep_1==1 && rep_2==1 && tdrift>30) h2_drift_tot->Fill(tdrift,s2/(charge*corr_s1)); //configurazione LNS Dic.19 run_1396
		
		//if(Ed==200 && f90_1>0.2 && f90_2<0.2 && s2>5000 && s2<14000 && rep_1==1 && rep_2==1) h2_drift_tot->Fill(tdrift,s2/(charge*corr_s1)); //configurazione LNS Dic.19 campi Vlad
		
		//if(Ed==400 && f90_1>0.2 && f90_2<0.2 && s2>15100 && s2<21800 && rep_1==1 && rep_2==1) h2_drift_tot->Fill(tdrift,s2/(charge*corr_s1));
		
		if(Ed==400 && f90_1>0.2 && f90_2<0.2 && s2>8000 && s2<17500 && rep_1==1 && rep_2==1) h2_drift_tot->Fill(tdrift,s2/(charge*corr_s1)); //LNS Dic. 19
		
		if(Ed==700 && f90_1>0.2 && f90_2<0.2 && s2>21000 && s2<28000 && rep_1==1 && rep_2==1) h2_drift_tot->Fill(tdrift,s2/(charge*corr_s1));
		
		if(Ed==1000 && f90_1>0.2 && f90_2<0.2 && s2>24000 && s2<33500 && rep_1==1 && rep_2==1) h2_drift_tot->Fill(tdrift,s2/(charge*corr_s1));
		
		
	}
  }		
  
  gStyle->SetOptFit(111);
  gStyle->SetOptStat("");     
  TF1 *fit2_tot = new TF1("fit2_tot",fun2,0,80,2);  
  //fit2_tot->SetParNames("Constant","#tau");
  //fit2_tot->SetParameters(2.95468e+00,-3.63235e-04);
     
     //fit
  c1->cd(2);
  gStyle->SetOptFit(111);
  gStyle->SetOptStat("");
  TProfile *p2 = h2_drift_tot->ProfileX();
  p2->SetMarkerStyle(21);
  p2->SetMarkerSize(0.5);
    
  fit2_tot->SetLineColor(kRed);
  p2->Fit("fit2_tot","EM","",pmin2,pmax2);
  gStyle->SetOptFit(111);
  gStyle->SetOptStat("");
/*  
  h2_drift_tot->GetXaxis()->SetTitle("t_{drift} (#mus)");
  h2_drift_tot->GetXaxis()->SetTitleOffset(1.16);
  h2_drift_tot->GetYaxis()->SetTitle("S2/S1");
  h2_drift_tot->GetYaxis()->SetTitleOffset(1.16);
  h2_drift_tot->GetXaxis()->SetRangeUser(r_min,r_max); //(12,26)
  //h2_drift_tot->GetYaxis()->SetRangeUser(0,500);
     
  //h2_drift_tot->Draw("colz");
*/
  p2->GetXaxis()->SetTitle("t_{drift} (#mus)");
  p2->GetXaxis()->SetTitleOffset(1.16);
  p2->GetYaxis()->SetTitle("S2/S1");
  p2->GetYaxis()->SetTitleOffset(1.16);
  p2->GetXaxis()->SetRangeUser(r_min,r_max); //(12,26)
  p2->GetYaxis()->SetRangeUser(s2s1tmin,s2s1tmax);
  
  p2->Draw("same");
  
  //calcolo le correzioni S2;
  double corr_s2 = 0;
  double tau = -1./fit2_tot->GetParameter(1);
  
  cout << "tau: " << tau << " us" << endl;
  
  for (int i=0;i<reco->GetEntries();i++){
      	
	reco->GetEntry(i); 
      
      	vector<double> f90 = evReco->GetF90();
      	vector<double> start_time = evReco->GetStartTime();
      	vector<double> charge = evReco->GetCharge();
      	vector<RDCluster *> clusters = evReco->GetClusters();
      	int nclusters = evReco->GetNClusters();
      
      	if (nclusters==2){
      
        	double f90_1 = clusters.at(0)->f90;
        	double f90_2 = clusters.at(1)->f90;
      
      		double charge_bottom = clusters.at(0)->tot_charge_bottom;
		double charge_top = clusters.at(0)->tot_charge_top;
		double charge = clusters.at(0)->charge;
		
		double s2_bottom = clusters.at(1)->tot_charge_bottom;
		double s2_top = clusters.at(1)->tot_charge_top;
		double s2 = clusters.at(1)->charge;
		
      		double tdrift = (clusters.at(1)->fixed_time-clusters.at(0)->start_time)*2./1000.;
		
		double rep_1 = clusters.at(0)->rep;
		double rep_2 = clusters.at(1)->rep;
		
		//s2
		corr_s2 = TMath::Exp(-tdrift/tau);
		
		if(f90_1>0.2 && f90_2<0.2 && rep_1==1 && rep_2==1){ 
		
			h2->Fill(s2);
			//h2_corr->Fill(s2/corr_s2); //s2 total
			
			 if(Ed==100 && charge>600 && charge<800) h2_corr->Fill(s2/corr_s2); 
			 //if(Ed==200 && charge>550 && charge<750) h2_corr->Fill(s2/corr_s2); 
			 if(Ed==200 && charge>550 && charge<750) h2_corr->Fill(s2/corr_s2); //config. LNS Dec.19
			 if(Ed==400 && charge>470 && charge<680) h2_corr->Fill(s2/corr_s2);  
			 if(Ed==700 && charge>400 && charge<610) h2_corr->Fill(s2/corr_s2); 
			 if(Ed==1000 && charge>354 && charge<556) h2_corr->Fill(s2/corr_s2);
		
		   }
		   
		//spettri di distribuzione delle correzioni S2
		h2_dist_tot->Fill(corr_s2);
		
		}
		
	}
	
  //disegno distribuzione correzioni
  c1->cd(3);
  
  h_dist_tot->SetTitle("Correction distribution S1");
  h_dist_tot->GetXaxis()->SetTitle("Correzione");
  h_dist_tot->GetXaxis()->SetTitleOffset(1.16);
  h_dist_tot->GetYaxis()->SetTitle("Counts (arb.)");
  h_dist_tot->GetYaxis()->SetTitleOffset(1.16);
  
  h_dist_tot->Draw();
  
  c1->cd(4);
  
  h2_dist_tot->SetTitle("Correction distribution S2");
  h2_dist_tot->GetXaxis()->SetTitle("Correzione");
  h2_dist_tot->GetXaxis()->SetTitleOffset(1.16);
  h2_dist_tot->GetYaxis()->SetTitle("Counts (arb.)");
  h2_dist_tot->GetYaxis()->SetTitleOffset(1.16);
  
  h2_dist_tot->Draw();
  
  c1->cd(5);
  
  h_time->SetTitle("Drift Time Spectrum");
  h_time->GetXaxis()->SetTitle("t_{drift} (#mus)");
  h_time->GetYaxis()->SetTitle("Counts (arb.)");
  h_time->Draw();
				

  //fitto e disegno S1 corretto
  
  TCanvas *c2 = new TCanvas("c2","c2",10,10,900,600);

/* 
  TF1 *fit_gaus = new TF1("fit_gaus","gaus(0)",mins1,maxs1); // fitting range 620-800 singola fase; 400-600 doppia fase
        
  fit_gaus->SetLineWidth(3);
  fit_gaus->SetLineColor(kRed);
    
  h_corr->Fit("fit_gaus","EMRQ");
*/

////////// FITTER DAVIDE ///////////////

  TFile *fmc = TFile::Open("am241.root");
  ham = (TH1F*) fmc->Get("ham");
  TF1 *fun3 = new TF1("fun3",fam241,0,1000,4);
  fun3->SetParameters(4,5,2.,1e4);
  fun3->SetLineColor(64);
  fun3->Draw();
  fun3->SetParNames("LY","#sigma_0","Fano","A");
  fun3->FixParameter(1,0);
  fun3->FixParameter(2,2);
 
  TH1F *hd = h_corr ;
  hd->GetXaxis()->SetTitle("Charge (PE)");
 
  hd->Fit("fun3","R","",mins1, maxs1);
  hd->Fit("fun3","R","",mins1, maxs1);
  hd->Fit("fun3","R","",mins1, maxs1);
  fun3->ReleaseParameter(2);
  hd->Fit("fun3","R","",mins1, maxs1);
  hd->Fit("fun3","R","",mins1, maxs1);
  hd->Fit("fun3","R","",mins1, maxs1);
    
  hd->GetXaxis()->SetTitle("Charge (PE)");
  hd->GetXaxis()->SetTitleOffset(1.16);
  hd->GetYaxis()->SetTitle("Counts (arb.)");
  hd->GetYaxis()->SetTitleOffset(1.16);
  hd->GetXaxis()->SetRangeUser(0,1000); //(12,26)
  
  //hd->SetFillColor(kGreen);
  hd->Draw("PE");
  
  //h->SetLineColor(kRed);
  //h->Draw("same");
  
  
  //fitto e disegno S2 corretto
  
  TCanvas *c3 = new TCanvas("c3","c3",10,10,900,600);
/* 
  TF1 *mygaus_s2 = new TF1("mygaus_s2","gaus(0)",mins2,maxs2); // fitting range 620-800 singola fase; 400-600 doppia fase
        
  mygaus_s2->SetLineWidth(3);
  mygaus_s2->SetLineColor(kRed);
    
  h2_corr->Fit("mygaus_s2","EMRQ");
  
  h2_corr->GetXaxis()->SetTitle("Charge (PE)");
  h2_corr->GetXaxis()->SetTitleOffset(1.16);
  h2_corr->GetYaxis()->SetTitle("Counts (arb.)");
  h2_corr->GetYaxis()->SetTitleOffset(1.16);
  
  //h2_corr->SetFillColor(kBlue);
  h2_corr->Draw("PE");
*/

  TF1 *fun32 = new TF1("fun32",fam241,0,1000,4);
  fun32->SetParameters(200,2000,3.,2e4);
  fun32->SetLineColor(64);
  fun32->Draw();
  fun32->SetParNames("LY","#sigma_0","Fano","A");
  fun32->FixParameter(1,0);
  fun32->FixParameter(2,2);

  TH1F *hfits2 = h2_corr ;
  hfits2->GetXaxis()->SetTitle("Charge (PE)");
 
  hfits2->Fit("fun32","R","",mins2, maxs2);
  hfits2->Fit("fun32","R","",mins2, maxs2);
  hfits2->Fit("fun32","R","",mins2, maxs2);
  fun32->ReleaseParameter(2);
  hfits2->Fit("fun32","R","",mins2, maxs2);
  hfits2->Fit("fun32","R","",mins2, maxs2);
  hfits2->Fit("fun32","R","",mins2, maxs2);
    
  hfits2->GetXaxis()->SetTitle("Charge (PE)");
  hfits2->GetXaxis()->SetTitleOffset(1.16);
  hfits2->GetYaxis()->SetTitle("Counts (arb.)");
  hfits2->GetYaxis()->SetTitleOffset(1.16);
  //hfits2->GetXaxis()->SetRangeUser(0,1000); //(12,26)
  
  hfits2->Draw("PE");
   
  //h2->SetLineColor(kRed);
  //h2->Draw("same");
  
  
/*  
  TCanvas *c5 = new TCanvas("c5","c5",10,10,900,600);
  c5->cd();
  
  h_time->SetTitle("Drift Time distribution");
  h_time->GetXaxis()->SetTitle("t_{drift} (ns)");
  h_time->GetXaxis()->SetTitleOffset(1.16);
  h_time->GetYaxis()->SetTitle("Counts (arb.)");
  h_time->GetYaxis()->SetTitleOffset(1.16);
  
  h_time->Draw();
 */

////////////// GROSS S1 LIGHT YIELD ////////////////////

    double energy = 59.5409; //g-line of 241Am (keV)
    double sigma = 0;
    
    sigma = sqrt(fun3->GetParameter(2)*(fun3->GetParameter(0)*energy));

    cout << "mu = " << fun3->GetParameter(0)*energy << " (PE)" << endl;
    
    cout << "S1LY = " << fun3->GetParameter(0) << " +/- " << fun3->GetParError(0) << " (PE/keV)" << endl;
    
    cout << "sigma = " << sigma << endl;
    
    cout << "Fano = " << fun3->GetParameter(2) << endl;
    
    cout << "sigma/mu = " << sigma/(fun3->GetParameter(0)*energy) << " " << "sigma/sqrt(mu) = " << sigma/sqrt((fun3->GetParameter(0))*energy) << endl;
    
    cout << "\n" << endl;
    
////////////// GROSS S2 LIGHT YIELD ////////////////////
/*
    double mu = mygaus_s2->GetParameter(1);   
    double S2LY = mu/energy; 
    
    cout << "S2 GROSS LY: " << S2LY <<  " +/- " << ((mygaus_s2->GetParError(1))/(mygaus_s2->GetParameter(1)))+(1/energy) << " [PE/keV]" << endl;
    
    cout << "sigma = " << mygaus_s2->GetParameter(2) << " " << "mu = " << mygaus_s2->GetParameter(1) << endl;
    
    cout << "sigma/mu = " << mygaus_s2->GetParameter(2)/mygaus_s2->GetParameter(1) << " " << "sigma/sqrt(mu) = " << mygaus_s2->GetParameter(2)/sqrt(mygaus_s2->GetParameter(1)) << endl;
    
    cout << "\n" << endl;
*/

    double sigma2 = 0;

    sigma2 = sqrt(fun32->GetParameter(2)*(fun32->GetParameter(0)*energy));

    cout << "mu = " << fun32->GetParameter(0)*energy << " (PE)" << endl;
    
    cout << "S2LY = " << fun32->GetParameter(0) << " +/- " << fun32->GetParError(0) << " (PE/keV)" << endl;
    
    cout << "sigma = " << sigma2 << endl;
    
    cout << "Fano = " << fun32->GetParameter(2) << endl;
    
    cout << "sigma/mu = " << sigma2/(fun32->GetParameter(0)*energy) << " " << "sigma/sqrt(mu) = " << sigma2/sqrt((fun32->GetParameter(0))*energy) << endl;
    
    cout << "\n" << endl;
    
////////////// SAVE ////////////////////

/*   
    ofstream ofile;
   
    ofile.open ("ly_top_bot_mod_corr.cfg", ofstream::out | ofstream::app);
    
    ofile << run << " " << S1LY <<  " " << ((fit_gaus->GetParError(1))/(fit_gaus->GetParameter(1)))+(1/energy) << " " << S1LY_top <<  " " << ((fit_gaus_top->GetParError(1))/(fit_gaus_top->GetParameter(1)))+(1/energy) << " " << S2LY <<  " " << ((fit_land->GetParError(1))/(fit_land->GetParameter(1)))+(1/energy) << " " << S2LY_top <<  " " << ((fit_land_top->GetParError(1))/(fit_land_top->GetParameter(1)))+(1/energy) << endl;
    
    ofile.close();
    
*/
/*   
    c1->cd(1);
    c2->SaveAs(Form("plot/naples/july19/ph2/s1_tdrift_%d.png", run));
    
    c1->cd(2);
    c2->SaveAs(Form("plot/naples/july19/ph2/s2_s1_tdrift_%d.png", run));   
   
    c2->cd();
    c2->SaveAs(Form("plot/naples/july19/ph2/s1_corr_%d.png", run));
    
    c3->cd();
    c3->SaveAs(Form("plot/naples/july19/ph2/s2_corr_%d.png", run)); 
*/
  
}
      
      
      
      
      
      
      
