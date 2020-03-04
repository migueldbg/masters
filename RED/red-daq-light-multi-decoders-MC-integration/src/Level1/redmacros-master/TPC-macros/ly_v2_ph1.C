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
  This macro is meant to extract informations about S1 for the single-phase 
  runs taken with Am241 source.
  
  It accepts the name of the root file as input together with some other
  parameters such as, in order: the run number (used in case of an autmatic
  save of spectra), the value of drift field and fitting ranges for plots.

  S1 is corrected starting from the s1 vs TBA distribution. An order two polinomial
  is used as fitting function, and then an evaluation of the TBA in the middle of
  the chamber is done. Finally, S1 is normalized by the ratio bewtween the above mentioned
  TBA over the total value. This is done per each value of drift field.  

  Finally, corrected S1 is fitted by using a Monte Carlo function provided by Davide Franco
  that accounts also for the gaussin smearing of the detector response.

  Error bars (when available) are purely statistical.
  
  It works under ROOT, as:
  .L ly_v2_ph1.C+
  ly_v2_ph1("root_filename",run_number,200,-0.38,-0.05,400,1000)
  in case of a 200 V/cm run taken in Catania.
  
  Last edit by Simone S. on 31 Jan. 2020

 */

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

void ly_v2_ph1(TString filename, int run, int Ed=200, 
	       double pmin=-0.35, double pmax=0.035, int mins1=400, 
	       int maxs1=950)
{

  TFile *f = new TFile(filename,"read");
    
//Titles
  TString title = "^{241}Am Spectrum corrected";
  TString title_tba_1 = "S1 vs TBA";
  TString title1 = "^{241}Am Spectrum";

//Canvases and Style
  TCanvas *c1 = new TCanvas("c1","c1",10,10,900,600);
  TCanvas *c2 = new TCanvas("c2","c2",10,10,900,600);
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(111);
 
  //let's start 
  TTree* reco = (TTree*) f->Get("reco");
  EvRec0* evReco = new EvRec0();
  reco->SetBranchAddress("recoevent",&evReco);
  cout << "reco entries: " << reco->GetEntries() << endl;
  
  //S1 vs TBA and S1_corr
  TH2F *h_tba_1 = new TH2F("h_tba_1",title_tba_1,300,-0.5,0.25,300,0,2000);
  TH2F *h_tba_2 = new TH2F("h_tba_2","S1_{corr} vs TBA",300,-0.5,0.25,300,0,2000);
  TH1F *h_s1 = new TH1F("h_s1",title1,300,0,2000);
  TH1F *h_s1_corr = new TH1F("h_s1_corr",title,300,0,2000); 
  
  //Correzione
  TH1F *h_corr = new TH1F("h_corr","",300,0,2);

  TH1D *s2s1source = new TH1D("s1s2source","S2/S1, source",150.,0,40);
   
  //S1 vs TBA
 
  for (int i=0;i<reco->GetEntries();i++){
      	
	reco->GetEntry(i); 
		
     
      	vector<double> f90 = evReco->GetF90();
      	vector<double> start_time = evReco->GetStartTime();
      	vector<double> charge = evReco->GetCharge();
      	vector<RDCluster *> clusters = evReco->GetClusters();
      	int nclusters = evReco->GetNClusters();
      
      	if (nclusters==1){
      
        	double f90_1 = clusters.at(0)->f90;
		
		double rep_1 = clusters.at(0)->rep;
      
		double charge_top = clusters.at(0)->tot_charge_top;
		double charge_bot = clusters.at(0)->tot_charge_bottom;
		double charge = clusters.at(0)->charge;
		
		double tba = (charge_top-charge_bot)/(charge_top+charge_bot);
		
		if(Ed==0 && f90_1>0.2 && charge>650 && charge<850 && rep_1==1) h_tba_1->Fill(tba,charge);
		
		if(Ed==50 && f90_1>0.2 && charge>600 && charge<800 && rep_1==1) h_tba_1->Fill(tba,charge);
		
		if(Ed==100 && f90_1>0.2 && charge>580 && charge<780 && rep_1==1) h_tba_1->Fill(tba,charge);

		if(Ed==200 && f90_1>0.2 && charge>540 && charge<740 && rep_1==1) h_tba_1->Fill(tba,charge);		

		if(Ed==300 && f90_1>0.2 && charge>500 && charge<700 && rep_1==1) h_tba_1->Fill(tba,charge);
	
		if(Ed==400 && f90_1>0.2 && charge>450 && charge<680 && rep_1==1) h_tba_1->Fill(tba,charge);
		
		if(Ed==500 && f90_1>0.2 && charge>450 && charge<660 && rep_1==1) h_tba_1->Fill(tba,charge);
      
		if(Ed==700 && f90_1>0.2 && charge>400 && charge<600 && rep_1==1) h_tba_1->Fill(tba,charge);
      
		if(Ed==1000 && f90_1>0.2 && charge>370 && charge<570 && rep_1==1) h_tba_1->Fill(tba,charge);
				        
      }
      
  }
  
  c1->Divide(1,2);
  c1->cd(1);
      
  TF1 *fit_tba_1 = new TF1("fit_tba_1","pol2",-0.5,0.5);   
     
  TProfile *p1 = h_tba_1->ProfileX();
  p1->SetMarkerStyle(21);
  p1->SetMarkerSize(0.5);
  fit_tba_1->SetLineColor(kRed);
  p1->Fit("fit_tba_1","EM","",pmin,pmax); //-0.35,0.035)
    
  p1->GetXaxis()->SetTitle("TBA");
  p1->GetXaxis()->SetTitleOffset(1.16);
  p1->GetYaxis()->SetTitle("Charge (PE)");
  p1->GetYaxis()->SetTitleOffset(1.16);
  //h_tba_1->GetXaxis()->SetRangeUser(r_min_bottom,r_max_bottom); //(12,26)
  
  if(Ed==400 || Ed==500 || Ed==700 || Ed==1000){
  	p1->GetYaxis()->SetRangeUser(400,600);
	
  }else{
	
	p1->GetYaxis()->SetRangeUser(500,800);
	
  }
    
  //h_tba_1->Draw("colz");
  p1->Draw();
  
  //calcolo correzioni S1 e riempio S1 corretto
  //double thalf_tba=0;
  double t_s1=0;
  double corr_tba=0;
  double s1_tba=0;
  
  for (int i=0;i<reco->GetEntries();i++){
      	
	reco->GetEntry(i); 
		
     
      	vector<double> f90 = evReco->GetF90();
      	vector<double> start_time = evReco->GetStartTime();
      	vector<double> charge = evReco->GetCharge();
      	vector<RDCluster *> clusters = evReco->GetClusters();
      	int nclusters = evReco->GetNClusters();
      
      	if (nclusters==2){
      
        	double f90_1 = clusters.at(0)->f90;
		
		double rep_1 = clusters.at(0)->rep;
      
		double charge_top = clusters.at(0)->tot_charge_top;
		double charge_bot = clusters.at(0)->tot_charge_bottom;
		double charge = clusters.at(0)->charge;
		
		double tba = (charge_top-charge_bot)/(charge_top+charge_bot);
		
		s1_tba = fit_tba_1->Eval(-0.15);

		t_s1 = fit_tba_1->Eval(tba);
		
		//cout << "io sono t_tba: " << t_s1 << endl;
		//cout << "io sono thalf_tba: " << thalf_tba << endl;
		
              
		corr_tba = t_s1/s1_tba; //correzione S1
		//corr_tba = 1.; //NO CORRECTION
		//cout << "io sono la correzione: " << corr_tba << endl;
		
		if(f90_1>0.2 && rep_1==1)
		  {
		    
		    Double_t S1 = charge/corr_tba;
		    
		    h_s1_corr->Fill(S1); 
		    h_s1->Fill(charge);
		    
		    if(Ed==200 && charge>540 && charge<740) 
		      h_tba_2->Fill(tba,S1);	 
		    if (clusters.at(1)->f90 < 0.2 && clusters.at(1)->rep == 1 &&
			S1 > 540 && S1 < 750)
		      {
			Double_t S2 = clusters.at(1)->charge;
			s2s1source->Fill(S2/S1);
		      }
		  }
		
		h_corr->Fill(corr_tba);
	
		    

      }
      
  }
  

  //disegno e fitto S1 corretto
  
  c1->cd(2);
  TProfile *p2 = h_tba_2->ProfileX();
  p2->SetMarkerStyle(21);
  p2->SetMarkerSize(0.5);
  //p2->SetMarkerColor(kRed);
  p2->GetXaxis()->SetTitle("TBA");
  p2->GetXaxis()->SetTitleOffset(1.16);
  p2->GetYaxis()->SetTitle("Charge_{corr} (PE)");
  p2->GetYaxis()->SetTitleOffset(1.16);
  p2->GetYaxis()->SetRangeUser(500,800);
  p2->DrawCopy();

  //return; 
  c2->cd();
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
 
  h_s1_corr->GetXaxis()->SetTitle("Charge (PE)");
 
  h_s1_corr->Fit("fun3","R","",mins1, maxs1);
  h_s1_corr->Fit("fun3","R","",mins1, maxs1);
  h_s1_corr->Fit("fun3","R","",mins1, maxs1);
  fun3->ReleaseParameter(2);
  h_s1_corr->Fit("fun3","R","",mins1, maxs1);
  h_s1_corr->Fit("fun3","R","",mins1, maxs1);
  h_s1_corr->Fit("fun3","R","",mins1, maxs1);
    
  h_s1_corr->GetXaxis()->SetTitle("Charge (PE)");
  h_s1_corr->GetXaxis()->SetTitleOffset(1.16);
  h_s1_corr->GetYaxis()->SetTitle("Counts (arb.)");
  h_s1_corr->GetYaxis()->SetTitleOffset(1.16);
  //h_s1_corr->GetXaxis()->SetRangeUser(0,1000); //(12,26)
  
  //h_s1_corr->SetFillColor(kGreen);
  h_s1_corr->DrawCopy();
  
  //h->SetLineColor(kRed);
  //h->Draw("same");
 
  
  TCanvas *c3 = new TCanvas("c3","c3",10,10,900,600);
  c3->cd();
  
  h_corr->SetTitle("TBA correction distribution");
  h_corr->GetXaxis()->SetTitle("Correction");
  h_corr->GetXaxis()->SetTitleOffset(1.16);
  h_corr->GetYaxis()->SetTitle("Counts (arb.)");
  h_corr->GetYaxis()->SetTitleOffset(1.16);
  
  h_corr->Draw();
 
 
  TCanvas *c4 = new TCanvas("c4","c4",10,10,900,600);
  c4->cd();
  
  h_s1->GetXaxis()->SetTitle("Charge (PE)");
  h_s1->GetXaxis()->SetTitleOffset(1.16);
  h_s1->GetYaxis()->SetTitle("Counts (arb.)");
  h_s1->GetYaxis()->SetTitleOffset(1.16);
  h_s1->SetLineColor(kRed);
  
  h_s1_corr->SetLineColor(kBlack);
  h_s1_corr->SetFillColor(0);
  
  h_s1->Draw();
  h_s1_corr->GetFunction("fun3")->SetBit(TF1::kNotDraw);
  h_s1_corr->Draw("same");
  
  TLegend *l1 = new TLegend(0.70,0.63,0.84,0.83,NULL,"brNDC");
  l1->AddEntry("h_s1","S1 not corrected","l");
  l1->AddEntry("h_s1_corr","Corrected S1","l");
  l1->SetLineWidth(0);
  l1->Draw();
  
  TCanvas *c5 = new TCanvas("c4","c4",10,10,900,600);
  c5->cd();
  
  s2s1source->GetXaxis()->SetTitle("S2/S1");
  s2s1source->GetXaxis()->SetTitleOffset(1.16);
  s2s1source->GetYaxis()->SetTitle("Counts (arb.)");
  s2s1source->GetYaxis()->SetTitleOffset(1.16);

  s2s1source->Draw("HIST");

  /*
  TCanvas *c5 = new TCanvas("c5","c5",10,10,900,600);
  c5->cd();
  
  h_s1->GetXaxis()->SetTitle("Charge (PE)");
  h_s1->GetXaxis()->SetTitleOffset(1.16);
  h_s1->GetYaxis()->SetTitle("Counts (arb.)");
  h_s1->GetYaxis()->SetTitleOffset(1.16);
  h_s1->SetLineColor(kRed);
  
  h_s1->GetXaxis()->SetTitle("Charge (PE)");
 
  h_s1->Fit("fun3","R","",mins1, maxs1);
  h_s1->Fit("fun3","R","",mins1, maxs1);
  h_s1->Fit("fun3","R","",mins1, maxs1);
  fun3->ReleaseParameter(2);
  h_s1->Fit("fun3","R","",mins1, maxs1);
  h_s1->Fit("fun3","R","",mins1, maxs1);
  h_s1->Fit("fun3","R","",mins1, maxs1);
    
  h_s1->GetXaxis()->SetTitle("Charge (PE)");
  h_s1->GetXaxis()->SetTitleOffset(1.16);
  h_s1->GetYaxis()->SetTitle("Counts (arb.)");
  h_s1->GetYaxis()->SetTitleOffset(1.16);
  //h_s1->GetXaxis()->SetRangeUser(0,1000); //(12,26)
  
  fun3->Draw();
  h_s1->Draw();
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

    
////////////// SAVE ////////////////////
/*   
    ofstream ofile;
   
    ofile.open ("ly_top_bot_mod_corr.cfg", ofstream::out | ofstream::app);
    
    ofile << run << " " << S1LY <<  " " << ((fit_gaus->GetParError(1))/(fit_gaus->GetParameter(1)))+(1/energy) << " " << S1LY_top <<  " " << ((fit_gaus_top->GetParError(1))/(fit_gaus_top->GetParameter(1)))+(1/energy) << " " << S2LY <<  " " << ((fit_land->GetParError(1))/(fit_land->GetParameter(1)))+(1/energy) << " " << S2LY_top <<  " " << ((fit_land_top->GetParError(1))/(fit_land_top->GetParameter(1)))+(1/energy) << endl;
    
    ofile.close();
*/
 
/*    
    c1->cd();
    c1->SaveAs(Form("plot/naples/tba/ph1/s1_tba_run_%d.png", run));
    
    c2->cd();
    c2->SaveAs(Form("plot/naples/tba/ph1/s1_corrtba_run_%d.png", run));
    
//    c4->cd();
//    c4->SaveAs(Form("plot/naples/tba/ph1/s1_s1corr_run_%d.png", run));
  
    c5->cd();
    c5->SaveAs(Form("plot/naples/tba/ph1/s1_raw_run_%d.png", run));
 
*/ 
 
    
}
      
      
      
      
      
      
      
