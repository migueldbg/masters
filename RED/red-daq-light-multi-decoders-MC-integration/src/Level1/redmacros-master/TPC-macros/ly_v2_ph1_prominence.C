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

void ly_v2_ph1_prominence(TString filename, int run, int Ed=200, 
			  double pmin=0.35, double pmax=0.55, int mins1=400, 
			  int maxs1=950)
{

  TFile *f = new TFile(filename,"read");
    
  //Titles
  TString title = "^{241}Am Spectrum corrected (P)";
  TString titleq = "^{241}Am Spectrum corrected (P)";
  TString title_tba_1 = "S1 vs TBA (prominence)";
  TString title_tba_q1 = "S1 vs TBA (charge)";
  TString title1 = "^{241}Am Spectrum";

  //Canvases and Style
  TCanvas *c1 = new TCanvas("c1","c1",10,10,900,600);
  TCanvas *c1q = new TCanvas("c1-charge","c1-charge",10,10,900,600);
  TCanvas *c2 = new TCanvas("c2","c2",10,10,900,600);
  TCanvas *c2q = new TCanvas("c2-charge","c2-charge",10,10,900,600);

  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(111);
 
  //Prominence Branches
  vector<Double_t> *clp_type = new vector<Double_t>;
  vector<Double_t> *clp_p = new vector<Double_t>;
  vector<Double_t> *clp_tba = new vector<Double_t>;
  vector<Double_t> *clp_f90 = new vector<Double_t>;
  vector<Double_t> *clp_startt = new vector<Double_t>;
  vector<Double_t> *clp_endt = new vector<Double_t>;
  Int_t nclp;

  //Charge branches
  vector<Double_t> *clq_type = new vector<Double_t>;
  vector<Double_t> *clq_q = new vector<Double_t>;
  vector<Double_t> *clq_tba = new vector<Double_t>;
  vector<Double_t> *clq_f90 = new vector<Double_t>;
  vector<Double_t> *clq_startt = new vector<Double_t>;
  vector<Double_t> *clq_endt = new vector<Double_t>;
  Int_t nclq;


  //let's start 
  TTree* reco = (TTree*) f->Get("dstree");
  reco->SetBranchAddress("nclp",&nclp);
  reco->SetBranchAddress("clp_type",&clp_type);
  reco->SetBranchAddress("clp_p",&clp_p);
  reco->SetBranchAddress("clp_tba",&clp_tba);
  reco->SetBranchAddress("clp_f90",&clp_f90);
  reco->SetBranchAddress("clp_startt",&clp_startt);
  reco->SetBranchAddress("clp_endt",&clp_endt);
  reco->SetBranchAddress("nclq",&nclq);
  reco->SetBranchAddress("clq_type",&clq_type);
  reco->SetBranchAddress("clq_q",&clq_q);
  reco->SetBranchAddress("clq_tba",&clq_tba);
  reco->SetBranchAddress("clq_f90",&clq_f90);
  reco->SetBranchAddress("clq_startt",&clq_startt);
  reco->SetBranchAddress("clq_endt",&clq_endt);
  cout << "reco entries: " << reco->GetEntries() << endl;
  
  //S1 vs TBA and S1_corr
  TH2F *h_tba_1 = new TH2F("h_tba_1",title_tba_1,300,0.2,0.7,300,0,2000);
  TH2F *h_tba_q1 = new TH2F("h_tba_q1",title_tba_q1,300,0.3,0.9,300,0,2000);
  TH2F *h_tba_2 = new TH2F("h_tba_2","S1_{corr} vs TBA",300,0.2,0.7,300,0,2000);
  TH2F *h_tba_q2 = new TH2F("h_tba_q2","S1_{corr} vs TBA (Q)",300,0.3,0.9,300,0,2000);
  TH2F *hcorrelation = new TH2F("hcorrelation","S1 (P) vs. S1 (Q)",300,0,2000,300,0,2000);
  TH1F *h_s1 = new TH1F("h_s1",title1,300,0,2000);
  TH1F *h_s1q = new TH1F("h_s1q","S1 (charge)",300,0,2000);
  TH1F *h_s1_corr = new TH1F("h_s1_corr",title,300,0,2000); 
  TH1F *h_s1q_corr = new TH1F("h_s1q_corr",titleq,300,0,2000); 
  
  //Correzione
  TH1F *h_corr = new TH1F("h_corr","",300,0,2);
   
  //S1 vs TBA
 
  for (int i=0;i<reco->GetEntries();i++){
      	
	reco->GetEntry(i); 
	
	//cout << "Sono qui: " << i << " " << endl;
	//cout << clp_p->size() << endl;

	if (clp_p->size() != (size_t)nclp)
	  cout << "WARNING" << clp_p->size() << " " << nclp <<  endl;

	//	int nclusters = nclp;

	//Prominence 
	for (Int_t icl=0;icl<nclp;icl++)
	  {
	    double f90_1 = clp_f90->at(icl);
	    double pcharge = clp_p->at(icl);		
	    double tba = clp_tba->at(icl);
	    double t0 = clp_startt->at(icl);
	    double t1 = clp_endt->at(icl);
	    
	    //cout << clp_type->at(icl) << " " << t0 << " " << t1 << " " << charge << " " << tba << endl;
	    	  
	    if (clp_type->at(icl) == 1 && pcharge>400 && pcharge<600) //check for S1
	      h_tba_1->Fill(tba,pcharge);
	  }	

	//cout << nclq << " " << clq_q->size() << endl;
	//Charge 
	for (Int_t icl=0;icl<nclq;icl++)
	  {
	    double f90_1 = clq_f90->at(icl);
	    double qcharge = clq_q->at(icl);		
	    double tba = clq_tba->at(icl);
	    double t0 = clq_startt->at(icl);
	    double t1 = clq_endt->at(icl);
	    
	    //cout << clp_type->at(icl) << " " << t0 << " " << t1 << " " << charge << " " << tba << endl;
	    if (clq_type->at(icl) == 1 && qcharge>650 && qcharge<850) //check for S1
	      h_tba_q1->Fill(tba,qcharge);

	  }
	

  }
  
  c1->Divide(1,2);
  c1->cd(1);
      
  //Fit prominence
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
  p1->GetYaxis()->SetRangeUser(400,600);       
  //h_tba_1->Draw("colz");
  p1->Draw();

  c1q->Divide(1,2);
  c1q->cd(1);

  //Fit charge
  TF1 *fit_tba_q1 = new TF1("fit_tba_q1","pol2",-0.5,0.5);   
     
  TProfile *pq1 = h_tba_q1->ProfileX();
  pq1->SetMarkerStyle(21);
  pq1->SetMarkerSize(0.5);
  fit_tba_q1->SetLineColor(kRed);
  pq1->Fit("fit_tba_q1","EM","",pmin+0.15,pmax+0.15); //-0.35,0.035)    
  pq1->GetXaxis()->SetTitle("TBA");
  pq1->GetXaxis()->SetTitleOffset(1.16);
  pq1->GetYaxis()->SetTitle("Charge (PE)");
  pq1->GetYaxis()->SetTitleOffset(1.16);
  //h_tba_1->GetXaxis()->SetRangeUser(r_min_bottom,r_max_bottom); //(12,26)
  pq1->GetYaxis()->SetRangeUser(600,800);       
  //h_tba_1->Draw("colz");
  pq1->Draw();

 
  //calcolo correzioni S1 e riempio S1 corretto
  //double thalf_tba=0;
  double t_s1=0;
  double corr_tba=0;
  double s1_tba=0;
  
  for (int i=0;i<reco->GetEntries();i++){
      	
	reco->GetEntry(i); 

	int nclusters = nclp;

	for (Int_t icl=0;icl<nclp;icl++)
	  {
	    double f90_1 = clp_f90->at(icl);
	    double charge = clp_p->at(icl);		
	    double tba = clp_tba->at(icl);
	    double t0 = clp_startt->at(icl);
	    double t1 = clp_endt->at(icl);
	    
	    if (clp_type->at(icl) == 1)
	      {
		s1_tba = fit_tba_1->Eval(0.45);	  
		t_s1 = fit_tba_1->Eval(tba);
	    
		corr_tba = t_s1/s1_tba; //correzione S1
		//corr_tba = 1.; //NO CORRECTION
		//cout << "io sono la correzione: " << corr_tba << endl;
		
		Double_t S1 = charge/corr_tba;
		
		h_s1_corr->Fill(S1); 
		h_s1->Fill(charge);
		
		//if(Ed==200 && charge>540 && charge<740) 
		
		if (charge > 400 && charge < 600)
		  h_tba_2->Fill(tba,S1);	 
		
		
		h_corr->Fill(corr_tba);
	      }	   
	  }

	//Charge 
	for (Int_t icl=0;icl<nclq;icl++)
	  {
	    double f90_1 = clq_f90->at(icl);
	    double qcharge = clq_q->at(icl);		
	    double tba = clq_tba->at(icl);
	    double t0 = clq_startt->at(icl);
	    double t1 = clq_endt->at(icl);
	    
	    if (clq_type->at(icl) == 1)
	      {
		s1_tba = fit_tba_q1->Eval(0.6);	  
		t_s1 = fit_tba_q1->Eval(tba);
	    
		corr_tba = t_s1/s1_tba; //correzione S1
		//corr_tba = 1.; //NO CORRECTION
		//cout << "io sono la correzione: " << corr_tba << endl;
		
		Double_t S1 = qcharge/corr_tba;
		
		h_s1q_corr->Fill(S1); 
		h_s1q->Fill(qcharge);
		
		//if(Ed==200 && charge>540 && charge<740) 
		
		if (qcharge > 650 && qcharge < 850)
		  h_tba_q2->Fill(tba,S1);	 			      	      
	      }
	  }
	
	if (nclq == 1 && nclp == 1 &&
	    clq_type->at(0) == 1 && 
	    clp_type->at(0) == 1)
	  hcorrelation->Fill(clq_q->at(0),clp_p->at(0));



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
  p2->GetYaxis()->SetRangeUser(400,600);
  p2->DrawCopy();

  c1q->cd(2);
  TProfile *pq2 = h_tba_q2->ProfileX();
  pq2->SetMarkerStyle(21);
  pq2->SetMarkerSize(0.5);
  //p2->SetMarkerColor(kRed);
  pq2->GetXaxis()->SetTitle("TBA");
  pq2->GetXaxis()->SetTitleOffset(1.16);
  pq2->GetYaxis()->SetTitle("Charge_{corr} (PE)");
  pq2->GetYaxis()->SetTitleOffset(1.16);
  pq2->GetYaxis()->SetRangeUser(650,850);
  pq2->DrawCopy();

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
  c2q->cd();
  h_s1q_corr->Fit("fun3","R","",mins1, maxs1);
  h_s1q_corr->Fit("fun3","R","",mins1, maxs1);
  h_s1q_corr->Fit("fun3","R","",mins1, maxs1);
  fun3->ReleaseParameter(2);
  h_s1q_corr->Fit("fun3","R","",mins1, maxs1);
  h_s1q_corr->Fit("fun3","R","",mins1, maxs1);
  h_s1q_corr->Fit("fun3","R","",mins1, maxs1);
    
  h_s1q_corr->GetXaxis()->SetTitle("Charge (PE)");
  h_s1q_corr->GetXaxis()->SetTitleOffset(1.16);
  h_s1q_corr->GetYaxis()->SetTitle("Counts (arb.)");
  h_s1q_corr->GetYaxis()->SetTitleOffset(1.16);
  //h_s1_corr->GetXaxis()->SetRangeUser(0,1000); //(12,26)
  
  //h_s1_corr->SetFillColor(kGreen);
  h_s1q_corr->DrawCopy();
  


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
  h_s1q->SetLineColor(kGreen);  
  h_s1q->Draw("same");
  h_s1q_corr->SetLineColor(kBlue);  
  h_s1q_corr->GetFunction("fun3")->SetBit(TF1::kNotDraw);
  h_s1q_corr->Draw("same");
 
  
  TLegend *l1 = new TLegend(0.70,0.63,0.84,0.83,NULL,"brNDC");
  l1->AddEntry("h_s1","S1 not corrected (P)","l");
  l1->AddEntry("h_s1_corr","Corrected S1 (P)","l");
  l1->AddEntry("h_s1q","S1 not corrected (Q)","l");
  l1->AddEntry("h_s1q_corr","Corrected S1 (Q)","l");
  l1->SetLineWidth(0);
  l1->Draw();

  TCanvas *c5 = new TCanvas("c5","c5",10,10,900,600);
  c5->cd();
  hcorrelation->GetXaxis()->SetTitle("S1 (charge)");
  hcorrelation->GetYaxis()->SetTitle("S1 (prominence)");
  hcorrelation->Draw("COLZ");

  
  TCanvas *c6 = new TCanvas("c6","c6",10,10,900,600);
  c6->cd();
  
  h_s1_corr->GetXaxis()->SetTitle("Charge (PE)");
  h_s1_corr->GetXaxis()->SetTitleOffset(1.16);
  h_s1_corr->GetYaxis()->SetTitle("Counts (arb.)");
  h_s1_corr->GetYaxis()->SetTitleOffset(1.16);
  h_s1_corr->SetLineColor(kBlack);
  h_s1_corr->SetFillColor(0);
  
  h_s1_corr->GetFunction("fun3")->SetBit(TF1::kNotDraw);
  h_s1_corr->Draw();
  h_s1q_corr->SetLineColor(kBlue);  
  h_s1q_corr->GetFunction("fun3")->SetBit(TF1::kNotDraw);
  h_s1q_corr->Draw("same");

  TLegend *l2 = new TLegend(0.70,0.63,0.84,0.83,NULL,"brNDC");
  l2->AddEntry("h_s1_corr","Prominence (filtered N_{phe})","l");
  l2->AddEntry("h_s1q_corr","Charge","l");
  l2->SetLineWidth(0);
  l2->Draw();

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
      
      
      
      
      
      
      
