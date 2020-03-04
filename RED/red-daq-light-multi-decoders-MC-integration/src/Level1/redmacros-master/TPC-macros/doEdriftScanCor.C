#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TROOT.h>
//#include <TPCAnalyzer.h>
#include <TCanvas.h>

void doEdriftScan () {
  TFile *_sf0    = new TFile("reco/naples/camp_VIII/run_1137.root");
  TFile *_sf50   = new TFile("reco/naples/camp_VIII/run_1141.root");
  TFile *_sf100  = new TFile("reco/naples/camp_VIII/run_1144.root");
  TFile *_sf200  = new TFile("reco/naples/camp_VIII/run_1145.root");
  TFile *_sf300  = new TFile("reco/naples/camp_VIII/run_1148.root");
  TFile *_sf400  = new TFile("reco/naples/camp_VIII/run_1151.root");
  TFile *_sf500  = new TFile("reco/naples/camp_VIII/run_1152.root");
  TFile *_sf700  = new TFile("reco/naples/camp_VIII/run_1155.root");
  TFile *_sf1000 = new TFile("reco/naples/camp_VIII/run_1156.root");
  //  .L ./TPCAnalyzer.C++  
  //  gSystem->Load( "./TPCAnalyzer.C+");
  //  gROOT->ProcessLine(".L ./TPCAnalyzer.C++");  
  FILE * flysp = fopen("REDLY_SP.txt","w+");
  TCanvas * c1 = new TCanvas();
  TCanvas * csp = new TCanvas("csp","Single Phase");

  //SP null field
  cout <<"Processing SP null field"<<endl;
  TPCAnalyzer Enull;
  _sf0->cd();
  TTree * tnull; _sf0->GetObject("reco",tnull);
  c1->cd();
  tnull->Process(&Enull);
  csp->cd();
  Enull.hS1cor->SetLineColor(1);
  TH1 * pippo = Enull.hS1cor->DrawNormalized();
  pippo->GetYaxis()->SetRangeUser(0.,0.03);

  Double_t LYSP0 = Enull.LYFit[0];
  Double_t eLYSP0 = Enull.LYFit[1];
  fprintf(flysp,"0 1.0  0. 0.\n"); 
  
  //SP 50 V/cm
  cout <<"Processing SP 50 V/cm field"<<endl;
  TPCAnalyzer E50;
  _sf50->cd();
  TTree * t50; _sf50->GetObject("reco",t50);
  c1->cd();
  t50->Process(&E50);
  csp->cd();
  E50.hS1cor->SetLineColor(2);
  E50.hS1cor->DrawNormalized("same");
  cout <<E50.LYFit[0]<<" "<<E50.FanoFit[0]<<endl;
  fprintf(flysp,"50 %f %f %f\n",E50.LYFit[0]/LYSP0, 50*0.01, 
  	  sqrt(pow(E50.LYFit[1]/E50.LYFit[0],2) +
  	       pow(eLYSP0/LYSP0,2))*E50.LYFit[0]/LYSP0); 

  //SP 100 V/cm
  cout <<"Processing SP 100 V/cm field"<<endl;
  TPCAnalyzer E100;
   _sf100->cd();
  TTree * t100; _sf100->GetObject("reco",t100);
  c1->cd();
  t100->Process(&E100);
  csp->cd();
  E100.hS1cor->SetLineColor(3);
  E100.hS1cor->DrawNormalized("same");
  fprintf(flysp,"100 %f %f %f\n",E100.LYFit[0]/LYSP0, 100*0.01, 
  	  sqrt(pow(E100.LYFit[1]/E100.LYFit[0],2) +
  	       pow(eLYSP0/LYSP0,2))*E100.LYFit[0]/LYSP0); 


  //SP 200 V/cm
  cout <<"Processing SP 200 V/cm field"<<endl;
  TPCAnalyzer E200;
   _sf200->cd();
  TTree * t200; _sf200->GetObject("reco",t200);
  c1->cd();
  t200->Process(&E200);
  csp->cd();
  E200.hS1cor->SetLineColor(4);
  E200.hS1cor->DrawNormalized("same");
  fprintf(flysp,"200 %f %f %f\n",E200.LYFit[0]/LYSP0, 200*0.01, 
  	  sqrt(pow(E200.LYFit[1]/E200.LYFit[0],2) +
  	       pow(eLYSP0/LYSP0,2))*E200.LYFit[0]/LYSP0); 


  //SP 300 V/cm
  cout <<"Processing SP 300 V/cm field"<<endl;
  TPCAnalyzer E300;
  _sf300->cd();
  TTree * t300; _sf300->GetObject("reco",t300);
  c1->cd();
  t300->Process(&E300);
  csp->cd();
  E300.hS1cor->SetLineColor(5);
  E300.hS1cor->DrawNormalized("same");
  fprintf(flysp,"300 %f %f %f\n",E300.LYFit[0]/LYSP0, 300*0.01, 
  	  sqrt(pow(E300.LYFit[1]/E300.LYFit[0],2) +
  	       pow(eLYSP0/LYSP0,2))*E300.LYFit[0]/LYSP0); 


  //SP 400 V/cm
  cout <<"Processing SP 400 V/cm field"<<endl;
  TPCAnalyzer E400;
  _sf400->cd();
  TTree * t400; _sf400->GetObject("reco",t400);
  c1->cd();
  t400->Process(&E400);
  csp->cd();
  E400.hS1cor->SetLineColor(6);
  E400.hS1cor->DrawNormalized("same");
  fprintf(flysp,"400 %f %f %f\n",E400.LYFit[0]/LYSP0, 400*0.01, 
  	  sqrt(pow(E400.LYFit[1]/E400.LYFit[0],2) +
  	       pow(eLYSP0/LYSP0,2))*E400.LYFit[0]/LYSP0); 


  //SP 500 V/cm
  cout <<"Processing SP 500 V/cm field"<<endl;
  TPCAnalyzer E500;
  _sf500->cd();
  TTree * t500; _sf500->GetObject("reco",t500);
  c1->cd();
  t500->Process(&E500);
  csp->cd();
  E500.hS1cor->SetLineColor(7);
  E500.hS1cor->DrawNormalized("same");
  fprintf(flysp,"500 %f %f %f\n",E500.LYFit[0]/LYSP0, 500*0.01, 
  	  sqrt(pow(E500.LYFit[1]/E500.LYFit[0],2) +
  	       pow(eLYSP0/LYSP0,2))*E500.LYFit[0]/LYSP0); 


  //SP 700 V/cm
  cout <<"Processing SP 700 V/cm field"<<endl;
  TPCAnalyzer E700;
  _sf700->cd();
  TTree * t700; _sf700->GetObject("reco",t700);
  c1->cd();
  t700->Process(&E700);
  csp->cd();
  E700.hS1cor->SetLineColor(8);
  E700.hS1cor->DrawNormalized("same");
  fprintf(flysp,"700 %f %f %f\n",E700.LYFit[0]/LYSP0, 700*0.01, 
  	  sqrt(pow(E700.LYFit[1]/E700.LYFit[0],2) +
  	       pow(eLYSP0/LYSP0,2))*E700.LYFit[0]/LYSP0); 


  //SP 1000 V/cm
  cout <<"Processing SP 1000 V/cm field"<<endl;
  TPCAnalyzer E1000;
  _sf1000->cd();
  TTree * t1000; _sf1000->GetObject("reco",t1000);
  c1->cd();
  t1000->Process(&E1000);
  csp->cd();
  E1000.hS1cor->SetLineColor(9);
  E1000.hS1cor->DrawNormalized("same");
  fprintf(flysp,"1000 %f %f %f\n",E1000.LYFit[0]/LYSP0, 1000*0.01, 
  	  sqrt(pow(E1000.LYFit[1]/E1000.LYFit[0],2) +
  	       pow(eLYSP0/LYSP0,2))*E1000.LYFit[0]/LYSP0); 

  fclose(flysp);


  TFile *_df0    = new TFile("reco/naples/camp_VIII/run_1159.root");
  TFile *_df100  = new TFile("reco/naples/camp_VIII/run_1164.root");
  TFile *_df200  = new TFile("reco/naples/camp_VIII/run_1135.root"); //1128 bkg //1115 Kr+Am
  TFile *_df400  = new TFile("reco/naples/camp_VIII/run_1165.root");
  TFile *_df700  = new TFile("reco/naples/camp_VIII/run_1160.root");
  TFile *_df1000 = new TFile("reco/naples/camp_VIII/run_1161.root");

  FILE * flydp = fopen("REDLY_DP.txt","w+");
  TCanvas * c2 = new TCanvas();
  TCanvas * cdp = new TCanvas("cdp","Double Phase");

  //DP null field
  cout <<"Processing DP null field"<<endl;
  TPCAnalyzer Enull_dp;
  _df0->cd();
  TTree * tnull_dp; _df0->GetObject("reco",tnull_dp);
  c2->cd();
  Enull_dp.isDoublePhase=0; Enull_dp.isDoublePhaseDAQ=1;
  tnull_dp->Process(&Enull_dp);
  cdp->cd();
  Enull_dp.hS1cor->SetLineColor(1);
  TH1 * pippo2 = Enull_dp.hS1cor->DrawNormalized();
  pippo2->GetYaxis()->SetRangeUser(0.,0.03);

  Double_t LYDP0 = Enull_dp.LYFit[0];
  Double_t eLYDP0 = Enull_dp.LYFit[1];
  fprintf(flydp,"0 1.0  0. 0.\n"); 
  
  //DP 100 V/cm
  cout <<"Processing DP 100 V/cm field"<<endl;
  TPCAnalyzer E100_dp;
  _df100->cd();
  TTree * t100_dp; _df100->GetObject("reco",t100_dp);
  c2->cd();
  E100_dp.isDoublePhase=1;  E100_dp.isDoublePhaseDAQ=2;
  t100_dp->Process(&E100_dp);
  cdp->cd();
  E100_dp.hS1cor->SetLineColor(3);
  E100_dp.hS1cor->DrawNormalized("same");
  cout <<E100_dp.LYFit[0]<<" "<<E100_dp.FanoFit[0]<<endl;
  fprintf(flydp,"100 %f %f %f\n",E100_dp.LYFit[0]/LYDP0, 100*0.01, 
	  sqrt(pow(E100_dp.LYFit[1]/E100_dp.LYFit[0],2) +
	       pow(eLYDP0/LYDP0,2))*E100_dp.LYFit[0]/LYDP0); 

  //DP 200 V/cm
  cout <<"Processing DP 200 V/cm field"<<endl;
  TPCAnalyzer E200_dp;
  _df200->cd();
  TTree * t200_dp; _df200->GetObject("reco",t200_dp);
  c2->cd();
  E200_dp.isDoublePhase=1;
  t200_dp->Process(&E200_dp);
  cdp->cd();
  E200_dp.hS1cor->SetLineColor(4);
  E200_dp.hS1cor->DrawNormalized("same");
  cout <<E200_dp.LYFit[0]<<" "<<E200_dp.FanoFit[0]<<endl;
  fprintf(flydp,"200 %f %f %f\n",E200_dp.LYFit[0]/LYDP0, 200*0.01, 
	  sqrt(pow(E200_dp.LYFit[1]/E200_dp.LYFit[0],2) +
	       pow(eLYDP0/LYDP0,2))*E200_dp.LYFit[0]/LYDP0); 

  //DP 400 V/cm
  cout <<"Processing DP 400 V/cm field"<<endl;
  TPCAnalyzer E400_dp;
  _df400->cd();
  TTree * t400_dp; _df400->GetObject("reco",t400_dp);
  c2->cd();
  E400_dp.isDoublePhase=1;
  t400_dp->Process(&E400_dp);
  cdp->cd();
  E400_dp.hS1cor->SetLineColor(6);
  E400_dp.hS1cor->DrawNormalized("same");
  cout <<E400_dp.LYFit[0]<<" "<<E400_dp.FanoFit[0]<<endl;
  fprintf(flydp,"400 %f %f %f\n",E400_dp.LYFit[0]/LYDP0, 400*0.01, 
	  sqrt(pow(E400_dp.LYFit[1]/E400_dp.LYFit[0],2) +
	       pow(eLYDP0/LYDP0,2))*E400_dp.LYFit[0]/LYDP0); 
 
  //DP 700 V/cm
  cout <<"Processing DP 700 V/cm field"<<endl;
  TPCAnalyzer E700_dp;
  _df700->cd();
  TTree * t700_dp; _df700->GetObject("reco",t700_dp);
  c2->cd();
  E700_dp.isDoublePhase=1;
  t700_dp->Process(&E700_dp);
  cdp->cd();
  E700_dp.hS1cor->SetLineColor(8);
  E700_dp.hS1cor->DrawNormalized("same");
  cout <<E700_dp.LYFit[0]<<" "<<E700_dp.FanoFit[0]<<endl;
  fprintf(flydp,"700 %f %f %f\n",E700_dp.LYFit[0]/LYDP0, 700*0.01, 
	  sqrt(pow(E700_dp.LYFit[1]/E700_dp.LYFit[0],2) +
	       pow(eLYDP0/LYDP0,2))*E700_dp.LYFit[0]/LYDP0); 

  //DP 1000 V/cm
  cout <<"Processing DP 1000 V/cm field"<<endl;
  TPCAnalyzer E1000_dp;
  _df1000->cd();
  TTree * t1000_dp; _df1000->GetObject("reco",t1000_dp);
  c2->cd();
  E1000_dp.isDoublePhase=1;
  t1000_dp->Process(&E1000_dp);
  cdp->cd();
  E1000_dp.hS1cor->SetLineColor(9);
  E1000_dp.hS1cor->DrawNormalized("same");
  cout <<E1000_dp.LYFit[0]<<" "<<E1000_dp.FanoFit[0]<<endl;
  fprintf(flydp,"1000 %f %f %f\n",E1000_dp.LYFit[0]/LYDP0, 1000*0.01, 
	  sqrt(pow(E1000_dp.LYFit[1]/E1000_dp.LYFit[0],2) +
	       pow(eLYDP0/LYDP0,2))*E1000_dp.LYFit[0]/LYDP0); 

  fclose(flydp);



}
