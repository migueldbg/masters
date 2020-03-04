#include "TH2Poly.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCut.h"
#include "TCutG.h"
#include "TChain.h"
#include "TLatex.h"
#include <iostream>
#include "red-daq/EvRec0.hh"

/*
  This macro scans for data of runs 121-122 ("triples" of May test beam) and 
  looks for double coincidences (Si+LSci4) and for triple coincidences with 
  the PMT0 at 45 deg and with the wheel at 33 deg.

  The following selection is performed:
  - energy in LSci4 at least 50 keV
  - DeltaE/E compatible with the 7Be blob (graphical cut read from the external 
    BeCut.root file)
  - f90 in LSci4 compatible with a neutron
  - tof between Si and LSci4 by no more than 15 ns
  - at least 10 keV in the "scattering" LSci

*/


double GetQuenched(Double_t energy,Double_t kb);
double GetEnergy(Double_t charge,TString id);

//Structure to save calibration coefficient. They are mapped 
//by using the detector name as the key
struct calibrationdata
{
  Double_t a;
  Double_t b;
  Double_t kb;
};

int plotTripleMay()
{
  //The nclocks is 0 for run121 and -2 for run 122.
  TChain* reco = new TChain("reco","SiAndPMT4");
  reco->AddFile("run_121.root");
  reco->AddFile("run_122.root");
  
  TH2D* fBanana = new TH2D("h0",
			   "Full banana",
			   400,1500,7000,400,600,2000);
  TH2D* f4Banana = new TH2D("h1",
			    "Coincident banana PMT4",400,1500,7000,400,600,2000);
  f4Banana->SetMarkerStyle(1);
  f4Banana->SetMarkerStyle(kRed);
  TH2D* f0Banana = new TH2D("h2",
			    "Coincident banana PMT0",400,1500,7000,400,600,2000);
  f0Banana->SetMarkerStyle(1);
  f0Banana->SetMarkerColor(kYellow);

  TH2D* fOthersBanana = new TH2D("h3",
				 "Coincident banana others",400,1500,7000,400,600,2000);
  fOthersBanana->SetMarkerStyle(1);
  fOthersBanana->SetMarkerColor(kGreen);


  TH1D* a1 = new TH1D("a1",
		      "All events in PMT4",
		      50,0,3000);
  TH1D* a2 = new TH1D("a2",
		      "Triples PMT4-PMT0",
		      25,0,3000);
  TH1D* a3 = new TH1D("a3",
		      "Triples PMT4-others",
		      50,0,3000);

 
  gStyle->SetOptStat(0);
 
 
  TH1D* dt1 = new TH1D("dt1","DeltaT PMT4-Si",300,-100,200);
  TH1D* dt2 = new TH1D("dt2","DeltaT PMT4-PMT0",300,-100,200);
  TH1D* dt3 = new TH1D("dt3","DeltaT PMT4-wheel",300,-100,200);
  TH1D* dt4 = new TH1D("dt4","DeltaT PMT4-wheel",300,-40,60);

  //
  TString cutname = "be"; //(runNb == 112 || runNb==115) ? 

  //Open cut file
  TFile* fcut = new TFile("BeCut.root");
  TCutG* cut = nullptr;
  if (fcut->IsOpen())    
    { 
      cut = static_cast<TCutG*>(fcut->Get(cutname));
      fcut->Close();
    }


  //Open data file
  TFile f("run_121.root");
  vector<TString> *chanIDs = new vector<TString>;
  //Read metadata
  TTree* metaevent = (TTree*) f.Get("metaevent");
  metaevent->SetBranchAddress("chanID",&chanIDs);
  metaevent->GetEntry(0);

  EvRec0* evReco = new EvRec0();
  reco->SetBranchAddress("recoevent",&evReco);
  
  Int_t becounts = 0;
  Int_t counter1 = 0;
  Int_t counter2 = 0;
  Int_t counter3 = 0;

  //Event loop
  for (Int_t iloop=0;iloop<reco->GetEntries();iloop++)
    { 
      Int_t nclocks = (iloop < 58599) ? 0 : -2; //0 for run 121 and -2 for run 122
      //if (!(iloop%1000000))
      //	cout << "Processing event #" << iloop << "/" << 
      //	  reco->GetEntries() << endl; 
      reco->GetEntry(iloop);
      vector<double> f90 = evReco->GetF90();
      vector<double> start_time = evReco->GetStartTime();
      vector<double> charge = evReco->GetCharge();
      vector<double> ymax = evReco->GetYmax();
      vector<double> basemean = evReco->GetBaseMean();
      EvHeader* evh = evReco->GetEvHeader();
      Long64_t b1 =  evh->GetBoardTime(1);
      Long64_t b0 =  evh->GetBoardTime(0);
      Double_t deltaE = ymax.at(18)-basemean.at(18);
      Double_t E = ymax.at(19)-basemean.at(19);
      fBanana->Fill(E,deltaE);
      Double_t startSi = 0.5*(start_time.at(16)+start_time.at(17));
      Long64_t boardCorrection = 
	8*(b1-b0)+ 8*nclocks;
      //cout << "Board correction: " << iloop << " " << boardCorrection << endl;
      //      boardCorrection = 0;
     //Fill bananas
      for (size_t ich=0;ich<8;ich++)
	{	 
	  Double_t energy = GetEnergy(charge.at(ich),chanIDs->at(ich));
	  Double_t deltaT = (start_time.at(ich)-startSi)*2.
	    	- (Double_t) boardCorrection; 
	  if (TMath::Abs(deltaT) < 15. && energy>20.)
	    {
	      if (ich==4)
		f4Banana->Fill(E,deltaE);
	      
	      if (ich==0 &&  GetEnergy(charge.at(4),chanIDs->at(4))<20. )
		f0Banana->Fill(E,deltaE);
	      	
	    }
	}
      
      //Apply cuts
      Bool_t isBe = false;
      if (cut) //The cut is found
	isBe = cut->IsInside(E,deltaE);
      else
	{
	  if (deltaE >1300 && deltaE<1600 && E>3500 && E<4400)
	    isBe = true;
	}

      if (!isBe)
	continue;     
      becounts++; //counts inside the Be window
      

      Int_t id=4;
      Double_t energy4 = GetEnergy(charge.at(id),chanIDs->at(id));
      //Look for signal in #4
      Double_t deltaT = (start_time.at(id)-startSi)*2. - 
	(Double_t) boardCorrection;     
      dt1->Fill(deltaT);
      if (TMath::Abs(deltaT) > 15.)
	continue;
      if (energy4 < 50.)
	continue;
      if (f90.at(4)<0.12)
	continue;
    
      a1->Fill(energy4);
      counter1++;

      for (size_t ich=0;ich<8;ich++)
	{
	  if (ich==4) 
	    continue;
	  Double_t energy = GetEnergy(charge.at(ich),chanIDs->at(ich));
	  

	  if (energy < 10.)
	    continue;
	  //cout << energy << " " << energy4 << " " << ich << endl;  
	  //In sync with #4
	  Double_t subtraction = (ich==0) ? 50. : 25;
	  Double_t threshold = (ich==0) ? 25. : 15;

	  Double_t delta2 = (start_time.at(ich)-start_time.at(4))*2.; //-subtraction;
	  if (ich == 0) 
	    {
	      dt2->Fill(delta2); 
	    }
	  else 
	    {	      
	      dt3->Fill(delta2);
	      dt4->Fill(delta2);
	    }
	  //events too far
	  if (TMath::Abs(delta2-subtraction)>threshold)
	    continue;
	  
	  /*	  if (ich == 0)
	  cout << "Found triple event: " << 
	    energy4 << " " << f90.at(4) << " " << deltaT << " / " <<
	    ich << " " << energy << " " << start_time.at(ich)-start_time.at(4) << endl;
	  */	
	  if (ich==0)
	    { 
	      a2->Fill(energy4);
	      counter2++;
	      //cout << "Sono qui: " << energy4 << " " << endl;
	    }
	  else
	    {
	      a3->Fill(energy4);
	      counter3++;
	    }	      
	}     
    }

  //Normalize according to events in the tail
  Int_t ibin = a1->FindBin(1860.);
  Double_t int1 = a1->Integral(ibin,50);
  Double_t int3 = a3->Integral(ibin,50);
  //cout << ibin << " " << a1->GetBinCenter(ibin) << " " << int1 << " " << int3 << endl;
  // Int_t ibin2 = a2->FindBin(1800.);
  //Double_t int2 = a2->Integral(ibin2,25);
  //cout << ibin2 << " " << a2->GetBinCenter(ibin2) << " " << int1 << " " << int2 << endl;
  
  cout << "Found " << counter1 << " events in PMT4 " << endl;
  cout << "Found " << counter3 << " events in PMT4 + Wheel " << endl;
  cout << "Found " << counter2 << " events in PMT4 + PMT0 " << endl;


  a1->Sumw2();
  a2->Sumw2();
  a3->Sumw2();

  //a2->Scale(int1/int2);
  a3->Scale(int1/int3);
  cout << "Scaling factor: " << int1/int3 << endl;
  a3->Rebin(2);
  a3->Scale(0.5);
  
  a2->Scale(int1/int3*counter3/(Double_t)counter2);
  a2->Scale(0.25);

  

  TCanvas* c1 = new TCanvas();
  c1->cd();
  // c1->SetLogz();
  fBanana->Draw("COLZ");
  f4Banana->Draw("Psame");
  //f0Banana->Draw("Psame");
  //fOthersBanana->Draw("Psame");
  if (cut)
    cut->Draw("same");
  TCanvas* c2 = new TCanvas();
  c2->Divide(1,3);
  TVirtualPad* p1 = c2->cd(1);
  p1->SetLogy();
  dt1->Draw();
  c2->cd(2);
  dt2->Draw();
  c2->cd(3);
  dt3->Draw();

  TCanvas* c3 = new TCanvas();
  c3->cd();
  a3->SetLineColor(kBlue);
  a3->SetMarkerColor(kBlue);
  a3->SetMarkerStyle(20);
  a2->SetLineColor(kRed);
  a2->SetMarkerColor(kRed);
  a2->SetMarkerStyle(21);
  a3->GetXaxis()->SetTitle("E_{PMT4} (keV_{ee})"); 
  a3->Draw("P");
  a1->SetLineColor(kBlack);

  a1->Draw("HISTsame");
  a2->Draw("same");
  
  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry(a1,Form("PMT4 and Si, %d events",counter1),"l");
  legend->AddEntry(a3,Form("Wheel PMT coincidences, %d events",counter3),"lep");
  legend->AddEntry(a2,Form("45 deg PMT coincidences, %d events",counter2),"lep");
  //legend->Draw("same");


  TFile of("triples.root","RECREATE");


  of.cd();
  fBanana->Write(fBanana->GetName());
  f4Banana->Write(f4Banana->GetName());
  dt1->Write(dt1->GetName());
  dt2->Write(dt2->GetName());
  dt3->Write(dt3->GetName());
  dt4->Write(dt4->GetName());

  of.Close();

  return 0;
  
}
double GetEnergy(Double_t charge,TString id)
{
  static bool readfile = false;
  static std::map<TString,calibrationdata> cali;
  if (!readfile)
    {
      ifstream iFile("fit_parameters_kb.cfg");
      TString name;
      Double_t a,b,c;
      //Load static data
      while (iFile >> name >> a >> b >> c) 
	{
	  calibrationdata entry;
	  entry.a = a;
	  entry.b = b;
	  entry.kb = c;
	  cali.insert(std::pair<TString,calibrationdata>(name,entry) );
	}	  
      iFile.close();
      cout << "Read file. Closing" << endl;
      readfile = true;
    }
  calibrationdata curve = cali.find(id)->second;

  Double_t energy0 = charge*curve.b + curve.a;
  Double_t energy =GetQuenched(energy0,curve.kb);

  //cout << charge << " " << id << " " << energy0 << " " << energy << endl;
  return energy;

}

//Quenching parametrization, written by Marco
Double_t GetQuenched(Double_t energy0, Double_t kB)
{
  if (energy0<=0)
    return 0.;

  static const Double_t GrauMalonda[21][7]={
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
  Double_t logE=log(energy0);

  Double_t num = GrauMalonda[ikb][1] + GrauMalonda[ikb][2]*logE + 
    GrauMalonda[ikb][3]*pow(logE,2) + GrauMalonda[ikb][4]*pow(logE,3);

  Double_t den = 1 + GrauMalonda[ikb][5]*logE + 
    GrauMalonda[ikb][6]*pow(logE,2) + GrauMalonda[ikb][4]*pow(logE,3);

  if (den!=0) 
    return energy0/(num/den);
  else 
    return 0;
}

