#include "TH2Poly.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCut.h"
#include "TCutG.h"
#include "TLatex.h"
#include <iostream>
#include "red-daq/EvRec0.hh"

/*
  This macro reads the output reco files from the LNS Test Beam 
  of May/July 2019. It is assumed that the "full" configuration is used 
  (i.e. LSci on board0 and Si on board1; trigger on Si only). 

  The argument of the macro is the run number. The root file is 
  expected to be found in the current directory.

  The flower counts are expresses as fraction (default) or as 
  counting rate, if useLiveTime is set to true.

  The macro select coincidence events in the flower requesting:
  - Be events in the Si
  - charge > 3000 ADC in the PMTs (~ 30 keVee)
  - PSD compatible with neutron (f90 > 0.12)
  - timing within 100 ns wits respect te the DeltaE-fast

  Inter-board timing corrections are performed according to the 
  recipe from the Cf252 calibration campaign. The input parameter 
  nclocks represents the number of DAQ clocks (8 ns each) 
  required to correct for the alignment of the boards.

  The Be7 is selected using a graphical cut read from BeCut.root. 
  Different cuts are used for different runs. 
  If the file is not available, a square box is used in DeltaE/E 
  (which is valid only for the May runs: gains are changed). 

  The macro produces: 2D histo with the full bananas (DeltaE/E) 
  superimposed with the events in coincidence; flower plot, normalized 
  to the total number of entries. All histograms are written on file.

*/

void Spot(int n, double x0, double y0, double r, 
	  double *px, double *py) 
{
  // Add points on a arc of circle from point 2 to n-2
  double da = 2*TMath::Pi()/n; // Angle delta
  for (int i = 0; i<n; i++) {
    Double_t a     = (i+0.5)*da;
    px[i] = r*TMath::Cos(a) + x0;
    py[i] = r*TMath::Sin(a) + y0;
  }
}

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

int plotFlower(Int_t runNb, Bool_t useLiveTime = false, Int_t nclocks = 0)
{
  TString angleString = " ";
  Double_t angle = 0.;
  Double_t livetime = 0;
  Double_t current = 0;
  if (runNb == 102)
    {
      angle = 5.0;
      angleString = "(5.0 #pm 0.5)";
      livetime = 3665.2;
      current = 1.43;
    }
  else if (runNb == 104)
    {
      angle = 5.2;
      angleString = "(5.2 #pm 0.2)" ;
      livetime = 4224.5;
      current = 1.05;
    }
  else if (runNb == 109)
    {
      angle = 5.26;
      angleString = "(5.26 #pm 0.13)" ;
      livetime = 10828.1;
      current = 2.13;
    }
  else if (runNb == 112)
    {
      angle = 4.61;
      angleString = "(4.61 #pm 0.13)" ;
      livetime = 5886.0;
      current = 1.33;
    }
  else if (runNb == 115)
    {
      angle = 4.61;
      angleString = "off - #phi" ;
      livetime = 4605.0;
      current = 1.92;
    }
  else if (runNb == 119)
    {
      angle = 5.13;
      angleString = "(5.13 #pm 0.13)" ;
      livetime = 1775.2;
      current = 13.2;
    }
  else if (runNb == 122)
    {
      angle = 5.13;
      angleString = "(5.13 #pm 0.13)" ;
      livetime = 265968.5;
      current = 13.2; //assumed
    }
  else if (runNb >= 162 && runNb <= 170)
    {
      angle = 5.1;
      angleString = "5.1";
    }
  else if (runNb == 171)
     {
      angle = 5.2;
      angleString = "5.2";
    }
  else if (runNb == 172)
   {
      angle = 5.0;
      angleString = "5.0";
    }
  else if (runNb == 173)
   {
      angle = 4.9;
      angleString = "4.9";
    }
  else if (runNb == 174)
   {
      angle = 4.95;
      angleString = "4.95";
    }
  else if (runNb == 175 || runNb==193)
   {
      angle = 5.05;
      angleString = "5.05";
    }  
  else if (runNb == 176)
    {
      angle = 5.1;
      angleString = "5.1";
    }
  else if (runNb == 184)
    {
      angle = 4.6;
      angleString = "4.6";
    }
  else if (runNb == 185)
    {
      angle = 5.45;
      angleString = "5.45";
    }
  
  gStyle->SetOptStat(0);
  
  const int NP = 50; // Number of point to build the current bin
  const int N = 5; // Number of point to build the current bin
  double px[NP];     // Bin's X positions
  double py[NP];     // Bin's Y positions

  Double_t radius = 3.81; //cm
  Double_t distance = 11.5; //cm, by Santi

  //Position of the flower 
  Double_t centerX[N] = {0.,-1*distance, distance, 0.,0.};
  Double_t centerY[N] = {0., 0.,0.,-1*distance,distance};
 			
  TString filename;
  filename.Form("run_%d.root",runNb);

  //May
  Double_t deltaEhigh = 2000;
  //July
  if (runNb > 126)
    deltaEhigh = 4000;

  TH2D* fBanana = new TH2D(Form("h0_%d",runNb),
			   Form("Full banana, run %d theta=%s",runNb,angleString.Data()),
			   400,1500,7000,400,600,deltaEhigh);
 
  TH2D* fCBanana = new TH2D(Form("h1_%d",runNb),
			    "Coincident banana",400,1500,7000,400,600,deltaEhigh);
  fCBanana->SetMarkerStyle(1);
  /*
  vector<TH1D*> timeDiff;
  for (size_t ich=0;ich<8;ich++)
    timeDiff.push_back(new TH1D(Form("time_%d",ich),Form("Time difference, Be events, ch %d",ich),
  				 200,-200,200));
  */
  TH1D* mul = new TH1D("m1","Multiplicity in LSci",5,-0.5,4.5);
  TH1D* h90 = new TH1D("h90","f90 for Be events",100,0.,1.);
  TH1D* dT = new TH1D("dT","DeltaT",400,-200,200);
  TH2D* dual = new TH2D("h90dt","f90 vs. deltaT",400,-200,200,100,0.,1.);
  TH1D* ene4 = new TH1D("ene4","Energy of petal 4",200,0.,3000.);

  TH2Poly *h2p = new TH2Poly();
  
  h2p->SetName(Form("h2p_%d",runNb));
  h2p->SetTitle(Form("Flower, run %d, theta %s",runNb, angleString.Data()));

  for (int i=0;i<N;i++)    
    {
      Spot(NP,centerX[i],centerY[i],radius,px,py);
      h2p->AddBin(NP,px,py);
    }      

  //Extra petal: runs 115 and 119
  Double_t xextra = distance + 13. + 2*radius;
  if (runNb >= 115 && runNb <= 119)
    {      
      Spot(NP,xextra,0.,radius,px,py);
      h2p->AddBin(NP,px,py);
    }



  //
  TString cutname = (runNb == 112 || runNb==115) ? 
    "betight" : "be";
  if (runNb >= 127) //May runs
    cutname = "be5deg-july";
  if (runNb == 184)
    cutname = "be46deg-july";


  //Open cut file
  TFile* fcut = new TFile("BeCut.root");
  TCutG* cut = nullptr;
  if (fcut->IsOpen())    
    { 
      cut = static_cast<TCutG*>(fcut->Get(cutname));
      fcut->Close();
    }
  if (cut)
    cout << " --> Using graphical cut called " << cutname << endl;


  //Open data file
  TFile* f = new TFile(filename);
  if (!(f->IsOpen()))
    {
      cout << "could not open file: " << filename << endl;
      return 0;
    }

  vector<TString> *chanIDs = new vector<TString>;
  //Read metadata
  TTree* metaevent = (TTree*) f->Get("metaevent");
  metaevent->SetBranchAddress("chanID",&chanIDs);
  metaevent->GetEntry(0);


  TTree* reco = (TTree*) f->Get("reco");
  EvRec0* evReco = new EvRec0();
  reco->SetBranchAddress("recoevent",&evReco);
  //cout << reco->GetEntries() << " entries" << endl;
  
  vector<Int_t> counter(8,0);
  vector<Int_t> counterBck(8,0);
  Int_t becounts = 0;

  //Event loop
  for (Int_t iloop=0;iloop<reco->GetEntries();iloop++)
    {            
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
      Double_t startSi = start_time.at(16);
      Double_t deltaE = ymax.at(18)-basemean.at(18);
      Double_t E = ymax.at(19)-basemean.at(19);
      fBanana->Fill(E,deltaE);
      
      Long64_t boardCorrection = 
	8*(b1-b0) + 
	8*nclocks;

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

      Int_t multiplicity= 0;
      //Loop over channels
      for (size_t ich=0;ich<8;ich++)
	{
	  //count multiplicity
	  Double_t energy = GetEnergy(charge.at(ich),chanIDs->at(ich));
	  //if (charge.at(ich)>3000)
	  //cout << charge.at(ich) << " " << chanIDs->at(ich) << " " << energy << " " << ich << endl;
	  if (energy > 20.)
	    multiplicity++;
	}
      mul->Fill(multiplicity);
      
      if (multiplicity < 1)
	continue;

      for (size_t ich=0;ich<8;ich++)
	{
	  Double_t energy = GetEnergy(charge.at(ich),chanIDs->at(ich));
	  if (energy < 20.)
	    continue;

	  Double_t tof = (start_time.at(ich)-startSi)*2. 
		- (Double_t) boardCorrection;

	  if (ich==4) //Check distributions on the central petal	    
	    ene4->Fill(energy);
	  h90->Fill(f90.at(ich));	      
	  dT->Fill(tof);
	  dual->Fill(tof,f90.at(ich));
	  //Keep neutrons only
	  if (f90.at(ich)<0.12)
	    continue;
	  
	  if (TMath::Abs(tof)>50 && 
	      TMath::Abs(tof)<200)
	    {
	      (counterBck.at(ich))++;
	    }
	  
	  if (TMath::Abs(tof)>50)
	    continue;

	  fCBanana->Fill(E,deltaE);
		
	  (counter.at(ich))++;
	}     
    }
  vector<Double_t> entries(8,0.);
  //Now subtract background and fill histos
  Double_t sum = 0;
  for (size_t ich=0;ich<8;ich++)
    {
      entries[ich] = counter.at(ich) - counterBck.at(ich)*0.25; //background: 300 ns, signal: 100 ns
      //cout << ich << " " << counter.at(ich) << " " << counterBck.at(ich) << " " << entries[ich] << endl;
      if (ich == 1) // PMT8
	h2p->Fill(-1*distance,0.,entries[ich]);
      else if (ich==2) 
	h2p->Fill(0.,-1*distance,entries[ich]);
      else if (ich==3) 
	h2p->Fill(0.,distance,entries[ich]);
      else if (ich==4) 
	h2p->Fill(0.,0.,entries[ich]);
      else if (ich==6) 
	h2p->Fill(distance,0.,entries[ich]);
      else if (ich == 5)
	{
	  if (runNb >= 115 && runNb <= 119)
	    h2p->Fill(xextra,0.,entries[ich]);
	}
      if (ich == 1 || ich == 2 || ich == 3 || ich == 4 || ich ==6 || (ich==5 && (runNb>=115 && runNb<=119)))
	{
	  sum += entries[ich];	
          cout << "Canale" << ich << " --> " << entries[ich] << " counts (bck subtracted)" << endl;
        }
   }
  
  Int_t totEvents = reco->GetEntries();
  f->Close();

  TFile of("flowers.root","UPDATE");
   

  /*cout << counter.at(1) << " " << counter.at(2) << " " << 
	counter.at(3) << " " << counter.at(4) << " " <<
	counter.at(6) << endl;
  */
  cout << "Angle: " << angleString << " deg" << endl;
  cout << "Total events: " << totEvents;
  if (livetime)
    cout << " --> " << totEvents/(livetime*current) << " Hz/nA ";
  cout << endl;
  cout << "Total events in Be band: " << becounts; 
  if (livetime)
    cout << " --> " << becounts/(livetime*current) << " Hz/nA ";
  cout << endl;
  cout << "Fraction of Be events with signal in LSci: " << (double)sum/becounts << " (" << 
    sum << "/" << becounts << ")" << endl;
  cout << "Fraction of Be events in central LSci: " << 
    entries[4]/becounts << 
    " (" << entries[4] << "/" << becounts << ")" << endl;

  TCanvas* c1 = new TCanvas();
  c1->cd();
  c1->SetLogz();
  fBanana->Draw("COLZ");
  fCBanana->Draw("Psame");
  if (cut)
    cut->Draw("same");
  
  TCanvas* c2 = new TCanvas("c2","Flower",4);
  c2->SetLogz();
  c2->cd();
  if (useLiveTime && livetime > 0)
    {
      h2p->Scale(1./(livetime*current));
      h2p->GetZaxis()->SetTitle("Hz/nA for 300 #mug/cm#^2#"); 
      h2p->GetZaxis()->SetRangeUser(1e-5,0.2);

    }
  else
    h2p->Scale(1./sum);  
  h2p->DrawCopy("COLZ TEXT");

  of.cd();
  fBanana->Write(fBanana->GetName());
  fCBanana->Write(fCBanana->GetName());
  h2p->Write(h2p->GetName());
  /*
    for (size_t ich=0;ich<timeDiff.size();ich++)
    timeDiff.at(ich)->Write(timeDiff.at(ich)->GetName());
  */
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

