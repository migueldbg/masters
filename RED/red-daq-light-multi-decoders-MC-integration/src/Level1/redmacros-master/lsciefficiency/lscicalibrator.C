#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include <iostream>
#include "red-daq/EvRec0.hh"
#include "red-daq/EvHeader.hh"

/*
  This macro reads the full reco file from red-daq-light and produced 
  a lighter root TTree which contains only coincidence events (thresholds 
  are 100 keV and 20 keV for near and far detector, respectively). The threshold 
  for the near detector can be changed.

  The output TTree contains: energies, f90s, tof, multiplicity and 
  a saturation flag.

  Energies are calibrated in keVee using the calibration curve read 
  from the file fit_parameters_kb.cfg (linear + Birks, taken from Marco).

  The input parameter nclocks represents the number of DAQ clocks (8 ns each) 
  required to correct for the alignment of the boards.
  
  Counting also the total number of fissions (tag-and-probe) from the near detector 
  above a given threshold and producing a histo with the full near-detector 
  spectrum (1D and 2D, with f90).

*/

Bool_t IsAGamma(Double_t energy,Double_t f90)
{
  //Flat below 1900 keVee and linear above.
  if (energy<1900)
    return (f90<0.13);
  else
    {
      Double_t thre = -6.2256e-2+energy*1.0119e-4;
      return (f90<thre);
    }
}

using namespace std;

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

int lscicalibrator(TString filename, bool isSingle, int nclocks=0, 
                   Double_t nearThreshold=100. /*keVee*/ )
{
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
  cout << reco->GetEntries() << " entries" << endl;
  //Fare un TTree con: tof, f90near, f90far, charge, charge2, nchan
  // multiplicity and saturation flag
  Double_t tof = 0;
  Double_t chargeFar = 0;
  Double_t chargeNear = 0;
  Double_t f90Far = 0;
  Double_t f90Near = 0;
  Int_t multiplicity = 0;
  Int_t issaturated = 0;
  Int_t nchan = 0;
  Int_t eventNb = 0;
  TFile* ofile = new TFile(filename + ".light","RECREATE");
  ofile->cd();
  TTree* oTree = new TTree("t1","tof output");
  oTree->Branch("eventNb",&eventNb,"eventNb/I");
  oTree->Branch("tof",&tof,"tof/D");
  oTree->Branch("chargeNear",&chargeNear,"chargeNear/D");
  oTree->Branch("chargeFar",&chargeFar,"chargeFar/D");
  oTree->Branch("f90Near",&f90Near,"f90Near/D");
  oTree->Branch("f90Far",&f90Far,"f90Far/D");
  oTree->Branch("nchan",&nchan,"nchan/I");
  oTree->Branch("multiplicity",&multiplicity,"multiplicity/I");
  oTree->Branch("issaturated",&issaturated,"issaturated/I");

  TH1D* nearSpe = new TH1D("h1","Spectrum of the near detector",1000,0.,5000.);
  TH2D* nearSpe2D = new TH2D("h2","Energy/f90 spectrum",1000,0.,5000.,200,0,1.);
  //Also count the total number of fissions, for runs in single-mode
  Int_t nFissions = 0;
  Int_t nFissionsGammas = 0; //Fissions where the near tag is from a gamma

  for (Int_t iloop=0;iloop<reco->GetEntries();iloop++)
    {            
      if (!(iloop%1000000))
	cout << "Processing event #" << iloop << "/" << 
	  reco->GetEntries() << endl; 
      reco->GetEntry(iloop);
      multiplicity=0;
      vector<double> f90 = evReco->GetF90();
      vector<double> start_time = evReco->GetStartTime();
      vector<double> base_mean = evReco->GetBaseMean();
      vector<double> charge = evReco->GetCharge();
      vector<double> ymin = evReco->GetYmin();
      EvHeader* evh = evReco->GetEvHeader();
      Long64_t b1 =  evh->GetBoardTime(1);
      Long64_t b0 =  evh->GetBoardTime(0);
      //Get info for near detector
      if (!charge.at(0))
	continue;
      Double_t startnear = start_time.at(0);
      chargeNear = GetEnergy(charge.at(0),chanIDs->at(0));
      f90Near = f90.at(0);
      eventNb = iloop;

      if (chargeNear > nearThreshold)
	nFissions++;
      if (chargeNear > nearThreshold && IsAGamma(chargeNear,f90Near))
        nFissionsGammas++;
      if (chargeNear > 0.1)
	{
	  nearSpe->Fill(chargeNear);
	  nearSpe2D->Fill(chargeNear,f90Near);
	}

      Long64_t boardCorrection = 
	8*(b1-b0) + 
	8*nclocks;
      
      vector<bool> isGood(charge.size(),false);

      //First loop on far channels: calculate multiplicity
      for (size_t ichan=1;ichan<charge.size();ichan++)
	{
	  if (!charge.at(ichan))
	    continue;
	  
	  //Selection here!
	  //Check for valid signal (in energy): above 100 keV in near and
	  //above 20 keV in far.
	  bool keepIt = (chargeNear > nearThreshold) & (GetEnergy(charge.at(ichan),chanIDs->at(ichan)) > 20.); 
	  //cout << charge.at(ichan) << " " << chanIDs->at(ichan) << " " << GetEnergy(charge.at(ichan),chanIDs->at(ichan)) << endl;
	  if (isSingle) //Apply extra selection 
	    keepIt &= (TMath::Abs(startnear-500)<50);

	  if (!keepIt) 
	    continue;
	  isGood.at(ichan) = true;
	  multiplicity++;
	}

      if (multiplicity == 0)
	continue;

      //Second loop: fill TTree
      for (size_t ichan=1;ichan<charge.size();ichan++)
	{	  	  
	  if (isGood.at(ichan))
	    {
	      issaturated = 0;
	      tof = (start_time.at(ichan)-startnear)*2. + (Double_t) boardCorrection;
	      /*	  cout << iloop << " " << start_time.at(ichan) << " " << 
			  startnear << " " <<  boardCorrection << " -->" << tof <<  endl;
			  cout << evh->GetBoardTime(1) << " " << evh->GetBoardTime(0) << " " << 
			  nclocks << endl;
	      */ 
	      chargeFar =  GetEnergy(charge.at(ichan),chanIDs->at(ichan));	      
	    
	      f90Far = f90.at(ichan);
	      if (ymin.at(ichan)< 10)
		issaturated = 1;
	      nchan = ichan; 	      
	      oTree->Fill();
	    }
	}  
    }
 
 

  oTree->Write(oTree->GetName());
  nearSpe->GetXaxis()->SetTitle("Energy (keV)");
  nearSpe->Write(nearSpe->GetName());
  nearSpe2D->GetXaxis()->SetTitle("Energy (keV)");
  nearSpe2D->Write(nearSpe2D->GetName());


  /*
  TCanvas* c1 = new TCanvas();
  c1->Divide(4,4);
 
  TCanvas* c2 = new TCanvas();
  c2->Divide(4,4);
  */
  ofile->Write();
  ofile->Close();
  
  if (isSingle)
    cout << "File: " << filename << " has " << nFissions << 
      " fissions tagged in the near detector above " << nearThreshold << " keV " << endl; 
    cout <<  " --> " << nFissionsGammas << " come from gamma-rays" << endl;   


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

