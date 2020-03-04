#include "TTree.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TVirtualPad.h"
#include "TChain.h"
#include "TMath.h"
#include <iostream>
#include "Math/Math.h"
#include "Math/Error.h"
#include "Math/ProbFuncMathCore.h"
#include "Math/SpecFuncMathCore.h"

using namespace std;

Double_t GetEnergy(Double_t val,Int_t ichan=0)
{
  return val*60./7500.;
}

int readlight(TString filename,Double_t nFissions=1, Bool_t background=false, 
	      Double_t offset=0)
{
  const Double_t mN = 939.5653; //MeV
  const Double_t c_light = 30.0; //cm/ns
  vector<Double_t> distances{99.1,99.5,99.6,99.5,99.7,99.5,99.1,98.9};
  Double_t distanceNear = 3.23; //cm
  
  Double_t timingOffset = 4.55; //ns

  Double_t threshold = 100.; //threshold on far detectors

  TFile* ifile = new TFile(filename);
  TTree* tree = (TTree*) ifile->Get("t1");
  
  // TF1* cffun = new TF1("Cf252","[0]*TMath::Exp(-0.88*x)*TMath::SinH(TMath::Sqrt(2.*x))",0,12);
  TF1* cffun = new TF1("Cf252-watt",
		       "[0]*TMath::Sqrt(x)*TMath::Power([1],1.5)*TMath::Exp(-x/[1])",
		       0,12);
  cffun->SetParameter(1,1.406);
  
  cffun->SetParameter(0,1);

  /*
  TChain* tree = new TChain("t1");
  tree->AddFile("run_75_int.root.light");
  tree->AddFile("run_77_int.root.light");
  tree->AddFile("run_78_int.root.light");
  tree->AddFile("run_79_int.root.light");
  */
  Double_t tof = 0;
  Double_t chargeFar = 0;
  Double_t chargeNear = 0;
  Double_t f90Far = 0;
  Double_t f90Near = 0;
  Int_t multiplicity = 0;
  Int_t nchan = 0;
  Int_t eventNb = 0;
  tree->SetBranchAddress("eventNb",&eventNb);
  tree->SetBranchAddress("tof",&tof);
  tree->SetBranchAddress("chargeNear",&chargeNear);
  tree->SetBranchAddress("chargeFar",&chargeFar);
  tree->SetBranchAddress("f90Near",&f90Near);
  tree->SetBranchAddress("f90Far",&f90Far);
  tree->SetBranchAddress("nchan",&nchan);
  tree->SetBranchAddress("multiplicity",&multiplicity);

  vector<TH1D*> spectra;
  vector<TH1D*> tofs;
  vector<TH1D*> f90s;
  vector<TH2D*> f90vsE;

  const Int_t nE = 70;
  Double_t tlow = 24;

  //tof between 20 and 200, 2 ns steps
  Double_t xbins[nE+1];
  for (size_t i=0;i<nE+1;i++)    
    {
      Double_t atof = tlow + i*2;
      Double_t velocity = (distances.at(0)+offset)/atof; //cm/ns
      Double_t beta = velocity/30.; //30 cm/ns = c
      Double_t gamma = 1./TMath::Sqrt(1-beta*beta);
      xbins[nE-i] = mN*(gamma-1);     
      //cout << i << " " << mN*(gamma-1) << " " << atof << endl;
    }
  TH1D* h1 = new TH1D("h1","Energy spectrum",nE,xbins);
  Int_t nslices=15; 
  Double_t emax = 200000.;
  for (size_t i=16;i<24;i++)
    {
      spectra.push_back(new TH1D(Form("h%d",i),
				 Form("Charge, far chan %d",i),
				 500,0,80000));
      tofs.push_back(new TH1D(Form("tof%d",i),
			      Form("Tof, far chan %d",i),
			      1600,-200,1400));
      f90s.push_back(new TH1D(Form("psd%d",i),
			      Form("PSD, far chan %d",i),
			      200,0.,0.5));      

      f90vsE.push_back(new TH2D(Form("psdvsE%d",i),
				Form("PSD vs. charge, chan %d",i),
				nslices,0,emax,
				100,0,0.5));
    }

  

  for (Int_t iev=0;iev<tree->GetEntries();iev++)
    {
      tree->GetEntry(iev);
      //Remove the offset.
      tof -= timingOffset;
      if (multiplicity==1)
	{
	  if (f90Far>0.13)
	    {
	      spectra.at(nchan-16)->Fill(chargeFar);
	      if (chargeFar> (threshold*125))
		tofs.at(nchan-16)->Fill(tof);
	    }
	  f90s.at(nchan-16)->Fill(f90Far);
	  f90vsE.at(nchan-16)->Fill(chargeFar,f90Far);

	  //Calculate neutron velocity
	  if (tof > 20 && f90Far>0.13 && chargeFar> (threshold*125))
	    {
	      //Double_t velocity = (distances.at(nchan-16)+offset)/tof; //cm/ns
	      //Double_t beta = velocity/c_light; //30 cm/ns = c

	      //Formula: 
	      // beta = x2/(c*tof-x1), with x2=source-far, x1=source-near
	      Double_t beta = (distances.at(nchan-16)+offset)/(c_light*tof - (distanceNear-offset));	     

	      Double_t gamma = 1./TMath::Sqrt(1-beta*beta);
	      Double_t kinenergy = mN*(gamma-1);
	      //cout << tof << " " << velocity << " " << beta << " " << gamma << " " << kinenergy << endl;
	      h1->Fill(kinenergy);
	    }
	}


    }

  //Find minimum of PSD spectra:
  TF1* gaus = new TF1("gaus","gaus(0)",0.,0.3);
  TF1* fun = new TF1("fun","gaus(0)+gaus(3)",0.,0.3);
  Double_t pars[6];
  Double_t xlow = 0.055;
  Double_t xmed = 0.12;
  Double_t xhigh = 0.23;    
  for (size_t i=0;i<f90s.size();i++)
    {
      f90s.at(i)->Fit(gaus,"Q","",xlow,xmed);
      gaus->GetParameters(&pars[0]);
      f90s.at(i)->Fit(gaus,"Q","",xmed,xhigh);
      gaus->GetParameters(&pars[3]);
      //for (size_t j=0;j<6;j++)
      //cout << "par[" << j << "] = " << pars[j] << endl;
      fun->SetParameters(pars);
      fun->SetLineColor(kBlue);
      f90s.at(i)->Fit(fun,"Q","",xlow,xhigh);
      //Write results
      cout << "Chan: " << i+16 << " gamma:   " << fun->GetParameter(1) 
	   << " +/-" << fun->GetParameter(2) << 
	endl; 
      cout << "         neutron: " << fun->GetParameter(4) 
	   << " +/-" << fun->GetParameter(5) << 
	endl; 

    }

  Double_t nMultiplicity = 3.76;

  TCanvas* c1 = new TCanvas();
  c1->Divide(3,3);
  for (size_t i=0;i<spectra.size();i++)
    {
      c1->cd(i+1);
      spectra.at(i)->DrawCopy();
      Double_t integral = spectra.at(i)->Integral();
      cout << "Chan: " << i+16 << ", integral in all spectrum --> " << 
	integral << endl;
    }

  TCanvas* c1b = new TCanvas();  
  c1b->cd();
  spectra.at(0)->DrawCopy("HIST");  
  spectra.at(4)->SetLineColor(kRed);
  spectra.at(4)->DrawCopy("same");
  /*
  for (size_t i=0;i<spectra.size();i++)
    {
      spectra.at(i)->SetLineColor(1+i);
      if (!i)
	spectra.at(i)->DrawCopy("HIST");
      else
	spectra.at(i)->DrawCopy("same");
    }
  */

  Double_t totalSolidAngle = 0.;
  TCanvas* c2 = new TCanvas();
  c2->Divide(3,3);
  c2->SetLogy();
  Double_t totalBackgroundCounts = 0.;
  for (size_t i=0;i<tofs.size();i++)
    {
      TVirtualPad* p = c2->cd(i+1);
      if (tofs.at(i)->GetEntries() < 10)
	continue;
      p->SetLogy();
      tofs.at(i)->Fit(gaus,"Q","",-1,12);
      cout << "Gamma tof (ch" << i+16 << ") = " << gaus->GetParameter(1) << " ns, sigma = " << gaus->GetParameter(2) << endl;

      tofs.at(i)->DrawCopy();
      Int_t bin1 = tofs.at(i)->FindBin(20.);
      Int_t bin2 = tofs.at(i)->FindBin(200.);
      Double_t integral = tofs.at(i)->Integral(bin1,bin2);
      // cout << "Chan: " << i+16 << ", integral in [20,200] ns --> " << 
      //integral << endl;
      bin1 = tofs.at(i)->FindBin(800.);
      bin2 = tofs.at(i)->FindBin(1200.);
      Double_t integral2 = tofs.at(i)->Integral(bin1,bin2);
      totalBackgroundCounts += integral2;
      //cout << "Chan: " << i+16 << ", integral in [800,1200] ns --> " << 
      //integral2 << endl;
      Double_t net = (background) ? integral : integral-(integral2*180./400.);
      if (background)
	cout << "Chan: " << i+16 << ", integral in [20,200] ns --> " << 
	  net << endl;
      else
	{
	cout << "Chan: " << i+16 << ", integral in [20,200] ns w/o bck --> " << 
	  net << endl;
	Double_t solidAngle = TMath::Pi()*3.81*3.81/((distances.at(i)+offset)*(distances.at(i)+offset));
	
	totalSolidAngle += solidAngle;
	cout << " Estimated efficiency:" << 
	  net/((solidAngle/(4.*TMath::Pi()))*nMultiplicity*nFissions) << " (" << solidAngle << " sr)" << endl;
	}
    }
  //Total background counts: calculated in 800,1200. Normalize it per ns
  totalBackgroundCounts /= 400;
  if (!background)
    cout << "Total background:"  << totalBackgroundCounts << " counts/ns" << endl;


  TCanvas* c3 = new TCanvas();
  c3->Divide(3,3);
  for (size_t i=0;i<f90s.size();i++)
    {
      c3->cd(i+1);
      f90s.at(i)->DrawCopy();
    }
  
  TCanvas* c3b = new TCanvas();  
  c3b->cd();
  tofs.at(0)->GetXaxis()->SetRangeUser(-5,20);
  tofs.at(0)->DrawCopy("HIST");
  for (size_t i=1;i<tofs.size();i++)
    {
      tofs.at(i)->SetLineColor(i+2);
      tofs.at(i)->DrawCopy("same");
    }
  //  tofs.at(0)->DrawCopy("HIST");
  //tofs.at(4)->SetLineColor(kRed);
  //tofs.at(4)->DrawCopy("same");
  

  TCanvas* c6 = new TCanvas();
  c6->Divide(4,4);
  xlow = 0.035;
  xmed = 0.14;
  xhigh = 0.28;   
  TGraphErrors* gneutrons = new TGraphErrors(nslices);
  TGraphErrors* ggammas = new TGraphErrors(nslices);
  TGraphErrors* gnleakage = new TGraphErrors(nslices);
  for (int i=1;i<=nslices;i++)
    {
      c6->cd(i);
      //      cout << "Processing slice: " << i << endl;
      TH1D* htemp = f90vsE.at(0)->ProjectionY("",i,i+1);
      htemp->Fit(gaus,"Q","",xlow,xmed);
      gaus->GetParameters(&pars[0]);
      htemp->Fit(gaus,"Q","",xmed,xhigh);
      gaus->GetParameters(&pars[3]);
      //for (size_t j=0;j<6;j++)
      //cout << "par[" << j << "] = " << pars[j] << endl;
      fun->SetParameters(pars);
      fun->SetLineColor(kBlue);
      htemp->Fit(fun,"QL","",xlow,xhigh);
      htemp->DrawCopy();
      cout << "Slice: " << i << " gamma:   " << fun->GetParameter(1) 
	   << " +/-" << fun->GetParameter(2) << 
	endl; 
      cout << "         neutron: " << fun->GetParameter(4) 
	   << " +/-" << fun->GetParameter(5) << " cut is " << (fun->GetParameter(4)-0.13)/fun->GetParameter(5) << " sigmas away" << 
	endl;
      gneutrons->SetPoint(i-1,i,fun->GetParameter(4));
      gneutrons->SetPointError(i-1,0.5,fun->GetParameter(5));

      ggammas->SetPoint(i-1,i,fun->GetParameter(1));
      ggammas->SetPointError(i-1,0.5,fun->GetParameter(2));

      Double_t leakage = 1. - ROOT::Math::normal_cdf((fun->GetParameter(4)-0.13),fun->GetParameter(5));
	
//TMath::NormQuantile((fun->GetParameter(4)-0.13)/fun->GetParameter(5));	
      gnleakage->SetPoint(i-1,GetEnergy((i+0.5)*emax/nslices),leakage);
    }
  TCanvas* c7 = new TCanvas();
  c7->Divide(1,2);
  TH2D* null = new TH2D("hnull","",100,0,nslices+2,100,0.,0.3);
  c7->cd(1);
  null->DrawCopy();
  gneutrons->Draw("ZP");
  ggammas->Draw("ZP");
  c7->cd(2);
  gnleakage->Draw("AZP");
  

  TCanvas* c4b = new TCanvas();  
  c4b->Divide(1,2);
  c4b->cd(1);
 
  
  //Normalize the spectrum to the same number of fissions
  Double_t int1 = cffun->Integral(0.,12.);
  Double_t param = int1/((totalSolidAngle/(4.*TMath::Pi()))*nMultiplicity*nFissions);
  cffun->SetParameter(0,1./param);
  cout << "Integral of Cf252 after renormalization: " << cffun->Integral(0.,12.)  << endl;
  cout << "Total solid angle: " << totalSolidAngle << " sr" << endl;
  Double_t  totalDataCounts = 0;

  TH1D* h1b = new TH1D("h1b","Energy spectrum, theoretical",nE,xbins);
  TH1D* h1c = new TH1D("h1c","Ratio",nE,xbins);
  for (int i=1;i<= h1->GetNbinsX();i++)
    {
      Double_t x = h1->GetBinCenter(i);
      Double_t y1 = h1->GetBinContent(i);
      Double_t x1 = h1->GetBinLowEdge(i);
      Double_t DeltaX = h1->GetBinWidth(i);
      //Remove also background. Each bin corresponds to 2 ns tof     
      //h1->SetBinContent(i,(y1-totalBackgroundCounts*2.)/DeltaX);
      h1->SetBinContent(i,(y1-totalBackgroundCounts*2.));
      totalDataCounts += y1-totalBackgroundCounts*2.;

      Double_t y= cffun->Integral(x1,x1+DeltaX);     
      h1b->SetBinContent(i,y);
      //h1b->SetBinContent(i,y/DeltaX);
      h1c->SetBinContent(i,h1->GetBinContent(i)/h1b->GetBinContent(i));
    }
  cout << "Integral of the data histogram: " << totalDataCounts << endl;
  //h1b->GetYaxis()->SetTitle("counts/MeV");
  h1b->SetLineColor(kGreen);
  h1b->Draw();
//cffun->Draw("same");
  h1->Draw("same");
  c4b->cd(2);
  TFile gg("efficiencies.root");
  TH1D* r3 = (TH1D*) gg.Get("r3");
  r3->DrawCopy();
  h1c->DrawCopy("same");
 

  return 0;
}
