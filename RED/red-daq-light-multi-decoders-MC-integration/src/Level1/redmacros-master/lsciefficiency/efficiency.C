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

/*
  This is the main macro for the calculation of the efficiency profile.
  
  Tof spectra are read from "light" files (by lscicalibrator()) and 
  corrected channel-by-channel for the time offset, derived from the 
  fit of timing_v2(). Distances between the near and far detectors are 
  taken from the output of timing_v2(): the near distance is the weighted 
  average from the output of the fit.

  Neutron-like events are selected and the 2-ns-binned tof distribution 
  is converted into a variabile-bin energy distribution. Accidental 
  background is subtracted. The threshold on the far detectors is an 
  input parameter (default = 100 keVee). 

  The experimental energy spectrum is compared against the theoretical 
  spectrum (Maxwell approximation). Absolute normalization is calculated 
  according to the total number of fissions, which has to be passed as an 
  input parameter (with background already subtracted).

  Efficiency profiles are superimposed with the MC results, if available. Plots and 
  results are saved as ROOT files and pdf.

  Error bars are purely statistical.

 */

Bool_t IsAGamma(Double_t energy,Double_t f90)
{
  return (f90<0.13);
  /*
  //Flat below 1900 keVee and linear above.
  if (energy<1900)
    return (f90<0.13);
  else
    {
      Double_t thre = -6.2256e-2+energy*1.0119e-4;
      return (f90<thre);
    }
  */
}

Bool_t IsANeutron(Double_t energy,Double_t f90)
{
  return (!IsAGamma(energy,f90));
}

using namespace std;

int efficiency(TString filename,
	       Double_t nFissions, //number of fissions from tag-and-probe
	       Double_t threshold = 100, //threshold on far detectors, keV
	       TString ofilename="output.root",
               Bool_t onlyNearGamma = false)

{
  //Constants and parameters
  const Double_t mN = 939.5653; //MeV
  const Double_t c_light = 30.0; //cm/ns
  const Double_t nMultiplicity = 3.76;
  vector<Double_t> distances(9,99.5); // cm; //{99.1,99.5,99.6,99.5,99.7,99.5,99.1,98.9};
  Double_t distanceNear = 0.; //cm
  vector<Double_t> timingOffset(9,4.55); //ns;

  Double_t extraTicks = 0;

  //Theoretical spectrum, with initial dummy normalization
  TF1* cffun = new TF1("Cf252-Maxwell",
		       "[0]*TMath::Sqrt(x)*TMath::Power([1],1.5)*TMath::Exp(-x/[1])",
		       0,12);
  cffun->SetParameter(1,1.406);  
  cffun->SetParameter(0,1);

  //Prepare energy spectrum with variable bin size
  //tof between 24 and 164, 2 ns steps
  const Int_t nE = 70; //number of bis
  Double_t tlow = 24; //left edge in tof
  const Double_t binSize = 2.;
  Double_t xbins[nE+1];
  for (size_t i=0;i<nE+1;i++)    
    {
      Double_t atof = tlow + i*binSize;
      Double_t velocity = distances.at(0)/atof; //cm/ns
      Double_t beta = velocity/30.; //30 cm/ns = c
      Double_t gamma = 1./TMath::Sqrt(1-beta*beta);
      xbins[nE-i] = mN*(gamma-1);     
      //cout << i << " " << mN*(gamma-1) << " " << atof << endl;
    }
  TH1D* h1 = new TH1D("h1","Energy spectrum: data",nE,xbins);
  TH1D* h2 = new TH1D("h2","Energy spectrum: theory",nE,xbins);
  TH1D* h3 = new TH1D("h3","Efficiency profile",nE,xbins);

  vector<TH1D*> h1det; //energy spectra for individual detectors
  vector<TH1D*> h3det; //efficiency profile for individual detectors
  vector<TH1D*> tofs; //vector of tof spectra
  vector<TH1D*> tofsgammas; //vector of tof spectra
  for (size_t i=16;i<24;i++)
    {
      tofs.push_back(new TH1D(Form("tof%d",i),
			      Form("Tof, far chan %d",i),
			      1600,-200,1400));
      h1det.push_back(new TH1D(Form("h1det_%d",i),
			       Form("Energy spectrum, data, chan %d",i),
			       nE,xbins));
      h3det.push_back(new TH1D(Form("h3det_%d",i),
			       Form("Efficiency profile, chan %d",i),
			       nE,xbins));
      tofsgammas.push_back(new TH1D(Form("tofg%d",i),
				    Form("Tof, gammas, far chan %d",i),
				    400,-40,40));
    }  


  //Read data from file!
  TChain* tree = new TChain("t1");
  
  if (filename.Contains("48")  || 
      filename.Contains("58"))
    {
      cout << "using run 48-58" << endl;
      tree->AddFile(filename);
      //Taken from run58
      timingOffset = {-2.65,-4.34,-2.98,-0.99,-4.45,-4.64,-2.30,-4.44};
      distances = {102.4,106.6,108.1,105.0,108.2,108.2,105.1,105.7};
      distanceNear = 8.06; //weighted average
      for (size_t ij=0;ij<timingOffset.size();ij++)
	timingOffset.at(ij) += 8.*extraTicks; //add ticks
    }      
  else if (filename.Contains("74"))
    {
      cout << "Using run 74" << endl;
      tree->AddFile(filename);
      timingOffset = {0.46,1.14};
      distances = {105.5,105.0};
      distanceNear = 8.26; //weighted average
      for (size_t ij=0;ij<timingOffset.size();ij++)
	timingOffset.at(ij) += 8.*extraTicks; //add ticks
    }
  else if (filename == " ")
    {
      tree->AddFile("run_75_int.root.light");
      tree->AddFile("run_77_int.root.light");
      tree->AddFile("run_78_int.root.light");
      tree->AddFile("run_79_int.root.light");
      timingOffset = {0.46,1.14};
      distances = {105.5,105.0};
      distanceNear = 8.26; //weighted average
      for (size_t ij=0;ij<timingOffset.size();ij++)
	timingOffset.at(ij) += 8.*extraTicks; //add ticks
    }
  else
    tree->AddFile(filename);

  cout << "Opening data files " << filename << ". Total entries: " << tree->GetEntries() << endl;

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

  

  // Event loop
  for (Int_t iev=0;iev<tree->GetEntries();iev++)
    {
      tree->GetEntry(iev);
      //Remove the offset.
      tof -= timingOffset.at(nchan-16);
      
      //Check if one wants only gammas or gammas+neutrons in the near detector
      Bool_t gammaSelectionNear = (onlyNearGamma) ? IsAGamma(chargeNear,f90Near) : true;      
      if (multiplicity==1 && IsANeutron(chargeFar,f90Far) && chargeFar> threshold && gammaSelectionNear)
	{
	  // Fill tof spectrum
	  tofs.at(nchan-16)->Fill(tof); 
	  if (tof> 20)
	    {
	      //Calculate neutron velocity
	      //Formula: 
	      // beta = (1/c)*x2/(tof-x1/c)
	      Double_t beta = (1./c_light)*(distances.at(nchan-16))/(tof-(distanceNear/c_light));

	      Double_t gamma = 1./TMath::Sqrt(1-beta*beta);
	      Double_t kinenergy = mN*(gamma-1);
	      //cout << tof << " " << velocity << " " << beta << " " << gamma << " " << kinenergy << endl;
	      h1->Fill(kinenergy);
	      h1det.at(nchan-16)->Fill(kinenergy);
	    }
	}	
      if (IsAGamma(chargeFar,f90Far) && IsAGamma(chargeNear,f90Near) && chargeFar>threshold)
	tofsgammas.at(nchan-16)->Fill(tof); 
    }
  //Done: data read.

  //Calculate solid angle and background subtraction  
  Double_t totalBackgroundCounts = 0.;
  Double_t totalSolidAngle = 0.;
  vector<Double_t> detSolidAngle; //solid angle of each detector
  vector<Double_t> detBackgroundCounts; //counts of each detector
  for (size_t i=0;i<tofs.size();i++)
    {
      //Background calculated in [800,1200] ns
      Int_t bin1 = tofs.at(i)->FindBin(800.);
      Int_t bin2 = tofs.at(i)->FindBin(1200.);
      totalBackgroundCounts += tofs.at(i)->Integral(bin1,bin2);
      detBackgroundCounts.push_back(tofs.at(i)->Integral(bin1,bin2));
      Double_t solidAngle = (tofs.at(i)->GetEntries() > 0) ? 
	TMath::Pi()*3.81*3.81/(distances.at(i)*distances.at(i)) : 0.;
      detSolidAngle.push_back(solidAngle); 
      totalSolidAngle += solidAngle;	
    }
  //Total background counts: calculated in 800,1200. Normalize it per ns
  totalBackgroundCounts /= 400;
  for (size_t i=0;i<detBackgroundCounts.size();i++)
    detBackgroundCounts.at(i) /= 400;
  cout << "Total background: "  << totalBackgroundCounts << " counts/ns" << endl;
  cout << "Solid angle: " << totalSolidAngle << " sr " << endl;
 
  //Now start the real thing.
  //Normalize the theoretical spectrum
  Double_t int0 = cffun->Integral(0.,12.);
  Double_t param = int0/((totalSolidAngle/(4.*TMath::Pi()))*nMultiplicity*nFissions);
  cffun->SetParameter(0,1./param);
  Double_t totalCfIntegral = cffun->Integral(0.,12.);
  cout << "Integral of Cf252 after renormalization: " << totalCfIntegral  << endl;
 
  Double_t  totalDataCounts = 0;
  //Handle the bin-per-bin calculations and ratio
  for (int i=1;i<= h1->GetNbinsX();i++)
    {
      //Data: subtract background
      Double_t x = h1->GetBinCenter(i);
      Double_t y1 = h1->GetBinContent(i);
      Double_t x1 = h1->GetBinLowEdge(i);
      Double_t DeltaX = h1->GetBinWidth(i);
      //Remove also background. Each bin corresponds to 2 ns tof     
      //h1->SetBinContent(i,(y1-totalBackgroundCounts*2.)/DeltaX);
      h1->SetBinContent(i,(y1-totalBackgroundCounts*binSize));
      totalDataCounts += y1-totalBackgroundCounts*2.;

      //Theory: calculate bin-per-bin integral
      Double_t y= cffun->Integral(x1,x1+DeltaX);     
      h2->SetBinContent(i,y);     

      //Ratio
      h3->SetBinContent(i,h1->GetBinContent(i)/h2->GetBinContent(i));
    }
  cout << "Integral of the data histogram: " << totalDataCounts << endl;
  cout << "Spectrum-averaged efficiency: " << totalDataCounts/totalCfIntegral << " (for " << threshold << " keV threshold)" << endl;

  //Now do channel by channel
  for (size_t idet = 0;idet<h1det.size();idet++)
    {     
      TH1D* h11 = h1det.at(idet);
      if (h11->GetEntries() == 0)
	continue;
      for (int i=1;i<= h11->GetNbinsX();i++)
	{
	  //Data: subtract background
	  Double_t x = h11->GetBinCenter(i);
	  Double_t y1 = h11->GetBinContent(i);
	  Double_t x1 = h11->GetBinLowEdge(i);
	  Double_t DeltaX = h11->GetBinWidth(i);
	  //Remove also background. Each bin corresponds to 2 ns tof     
	  
	  h11->SetBinContent(i,(y1-detBackgroundCounts.at(idet)*binSize));
	  h11->SetBinError(i,TMath::Sqrt(y1));
	  //cout << "Bin: " << i << " contenuto: " << y1 << " errore " << h11->GetBinError(i) << endl;

	  //Theory: will scale only according to solid angle
	  Double_t val = h2->GetBinContent(i)*detSolidAngle.at(idet)/totalSolidAngle;
	  
	  //Ratio
	  h3det.at(idet)->SetBinContent(i,h11->GetBinContent(i)/val);
	  h3det.at(idet)->SetBinError(i,h11->GetBinError(i)/val);
	}
    }

  //Now plot
  TCanvas* c1 = new TCanvas();
  c1->Divide(1,2);
  c1->cd(1);

  //Theory and data
  h2->SetLineColor(kGreen);
  h2->Draw();
  //cffun->Draw("same");
  h1->Draw("same");

  //Plot individual channels. Separate parts
  TFile* mcfile = new TFile("efficiencies_Inelastic_HP.root");
  TH1D* mcprofile = nullptr;
  if (mcfile)
    {
      mcprofile = (TH1D*) mcfile->Get("r3");
      cout << "Found histogram: " << mcprofile << endl;
      mcprofile->SetMarkerStyle(21);
      //mcprofile->SetMarkerSize(0.3);
      mcprofile->Rebin(5);
      mcprofile->Scale(0.2);
    }
  //ratio
  c1->cd(2); 
  h3->DrawCopy();
  if (mcprofile)
    mcprofile->DrawCopy("Psame");
 
  TCanvas* c2 = new TCanvas();
  c2->Divide(3,3);
  for (size_t i=0;i<h3det.size();i++)
    {
      if (h3det.at(i)->GetEntries())
	{
	  c2->cd(i+1);
	  h3det.at(i)->GetXaxis()->SetRangeUser(0.5,9.);
	  h3det.at(i)->DrawCopy();
	  if (mcprofile)
	    mcprofile->DrawCopy("Psame");
	}
    }

  //.. and all together
  TCanvas* c3 = new TCanvas();
  c3->cd();
 
  h3det.at(0)->DrawCopy();
  for (size_t i=1;i<h3det.size();i++)
    {
      if (h3det.at(i)->GetEntries())
	{
	  h3det.at(i)->SetLineColor(i+2);
	  h3det.at(i)->DrawCopy("same");
	}
      if (mcprofile)
	mcprofile->DrawCopy("Psame");
    }
  TString pdfname = ofilename;
  pdfname.ReplaceAll("root","pdf");
  //Save everything on file!
  TFile* output = new TFile(ofilename,"RECREATE");
  output->cd();
  for (size_t i=0;i<h3det.size();i++)
    h3det.at(i)->Write(h3det.at(i)->GetName());
  for (size_t i=0;i<h1det.size();i++)
    h1det.at(i)->Write(h1det.at(i)->GetName());
  for (size_t i=0;i<tofs.size();i++)
    tofs.at(i)->Write(tofs.at(i)->GetName());
  for (size_t i=0;i<tofsgammas.size();i++)
    tofsgammas.at(i)->Write(tofsgammas.at(i)->GetName());
  h1->Write(h1->GetName());
  h2->Write(h2->GetName());
  h3->Write(h3->GetName());
  output->Write();
  output->Close();

  //SAVE PDF FILES
  c1->Print(pdfname+"(","pdf");
  c2->Print(pdfname);
  c3->Print(pdfname+")","pdf");

  return 0;
 
}
