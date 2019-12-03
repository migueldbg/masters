#include "TROOT.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include <iostream>
#include "red-daq/EvRec0.hh"
#include "red-daq/EvRaw0.hh"

using namespace std;

int average(TString filename, int run)
{
  TFile* f = new TFile(filename);
  if (!(f->IsOpen()))
    {
      cout << "could not open file: " << filename << endl;
      return 0;
    }
  TTree* reco = (TTree*) f->Get("reco");
  TTree* raw = (TTree*) f->Get("raw");
  if (!raw)
    {
      cout << "TTree of rawdata not found! " << endl;
      cout << "Please re-process with RedLevel1 -w " << endl;
      return 0;
    }
  EvRec0* evReco = new EvRec0();
  EvRaw0* evRaw = new EvRaw0();
  reco->SetBranchAddress("recoevent",&evReco);
  raw->SetBranchAddress("rawevent",&evRaw);
  cout << reco->GetEntries() << " " << raw->GetEntries() << 
    " entries" << endl;
 
  raw->GetEntry(0);
  //Check what's the last channel:
  int imax = evRaw->GetWFs().rbegin()->first;	  
  cout << "Max channel is: " << imax << endl;
  vector<TH1D*> fAverages(imax+1,nullptr);
  vector<Int_t> fMean(imax+1,0);

  for (Int_t iloop=0;iloop<reco->GetEntries();iloop++)
    {      
      cout << iloop << endl;
      reco->GetEntry(iloop);
      raw->GetEntry(iloop);
      map<int,vector<double>* > wfs = evRaw->GetWFs();
      for (map<int,vector<double>* >::iterator it=wfs.begin(); 
       it!=wfs.end(); ++it)
	{
	  int i = it->first; //channel number  
	  vector<double>* wf = evRaw->GetWF(i);	  
	  //cout << "Leggo la wf: " << wf << endl;
	  int nchannels= wf->size();
	  //cout << "Channels: " << nchannels << endl;
	  if (fAverages.at(i) == nullptr)
	  {
	    fAverages.at(i) = new TH1D(Form("h%d",i),Form("h%d",i),
				       nchannels,0,nchannels-1);  
	    cout << "Registered histo " << i << endl;
	  }
	  Double_t aCharge = evReco->GetChargeTot();
	  //cout << iloop << " " << i << " " << aCharge << endl;
	  if (aCharge>400 && aCharge<800)
	    {
	      for (Int_t ch=1;ch<=fAverages.at(i)->GetNbinsX();ch++)
		{
		  Double_t val = fAverages.at(i)->GetBinContent(ch);
		  val += wf->at(ch-1)/aCharge;
		  fAverages.at(i)->SetBinContent(ch,val);
		}
	      fMean.at(i)++;
	    }	  
	}  
    }
  TFile* ofile = new TFile(Form("averages_run_%d.root",run),"RECREATE");
  ofile->cd();
  for (size_t i=0;i<fAverages.size();i++)
    {
      if (fMean.at(i))
	fAverages.at(i)->Scale(1./fMean.at(i));
      TString title;
      title.Form("ch%d, %d events",i,fMean.at(i));
      if (fAverages.at(i))
        {
      	  fAverages.at(i)->SetTitle(title);
          fAverages.at(i)->Write(fAverages.at(i)->GetName());
        }
    } 
  
  TCanvas* c1 = new TCanvas();
  c1->Divide(4,4);
  for (size_t i=0;i<fAverages.size();i++)
    {
      c1->cd(i+1);
      if (fAverages.at(i))
       {
         //fAverages.at(i)->GetXaxis()->SetRangeUser(2600,3800);
         //fAverages.at(i)->Fit("expo","","",3000,3800);
         fAverages.at(i)->DrawCopy();
         //cout << "Slope: " << i << " " << 
         //fAverages.at(i)->GetFunction("expo")->GetParameter(1) << endl;
       }  
  }
  TCanvas* c2 = new TCanvas();
  c2->Divide(4,4);
  for (size_t i=0;i<fAverages.size();i++)
    {
      c2->cd(i+1);
      if (fAverages.at(i))
	{
         fAverages.at(i)->GetXaxis()->SetRangeUser(2900,3000);
         fAverages.at(i)->DrawCopy();
        }
    }
  
  //Calculate RT
  /*
  vector<Double_t> theMax(fAverages.size(),0);
  vector<Double_t> theMin(fAverages.size(),0);
  for (size_t i=0;i<fAverages.size();i++)
    {
      Double_t amax = fAverages.at(i)->GetBinContent(1);
      Double_t amin = amax;
      for (Int_t ch=2;ch<=fAverages.at(i)->GetNbinsX();ch++)
	{
	  Double_t val = fAverages.at(i)->GetBinContent(ch);
	  if (val > amax)
	    amax = val;
	  if (val < amin)
	    amin = val;
	}
      theMax.at(i) = amax;
      theMin.at(i) = amin;
    }

  //10-90%
  for (size_t i=0;i<fAverages.size();i++)
    {
      Double_t threshold1 = theMax.at(i)-(theMax.at(i)-theMin.at(i))*0.1;
      Double_t threshold2 = theMax.at(i)-(theMax.at(i)-theMin.at(i))*0.9;
      cout << i << " " << theMin.at(i) << " " << theMax.at(i) << " " << threshold1 << " " << 
	threshold2 << endl;
      Int_t i1 = 0;
      Int_t i2 = 0;
      Int_t i3 = 0;
      Int_t i4 = 0;
      for (Int_t ch=2;ch<=fAverages.at(i)->GetNbinsX();ch++)
	{
	  Double_t oldval = fAverages.at(i)->GetBinContent(ch-1);
	  Double_t val = fAverages.at(i)->GetBinContent(ch);
	  if (oldval > threshold1 && val < threshold1 && !i1)
	    i1 = ch;
	  if (oldval > threshold2 && val < threshold2 && !i2)
	    i2 = ch;
	  if (oldval < threshold1 && val > threshold1 && !i3)
	    i3 = ch;
	  if (oldval < threshold2 && val > threshold2 && !i4)
	    i4 = ch;	  
	}
      cout << i1 << " " << i2 << " " << i3 << " " << i4 << endl;
      cout << "Rise time ch:" << i << " = " << (i2-i1) << " samples " << endl;
      cout << "Fall time ch:" << i << " = " << (i4-i3) << " samples " << endl;

    }
  */
  ofile->Write();
  ofile->Close();
  return 0;
}
