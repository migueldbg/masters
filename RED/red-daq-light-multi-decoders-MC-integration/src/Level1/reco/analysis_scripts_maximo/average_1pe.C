#include "TROOT.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include <iostream>
#include "red-daq/EvRec0.hh"
#include "red-daq/EvRaw0.hh"
#include "../Root_Plot.cc"
using namespace std;

TGraph *tSer;
std::vector<int> sipm;

void SetSer(int run)
{
  tSer=new TGraph();
  std::string SerName=Form("ser_%d.cfg",run);
  std::ifstream Ser(SerName);
  int isipm=0;
  for (int i=0;i<48;i++)
    {
      double adc;
      Ser>> adc;
      tSer->SetPoint(tSer->GetN(),double(i),adc);

      sipm.push_back(adc>0?isipm:-100);
      std::cout << adc << "  " << isipm <<"  " << i << std::endl;
      if ( adc>0 ) isipm++;
    }
  Ser.close();
}




int average(TString filename, int run)
{
  SetSer(run);
  SetStyle();

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
  vector<TProfile*> fAverages(imax+1,nullptr);
  vector<TProfile*> fCumulative(imax+1,nullptr);
  vector<TH1D*> fCharge(imax+1,nullptr);
  vector<Int_t> fN(imax+1,0);


  for (Int_t iloop=0;iloop<reco->GetEntries();iloop++)
    {      
      //if (iloop>100) break;
      reco->GetEntry(iloop);
      raw->GetEntry(iloop);
      map<int,vector<double>* > wfs = evRaw->GetWFs();
      for (map<int,vector<double>* >::iterator it=wfs.begin(); 
       it!=wfs.end(); ++it)
	{
	  int i = it->first; //channel number  
	  vector<double>* wf = evRaw->GetWF(i);	  

	  // Select SiPMs
	  int isipm=sipm[i];
	  int nsa= wf->size();

	  Double_t bm = evReco->GetBaseMean().at(i);
	  Double_t br = evReco->GetBaseRMS().at(i);
	  Double_t chargeADC=0.;
	  Double_t peakADC=0;
	  for (Int_t ch=2950;ch<4050;ch++) 
	    {
	      double val=(bm-wf->at(ch));
	      chargeADC+=val;
	      if ( val>peakADC ) peakADC=val;
	    }
	  double charge=chargeADC/tSer->Eval(i);
	  
	  if ( isipm<0 ) continue;
	  if (isipm<4 && br>6 ) continue;
	  if (isipm>3 && br>4 ) continue;
	  if (fCharge.at(isipm) == nullptr)
	    fCharge.at(isipm) = new TH1D(Form("hSer%d",isipm),Form("hSer%d",isipm),
					    300,0,10. );	    
	  fCharge.at(isipm)->Fill(charge);
	  if (charge<0.4 || charge>1.4) continue;
	  printf("iloop=%d isipm=%d (%d)  nsa=%d  (bm,br)=(%3.2lf,%3.2lf)  Charge=%3.2lf  Ser=%3.2lf\n",iloop,isipm,i,nsa,bm,br,charge,tSer->Eval(i));

	  if (fAverages.at(isipm) == nullptr)
	    {
	      fAverages.at(isipm) = new TProfile(Form("h%d",isipm),Form("h%d",isipm),
					     nsa,0,nsa-1);  
	      fCumulative.at(isipm) = new TProfile(Form("hC%d",isipm),Form("hC%d",isipm),
					      nsa,0,nsa-1);  
	    }

	  std::vector<double> *wfc = new vector<double>(wf->size());
	  wfc->at(0) = (bm-wf->at(0))/chargeADC;
	  for (int i = 1; i < wf->size(); i++ ) 
	    {
	      double val=(bm-wf->at(i))/chargeADC;
	      wfc->at(i) = wfc->at(i - 1) + val;
	    }

	  fN.at(isipm)++;
	  for (Int_t ch=0;ch<wf->size();ch++)
	    {
	      fAverages.at(isipm)->Fill(ch,  (wf->at(ch)-bm)/charge );	    
	      fCumulative.at(isipm)->Fill(ch,  wfc->at(ch) );	    
	    }
	}  
    }

  //----------------------------------------------------------------------
  TGraph      *tPETraceMaxChan[28];  // CDF
  for (int i=0;i<28;i++)
    {
      tPETraceMaxChan[i]=new TGraph();
      TGraph *t0=tPETraceMaxChan[i];
      std::string nameT=Form("tPETraceMaxChan_%d_%d",i,run);	          
      t0->SetName(nameT.c_str());

      if (fCumulative.at(i))
	{
	  TProfile *h0=fCumulative.at(i);
	  for (int ix=0;ix<h0->GetNbinsX();ix++)
	    {
	      double y=h0->GetBinContent(ix);
	      double x=h0->GetBinCenter(ix);
	      if ( x-2950 >2000 ) continue;
	      t0->SetPoint(t0->GetN(),(x-2950),y);
	    }
	}
    }

  TFile* ofile = new TFile(Form("%d.results.root",run),"RECREATE");
  ofile->cd();
  for (size_t i=0;i<28;i++)
    tPETraceMaxChan[i]->Write();
  ofile->Close();

  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  //----------------------------------------------------------------------
  int color[4]={ kBlue,kRed,kGreen+2,kOrange};


  for (int iPlot=0;iPlot<3;iPlot++)
    {
      int isipm=0;
      for (int top=0;top<2;top++)
	{
	  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,1.5*wh);
	  if ( top==1 ) cc[canvas]->Divide(3,2);
	  cc[canvas]->Draw();  

	  for (int ipad=0;ipad<(top==0?1:6);ipad++)
	    for (int ig=0;ig<4;ig++)
	      {
		cc[canvas]->cd(ipad+1);

		if (iPlot<2 )
		  {
		    TProfile *h0=(iPlot==0?fAverages.at(isipm):fCumulative.at(isipm));
		    if (h0)
		      {
			h0->SetLineColor(color[ig]);
			if ( iPlot==0 ) h0->SetMinimum(-30);
			else            h0->SetMaximum(1.5);
			h0->GetXaxis()->SetRange(2800,6000);
			h0->GetXaxis()->SetTitle("Samples");
			h0->GetYaxis()->SetTitle(iPlot==0?"PDF":"CDF");

			h0->DrawCopy((ig==0?"HIST":"HISTSAME"));
		      }
		  }
		else
		  {
		    TH1D *h0=fCharge.at(isipm);
		    if ( h0 )
		      {
			h0->SetLineColor(color[ig]);
			h0->GetXaxis()->SetTitle("Charge [PE]");
			h0->GetYaxis()->SetTitle("A.U.");
			h0->DrawCopy((ig==0?"HIST":"HISTSAME"));
		      }
		  }
		isipm++;
	      }//graph
	  canvas++;
	}//top
    }//  PDF/CDF/Charge
	

  return 0;
}
