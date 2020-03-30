#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "TROOT.h"
#include "red-daq/EvRec0.hh"
#include "Root_Plot.cc"  





TTree* t ;
void Load(int run_in, bool UseStart )
{
  SetStyle();


  std::string FileName=Form("run_%d.root",run_in);
      
  TFile *f=new TFile(FileName.c_str());
  t = (TTree*) f->Get("reco");

  EvRec0* e=new EvRec0();
  t->SetBranchAddress("recoevent",&e);

  TGraph tlsci[2][2];
  TGraph tbanana[2];
  TH1F *hdt=new TH1F("hdt","", 130,  -50, 15.);
  TH1F *hcharge=new TH1F("hcharge","", 100,  0, 500000. );

  TH1F *hpsd[2];
  for (int ipsd=0;ipsd<2;ipsd++)
    hpsd[ipsd]=new TH1F(Form("hpsd%d",ipsd),"", 30,  0.15, 0.45 );
  
  double chargeMin[2]={100.e3,200.e3},chargeMax[2]={200.e3,300.e3};
  int N[3]={0,0,0};
  for (int i=0;i<t->GetEntries();i++)
    {
      t->GetEntry(i);
      // Select Low Blob Be7 events
      RDChannel* SiE=e->GetSiDet(0); // slow
      RDChannel* SidE=e->GetSiDet(2); // slow
      RDChannel* SiE_Fast=e->GetSiDet(1); // fast
      RDChannel* SidE_Fast=e->GetSiDet(3); // fast

      if ( SiE->npeaks>2 ||  SidE->npeaks>2  || SiE_Fast->npeaks>2 ||  SidE_Fast->npeaks>2 ) continue;

      bool IsBe7= ( SidE->ymax[1] > 7000 &&  SidE->ymax[1]<8000 && SiE->ymax[1]>4800. &&  SiE->ymax[1]<5600  );
      if ( SiE->ymax[1] >3000. && SiE->ymax[1]<9000. && SidE->ymax[1]>3000. && SidE->ymax[1]<9000.)
	tbanana[0].SetPoint(tbanana[0].GetN(),SiE->ymax[1], SidE->ymax[1]);
      if ( !IsBe7 ) continue;
      tbanana[1].SetPoint(tbanana[0].GetN(),SiE->ymax[1], SidE->ymax[1]);

      N[0]++;

  
      // Select Time Coincidences
      double tSiTel=(UseStart?SiE_Fast->start[1]: SiE_Fast->xmin[1]); // Fast
      RDChannel* lsci=e->GetLSci(0); 

      bool ok=false;
      double charge=-100.,fprompt[2]={-100,-100},isSat=false;;
      for (int ipeak=1;ipeak<lsci->npeaks;ipeak++)
	{
	  double tlsci=(UseStart?double(lsci->start[ipeak]):double(lsci->xmin[ipeak]));
	  double dt=tlsci-tSiTel;
	  if (1.-lsci->fprompt[ipeak]>0.2 ) hdt->Fill(dt);
	  if ( dt>-50 && dt< 15.) 
	    {
	      if ( lsci->charge[ipeak] > charge ) 
		{
		  charge=lsci->charge[ipeak];
		  fprompt[0]=1.-lsci->fprompt[ipeak];
		  //if ( ipeak==lsci->ipeakPSD ) fprompt[1]=1.-lsci->fxx[5];
		  isSat=lsci->isSat[ipeak];
		}
	      //printf(" E %5.2lf %5.2lf | dE %5.2lf %5.2lf | Lsci=%5.2lf \n",SiE_Fast->xmin[1],SiE->xmax[1],SidE_Fast->xmin[1],SidE->xmax[1],tlsci);
	      ok=true;
	    }
	}
      if ( ok )
	{
	  if ( charge> 2000.) 
	    {
	      for (int ifxx=0;ifxx<2;ifxx++)
		{
		  tlsci[ifxx][0].SetPoint(tlsci[ifxx][0].GetN(),charge,fprompt[ifxx]);
		  if ( isSat ) tlsci[ifxx][1].SetPoint(tlsci[ifxx][1].GetN(),charge,fprompt[ifxx]);
		}
	    }
	  hcharge->Fill(charge);
	  for (int ipsd=0;ipsd<2;ipsd++)
	    if ( charge>chargeMin[ipsd] && charge<chargeMax[ipsd]) hpsd[ipsd]->Fill(fprompt[0]);
	  N[1]++;
	  if ( charge> 2000.  ) N[2]++;
	}
    }
  printf("N=%d %d %d | Eff=%5.2lf %5.2lf \n",N[0],N[1],N[2],double(N[2])/double(N[0]),double(N[1])/double(N[0]));

  cc[canvas]=new TCanvas(Form("c%d",canvas),"",1.5*ww,1.5*wh);
  cc[canvas]->Divide(2,2);
  cc[canvas]->Draw();  

  cc[canvas]->cd(1);
  tbanana[0].SetMarkerStyle(10);
  tbanana[0].GetYaxis()->SetTitle("SidE");
  tbanana[0].GetXaxis()->SetTitle("SiE");
  tbanana[0].DrawClone("AP");
  tbanana[1].SetMarkerStyle(10);
  tbanana[1].SetMarkerColor(kRed);
  tbanana[1].DrawClone("P");
  
  cc[canvas]->cd(2);
  hdt->GetYaxis()->SetTitle("A.U.");
  hdt->GetXaxis()->SetTitle("t_{lsci}-t_{SiTel} [sa]");
  hdt->DrawCopy();

  gPad->SetLogy();
  //hdt->Fit("gaus");

  TH1F *hdt2=(TH1F*) hdt->Clone();

  double xmin=hdt->GetMean()-3.*hdt->GetRMS();
  double xmax=hdt->GetMean()+3.*hdt->GetRMS();
  int imin=hdt->FindBin(xmin),imax=hdt->FindBin(xmax);
  hdt2->GetXaxis()->SetRange(imin,imax);
  l.DrawLatex(0.2,0.85,Form("rms=%5.2lf (sa=%d to %d) [sa]",hdt2->GetRMS(),int(xmin),int(xmax)));

  cc[canvas]->cd(3);

  int color[2][2]={ {kBlue,kRed}, {kBlue+2,kRed+2} };
  for (int ifxx=0;ifxx<2;ifxx++)
    for (int isat=0;isat<2;isat++)
      {
	tlsci[ifxx][isat].SetMaximum(0.35);
	tlsci[ifxx][isat].SetMinimum(0.);
	tlsci[ifxx][isat].SetMarkerSize(0.5);
	tlsci[ifxx][isat].SetMarkerColor(color[ifxx][isat]);
	tlsci[ifxx][isat].GetXaxis()->SetTitle("LSci Charge");
	tlsci[ifxx][isat].GetYaxis()->SetTitle("LSci PSD");
	tlsci[ifxx][isat].DrawClone((ifxx==0&&isat==0?"AP":"P"));
      }

  cc[canvas]->cd(4);
  for (int ipsd=0;ipsd<2;ipsd++)
    {
      hpsd[ipsd]->Scale(1./hpsd[ipsd]->Integral());      
      hpsd[ipsd]->SetMaximum(1.0);
      hpsd[ipsd]->SetLineColor(ipsd==0?kBlue:kRed); hpsd[ipsd]->SetMarkerColor(ipsd==0?kBlue:kRed);
      hpsd[ipsd]->GetYaxis()->SetTitle("A.U.");
      hpsd[ipsd]->GetXaxis()->SetTitle("LSci PSD");
      
      hpsd[ipsd]->DrawCopy((ipsd==0?"E1":"E1SAME"));
      
      hpsd[ipsd]->Fit("gaus","Q0");
      hpsd[ipsd]->GetFunction("gaus")->SetLineColor(ipsd==0?kBlue:kRed);
      hpsd[ipsd]->GetFunction("gaus")->DrawCopy("SAME");
      l.SetTextColor(ipsd==0?kBlue:kRed);
      l.DrawLatex(0.2,0.9-double(ipsd)*0.05,Form("#sigma=%5.3lf (%3.0lf-%3.0lf e3 ADC)",hpsd[ipsd]->GetFunction("gaus")->GetParameter(2), chargeMin[ipsd]/1.e3,
						 chargeMax[ipsd]/1.e3 ) );
    }

  gPad->SetLogy();
  canvas++;
}

