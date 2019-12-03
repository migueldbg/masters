#include "../../DarkSide/Project/src/Root_Plot.cc"
#include <TTree.h>
#include <TFile.h>

void PlotHisto()
{

  SetStyle();
  /*  TPC bipolar noise */
  const int nF=5;
  std::string FileName[nF]={ 
    "run_588.root",
    "run_589.root",
    "run_590.root",
    "run_591.root",
    "run_592.root" };


  int Opt=0;
  double xmin[3]={ 0.,-200} ,xmax[3]={60.,200.};
  int nbins=200;
  int meanB=14500;


  TTree *t[nF];
  TFile *f[nF];
  TH1F *h[nF];
  int color[5]={ kRed,kBlue,kGreen,kOrange,kBlack};
  
  int chan=16;
  for (int i=0;i<nF;i++)
    {
      f[i]=new TFile(FileName[i].c_str());
      f[i]->GetObject("reco",t[i]);
      t[i]->SetLineColor(i==0?kBlue:kRed);   t[i]->SetMarkerColor(i==0?kBlue:kRed); t[i]->SetMarkerStyle(i==0?24:20);

      const char *hname=Form("h%d",i);
      h[i]=new TH1F(hname,hname,nbins, xmin[Opt],xmax[Opt] );

      std::string x;
      if (Opt==0 )
	x=Form("baseline_mean[%d]-ymin[%d]",chan,chan);
      else if ( Opt==1 )
	x=Form("%d-ymin[%d]",meanB,chan);
      
      t[i]->Project(hname,x.c_str());
      h[i]->SetDirectory(0);

      h[i]->Scale(1./h[i]->Integral());
    }

  cc[canvas]=new TCanvas(Form("c%d",canvas),"",1.5*ww,1.*wh);
  cc[canvas]->Draw();  

  for (int i=0;i<nF;i++)
    {
      //if (i>1 ) continue;
      h[i]->SetLineColor(color[i]);
      h[i]->GetXaxis()->SetTitle((Opt==0?"baseline_mean-ymin":Form("%d-ymin",meanB)));
      h[i]->GetYaxis()->SetTitle("A.U.");
      h[i]->SetMaximum((Opt==0?0.15:0.3));
      h[i]->DrawCopy(i==0?"HIST":"HISTSAME");
      l.SetTextColor(color[i]);
      l.DrawLatex(0.6,0.85-0.05*double(i),FileName[i].c_str());
    }
  l.SetTextColor(kBlack);
  l.DrawLatex(0.2,0.2,Form("chan%d",chan));
  gPad->SetLogy();
  canvas++;
  

  

}
