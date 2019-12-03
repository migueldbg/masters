#include "../../DarkSide/Project/src/Root_Plot.cc"
#include <TTree.h>
#include <TFile.h>



// Alpha peaks
double EnergyA[7] = {    5144.3, 5156.59, 5388.2,   5442.8, 5485.56, 5762.7, 5804.82};
double intensityA[7]={ 0.20  ,    1.00,  0.014,     0.15,   1.0 ,  0.31,  1. };
int group[7]        ={  0, 0, 1,1,1,  2, 2 };

double Func(double *v, double *par)
{
  double x=v[0];
  double a=par[0];
  double b=par[1];
  double sigma=par[2];
  double norm[3]={par[3],par[4],par[5]};

 
  double val=0.;
  for (int i=0;i<7;i++)
    {
      double mean=(EnergyA[i]-a)/b;
      val+=( TMath::Gaus(x,mean,sigma)*intensityA[i]*norm[group[i]]);
    }
  return val;
  
}

void Alpha()
{

  SetStyle();

  const int nF=2;
  std::string FileName[nF]={ 
    "reco/run_565.root",   // E
    "reco/run_578.root"};  // Monitor
  int chan[nF]={ 11, 13 };


  double xmin=1800.,xmax=2700.;
  int nbins=1000;


  TTree *t[nF];
  TFile *f[nF];
  TH1F *h[nF];
  int color[5]={ kRed,kBlue,kGreen,kOrange,kBlack};
  
  for (int i=0;i<nF;i++)
    {

      f[i]=new TFile(FileName[i].c_str());
      f[i]->GetObject("reco",t[i]);
      t[i]->SetLineColor(i==0?kBlue:kRed);   t[i]->SetMarkerColor(i==0?kBlue:kRed); t[i]->SetMarkerStyle(i==0?24:20);

      const char *hname=Form("h%d",i);
      h[i]=new TH1F(hname,hname,nbins, xmin,xmax );
      h[i]->Sumw2();

      std::string x;
      if ( chan[i]==13) x=Form("baseline_mean[%d]-ymin[%d]",chan[i],chan[i]);
      else              x=Form("ymax[%d]-baseline_mean[%d]",chan[i],chan[i]);

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
      h[i]->GetXaxis()->SetTitle("baseline_mean-ymin");
      h[i]->GetYaxis()->SetTitle("A.U.");
      h[i]->SetMaximum(1.);
      h[i]->SetMinimum(1.e-6);

      h[i]->DrawCopy(i==0?"HIST":"HISTSAME");
      l.SetTextColor(color[i]);
      l.DrawLatex(0.6,0.85-0.05*double(i),FileName[i].c_str());

      TF1 *f=new TF1("f",Func,xmin,xmax,6);
      f->SetParameter(0,0.);
      f->SetParameter(1,(i==0?2.57:2.33));
      f->SetParameter(2,5.);
      for (int ii=0;ii<3;ii++) f->SetParameter(ii+3,1./30.);
      h[i]->Fit(f,"QN","",xmin,xmax);
      f->ReleaseParameter(2);
      h[i]->Fit(f,"QN","",xmin,xmax);

      f->SetLineColor(color[i]);
      f->DrawCopy("SAME");
      
      printf("Calibration (%s) : a=%5.4lf b=%5.4lf \n",(i==0?"Si E":"Si Monitor"),f->GetParameter(0),f->GetParameter(1));
    }
  l.SetTextColor(kBlack);

  gPad->SetLogy();



  canvas++;
  

  

}
