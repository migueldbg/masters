#include "Utils.C"



TTree *t[2];


void Load(int run, double time=1000 )
{
  SetStyle();
  int iF=0;
  frun_i[iF]=run; ftime_i[iF]=10.;  IsTPC=false;
  iF_plot.push_back(iF);

  std::string FileName=Form("reco/run_%d.root",run);
  TFile *f=new TFile(FileName.c_str());
  f->GetObject("reco",t[iF]);
  Summary(t[iF],iF);
  //SetHistos(t[iF],iF);
  //PlotHistos();
}
