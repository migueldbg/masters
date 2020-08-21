
#include <iostream>
#include <string>
#include <sstream>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TPad.h>
#include <TH2D.h>


using namespace std;
int draw_netrons_fromTPC(int run_number,string suffix="",string filepath="/home/pandola/Level1/");

int draw_netrons_fromTPC(int run_number,string suffix,string filepath){

  stringstream file_name;
  if (suffix.size()==0)
    file_name<<filepath.c_str()<<"/run_"<<run_number<<".root";
  else
    file_name<<filepath.c_str()<<"/run_"<<run_number<<"_"<<suffix.c_str()<<".root";

  cout<<" file_name "<<file_name.str()<<endl;

  TFile *_file0 = TFile::Open(file_name.str().c_str());
  if(! _file0) return -1;
  TTree *_reco= (TTree*) _file0->Get("reco");
  if(! _reco) return -2;
  cout<<" Entries "<<_reco->GetEntries()<<endl;

  TCut gf90  = TCut("clusters[0].f90>0. && clusters[0].f90<1."); 
  TCut gf9045= TCut("clusters[0].f90>0.45 && clusters[0].f90<1. && number_of_clusters==1 ");
  TCut gNcl1 = TCut("number_of_clusters==1");
  int n_entries=_reco->GetEntries();
  TCanvas *c1= new TCanvas("first_plots");
  _reco->SetMarkerColor(kBlue);
  c1->Divide(2,2);
  TPad *_gpad= (TPad *)  c1->cd(1);
  _gpad->SetLogy();
  _gpad->SetGridx();
  _gpad->SetGridy();

  int n_gNcl1=_reco->GetEntries(gNcl1);
  int n_gNcl1_gf90= _reco->Draw("clusters[0].f90",gNcl1&&gf90);
  
  _gpad= (TPad *)  c1->cd(2);
  _gpad->SetGridx();
  _gpad->SetGridy();
  _reco->Draw("clusters[0].f90:clusters[0].charge",gNcl1&&gf90);
  _reco->SetMarkerColor(kRed);
  int n_gNcl1_gf90_gf9045= _reco->Draw("clusters[0].f90:clusters[0].charge",gNcl1&&gf90&&gf9045,"same");
  _reco->SetMarkerColor(kBlue);


  _gpad= (TPad *)  c1->cd(3);
  _gpad->SetGridx();
  _gpad->SetGridy();
  TH2D *_h2=new TH2D("cdf_time_xmin","cdf_time_xmin",100,2000.,9000.,100,4800., 5000. );
  for(int i_ch=1; i_ch<9; i_ch++){
    _reco->SetLineColor(kRed+i_ch);
    if(i_ch==1) 
      _reco->Draw(Form("clusters[0].cdf_time:xmin[%d]>>cdf_time_xmin",i_ch),gNcl1&&gf90&&gf9045,"L");
    else 
      _reco->Draw(Form("clusters[0].cdf_time:xmin[%d]>>cdf_time_xmin",i_ch),gNcl1&&gf90&&gf9045,"Lsame");
  }
  _h2->Draw("colz");


  _gpad= (TPad *)  c1->cd(4);
  _gpad->SetGridx();
  _gpad->SetGridy();
  _reco->SetMarkerColor(kBlue);
  _reco->Draw("ymax[9]-baseline_mean[9]:ymax[10]-baseline_mean[10]",gNcl1&&gf90&&gf9045);
  cout<<" Entries :           "<<n_entries<<"\n"
      <<" gNcl1               "<<n_gNcl1<<" "<< n_gNcl1/float(n_entries)<<"\n"
      <<" n_gNcl1_gf90        "<<n_gNcl1_gf90<<" "<< n_gNcl1_gf90/float(n_entries)<<"\n"
      <<" n_gNcl1_gf90_gf9045 "<<n_gNcl1_gf90_gf9045<<" "<<n_gNcl1_gf90_gf9045/float(n_entries)
      <<endl; 


  return 0;
}



