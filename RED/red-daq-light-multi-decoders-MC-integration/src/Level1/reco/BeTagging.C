#include "Utils.C"
#include "Math/SpecFunc.h"
#include "Math/Functor.h"
#include "Fit/Fitter.h"
#include "TF2.h"

#define UNDEF -100

//  E = 28 MeV

//  Goal : 1) Determine which Be are giving coincidences in the TPC and PMT for different positions
//         2) Benchmark the LowBe blob energy cut
//
//
// Banana runs with High Threshold
//          Coll   Run      Coin
//          ----------------------
//           0     645      PMT0 close central
//           1     715      TPC nominal
//           1     717      TPC 4 spaces Up
//           1     719      TPC 8 spaces Up

//  Require a neutron in the TPC/PMT in time coincidence with a Be in SiTel


void BeTaggingSlave()
{
  SetStyle();
  std::string NameSet[3]={   "TPC 4 spaces up C1 ", "PMT0 close central C0" };
  const int nF_in=2;
  const int run_i[2]={717, 645};

  TH1F   *hTime[nF_in], *hSiE[2][nF_in], *hSiTel[2][nF_in];

  int color[4]={kBlue, kRed, kGreen+2, kOrange };

  nF=nF_in;   // nF -> Defined in "Utils.C", indicates the actual number of runs.
  for (int iF=0;iF<nF;iF++)
  {
    iF_plot.push_back(iF);
    frun_i[iF] = run_i[iF];   ftime_i[iF] = 10.;  IsTPC = (iF==1?false:true);

    TTree *t;
    std::string FileName = Form("runs/run_%d.root", run_i[iF]);
    TFile *f = new TFile(FileName.c_str());
    f -> GetObject("reco", t);
    Summary(t,iF);

    //-------------------------------------------------------------------------------
    //   SiE and SiTel energies for Be Band events  (all and with a neutron tagged by TPC/PMT in time coincidence)
    //-------------------------------------------------------------------------------
    for (int ii = 0; ii<2; ii++)
	  {
	    std::string ok = "1";
	    if ( ii == 1 )
	    {
	      if ( IsTPC ) ok = Form("number_of_clusters == 1 && %s && %s ", IsNeutron_TPC.c_str(), IsTime_TPC.c_str()  );
	      else         ok = Form("%s && %s " , IsNeutron_LSci.c_str(), IsTime_LSci.c_str() );
	    }

	    const char *hnameE = Form("hSiE_%d_%d", ii, iF);
	    hSiE[ii][iF] = new TH1F(hnameE, hnameE, 50, SiE_Range[0], SiE_Range[1]);
	    hSiE[ii][iF] -> Sumw2();
	    t -> Project( hnameE, SiE_dep.c_str(), Form(" %s < 100. && %s && %s", monit.c_str(), IsBeBand.c_str(), ok.c_str()) );
	    hSiE[ii][iF] -> SetDirectory(0);

	    const char *hnameE_C = Form("hSiTel_%d_%d", ii, iF);
	    hSiTel[ii][iF] = new TH1F(hnameE_C, hnameE_C, 100, 15., 25.);
	    hSiTel[ii][iF] -> Sumw2();
	    t -> Project( hnameE_C, SiEnergy.c_str(), Form(" %s < 100. && %s  && %s", monit.c_str(), IsBeBand.c_str(), ok.c_str()) );
	    hSiTel[ii][iF] -> SetDirectory(0);
	  }


    //------------------------------------------------------------------------------
    //----Time coincidence when SiTel is in BeBand and PMT/TPC detected a neutron
    //------------------------------------------------------------------------------
    std::string IsOK_Si = Form("%s <100. && %s && %s >17. && %s <26.", monit.c_str(), IsBeBand.c_str(), SiEnergy.c_str(), SiEnergy.c_str() );


    const char *hname = Form("hTime%d", iF);
    hTime[iF] = new TH1F(hname, hname, 30, 0., 50.);
    hTime[iF] -> Sumw2();

      if ( IsTPC )
	    t -> Project( hname, dt_tpc.c_str(), Form("%s && %s", IsOK_Si.c_str(), IsNeutron_TPC.c_str()) );
      else
	    t -> Project( hname, dt_lsci.c_str(), Form("%s && %s", IsOK_Si.c_str(), IsNeutron_LSci.c_str() ) );

    hTime[iF]->SetDirectory(0);

    f->Close();
  }


  //----------------------------------------
  //----------------------------------------


  cc[canvas] = new TCanvas(Form("c%d",canvas),"", 2.5*ww, 1.*wh);
  cc[canvas] -> Divide(3,1);
  cc[canvas] -> Draw();


  for (int ipad = 0; ipad < 3; ipad++)
  {
    cc[canvas] -> cd(ipad+1);

    std::string xaxis[3] = {"SiE [MeV]", "SiTel [MeV] (E+dE)", "#Delta t [samples]"};

    double xmin[3] = { 11.5, 18., 0. }, xmax[3] = { 20., 25., 50 };
    TH2D *hFrame = new TH2D(Form("hX%d",ipad), "", 1000, xmin[ipad], xmax[ipad] ,1000, 1.e-4, 1.);
    hFrame -> GetXaxis() -> SetTitle(xaxis[ipad].c_str());
    hFrame -> GetYaxis() -> SetTitle("A.U");
    hFrame -> DrawCopy();
    gPad -> SetLogy();

    for (int iF = 0; iF < nF; iF++)
    {
      //----------------------------------

      double Rate[2];
      for (int ii = 0; ii < 2; ii++)
      {
        if ( ipad == 2 && ii == 1 ) continue;
        TH1F *h0;
        if (ipad == 0)      h0 = hSiE[ii][iF];      // ii = 0 corresponds to all events in the BE band.
        else if (ipad == 1) h0 = hSiTel[ii][iF];    // ii = 1 corresponds to events in the BE band with a PMT/TPC detected neutron.
        else if (ipad == 2) h0 = hTime[iF];
        Rate[ii] = h0->Integral();
      }

      //--------------------------------

      for (int ii=0;ii<2;ii++)
      {
        if ( ipad==2 && ii==1 ) continue;
        TH1F *h0;
        if (ipad == 0)      h0 = hSiE[ii][iF];
        else if (ipad == 1) h0 = hSiTel[ii][iF];
        else if (ipad == 2) h0 = hTime[iF];

        int color_i = color[iF];
        double norm = (ii == 0? 1.:Rate[1]/Rate[0]);
        h0 -> Scale(norm/h0->Integral());

        h0 -> SetLineColor(color_i);
        h0 -> SetLineStyle(2-ii);
        h0 -> DrawCopy("HISTSAME");
      }

      l.SetTextColor(color[iF]);
      l.SetTextSize(0.03);
      l.DrawLatex(0.2, 0.85-double(iF)*0.05 ,NameSet[iF].c_str());
    }
  }//ipad

}


//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------




//  E=28 MeV

//  Goal : Check if the Be spectrum for coincidences in the Vertical/Horizontal Scan with PMT0 at LSci ring
//
// Banana runs with High Threshold
//          Coll   Run      Coin
//          ----------------------
//           1     679      PMT0 at LSci ring  ( dz=8 cm    dx=0. in Vertical bar   )
//           1     713      PMT0 at LSci ring  ( dz=6-8 cm  dx=0. in Horizontal bar )

//  Require a neutron in the PMT in time coincidence with a Be in SiTel

void BeTaggingCoincidence()
{
  SetStyle();
  std::string NameSet[2]={   "PMT0 at LSci ring C1 (VS) ", "PMT0 at LSci ring C1 (HS) "};
  const int nF_in=2;
  const int run_i[3]={  679, 713 };

  TH1F   *hTime[nF_in], *hSiE[2][nF_in], *hSiTel[2][nF_in];


  //--------------------------------------------------------------------
  TFile *f1 = new TFile("bananas.root", "READ" );
  for (int iF=0;iF<nF_in;iF++)
    {
      int iBand=1; //Be band
      int iRef=(run_i[iF]>670?1:0);
      {
	std::string name=Form("h1D_%d%d",iRef,iBand);
	hSiTel[0][iF]=(TH1F*) f1->Get(name.c_str())->Clone(Form("hSiTel_%d_%d",0,iF));
	hSiTel[0][iF]->SetDirectory(0);
      }
      {
	std::string name=Form("h1D_SiE_%d%d",iRef,iBand);
	hSiE[0][iF]=(TH1F*) f1->Get(name.c_str())->Clone(Form("hSiE_%d_%d",0,iF));
	hSiE[0][iF]->SetDirectory(0);
      }
    }
  f1->Close();
  //--------------------------------------------------------------------


  int color[4]={kBlue, kRed, kGreen+2, kOrange };

  nF=nF_in;
  for (int iF=0;iF<nF;iF++)
    {
      iF_plot.push_back(iF);
      frun_i[iF]=run_i[iF]; ftime_i[iF]=10.;  IsTPC=false;

      TTree *t;
      std::string FileName=Form("runs/run_%d.root",run_i[iF]);
      TFile *f=new TFile(FileName.c_str());
      f->GetObject("reco",t);
      Summary(t,iF);

      //-------------------------------------------------------------------------------
      //   SiE and SiTel energies for Be Band events  (all and with a neutron tagged by TPC/PMT in time coincidence)
      //-------------------------------------------------------------------------------
      int ii=1;
      std::string ok;
      if ( IsTPC ) ok=Form("number_of_clusters==1 && %s && %s " ,IsNeutron_TPC.c_str() ,IsTime_TPC.c_str()  );
      else         ok=Form("%s && %s " ,IsNeutron_LSci.c_str(),IsTime_LSci.c_str() );

      const char *hnameE=Form("hSiE_%d_%d",ii,iF);
      hSiE[ii][iF]=new TH1F(hnameE,hnameE,50,SiE_Range[0],SiE_Range[1]);
      hSiE[ii][iF]->Sumw2();
      t->Project(hnameE,SiE_dep.c_str(),Form(" %s <100. && %s && %s ",monit.c_str(),IsBeBand.c_str(),ok.c_str()  ));
      hSiE[ii][iF]->SetDirectory(0);

      const char *hnameE_C=Form("hSiTel_%d_%d",ii,iF);
      hSiTel[ii][iF]=new TH1F(hnameE_C,hnameE_C,100,15.,25.);
      hSiTel[ii][iF]->Sumw2();
      t->Project(hnameE_C,SiEnergy.c_str(),Form(" %s <100. && %s  && %s ",monit.c_str(),IsBeBand.c_str(),ok.c_str() ) );
      hSiTel[ii][iF]->SetDirectory(0);


      //------------------------------------------------------------------------------
      //----Time coincidence when SiTel is in BeBand and PMT/TPC detected a neutron
      //------------------------------------------------------------------------------
      std::string IsOK_Si=Form("%s <100. && %s && %s >17. && %s <26.",monit.c_str(), IsBeBand.c_str(), SiEnergy.c_str(), SiEnergy.c_str() );


      const char *hname=Form("hTime%d",iF);
      hTime[iF]=new TH1F(hname,hname,30,0.,50.);
      hTime[iF]->Sumw2();
      if ( IsTPC )
	t->Project(hname,dt_tpc.c_str(),Form("%s && %s " ,IsOK_Si.c_str(), IsNeutron_TPC.c_str() ) );
      else
	t->Project(hname,dt_lsci.c_str(),Form("%s && %s ",IsOK_Si.c_str(), IsNeutron_LSci.c_str() ) );
      hTime[iF]->SetDirectory(0);

      f->Close();
    }


  //----------------------------------------
  //----------------------------------------


  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,1.*wh);
  cc[canvas]->Divide(3,1);
  cc[canvas]->Draw();


  for (int ipad=0;ipad<3;ipad++)
    {
      cc[canvas]->cd(ipad+1);

      std::string xaxis[3]={ "SiE [MeV]", "SiTel [MeV] (E+dE)","#Delta t [samples]" };

      double xmin[3]={ 11.5,17.5,0. }, xmax[3]={ 20., 25., 50};
      TH2D *hFrame=new TH2D(Form("hX%d",ipad),"",1000,xmin[ipad],xmax[ipad] ,1000,1.e-2,1.);
      hFrame->GetXaxis()->SetTitle(xaxis[ipad].c_str());
      hFrame->GetYaxis()->SetTitle("A.U");
      hFrame->DrawCopy();
      gPad->SetLogy();

      for (int iF=0;iF<nF;iF++)
	{
	  for (int ii=0;ii<2;ii++)
	    {
	      if ( ipad==2 && ii==1 ) continue;
	      TH1F *h0;
	      if (ipad==0 )      h0=hSiE[ii][iF];
	      else if (ipad==1 ) h0=hSiTel[ii][iF];
	      else if (ipad==2 ) h0=hTime[iF];

	      int color_i=color[iF];
	      h0->Scale(1./h0->Integral() );
	      if ( ii==0 && ipad<2  ) h0->Scale(5.);

	      h0->SetLineColor((ii==0&&ipad<2?kBlack:color_i) );
	      h0->DrawCopy("HISTSAME");
	    }

	  l.SetTextColor(color[iF]);
	  l.SetTextSize(0.03);
	  l.DrawLatex(0.2,0.85-double(iF)*0.05,NameSet[iF].c_str());
	}
    }//ipad
  canvas++;
}


//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
