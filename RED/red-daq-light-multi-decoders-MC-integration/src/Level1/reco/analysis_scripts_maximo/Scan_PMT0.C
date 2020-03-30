#include "Utils.C"

//NOTE:  Coin/Be in PMT0 or TPC slave mode monitor changes in the 
//           convolution of detection efficiency + cone intersection.
//       Coin/Monitor as well as long as we stay at the same energy.
//
//  The neutron defficiency problem is a STRONG problem. The neutron beam has not 
//  been characterized in this shift. And without that, High/Low mass results are
//  meaningless. 
//  Take a look a little bit more, but otherwise do MC for prospects.

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//   E=28 MeV         
//        PMT0 placed between TPC and scattering chamber 
//           640/654     -> Trigger :  ( SiTel AND PMT0 ) OR SiMon 
//           645         -> Trigger :  ( SiTel OR SiMon )
// 
//   E=27 MeV     
//        PMT0 placed in the LSci ring  delta_Z=6-8 cm ( on Horizontal Bar)
//           737             -> Trigger :  ( SiTel AND PMT0 ) OR SiMon 
//           735/736         -> Trigger :  ( SiTel OR SiMon )
//---------------------------------------------------------------------
//---------------------------------------------------------------------


// TBD:  1) Calibrate the fastE Be blob
//       2) Repeat all calculations in the macros with fastE
//       3) Get the latest data with runs at 27/28/29 MeV and trigger on (Si OR Monitor) / (Si And PMT ) OR Monitor to cross check.
//       4) Introduce a consistency check of ESlow/EFast to warn when the Amplifier is dying so the run cannot be taken.
//       5) Get the scaler numbers from log book to discriminate between runs with high dead time or not.


void SlaveModePMT()
{
  double E0=28.;
  SetStyle();

  const int nF_in=3;
  int    scaler[nF_in]={  1181,  8450,  31574  };//, 230064,  100646, 1845};
  int    nevt_log[nF_in]={1079,  7761  , 31564 };//, 229617,  100855, 1882  };
  double time[nF_in]=  {  1481 , 10000 , 2015  };//, 2923., 1357.    , 10000. };
  int    run_i[nF_in]= { 654  ,  640  ,  645  };//, 735  , 736      , 737  };
  int    plot_i[nF_in]={   1,   0,   1  };//,  0 ,  0 , 0 };

  nF=nF_in; IsTPC=false; iF_plot.clear();
  for (int iF=0;iF<nF;iF++) { frun_i[iF]=run_i[iF]; ftime_i[iF]=time[iF]; if ( plot_i[iF]) iF_plot.push_back(iF); }
  Analyze();
  
}



//---------------------------------------------------------------------
//---------------------------------------------------------------------
//   E=34 MeV       Trigger :  ( SiTel AND PMT0 ) OR SiMon  ,       except run 666 ->  SiTel OR SiMon  (TPC slave)
//           663-666 -> PMT0 placed behind TPC (between TPC and Lsci)  
//           669     -> PMT0 placed in the ring  
//                                                             
//---------------------------------------------------------------------
//---------------------------------------------------------------------

void Ene34MeV() // IsBe/IsMonitor (needs to be adjusted)
{
  double E0=34.;
  SetStyle();

  const int nF_in=4;
  double time[nF_in]=  { 1590., 9634., 7417. , 3033. };
  int    run_i[nF_in]= { 666, 663, 664, 669 };
  int    plot_i[nF_in]={   1,   0,   0,   1 };

  nF=nF_in; IsTPC=false; iF_plot.clear();
  for (int iF=0;iF<nF;iF++) { frun_i[iF]=run_i[iF]; ftime_i[iF]=time[iF]; if ( plot_i[iF]) iF_plot.push_back(iF); }
  Analyze();
  
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//   E=28 MeV       Trigger :  ( SiTel AND PMT0 ) OR SiMon |     PMT0 placed in different places along the beam line 
//                                                             
//---------------------------------------------------------------------
//---------------------------------------------------------------------
void ScanBeamLine()
{
  SetStyle();
  const double E0=28.;

  //   E=28 MeV       Trigger :  ( SiTel AND PMT0 ) OR SiMon
  //
  //           654     ->   PMT0 placed between TPC and Scattering chamber (central position)
  //           662     ->   PMT0 placed between TPC and LSci ring
  //           649     ->   PMT0 placed in the center of LSci ring 
  //           674     ->    ( New SiTel collimator : hole size but different hole position) 
  //           677     ->   PMT0 placed in the center of LSci ring, 8 cm vertical shift, down (  New SiTel collimator )
  //           679     ->   PMT0 placed in the center of LSci ring, 8 cm vertical shift, Up   (  New SiTel collimator )  --> Differences in <SiE> are taken into account but not plotted

 
  const int nF_in=6;
  double time[nF_in]={ 1481.,   3693., 1975., 3283., 1849. ,1495. };
  int    run_i[nF_in]={ 654, 662, 649, 674, 677, 679 };
  int    plot_i[nF_in]={   1,   0,  0,   0,   0,   1 };
  
  nF=nF_in; IsTPC=false; iF_plot.clear();
  for (int iF=0;iF<nF;iF++) { frun_i[iF]=run_i[iF]; ftime_i[iF]=time[iF]; if ( plot_i[iF]) iF_plot.push_back(iF); }
  Analyze();


}

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//   E=28 MeV       Trigger :  ( SiTel AND PMT0 ) OR SiMon |     PMT0 placed in front of Scattering chamber   | Collimator 0
//                                                              Run 662 ( between TPC and LSci Ring )
//---------------------------------------------------------------------
//---------------------------------------------------------------------
void ScanHorizontal_ScatteringChamber()
{
  SetStyle();

  const int nF_in=7;
  double time[nF_in]    ={ 1481., 1933., 1356.+864.,1815.,1732., 3693.,1015.};
  int run_i[nF_in]      ={ 654  , 655  , 656       , 659 , 660 , 662  ,645  };
  int    plot_i[nF_in]  ={   1,   0    ,   0       ,   0 ,   0,   0   , 1   };
  double delta_cm[nF_in]={   0,   5.   ,  10.      ,  -5., -10.,  0.  , 0   }; 


  nF=nF_in;IsTPC=false; iF_plot.clear();
  for (int iF=0;iF<nF;iF++) { frun_i[iF]=run_i[iF]; ftime_i[iF]=time[iF]; if ( plot_i[iF]) iF_plot.push_back(iF); }
  Analyze();


  TGraphErrors *tg=new TGraphErrors();
  TGraphErrors *tgB1=new TGraphErrors();
  TGraphErrors *tgB2=new TGraphErrors();
  for (int iF=0;iF<nF;iF++)
    {                    
      TGraphErrors *t0;
      if ( iF<5   )    t0=tg;
      else if ( iF==5) t0=tgB1;
      else if ( iF==6) t0=tgB2;
      int ipoint=t0->GetN();
      t0->SetPoint(ipoint,delta_cm[iF],fRateNorm[iF]);
      t0->SetPointError(ipoint,0.,efRateNorm[iF]);
    }
  
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",1.*ww,1.*wh);
  cc[canvas]->Draw();  

  tg->GetXaxis()->SetTitle("Lateral displacement [cm]");
  tg->GetYaxis()->SetTitle("Normalized Rate [Hz /Hz^{mon}] ");
  tg->SetMaximum(2.5);
  tg->DrawClone("AP");
  
  TF1 *f=new TF1("f","[0]*exp(-pow(x-[1],2.)/2./[2]/[2])",-15.,15.);
  double par[3]={1.,0.,1.};
  f->SetParameters(par);
  f->FixParameter(1,0.);
  tg->Fit("f","QN",0,12.);
  f->DrawCopy("SAME");
  
  double mean=f->GetParameter(1);
  double rms=f->GetParameter(2);
  double norm=f->GetParameter(0);
  l.SetTextSize(0.03);
  l.DrawLatex(0.2,0.85,Form("rms=%5.2lf cm",rms));
  l.DrawLatex(0.2,0.80,Form("Rate=%5.2lf [Hz/Hz^{mon}]",IntegrateRate(rms,norm) ));
	 
  l.SetTextColor(kOrange);
  l.DrawLatex(0.2,0.65,"Collimator 0");	  
  
  l.SetTextColor(kBlack);
  l.DrawLatex(0.5,0.85,"PMT0 close to Scatt. Chamber");
  l.SetTextColor(kGreen+2);
  l.DrawLatex(0.5,0.80,"PMT0 close to Scatt. Chamber (Slave)");
  l.SetTextColor(kRed);
  l.DrawLatex(0.55,0.75,"Just Behind  TPC");
  
  
  tgB1->SetLineColor(kRed);
  tgB1->SetMarkerColor(kRed);
  tgB1->DrawClone("P");
  tgB2->SetLineColor(kGreen+2);
  tgB2->SetMarkerColor(kGreen+2);
  tgB2->DrawClone("P");

  canvas++;

}

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//   E=28 MeV       Trigger :  ( SiTel AND PMT0 ) OR SiMon |     PMT0 placed in the LSci ring   | Collimator 1
//---------------------------------------------------------------------
//---------------------------------------------------------------------
// Vertical Scan 
void ScanVertical()
{
  SetStyle();

  //           674     ->   Central Position
  //           675     ->   Central Position
  //           677     ->   8 cm vertical shift, down 
  //           679     ->   8 cm vertical shift, Up   
  //           687     ->   8 cm vertical shift, Up   
  //           688     ->   4 cm vertical shift, Up
  //           689     ->  12 cm vertical shift, Up
  //           690     ->   6 cm vertical shift, Up

  const int nF_in=7;
  int    run_i[nF_in]   ={    674,  677  , 679 ,  687  ,  688  ,  689 , 690   };
  double time[nF_in]    ={  3283., 1849. ,1495., 2047. , 1854. , 1840., 1892. };
  int    plot_i[nF_in]  ={   1,   0    ,   0   ,   0   ,   0   ,    0 ,   1   };
  double delta_Z[nF_in] ={  0.,  -8.  ,    8.  ,   8.  ,   4.  ,   12.,   6.  };

  nF=nF_in; IsTPC=false; iF_plot.clear();
  for (int iF=0;iF<nF;iF++) { frun_i[iF]=run_i[iF]; ftime_i[iF]=time[iF]; if ( plot_i[iF]) iF_plot.push_back(iF); }
  Analyze();


  TGraphErrors *t0=new TGraphErrors();
  for (int iF=0;iF<nF;iF++)
    {        
      int ipoint=t0->GetN();
      t0->SetPoint(ipoint,delta_Z[iF],fRateNorm[iF]);
      t0->SetPointError(ipoint,0.,efRateNorm[iF]);
    }

  cc[canvas]=new TCanvas(Form("c%d",canvas),"",1.*ww,1.*wh);
  cc[canvas]->Draw();  

  t0->GetXaxis()->SetTitle("Vertical displacement [cm]");
  t0->GetYaxis()->SetTitle("Normalized Rate [Hz/Hz^{mon} ");
  t0->DrawClone("AP");
  
  TF1 *f=new TF1("f","[0]*exp(-pow(x-[1],2.)/2./[2]/[2])",-15.,15.);
  double par[3]={1.,6.,1.};
  f->SetParameters(par);
  
  t0->Fit("f","QN",0,12.);
  f->DrawCopy("SAME");
  
  double mean=f->GetParameter(1);
  double rms=f->GetParameter(2);
  double norm=f->GetParameter(0);
  
  l.SetTextSize(0.03);
  l.DrawLatex(0.2,0.85,Form("mean=%5.2lf",mean));
  l.DrawLatex(0.2,0.80,Form("rms=%5.2lf",rms));
  l.DrawLatex(0.2,0.75,Form("Total Rate=%5.2lf [Hz/Hz^{mon}]",IntegrateRate(rms,norm) ));

  l.SetTextColor(kRed);
  l.DrawLatex(0.2,0.40,"In LSci Ring");
  l.DrawLatex(0.2,0.35,"Collimator 1");

  canvas++;



}

// Horizontal Scan 
void ScanHorizontal()
{
  SetStyle();

  const double E0=28.;

  
  // 711,712 SlowE amplifier is malfunctioning
  const int nF_in=10;
  const int    run_i[nF_in]={  687 , 706   ,  707 , 708 , 709 , 710, 711 ,712  ,713   ,714 };
  const double time[nF_in] ={  2047., 926. , 1793.,1099.,1710.,1225,1234.,1047.,1289.,1120.};
  const int    plot_i[nF_in]={     1,   1  ,   0  ,   0 ,   0 ,  0 ,    0,    0,    0,   1 };
  
  const double delta_Z[nF_in]={ 8.,8. ,8.  ,8. ,8., 8.,8.,8.,8.,8.};
  const double delta_X[nF_in]={ 0.,0.,-8. , 8.,-4.,+4.,0.,0.,0.,0.};

  nF=nF_in; IsTPC=false; iF_plot.clear();
  for (int iF=0;iF<nF;iF++) { frun_i[iF]=run_i[iF]; ftime_i[iF]=time[iF]; if ( plot_i[iF]) iF_plot.push_back(iF); }
  Analyze();

  TGraphErrors *tg=new TGraphErrors();
  TGraphErrors *tgB=new TGraphErrors();
  

  for (int iF=0;iF<nF;iF++)
    {              
      TGraphErrors *t0=(iF==0?tgB:tg);
      int ipoint=t0->GetN();
      t0->SetPoint(ipoint,delta_X[iF],fRateNorm[iF]);
      t0->SetPointError(ipoint,0.,efRateNorm[iF]);
    }

  cc[canvas]=new TCanvas(Form("c%d",canvas),"",1.*ww,1.*wh);
  cc[canvas]->Draw();  

  tg->SetMaximum(0.18);	
  tg->GetXaxis()->SetTitle("Horizontal displacement [cm]");
  tg->GetYaxis()->SetTitle("Normalized Rate [Hz/Hz^{mon}] ");
  tg->DrawClone("AP");
  
  tgB->SetMarkerColor(kOrange);
  tgB->DrawClone("P");
	  
	  
  l.SetTextColor(kRed);
  l.DrawLatex(0.2,0.85,"In LSci Ring");
  l.DrawLatex(0.2,0.80,"Collimator 1, Z=6-8 cm");
  l.SetTextColor(kOrange);
  l.DrawLatex(0.2,0.75,"Vertical Scan maximum");

  canvas++;


}
