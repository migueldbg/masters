#include "Root_Plot.cc"
#include <TFile.h>

#include <stdio.h>
#include <math.h>

#include <iostream>     // std::cout
#include <fstream>      // std::ifstream

#include <TFitResult.h>

#define UNDEF -100

const double FADCWidth=2.;
bool debug=false;


namespace LArProperties {
  
  static const double rhoLAr=1400.; //[kg/m3]
  static const double g=9.81;        //[m/s2]
  static const double dH=0.3;        // [m]  Depth of Gas Pocket
  static const double pLAr_ReD=1000.e2; //1000.e2; LUCIANO                // Pa 
  static const double pGAr_ReD=pLAr_ReD+  rhoLAr*g*dH; // Pa 
  
  //-------------------------------------------
  static const double Tc=150.687, pc=4.863 /*MPa*/;
  static const double Tt=83.8058, pt=0.068891/*MPa*/;

  double GetBoilingLine( double T /*K*/)  // Boiling line   [Pa]
  {
    double phi=(1.-T/Tc);
    double a1=-5.9409785, a2=1.3553888, a3=-0.46497607, a4=-1.5399043;
    double p=pc*exp(Tc/T  *  (a1*phi+a2*pow(phi,1.5)+a3*pow(phi,2.)+a4*pow(phi,4.5)));
    return p*1.e6;
  }

  double GetMeltingLine( double T /*K*/)  // Boiling line   [Pa]
  {

    double a1=-7476.2665;
    double a2=9959.0613;
    double phi=T/Tt;
    double p=pt*(1.+a1*(pow(phi,1.05)-1.) + a2*(pow(phi,1.275)-1.));
    return p*1.e6;
  }

  double GetSublimationLine( double T /*K*/)  // Boiling line   [Pa]
  {
    double a1=-11.391604;
    double a2=-0.39513431;
    double phi=T/Tt;
    double p=pt*exp( phi* ( a1*(1.-pow(phi,1.)) + a2*(1.-pow(phi,2.7)) ) );
    return p*1.e6;
  }

  
  double GetBoilingTemperature( double p ) // [Pa] --> [K]
   {
    double Tmin=87.,Tmax=150.;
    
    for (int i=0;i<2000;i++)
      {
	double Ti=Tmin+(Tmax-Tmin)/2000.*double(i);
	double pb=GetBoilingLine(Ti);
	if ( pb>p ) return Ti;
      }    
    return -100;
  }

  //------------------------------------------------------------------
  double GetLArMobility( double T /*K*/, double E /*kV/cm*/ ) // [cm2/V/s]
  {
    //  Global Fit in --> Astro-ph 1508.07059
    double a[6]={ 551.6, 7953.7, 4440.43, 4.29, 43.63, 0.2053 };
    double T0=89.;
    double val=a[0]+a[1]*E+a[2]*pow(E,3./2.)+a[3]*pow(E,5./2);
    val/=(1.+a[1]/a[0]*E+a[4]*E*E+a[5]*E*E*E);
    return val*pow(T/T0,-3./2.);
  }

  double GetLArEpsL( double T /*K*/, double E /*kV/cm*/ )  // [eV]
  {
    //  Global Fit in --> Astro-ph 1508.07059
    double b[4]={ 0.0075, 742.9, 3269.6 , 31678.2 };
    double T0=87.;
    double val=b[0]+b[1]*E+b[2]*E*E;
    val/=(1.+b[1]/b[0]*E+b[3]*E*E);
    return val*(T/T0);    
  }

  double GetLArDL( double T/*K*/, double E /*kV/cm*/)   // [mm2/us]
  {
    //  Global Fit in --> Astro-ph 1508.07059
    return GetLArMobility(T,E)*GetLArEpsL(T,E)*1.e-2;
  }

  double GetLArDriftSpeed(double T /*K*/,double E /*kV/cm*/ ) // [mm/us]
  {
    return GetLArMobility(T,E)*E*1.e-2;
  }

  //------------------------------------------------------------------
  
  //  MC --> Astro-ph 1803.05329   cm/us
  static const double Td_mc[19]={ 0.142, 0.414, 0.851, 1.288, 1.868, 2.506, 3.262, 3.877, 4.563, 5.130, 5.650, 6.170, 6.738, 7.317, 7.754, 8.262, 8.853, 9.480, 9.988};
  static const double vGAr_mc[19]={ 0.174,  0.234,  0.286,  0.308, 0.339, 0.365,  0.390,  0.410,  0.447, 0.486,0.529,  0.580,  0.634,  0.694, 0.737,0.794,  0.850,  0.916, 0.964};
  TGraph *tGArDriftSpeed_mc  =new TGraph(19,Td_mc,vGAr_mc);

 
  //  Data --> Astro-ph 1508.07059 cm/us
  static const double Td_data[12]={ 0.027, 0.054, 0.110, 0.203, 0.390, 0.654, 1.098, 1.809, 2.923, 4.292,  5.609 , 8.0};
  static const double vGAr_data[12] ={ 0.137, 0.150, 0.173, 0.200, 0.230, 0.261, 0.296, 0.321, 0.355, 0.377,  0.423,0.55}; // cm/us
  TGraph *tGArDriftSpeed_data=new TGraph(12,Td_data,vGAr_data);


  double GetGArN()
  {
    const double pGAr=pGAr_ReD;//Pa
    const double TGAr=GetBoilingTemperature(pGAr);

    const double Na=6.022e23;
    const double A=39.948;
    double N=1.723e-3*(pGAr/1013e2)/(TGAr/293.15)*Na/A;   /*atoms/cm3 */
    //printf("N=%e pGAr=%e  TGar=%e  \n",N,pGAr,TGAr);
    return N;
  }

  double GetGArTd(double E /*kV/cm*/)
  {
    double N=GetGArN();
    double Td=E*1.e3/N/1.e-17;    
    return Td;
  }

  double GetGArDriftSpeed( double Td, bool UseMC) //[cm/us]
  {
    if ( UseMC ) return tGArDriftSpeed_mc->Eval(Td);
    else         return (Td<9.?tGArDriftSpeed_data->Eval(Td):0.);
  }

  double GetGArDriftSpeed_kV(double E /*kV/cm*/ ) // [cm/us]
  {
    double Td=GetGArTd(E);
    return GetGArDriftSpeed(Td,true);
  }

  //------------------------------------------------------------------

  static const double Td_yield[20]={ 4.523, 5.505, 6.360, 7.172, 8.155, 9.351, 10.761, 12.171, 13.324, 14.606, 15.760, 17.043, 17.814, 18.713, 19.185, 19.658, 20.044, 20.389, 20.734, 21.037};
  static const double ELYield[20] ={ 0.042, 0.127, 0.223, 0.325, 0.437, 0.570, 0.714,  0.848,  0.960,  1.104,   1.253,  1.467,  1.659,  1.926,  2.124,  2.317,  2.509,  2.766,  3.022,  3.338};
  TGraph *tELYield_mc=new TGraph(20,Td_yield,ELYield);
  
  double GetReducedLuminiscenceYield(double Td)// Y/N  -->  [10^{-17} ph/ electron cm^{2} atom^{-1} ]
  {
    return tELYield_mc->Eval(Td);
  }
  double GetLuminiscenceYield_kV(double E/*kV/cm*/)  //--> ph/cm 
  {
    double N=GetGArN();
    double Td=GetGArTd(E);
    return tELYield_mc->Eval(Td)*N*1.e-17;
  }
  double GetLuminiscenceYield(double Td)  //--> ph/cm 
  {
    double N=GetGArN();
    return tELYield_mc->Eval(Td)*N*1.e-17;
  }

}



//--------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------


double DLFunc(double *v,double *par)
{
  double T=par[0];
  double E=v[0];
  return LArProperties::GetLArDL(T,E);
}
double EpsLFunc(double *v,double *par)
{
  double T=par[0];
  double E=v[0];
  return LArProperties::GetLArEpsL(T,E);
}
double SpeedFunc(double *v,double *par)
{
  double T=par[0];
  double E=v[0];
  return LArProperties::GetLArDriftSpeed(T,E);
}

double SpeedGArFunc(double *v,double *par)
{
  double Td=v[0];
  bool UseMC=int(par[0]);
  return LArProperties::GetGArDriftSpeed(Td,UseMC);
}


double FuncRatioField(double *v, double *par)
{
  double F=v[0];
  
  const double r0=0.6409,r1=518.;
  const double f=0.755;

  double r    =r0*exp(-F/r1);
  double r200 =r0*exp(-200./r1);


  double s2Overs1    = (1.-r)/(r+f);
  double s2Overs1_200= (1.-r200)/(r200+f);

  return log10(s2Overs1/s2Overs1_200);
}




void Plot()
{


  SetStyle();

  //--------------------------------------------------------------------------------------------
  //   S1 field quenching from ReD
  double F[6]={50.,100.,200,400.,700.,1000.};
  double eF[6]={0,0,0, 0,0,0};
  double S1[6]={ 0.944, 0.916, 0.847, 0.760, 0.667, 0.601 };
  double eS1[6]={ 0.005, 0.005,0.005, 0.005,0.005,0.005};
  TGraphErrors tQ(6,F,S1,eF,eS1);
  
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.0*ww,1.2*wh);
  cc[canvas]->Divide(2,1);
  cc[canvas]->Draw();  
  
  cc[canvas]->cd(1);

  TH2F *hFrameQ=new TH2F("hFrameQ","",1000,0.,3500.,1000,0.5,1.);
  hFrameQ->GetXaxis()->SetTitle("Field [V/cm]");
  hFrameQ->GetYaxis()->SetTitle("S1/S1^{0}");
  hFrameQ->DrawCopy();

  tQ.DrawClone("P");

  TF1 *fR=new TF1("f","([0]*exp(-x/[1])+[2])/([0]+[2])",0.,10000.);
  fR->FixParameter(2,0.755);
  fR->SetParameter(0,0.6);
  fR->SetParameter(1,200.);
  tQ.Fit(fR,"QN","",0.,10000.);
  fR->DrawCopy("SAME");

  printf("Field quenching parameters r0=%5.4lf r1=%5.4lf | %5.4lf \n",fR->GetParameter(0),fR->GetParameter(1),log10(fR->Eval(6000.))-log10(fR->Eval(200.)));


  cc[canvas]->cd(2);

  TH2F *hFrameQ2=new TH2F("hFrameQ","",1000,200.,5000.,1000,0.,0.5);
  hFrameQ2->GetXaxis()->SetTitle("Field [V/cm]");
  hFrameQ2->GetYaxis()->SetTitle("log_{10} S2/S1 - Val at (200 V/cm)");
  hFrameQ2->DrawCopy();
  TF1 *fQ2=new TF1("fQ2",FuncRatioField,200.,6000.);
  fQ2->DrawCopy("SAME");

  canvas++;
  //--------------------------------------------------------------------------------------------


  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.0*ww,1.2*wh);
  cc[canvas]->Divide(2,1);
  cc[canvas]->Draw();  
  
  cc[canvas]->cd(1);

  TH2F *hFrame=new TH2F("hFramePhase","",1000,70.,150.,1000,0.1,100.);
  hFrame->GetXaxis()->SetTitle("T [k]");
  hFrame->GetYaxis()->SetTitle("Pressure [bar]");
  hFrame->DrawCopy();
  gPad->SetLogy();

  double pRed[2]={1.,1.1};
  double TRed[2][2]={{UNDEF,UNDEF},{UNDEF,UNDEF}};

  for (int iline=0;iline<3;iline++)
    {
      int color[3]={kRed, kBlue,kGreen+2}; 
      double Tstart=(iline==2?70.:LArProperties::Tt);
      double Tend  =(iline==2?LArProperties::Tt :LArProperties::Tc);
      TGraph *t=new TGraph();
      t->SetLineColor(color[iline]);
      for (int ii=0;ii<1000;ii++)
	{
	  double Ti=Tstart+(Tend-Tstart)/1000.*double(ii); 
	  double pi;
	  if ( iline==0 ) pi=LArProperties::GetBoilingLine(Ti);
	  if ( iline==1 ) pi=LArProperties::GetMeltingLine(Ti);
	  if ( iline==2 ) pi=LArProperties::GetSublimationLine(Ti);
	  
	  t->SetPoint(ii,Ti,pi*1.e-5);
	  for (int jj=0;jj<2;jj++)
	    {
	      if ( iline==1 && pi*1.e-5>pRed[jj] && TRed[jj][0]==UNDEF ) TRed[jj][0]=Ti;
	      if ( iline==0 && pi*1.e-5>pRed[jj] && TRed[jj][1]==UNDEF ) TRed[jj][1]=Ti;
	    }
	  //printf("line=%d  Ti %e  pi %e | %e %e \n",iline,Ti,pi/10.,TRed[0][0],TRed[0][1]);	  
	}
      t->SetLineWidth(2.);
      t->DrawClone("L");
    }

  for (int i=0;i<2;i++)
    {
      TF1 *f=new TF1("f","[0]",TRed[i][0],TRed[i][1]);
      f->SetParameter(0,pRed[i]);
      f->DrawCopy("SAME");
    }
  l.DrawLatex(0.5,0.85,Form("T=%3.2lf-%3.2lf",TRed[0][0],TRed[0][1]));
  l.DrawLatex(0.5,0.80,Form("T=%3.2lf-%3.2lf",TRed[1][0],TRed[1][1]));


  //---
  cc[canvas]->cd(2);
  TGraph *t0=LArProperties::tELYield_mc;
  t0->GetXaxis()->SetTitle("E/N [Td]");
  t0->GetYaxis()->SetTitle("Y/N [10^{-17} ph/ electron cm^{2} atom^{-1} ]");
  t0->DrawClone("ALP");
  l.DrawLatex(0.2,0.85,Form("Yield=%5.2lf ph/cm at 7 Td",LArProperties::GetLuminiscenceYield(7.)));
  canvas++;

  //--------------------------------------------------------------------------------------------

  double T[2]={ 85.,88. };
  double xmin=1.e-2,xmax=10;

  TF1 *fspeed = new TF1("fspeed", SpeedFunc, xmin,xmax,1);
  fspeed->SetNpx(10000);
  
  TF1 *fspeedGAr = new TF1("fspeedGar", SpeedGArFunc,xmin,xmax,1);
  fspeedGAr->SetNpx(10000);

  TF1 *fEpsL = new TF1("fepsl", EpsLFunc, xmin,xmax,1);
  fEpsL->SetNpx(10000);

  TF1 *fDL = new TF1("fdl", DLFunc, xmin,xmax,1);
  fDL->SetNpx(10000);


  
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",1.5*ww,1.5*wh);
  cc[canvas]->Divide(1,1);
  cc[canvas]->Draw();  


  std::string nameV[4]={ "v (LAr) [mm/#mu s]","#eps_{L} (LAr) [eV]",   "D_{L} (LAr) [mm^{2}/#mu s]",  "v (GAr) [cm/#mu s]" };

  int color[2]={ kRed,kBlue};

  for (int iloop=0;iloop<1;iloop++)
    {
      int ipad = iloop+3;
      cc[canvas]->cd(iloop+1);
      
      TF1 *f;
      if ( ipad==0 ) f=fspeed;
      else if (ipad==1 ) f=fEpsL;
      else if ( ipad==2 ) f=fDL;
      else if ( ipad==3 ) f=fspeedGAr;

      if ( ipad==3 ) f->SetParameter(0,double(1));
      else           f->SetParameter(0,T[0]);

      double min=f->Eval(xmin);
      double max=f->Eval((ipad==2?1.:xmax));
      printf("min=%5.2lf max=%5.2lf \n",min,max);

      TH2F *h1=new TH2F(Form("h1%d",ipad),"",1000,xmin,xmax,1000, TMath::Min(min,max), TMath::Max(min,max) );
      h1->GetYaxis()->SetTitle(nameV[ipad].c_str());
      h1->GetXaxis()->SetTitle((ipad==3?"E/N [Td]":"E [kV/cm]"));
      h1->DrawCopy();

      if ( ipad!=3 )
	{
	  gPad->SetLogx();      
	  gPad->SetLogy();      
	}
      for (int i=0;i<2;i++)
	{
	  if ( ipad==3 ) f->SetParameter(0,double(i));
	  else           f->SetParameter(0,T[i]);
	  f->SetLineColor(color[i]);
	  f->DrawCopy("SAME");	 
	  f->SetRange(xmin,xmax); 

	  l.SetTextColor(color[i]);
	  if ( ipad==3 )
	    {
	      const double ratioEps=1.5;
	      const double h_g=0.8,h_e=0.2;
	      double VAnode=5.211;
	      double E=VAnode/(h_g+h_e/ratioEps); 	      
	      double Td=LArProperties::GetGArTd(E);
	      
	      l.DrawLatex(0.2,0.85-0.05*double(i),Form("%5.2e at %5.2lf Td ",f->Eval(Td),Td));
	    }
	  else           
	    {
	      l.DrawLatex(0.2,0.85-0.05*double(i),Form("%5.2e at 200 V/cm",f->Eval(0.2)));
	      l.DrawLatex(0.5,0.25-0.05*double(i),Form("T=%5.2lf [K]",T[i]));
	    }
	}
    }


  canvas++;
}


void Dump()
{
  const double ratioEps=1.5;
  const double h_g=0.8,h_e=0.2;
  double VAnode=5.211;
  double E=VAnode/(h_g+h_e/ratioEps); 	   
  cout << "Pressure: " << LArProperties::pGAr_ReD/100. << " mbar" << endl;
  double Td=7.;
  //cout << "EL field: " << E << " kV/cm" << endl;
  //cout << "Townsend: " <<  LArProperties::GetGArTd(E) << " Td" << endl;
  cout << "EL yield: " <<  LArProperties::GetLuminiscenceYield(Td) << " ph/cm" << endl;
}
