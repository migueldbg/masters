#include "Root_Plot.cc"
#include <TTree.h>
#include <TFile.h>

// Change of Collimators:  Setting 0:  untill 670
//                         Setting 1:  after 670

// Channel mapping settings :  Setting 0 > run=555-637
//                             Setting 1 > run=637-680  run=685-714  run>=730
//                             Setting 2 > run=680-684  run=715-729

// SER -> Run  670

const double u    = 931.4940954;        // Atomic Mass Unit (in MeV).
const double m_e  = 0.000548579909 *u;
const double m_n  = 1.00866491588  *u;
const double m_p  = 1.007276466879 *u;

// Nuclear masses (in MeV)
const Double_t m_Li   = 7.0160034366*u-3.*m_e;        // 7Li
const double   m_C    = 12.000000000*u-12*m_e;        // 12C
const double   m_Au   = 196.96656879*u-79*m_e;        // 197Au
const double   m_Be   = 7.016928717*u-4.*m_e ;        // 7Be
const double   m_BeEx = m_Be+0.454 ;                  // 7Be (excited)

const int nChannels = 5;
std::string ChannelName[nChannels] = { "Li+C ->Li+C" ,"Li+Au ->Li+Au", "Li+p ->Be+n","Li+p ->Li+p", "Li+p ->Be^{ex}+n"};

//  Simone 2D fits   E0=28 MeV  (I believe these are bi-gaussian fits to the deltaE vs E plot)
/*
Energy   <SiE>      err      <SidE>       err        rho      Blob

  28    5329.16   6.36838    7433.82    5.62415   -0.726917   hiBe    -> Low energy Be
  28    6911.34   5.32618    6487.78    4.36375   -0.701066   lowBe   -> High energy Be
  28    9128.59   3.79272    3183.8     5.23023   -0.291194   lip2    -> Elastic Li+p
  28    9882.1    0.281556   3017.61    0.767732  -0.499893   liC     -> Elastic Li+C
  26    4642.74   24.2769    7882.29    27.3499   -0.735457   hiBe
  26    6169.54   6.97456    6902.27    8.69486   -0.703987   lowBe
  26    8407.22   2.9261     3369.7     5.17021   -0.382229   lip2
  26    9053.87   0.43706    3206.58    1.05283   -0.494892   liC
  18    5573.14   0.66262    4335.97    1.32331   -0.545126   liC
*/

int iFastdE, iFastE, iSlowE, iSlowdE, iMonitor;

std::string IsTime_TPC,   IsTime_LSci;
std::string IsBe,IsLiBand,IsBeBand;
std::string IsNeutron_TPC,IsNeutron_LSci;
std::string IsMonitor;
std::string IsGoodSiTel;

std::string monit,banana,banana_dep,fastslowE,fastslowdE;
std::string dt_lsci,psd_lsci;
std::string dt_tpc,psd_tpc;
std::string SiE,SidE;
std::string SiE_dep,SidE_dep,SiEnergy,MonitEnergy;
std::string SiEFast,SiESlow,SidEFast,SidESlow;


double SiTel_Range[2],SiMon_Range[2],SiE_Range[2];
const double tmin_c=0., tmax_c=50.;

const int nF_max=20;
TH1F   *hT[nF_max], *hE[nF_max], *hE_C[nF_max], *hMon[nF_max];
TH2F   *hPSD[nF_max],*hBanana[nF_max];
TGraph *tPSD[nF_max],*tBanana[nF_max];
int frun_i[nF_max];
double  ftime_i[nF_max], fRateNorm[nF_max], efRateNorm[nF_max];
std::vector<int> iF_plot;
int nF; // Actual number of runs
bool IsTPC;


// MeV/ADC  conversion -->  [iColl] Coll 0,1 [iSi] SiDE, SiE, Monitor
// MeV = A_Conv + B_Conv *adc
const double A_Conv[2][3]={ { 0., 0.,  4.312038e-01 }, { 0., 0.,  4.312038e-01 } };
const double B_Conv[2][3]={ { 0.805202e-3 , 2.56047e-3, 2.335823e-03 }, { 0.805202e-3*1.01735 , 2.56047e-3* 1.02692, 2.335823e-03 } };


// Selection of Li/Be Band in SiTel
const double pLiBand[2]={ 4.18,-4.37 }; // Warning: if you change (A_Conv,B_Conv) this has to be retuned.
const double pBeBand[2]={ 6.84,-6.55 };



//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

void GetKinematics( double *masses, double BeamEnergy,  std::vector<double> &fE3,  std::vector<double> &fTheta3,  std::vector<double> &fE4,  std::vector<double> &fTheta4, bool SelForward)
{

  fE3.clear(); fTheta3.clear(); fE4.clear(); fTheta4.clear();
  //  ------------
  Double_t m1=masses[0];         // 7Li mass
  Double_t m2=masses[1];         // Proton mass (can also be Au or C mass, depending on the channel being considered)

  double E1=BeamEnergy+masses[0];
  double p1=sqrt(E1*E1-m1*m1);

  double s=m1*m1+m2*m2+2*E1*m2;
  double gamma_cm=(E1+m2)/sqrt(s);
  double beta_cm =sqrt( 1.- 1/gamma_cm/gamma_cm);


  //  ------------
  Double_t m3=masses[2];          // 7Be mass (can also be Li depending on the channel being considered)
  Double_t m4=masses[3];          // Neutron mass (can also be Au or C mass, depending on the channel being considered)

  double E3_p=(m3*m3-m4*m4+s)/2./sqrt(s);
  double p3_p=sqrt(E3_p*E3_p-m3*m3);

  double E4_p=sqrt(s)-E3_p;
  double p4_p=sqrt(E4_p*E4_p-m4*m4);

  const double ELab=E1+masses[1];
  const double pLab=p1;

  const int nSteps=1000;
  for (int i=0;i<nSteps;i++)
    {
      double theta_cm=double(i)/nSteps*TMath::Pi();

      double E3   = gamma_cm * (  E3_p               +  beta_cm*p3_p*cos(theta_cm));
      double p3_x = gamma_cm * (  p3_p*cos(theta_cm) +  beta_cm*E3_p              );
      double p3_y = p3_p*sin(theta_cm);
      double p3_z = 0.;
      double theta3_lab = atan2( p3_y,p3_x ); // Berillium angle

      bool Forward=(theta_cm > -TMath::Pi()/2. && theta_cm < TMath::Pi()/2. ? true:false );
      if (  Forward && !SelForward ) continue;
      if ( !Forward &&  SelForward ) continue;

      fTheta3.push_back(theta3_lab*TMath::RadToDeg());
      fE3.push_back(E3-m3);

      double E4   = gamma_cm * (  E4_p                            +  beta_cm  *  p4_p  *  cos(theta_cm+TMath::Pi() )  );
      double p4_x = gamma_cm * (  p4_p*cos(theta_cm+TMath::Pi() ) +  beta_cm  *  E4_p );
      double p4_y = p4_p*sin(theta_cm+TMath::Pi());
      double p4_z = 0.;
      double theta4_lab = fabs(atan2( p4_y,p4_x )); // Neutron angle

      fTheta4.push_back(theta4_lab*TMath::RadToDeg());
      fE4.push_back(E4-m4);
    }

}

void GetThetaNeutron(  double beamEne,  int channel, bool SelForward, TGraph *tg, double thetaMax=25. )
{
  double masses[5][4]={
    { m_Li, m_C,  m_Li,   m_C  },
    { m_Li, m_Au, m_Li,   m_Au },
    { m_Li, m_p,  m_Be,   m_n  },
    { m_Li, m_p,  m_Li,   m_p  },
    { m_Li, m_p,  m_BeEx, m_n  }};


  std::vector<double> fE3,fTheta3,fE4,fTheta4;
  GetKinematics(masses[channel],beamEne,fE3,fTheta3,fE4,fTheta4, SelForward);
  for (int i=0;i<int(fTheta4.size());i++)
    {
      if ( fTheta3[i]>thetaMax )continue; // Breaks the interation and continues into the next one.
      tg->SetPoint(i,fTheta3[i],fTheta4[i]);
    }
}

void GetEnergy(  double beamEne,  int channel, bool SelForward, TGraph *tg, double thetaMax=25. )
{
  double masses[5][4]={
    { m_Li, m_C,  m_Li,   m_C  },   // Li + C  -> Li + C
    { m_Li, m_Au, m_Li,   m_Au },   // Li + Au -> Li + Au
    { m_Li, m_p,  m_Be,   m_n  },   // Li + p  -> Be + n
    { m_Li, m_p,  m_Li,   m_p  },   // Li + p  -> Li + p
    { m_Li, m_p,  m_BeEx, m_n  }};  // Li + p  -> Be* + n


  std::vector<double> fE3,fTheta3,fE4,fTheta4;
  GetKinematics(masses[channel],beamEne,fE3,fTheta3,fE4,fTheta4, SelForward);
  for (int i=0;i<int(fTheta4.size());i++)
    {
      if ( fTheta3[i]>thetaMax )continue;
      tg->SetPoint(i,fTheta3[i],fE3[i]);
    }
}


double GetEnergyAtSi(  double beamEne,  int channel, bool SelForward, double thetaSi )
{
  TGraph *tg=new TGraph();
  GetEnergy(beamEne,channel,SelForward,tg);
  return tg->Eval(thetaSi);
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

double Get_nA(double Hz, double E0, double rho=292.)
{
  if ( E0 == 28. )     return Hz/(rho/292.)/3.75;
  else if ( E0 == 34.) return Hz/(rho/292.)/2.54;
  else return 0.;
}


double GetSiE_Energy(double peak)
{
  return (peak+125.)/0.4112;
}


double GetSiMon_Energy(double peak)
{
  return (peak+28.77)/0.4227;
}

double IntegrateRate(double sigma_cm, double Hz)
{
  TF1 *f=new TF1("f","[0]*exp(-pow(x-[1],2.)/2./[2]/[2])*TMath::Pi()*2*x",-15.,15.);
  double AreaLSci=TMath::Pi()*4.*4.;
  f->SetParameter(1,0.);
  f->SetParameter(0,Hz/AreaLSci);
  f->SetParameter(2,sigma_cm);
  return f->Integral(0.,15.);

}


//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------


void Summary(TTree *t,  int iF )
{
  int     run = frun_i[iF];
  double time = ftime_i[iF];

  double E0 = 28.;
  int iColl = (run > 670 ? 1:0);

  //-----------------------------------------------
  //   Collimator/Beam energy setting (Calibration runs at 18,26 MeV and the last night energy scan 27,29 MeV are not considered)
  //-----------------------------------------------

  int iSetting=0;  //[0] 28 MeV Collimator 0  [1] 34 MeV Collimator 0  [2] 28 MeV Collimator 1
  if ( run>=663 && run <=669 )  { E0=34.;   iSetting=1; }
  else if ( run > 670 )  iSetting=2;

  // [isetting]
  double SiTel_Range_i[3][2]={ { 18.5, 20.5 }, { -100, -100 }, { 17.5, 19.5 } }; // Low energy Berillium location (LowBe blob).

  SiTel_Range[0]=SiTel_Range_i[iSetting][0];
  SiTel_Range[1]=SiTel_Range_i[iSetting][1];

  double ELi=GetEnergyAtSi(E0,0,true,23.); // Li-C elastic
  SiMon_Range[0]=20.;
  SiMon_Range[1]=26.;
  SiE_Range[0]=10.;  SiE_Range[1]=20;    // SiE_dep  histogram in the Be Band


  //-----------------------------------------------
  //  Channel mapping setting
  //-----------------------------------------------

  // LNS-July
  if ( run>=393 && run<=555 )
  {

  }
  // LNS-September
  //   ch 0-8                 LSci
  //   ch 16-18               Bot TPC
  //   ch 20-24, 26-28,32-46  Top TPC
  //   ch 9                   SidE Slow
  //   ch 11                  SiE Slow
  else if ( run >555 )
  {
    iSlowE=11; iSlowdE=9;
    // Setting 0
    if ( run>555 && run<637 )
    {
      iMonitor=13;
      iFastdE=10;
      iFastE=12;
    }
    // Setting 2
    else if ( (run>=715 &&run<=729) || ( run>=680 && run<=684 ) || run==705 )
    {
      iMonitor=13;
      iFastdE=30;
      iFastE=31;
    }
    // Setting 1
    else
	  {
	    iMonitor=31;
	    iFastdE=10;
	    iFastE=12;
	  }
  }

  //-----------------------------------------------
  //  Variables
  //-----------------------------------------------
  SiEFast = Form("(-ymin[%d]+baseline_mean[%d])", iFastE, iFastE);
  SiESlow = Form("(ymax[%d]-baseline_mean[%d])" , iSlowE, iSlowE);

  SidEFast = Form("(-ymin[%d]+baseline_mean[%d])", iFastdE, iFastdE);
  SidESlow = Form("(ymax[%d]-baseline_mean[%d])" , iSlowdE, iSlowdE);

  SiE  = SiESlow ;
  SidE = SidESlow;
  monit = Form("(baseline_mean[%d]-ymin[%d])", iMonitor, iMonitor);


  SiE_dep  = Form("((%s*%e) + %e)", SiESlow.c_str() , B_Conv[iColl][1] , A_Conv[iColl][1]);
  SidE_dep = Form("((%s*%e) + %e)", SidESlow.c_str(), B_Conv[iColl][0] , A_Conv[iColl][0]);
  SiEnergy = Form("(%s*%e + %s*%e + %e)", SidESlow.c_str(), B_Conv[iColl][0],  SiESlow.c_str(), B_Conv[iColl][1], A_Conv[iColl][0] + A_Conv[iColl][1] );

  IsLiBand = Form("%s < (%e + %e*log10(%s/10.))*1.15 && %s > (%e + %e*log10(%s/10.))*0.85 ", SidE_dep.c_str(), pLiBand[0], pLiBand[1], SiE_dep.c_str(),
		SidE_dep.c_str(), pLiBand[0], pLiBand[1], SiE_dep.c_str() );
  IsBeBand = Form("%s < (%e + %e*log10(%s/10.))*1.15 && %s > (%e + %e*log10(%s/10.))*0.85 ", SidE_dep.c_str(), pBeBand[0], pBeBand[1], SiE_dep.c_str(),
		SidE_dep.c_str(), pBeBand[0], pBeBand[1], SiE_dep.c_str() );

  cout << "********************************************" << endl;
  cout << "Li band: " << IsLiBand << endl << endl;
  cout << "Be band: " << IsBeBand << endl << endl;

  MonitEnergy = Form("((%s*%e) + %e)", monit.c_str(), B_Conv[iColl][2], A_Conv[iColl][2]);

  banana     = Form("%s:%s", SidE.c_str(), SiE.c_str());
  banana_dep = Form("%s:%s", SidE_dep.c_str(), SiE_dep.c_str());


  fastslowE  = Form("%s:%s", SiEFast.c_str(), SiESlow.c_str());
  fastslowdE = Form("%s:%s", SidEFast.c_str(), SidESlow.c_str());

  dt_lsci  = Form("xmin[%d] - xmin[0]", iFastE);
  psd_lsci = "f90[0]:charge[0]";

  dt_tpc = Form("xmin[%d] - clusters.min_x[0]", iFastE);
  psd_tpc = "clusters.f90[0]:clusters.charge[0]";


  //-----------------------------------------------
  //  Cuts (To be done, implement cuts to avoid pile ups in Silicon/LSci detectors)
  //-----------------------------------------------

  IsNeutron_LSci = "charge[0] > 20000.&& f90[0] > 0.15 && f90[0] < 0.5";
  IsNeutron_TPC  = "number_of_clusters ==1 && clusters.charge[0] > 0. && clusters.f90[0] < 1. && clusters.f90[0] > 0.45";

  IsTime_LSci = Form("charge[0] > 20000. && xmin[%d] - xmin[0] > %5.2lf && xmin[%d] - xmin[0]<%5.2lf", iFastE, tmin_c, iFastE, tmax_c);
  IsTime_TPC  = Form("number_of_clusters == 1 && clusters.charge[0] > 0. && clusters.f90[0] < 1. && xmin[%d] - clusters.min_x[0] > %5.2lf && xmin[%d]-clusters.min_x[0] < %5.2lf", iFastE, tmin_c, iFastE, tmax_c);

  IsBe = Form("%s && %s > %6.2lf &&  %s < %6.2lf ", IsBeBand.c_str(), SiEnergy.c_str(), SiTel_Range[0], SiEnergy.c_str(),  SiTel_Range[1] );
  IsMonitor = Form("%s>%6.2lf &&  %s < %6.2lf ", MonitEnergy.c_str(), SiMon_Range[0], MonitEnergy.c_str(),  SiMon_Range[1] );


  const char *hname_monit = Form("hMonit%d", run);
  TH1F *hMonit = new TH1F(hname_monit, hname_monit, 200, SiMon_Range[0], SiMon_Range[1]);
  hMonit -> Sumw2();
  t -> Project(hname_monit, MonitEnergy.c_str());

  TH1F *h0 = hMonit;
  double Epeak_Monit = h0 -> GetBinCenter(h0 -> GetMaximumBin()), dx = 0.5, Nmonit_peak = 0.;

  for (int ibin = 0; ibin < h0 -> GetNbinsX(); ibin++)
    {
      double x = h0 -> GetBinLowEdge(ibin);
      if ( x < Epeak_Monit - dx || x > Epeak_Monit + dx ) continue;
      double y = h0 -> GetBinContent(ibin);
      Nmonit_peak += y;
    }

  //-----------------------------------------------
  //  Print statistics
  //-----------------------------------------------

  double rate[5], erate[5];
  int N[5];
  N[0] = t -> Draw("event_number", "", "goff");
  N[1] = t -> Draw("event_number", Form("%s", IsBe.c_str()), "goff");

  if ( !IsTPC )
    {
      N[2] = t -> Draw("event_number", Form("%s && %s", IsNeutron_LSci.c_str(), IsBe.c_str()), "goff");
      N[3] = t -> Draw("event_number", Form("%s && %s && %s", IsTime_LSci.c_str(), IsBe.c_str(), IsNeutron_LSci.c_str()),"goff");
    }
  else
    {
      N[2] = t -> Draw("event_number", Form("%s && %s", IsNeutron_TPC.c_str(), IsBe.c_str()), "goff");
      N[3] = t -> Draw("event_number", Form("%s && %s && %s", IsTime_TPC.c_str(), IsBe.c_str(), IsNeutron_TPC.c_str()), "goff");
    }

  N[4] = int (Nmonit_peak);
  int Nmonit_range = t -> Draw("event_number", Form("%s", IsMonitor.c_str()), "goff");

  for (int ii=0;ii<5;ii++)
    {
      rate[ii] = double(N[ii])/time;
      double eN = (N[ii] == 0? 1.:sqrt( double( N[ii] ) ));
      erate[ii] = eN/time;
    }

  double nA      = Get_nA(rate[4],E0);
  fRateNorm[iF]  = rate[3]/rate[4];
  efRateNorm[iF] = erate[3]/rate[4];

  printf("Run = %d  N = %d/%5.0lf [s] (%5.2lf Hz)   |   %5d (Be)  %5d (Neutron)  %5d (dT)   |   Monitor (%5d %5d)  ELi = %5.2lf  Epeak = %5.2lf\n",
         run, N[0], time, rate[0], N[1], N[2], N[3], N[4], Nmonit_range, ELi, Epeak_Monit);
  printf("                           Rate  %5.3lf (C)  %5.3lf (M) [Hz]\n", rate[3], rate[4] );
  printf("                           RateNorm  %5.3lf  [Hz/Hz_mon]\n",fRateNorm[iF] );
  printf("                           Coin/Monitor  %5.3lf   |   Coin/Be %5.3lf \n", double(N[3])/double(N[4]), double(N[3])/double(N[1]) );
  printf("********************************************\n");



}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

void SetHistos(TTree *t,  int iF  )
{

  //----
  const char *hname=Form("hT%d",iF);
  hT[iF]=new TH1F(hname,hname,30,0.,50.);
  hT[iF]->Sumw2();
  if ( IsTPC )
    t->Project(hname,dt_tpc.c_str(),Form("%s <100. && %s && %s",monit.c_str(),IsNeutron_TPC.c_str(),IsBe.c_str()) );
  else
    t->Project(hname,dt_lsci.c_str(),Form("%s <100. && %s && %s",monit.c_str(),IsNeutron_LSci.c_str(),IsBe.c_str()) );
  hT[iF]->SetDirectory(0);


  const char *hnameE=Form("hE%d",iF);
  hE[iF]=new TH1F(hnameE,hnameE,50,SiE_Range[0],SiE_Range[1]);
  hE[iF]->Sumw2();
  t->Project(hnameE,SiE_dep.c_str(),Form(" %s <100. && %s",monit.c_str(),IsBeBand.c_str()));
  hE[iF]->SetDirectory(0);

  const char *hnameE_C=Form("hE%d_comb",iF);
  hE_C[iF]=new TH1F(hnameE_C,hnameE_C,100,15.,25.);
  hE_C[iF]->Sumw2();
  t->Project(hnameE_C,SiEnergy.c_str(),Form(" %s <100. && %s",monit.c_str(),IsBeBand.c_str()) );
  hE_C[iF]->SetDirectory(0);


  const char *hnameMon=Form("hMon%d",iF);
  hMon[iF]=new TH1F(hnameMon,hnameMon,200,SiMon_Range[0],SiMon_Range[1]);
  hMon[iF]->Sumw2();
  t->Project(hnameMon,MonitEnergy.c_str());
  hMon[iF]->SetDirectory(0);

  //----
  const char *hnamePSD=Form("hPSD%d",iF);
  hPSD[iF]=new TH2F(hnamePSD,hnamePSD,1000,0.,(IsTPC?2000:500.e3),1000,0.,1.);
  hPSD[iF]->Sumw2();
  int nn;
  if ( IsTPC )
    nn=t->Project(hnamePSD,psd_tpc.c_str(),"number_of_clusters==1&&clusters.charge[0]>0.&&clusters.charge[0]<2000.&&clusters.f90[0]<1.&&clusters.f90[0]>0.");
  else
    nn=t->Project(hnamePSD,psd_lsci.c_str(),"f90[0]>0&&f90[0]<1&&charge[0]>100&&charge[0]<500.e3");
  hPSD[iF]->SetDirectory(0);

  tPSD[iF]=new TGraph(nn,t->GetV2(),t->GetV1());

  //----
  const char *hnameBanana=Form("hBanana%d",iF);
  hBanana[iF]=new TH2F(hnameBanana,hnameBanana,300, 5., 28, 300, 1.5, 7.5 );
  hBanana[iF]->Sumw2();
  nn=t->Project(hnameBanana,banana_dep.c_str(),Form("%s>5.&&%s<28.&&%s>1.5&&%s<7.5",SiE_dep.c_str(),SiE_dep.c_str(),SidE_dep.c_str(),SidE_dep.c_str()));
  hBanana[iF]->SetDirectory(0);
  tBanana[iF]=new TGraph(nn,t->GetV2(),t->GetV1());

}


//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

void PlotHistos( int mstyle=10 )
{
  //----------------
  // Get references
  //-----------------

  const int nF_ref=2;
  TH1F   *hE_Ref[2][nF_ref];
  TH1F   *hSiE_Ref[2][nF_ref];
  TH2F   *hBananas_Ref[nF_ref];

  TFile *f1 = new TFile("bananas.root", "READ" );
  for (int iF=0;iF<nF_ref;iF++)
    {
      std::string name=Form("h2D%d",iF);
      if ( f1->Get(name.c_str())==NULL )    hBananas_Ref[iF]=new TH2F();
      else                                  hBananas_Ref[iF]= (TH2F*)   f1->Get(name.c_str());
      hBananas_Ref[iF]->SetDirectory(0);

      for (int iBand=0;iBand<2;iBand++)
	{
	  std::string name=Form("h1D_%d%d",iF,iBand);
	  if ( f1->Get(name.c_str())==NULL )   {printf("NULL\n"); hE_Ref[iBand][iF]=new TH1F();}
	  else                                  hE_Ref[iBand][iF]=(TH1F*) f1->Get(name.c_str());
	  hE_Ref[iBand][iF]->SetDirectory(0);
	  hE_Ref[iBand][iF]->SetLineColor(kBlack);
	  hE_Ref[iBand][iF]->SetLineStyle(iF+1);
	  hE_Ref[iBand][iF]->SetLineWidth(2);
	}
      for (int iBand=0;iBand<2;iBand++)
	{
	  std::string name=Form("h1D_SiE_%d%d",iF,iBand);
	  if ( f1->Get(name.c_str())==NULL )    hSiE_Ref[iBand][iF]=new TH1F();
	  else                                  hSiE_Ref[iBand][iF]=(TH1F*) f1->Get(name.c_str());
	  hSiE_Ref[iBand][iF]->SetDirectory(0);
	  hSiE_Ref[iBand][iF]->SetLineColor(kBlack);
	  hSiE_Ref[iBand][iF]->SetLineStyle(iF+1);
	  hSiE_Ref[iBand][iF]->SetLineWidth(2);
	}
    }
  f1->Close();


  //---------------
  //---------------
  int color[8]={ kRed,kBlue,kOrange+2,kGreen+2,kYellow,kMagenta,kBlack,kGreen+2};
  for (int iplot=0;iplot<2;iplot++)
    {
      cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.*ww,(iplot==0?2.:1.)*wh);
      cc[canvas]->Divide(2,(iplot==0?2:1));
      cc[canvas]->Draw();


      for (int ipad=0;ipad<(iplot==0?4:2);ipad++)
	{
	  cc[canvas]->cd(ipad+1);
	  for (int igraph=0;igraph<(int)iF_plot.size();igraph++)
	    {
	      int iF=iF_plot[igraph];
	      int color_i=color[igraph>7?7:igraph];

	      if ( iplot==0 )
		{
		  TH1F *h0;
		  std::string xaxis;
		  if ( ipad==0)      { h0=hT[iF];  xaxis=Form("t_{Si_E_Fast}-t_{%s} [samples]",(IsTPC?"TPC":"LSci")); }
		  else if (ipad==1 ) { h0=hE[iF];  xaxis="SiE [MeV]";}
		  else if (ipad==2 ) { h0=hE_C[iF]; xaxis="SiTel [MeV] (E+dE)";}
		  else if (ipad==3 ) { h0=hMon[iF]; xaxis="SiMonit [MeV]";}

		  h0->GetXaxis()->SetTitle(xaxis.c_str());
		  h0->GetYaxis()->SetTitle("A.U.");

		  h0->SetLineColor(color_i);
		  h0->Scale(1./h0->Integral());
		  h0->SetMaximum((ipad==0?0.8:0.5));
		  h0->DrawCopy(iF==0?"HIST":"HISTSAME");

		  // Plot Ref
		  if ( ipad==2 || ipad==1 )
		    {
		      for (int iref=1;iref<2;iref++)
			{
			  if (ipad==1 ) hSiE_Ref[1][iref]->DrawCopy("HISTSAME");
			  else          hE_Ref[1][iref]->DrawCopy("HISTSAME");
			}
		    }

		  if ( ipad>0 )
		    {
		      double peak=h0->GetBinCenter(h0->GetMaximumBin());

		      gPad->SetLogy();
		      h0->Fit("gaus","Q0","",peak-1.,peak+1.);
		      h0->GetFunction("gaus")->SetLineWidth(1);
		      h0->GetFunction("gaus")->SetLineColor(color_i);
		      //h0->GetFunction("gaus")->DrawCopy("SAME");
		      l.SetTextColor(color_i);
		      //l.DrawLatex(0.2,0.85-double(igraph)*0.05,Form("mean=%5.2lf rms=%5.2lf",h0->GetFunction("gaus")->GetParameter(1),h0->GetFunction("gaus")->GetParameter(2)));

		    }

		}
	      else
		{
		  TGraph *t0=(ipad==0?tPSD[iF]:tBanana[iF]);

		  t0->GetYaxis()->SetTitle(ipad==0?(IsTPC?"PSD TPC":"PSD_LSci"):"Si_dE [MeV]");
		  t0->GetXaxis()->SetTitle(ipad==0?(IsTPC?"Charge TPC":"Charge_LSci"):"Si_E [MeV]");
		  t0->SetMarkerStyle(mstyle);
		  if ( mstyle==20) t0->SetMarkerSize(0.5);

		  if ( ipad==0 ) t0->SetMinimum(0.);

		  t0->SetMarkerColor(color_i);
		  t0->DrawClone((iF==0?"AP":"P"));

		  if ( ipad==1 )
		    {
		      for (int iBand=0;iBand<2;iBand++)
			{
			  TF1 *f=new TF1("f","([0]+[1]*(log10(x/10.)))*[2]",0.,30.);
			  for (int ii=0;ii<3;ii++)
			    {
			      double par[3]={ (iBand==0?pLiBand[0]:pBeBand[0]), (iBand==0?pLiBand[1]:pBeBand[1]), (ii==0?1.:(ii==1?0.9:1.1)) };
			      f->SetParameters(par);
			      f->SetLineStyle(ii==0?1:2);
			      f->DrawCopy("SAME");
			    }
			}
		    }
		}
	      if ( ipad==0 && iplot==0 )
		{
		  l.SetTextColor(color_i);
		  l.DrawLatex(0.2,0.85-0.05*double(igraph),Form("run=%d",frun_i[iF]));
		}
	    }
	}// ipad
      canvas++;
    }//iplot
}



void Analyze()
{

  for (int iF=0;iF<nF;iF++)
    {
      std::string FileName=Form("reco/run_%d.root",frun_i[iF]);

      TTree *t;
      TFile *f=new TFile(FileName.c_str());

      f->GetObject("reco",t);

      Summary(t,iF);
      SetHistos(t,iF);

      f->Close();
    }
  PlotHistos();
}
