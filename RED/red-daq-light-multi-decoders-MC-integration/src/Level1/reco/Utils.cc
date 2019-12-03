#include "Root_Plot.cc"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "red-daq/EvRec0.hh"

#include <TFitResult.h>

#define UNDEF -100
static const double FADCWidth=2.;

//-----------------------------------------------------------
//  Channel names
//-----------------------------------------------------------
const bool PlotXY=true;

const std::string NameChannel[28]={ "F2", "F3", "F4", "F5",
			      "A1", "A2", "A3", "A4", "A5",
			      "B1", "B2", "B3", "B4", "B5",
			      "C1", "C2", "C3", "C4", "C5",
			      "D1", "D2", "D3", "D4", "D5",
				"E2", "E3", "E4", "E5" };

const double Top_XY[6][4] = {{0, 1, 2, 3},
			   {5, 6, 4, 8},
			   {10, 11, 7, 9},
			   {2+12, 3+12, 1+12, 0+12},
			   {7+12, 8+12, 6+12, 5+12},
			   {9+12, 10+12, 4+12, 11+12}};

const int XYBorder[16]={ 0, 1, 2, 3, 8, 9, 12, 17, 23, 16, 22, 21, 19,14, 10, 5 }; // s2 max chan top of border
const int XYCentral[8]={ 6,4, 11,7,15,13,20,18}; // s2 max chan top of central


const double s1min=100.,s1max=1000;                       // Histogramming and Extraction region
const double s1min_tba=350.,s1max_tba=550.,s2min_tba=1000.;              // TBA studies
const double tdrift_min_tba=23, tdrift_max_tba=33.;       

const double s2max=500000.;
const double s2minFits=5000.;
const double s2maxEcho=200.;

const bool IsAm=false;
// kr --> Kr selection (S1 Z corrections and S2 eLifetime)  --> Must contain a mono peak to work, otherwise no eLifetime
// s2 --> S2XY corrections, S2 resolution and Extaction region studies (Also test eLifetime corrections in other S1 range)
const double s1min_kr=(!IsAm?350.:500.),s1max_kr=(!IsAm?550:800.);    
const double s1min_S2=(!IsAm?550.:500.),s1max_S2=(!IsAm?1000.:800.);  

const double tdriftMax=100,tdriftEcho[2]={30,80.};

// Variables saved in info file
double tdrift0=UNDEF,tdrift1=UNDEF; 
double eLifetime=UNDEF;
double parS2S1[2][2]={{UNDEF,UNDEF},{UNDEF,UNDEF}};  // [S1/S1p][Intercept/Slope]


//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------


//-----------------
// tdrift/S1 distribution
//-----------------
TH1F *hdrift[3][2];     //[all/without echo/with echo][inner/outer]  --> s1min_kr<s1<s1max_kr
TH1F *hS1[2],*hS2[2];

//-------------------
//  S1/S2 vs tdrift   s1min_kr<s1<s1max_kr
//-------------------
TH2F *hS1Z[2],*hS2Z[2]; // [Inner/Outer]

//-----------
//  S2 vs S1  
//-------------
TH2F      *hS1S2[2],*hS1S2p[2];  // [Inner/Outer]
TGraph    *tS1S2[2][2];          // [region][S1/S1p]

TGraphErrors *tS1S2_param[2];                // [S1/S1p]  log10(S2) vs log10(S1) (Inner)
TGraphErrors *tS2_XY[2];                     // Mean/rms of residuals (w.r.t fit) vs max chan  
TGraphErrors *tS1_XY[2];                     // Mean/rms of S1 (s1min_kr<s1<s1max_kr) vs max chan  

//------------
//  F100
//------------
TH2F   *hS1F90;


//-------------
// TBA
//-------------
TH2F *hTBA_Z[3][2]; //[S1/S2/S3][region]   tba vs tdrift
TH1F *hTBA[3][2];   //[S1/S2/S3][region]   tba distribution

//-------------
// Echo
//-------------
TH1F *htEcho[2]; // tEcho(start)-tS2 (start)  , tEcho(start)-tS2 (50%)   
TH1F *hEcho;     // charge
TH2F *hEcho1,*hEcho2; // f100 vs charge , charge S2 vs charge S3  


//-------------
// S2 shape fits
//-------------
namespace S2Shape {

  // Parameters->  [0] tau1 (fixed)  [1]  tau2   [2]   p  [3]  T   [4] sigma  [5] norm [6] t0  [7] start cluster (fixed) 
  // Parameters saved by this macro: 
  //               [0] tau2 [1] p  [2] T [3] sigma 
  
  int FitS2Opt;

  std::string nameV[4]={ "#tau_{2} [#mu s]","p","T [#mu s]", "#sigma [#mu s]" };
  const double ymin[4]={ 3.2, 0.    , 0.5, 0.0 }; // ranges
  const double ymax[4]={ 3.8, 0.35  , 2.3, 0.5  }; // ranges
  const double minR=-0.3,maxR=0.3;

  const double s2min=10000.;
  const double xmin=0.,xmax=60.;
  const double tdrift_min=10.;  // To fit Sigma(tdrift) and calculate PDFs

  // Cuts 
  const double eT_max=0.2,ep_max=0.05, esigma_max=0.1, et0_max=0.15;
  const double chi2Max=3.;
  const double p_min=0.02;


  const int nRegions=3;
  double norm[2][4][nRegions],mean[2][4][nRegions],rms1[2][4][nRegions],rms2[2][4][nRegions];      // [driftLow/High][ipar][outer/inner/all]   Fit to asymmetric gaussian of hPar_PDF
  double parSigma[nRegions][3];                            // Fit to the tdrift dependence of hPar_Drift  Sigma Only
  TH1F     *hPar_PDF[2][nRegions][4];   // [high/low tdrift][inner/outer/all][ipar]    Distribution all tdrift
  TProfile *hPar_Drift[nRegions][4];    // [outer/inner/all][ipar]    tdrift dependence   
  TGraph   *tPar_Drift[nRegions][4];    // [outer/inner/all][ipar]    tdrift dependence   
  TGraphErrors *tPar_XY[2][4][2]; //[idrift][ipar][Mean/RMS]       Mean/RMS of residuals (w.r.t mean at all tdrift) vs XY
  TH2F     *ht0_td;   // t0 vs tdrift
  TH1F     *htdrift;  // tdrift distribution using ( tS1-t0 )
  
  struct FitResult {
    int type; //[0] Deconv trace [1] Raw trace
    int EvtNumber;
    double td;
    double s2,s1;
    int xy;
    int region;
    int rep;
    double tau1,tau2,p,T,sigma,norm,t0,start;
    double chi2,ndf;
    int status,covstatus;
    double etau1,etau2,ep,eT,esigma,enorm,et0,estart;
    double tauS1,etauS1;
  };
}
TGraph *tTauS1;

//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------

double FuncGausWithBG(double *v, double *par)
{
  double t=v[0];
  double norm=par[0];
  double mean=par[1];
  double sigma=par[2];
  double a_bg=par[3];
  double b_bg=par[4];

  double val=a_bg+b_bg*t+norm*exp(-pow(t-mean,2.)/2./sigma/sigma);

  return val;
}

double FuncAsymGaus(double *v, double *par)
{
  double t=v[0];
  double norm=par[0];
  double mean=par[1];
  double sigma1=par[2];
  double sigma2=par[3];
  if ( sigma2==UNDEF ) sigma2=sigma1;

  double sigma=(t<mean?sigma1:sigma2);
  double val=norm*exp(-pow(t-mean,2.)/2./sigma/sigma);

  return val;

}

double FuncDrift(double *v, double *par)
{
  double t=v[0];
  double width=par[0];
  double sigma=par[1];
  
  double A[2]={ t/(TMath::Sqrt2()*sigma), (t-width)/(TMath::Sqrt2()*sigma) };
  double val=  (TMath::Erf(A[0])- TMath::Erf(A[1]) )/2.0*width;
  return par[2]*val;
}

double FuncDrift0(double *v, double *par)
{
  double t=v[0];
  double A=par[0];
  double t0=par[1];
  double t1=par[2];
  double B =par[3];

  if ( t<t0 )      return A;
  else if ( t>t1 ) return B;
  else             return A-(t-t0)*(A-B)/(t1-t0);


}

double FuncEps2(double *v, double *par)
{
  double x=v[0];
  double eps2=par[0];
  double n[6];
  for (int i=0;i<6;i++) n[i]=par[1+i];
 
  double val=0.;
  for (int ie=0;ie<6;ie++)
    {
      double npe=eps2*double(ie+1);
      val+=(n[ie]*TMath::Gaus(x,npe,2.*sqrt(npe)));
    }
  return val;

}
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------


void GetQuantiles(TH2F *h0, TGraph *t0, double q)
{

  int nbinsX=h0->GetNbinsX();
  for (int index=0;index<nbinsX;index++)
    {
      TString name=Form("ProjY_%s_%d",h0->GetName(),index);
      TH1D *hY=h0->ProjectionY(name, index, index );
      double x=h0->GetXaxis()->GetBinCenter(index);
      if (hY->Integral()<20 ) continue;

      double xq[1]={q},yq[1];
      hY->GetQuantiles(1,yq,xq); 
      int ipoint=t0->GetN();
      t0->SetPoint(ipoint,x,yq[0]);
    }
}


void GetProfile(TH2F *h0, TGraphErrors *t0)
{

  int nbinsX=h0->GetNbinsX();
  for (int index=0;index<nbinsX;index++)
    {
      TString name=Form("ProjY_%s_%d",h0->GetName(),index);
      TH1D *hY=h0->ProjectionY(name, index, index );
      double x=h0->GetXaxis()->GetBinCenter(index);
      if (hY->Integral()<20 ) continue;
      double y=hY->GetMean();
      double N=double(hY->Integral());
      double ey=hY->GetRMS()/sqrt(N);
      int ipoint=t0->GetN();
      t0->SetPoint(ipoint,x,y);
      t0->SetPointError(ipoint,0.,ey);
    }
}

void GetFittedProfile(TH2F *h0, TGraphErrors *t0, double min, double max, bool FitBG )
{
  int nbinsX=h0->GetNbinsX();
  TF1 *f=new TF1("f",FuncGausWithBG,0.,1000.,5);

  for (int index=0;index<nbinsX;index++)
    {
      TString name=Form("ProjY_%s_%d",h0->GetName(),index);
      TH1D *hY=h0->ProjectionY(name, index, index );
      double x=h0->GetXaxis()->GetBinCenter(index);
      if (hY->Integral()<100 ) continue;

      double mean=(min+max)/2.,rms=mean*0.1;
      double par[5]={ hY->GetMaximum(), mean,rms, 0.,0. };
      f->SetParameters(par);

      f->FixParameter(4,0.);
      if ( !FitBG  ) f->FixParameter(3,0.);

      TFitResultPtr rp=hY->Fit(f,"QNS","",min,max);
      TFitResult *r=rp.Get();
      if ( !(r->Status()==0 && r->CovMatrixStatus()==3 ) ) continue;

      double y=f->GetParameter(1);
      double ey=f->GetParError(1);
      int ipoint=t0->GetN();
      t0->SetPoint(ipoint,x,y);
      t0->SetPointError(ipoint,0.,ey);
    }
}
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------

void GetXY(int maxchan, double &xcm, double &ycm)
{
  int nx=4;
  int ny=6;
  for (int iy = 0; iy < ny; iy++) 
    for (int ix = 0; ix < nx; ix++)   
      {
	double x=((0.5 + double(ix) )*5./double(nx) ); 
	double y=5.-((0.5 + double(iy) )*5./double(ny) );
	int XY=Top_XY[iy][ix];
	if ( XY==maxchan ) { xcm=x-2.5; ycm=y-2.5; return;}
      }
}

void GetMapXY(TGraphErrors *t0, TGraphErrors *tO, TGraphErrors *tI)
{
  for (int i=0;i<t0->GetN();i++)
    {
      double x,y;
      t0->GetPoint(i,x,y);
      double ey=t0->GetErrorY(i);
      bool IsInner=false;
      for (int ixy=0;ixy<8;ixy++)
	if ( XYCentral[ixy]==int(x) ) IsInner=true;
      TGraphErrors *t1=(IsInner?tI:tO);

      double xcm,ycm;
      GetXY(int(x),xcm,ycm);

      double phi=atan2(ycm,xcm)*TMath::RadToDeg();

      int ipoint=t1->GetN();
      t1->SetPoint(ipoint,phi,y);
      t1->SetPointError(ipoint,0.,ey);
    }

}

void GetDS50(int Opt, TGraphErrors *t0)
{
  double shiftS1=442./302.;
  double shiftS2=9438./7847./1.2;

  if ( Opt==2 )
    {
      double logS1[15*6]={  0.99, 1.01, 1.04, 1.07, 1.09, 1.12, 1.15, 1.18, 1.20, 1.23, 1.26, 1.28, 1.31, 1.34, 1.36,
			    1.39, 1.42, 1.45, 1.47, 1.50, 1.53, 1.55, 1.58, 1.61, 1.63, 1.66, 1.69, 1.72, 1.74, 1.77, 
			    1.80, 1.82, 1.85, 1.88, 1.90, 1.93, 1.96, 1.99, 2.01, 2.04, 2.07, 2.09, 2.12, 2.15, 2.17,
			    2.20, 2.23, 2.26, 2.28, 2.31, 2.34, 2.36, 2.39, 2.42, 2.44, 2.47, 2.50, 2.52, 2.55, 2.58,
			    2.61, 2.63, 2.66, 2.69, 2.71, 2.74, 2.77, 2.79, 2.82, 2.85, 2.88, 2.90, 2.93, 2.96, 2.98, 
			    3.01, 3.04, 3.06, 3.09, 3.12, 3.15, 3.17, 3.20, 3.23, 3.25, 3.28, 3.31, 3.33, 3.36, 3.39 };
      
      double logRat[15*6]={ 2.04, 1.98, 1.97, 1.94, 1.93, 1.91, 1.90, 1.87, 1.86, 1.84, 1.83, 1.81, 1.79, 1.78, 1.76,
			    1.74, 1.73, 1.71, 1.69, 1.67, 1.66, 1.64, 1.63, 1.61, 1.60, 1.59, 1.57, 1.56, 1.55, 1.53, 
			    1.52, 1.51, 1.50, 1.49, 1.48, 1.47, 1.46, 1.46, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 
			    1.46, 1.46, 1.47, 1.48, 1.49, 1.50, 1.50, 1.52, 1.53, 1.54, 1.55, 1.57, 1.58, 1.59, 1.61,
			    1.62, 1.64, 1.65, 1.67, 1.68, 1.70, 1.71, 1.72, 1.73, 1.75, 1.77, 1.78, 1.79, 1.80, 1.81, 
			    1.81, 1.82, 1.83, 1.83, 1.84, 1.86, 1.86, 1.86, 1.87, 1.87, 1.88, 1.88, 1.87, 1.88, 1.86 };
      
      for (int i=0;i<15*6;i++)
	t0->SetPoint(i,logS1[i]+log10(shiftS1),logRat[i]-log10(shiftS1)+log10(shiftS2) );
    }
  else
    { 
      double fprompt[90]={  0.375 ,  0.367 ,  0.358 ,  0.351 ,  0.345 ,  0.339 ,  0.334 ,  0.329 ,  0.325 ,  0.321 , 
			    0.318 ,  0.314 ,  0.312 ,  0.309 ,  0.307 ,  0.305 ,  0.303 ,  0.301 ,  0.300 ,  0.299 , 
			    0.298 ,  0.296 ,  0.295 ,  0.294 ,  0.294 ,  0.293 ,  0.292 ,  0.292 ,  0.291 ,  0.290 , 
			    0.290 ,  0.289 ,  0.289 ,  0.288 ,  0.288 ,  0.287 ,  0.287 ,  0.287 ,  0.286 ,  0.286 , 
			    0.286 ,  0.285 ,  0.285 ,  0.285 ,  0.285 ,  0.285 ,  0.285 ,  0.285 ,  0.284 ,  0.284 , 
			    0.284 ,  0.284 ,  0.284 ,  0.284 ,  0.284 ,  0.284 ,  0.283 ,  0.283 ,  0.283 ,  0.283 , 
			    0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 , 
			    0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 , 
			    0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 ,  0.283 };
      
      double fprompt_rms[90]={  0.094 ,  0.089 ,  0.084 ,  0.080 ,  0.076 ,  0.073 ,  0.070 ,  0.068 ,  0.065 ,  0.063 , 
				0.061 ,  0.059 ,  0.058 ,  0.056 ,  0.054 ,  0.053 ,  0.052 ,  0.051 ,  0.050 ,  0.049 , 
				0.048 ,  0.047 ,  0.046 ,  0.045 ,  0.045 ,  0.044 ,  0.043 ,  0.042 ,  0.042 ,  0.041 , 
				0.041 ,  0.040 ,  0.040 ,  0.039 ,  0.039 ,  0.038 ,  0.038 ,  0.037 ,  0.037 ,  0.036 , 
				0.036 ,  0.036 ,  0.035 ,  0.035 ,  0.035 ,  0.034 ,  0.034 ,  0.034 ,  0.033 ,  0.033 , 
				0.033 ,  0.033 ,  0.032 ,  0.032 ,  0.032 ,  0.032 ,  0.031 ,  0.031 ,  0.031 ,  0.031 , 
				0.030 ,  0.030 ,  0.030 ,  0.030 ,  0.030 ,  0.029 ,  0.029 ,  0.029 ,  0.029 ,  0.029 , 
				0.028 ,  0.028 ,  0.028 ,  0.028 ,  0.028 ,  0.027 ,  0.027 ,  0.027 ,  0.027 ,  0.027 , 
				0.027 ,  0.027 ,  0.026 ,  0.026 ,  0.026 ,  0.026 ,  0.026 ,  0.026 ,  0.026 ,  0.026 };
      
      double s1_fprompt[90]={  32.500 ,  37.500 ,  42.500 ,  47.500 ,  52.500 ,  57.500 ,  62.500 ,  67.500 ,  72.500 ,  77.500 , 
			       82.500 ,  87.500 ,  92.500 ,  97.500 ,  102.500 ,  107.500 ,  112.500 ,  117.500 ,  122.500 ,  127.500 , 
			       132.500 ,  137.500 ,  142.500 ,  147.500 ,  152.500 ,  157.500 ,  162.500 ,  167.500 ,  172.500 ,  177.500 , 
			       182.500 ,  187.500 ,  192.500 ,  197.500 ,  202.500 ,  207.500 ,  212.500 ,  217.500 ,  222.500 ,  227.500 , 
			       232.500 ,  237.500 ,  242.500 ,  247.500 ,  252.500 ,  257.500 ,  262.500 ,  267.500 ,  272.500 ,  277.500 , 
			       282.500 ,  287.500 ,  292.500 ,  297.500 ,  302.500 ,  307.500 ,  312.500 ,  317.500 ,  322.500 ,  327.500 , 
			       332.500 ,  337.500 ,  342.500 ,  347.500 ,  352.500 ,  357.500 ,  362.500 ,  367.500 ,  372.500 ,  377.500 , 
			       382.500 ,  387.500 ,  392.500 ,  397.500 ,  402.500 ,  407.500 ,  412.500 ,  417.500 ,  422.500 ,  427.500 , 
			       432.500 ,  437.500 ,  442.500 ,  442.500 ,  442.500 ,  442.500 ,  442.500 ,  442.500 ,  442.500 ,  442.500 };
      for (int i=0;i<90;i++)
	{
	  double x=s1_fprompt[i]*shiftS1;
	  double y=( Opt==0?fprompt[i]:fprompt_rms[i]);
	  t0->SetPoint(i,x,y);
	}
    }

}

//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
void ResetHisto()
{
  
  for (int type=0;type<3;type++)
    for (int region=0;region<2;region++)
      if ( hdrift[type][region]!=NULL ) hdrift[type][region]->Delete();

  for (int region=0;region<2;region++)
    {
      if ( hS1[region] !=NULL ) hS1[region]->Delete();
      if ( hS2[region] !=NULL ) hS2[region]->Delete();
      if ( hS1S2[region]!=NULL ) hS1S2[region]->Delete();
      if ( hS1Z[region] !=NULL ) hS1Z[region]->Delete();
      if ( hS2Z[region]!=NULL )  hS2Z[region]->Delete();
    }

  for (int region=0;region<2;region++)
    for (int type=0;type<2;type++)
      if ( tS1S2[region][type]!=NULL ) tS1S2[region][type]->Delete();
  
  if (hS1F90 !=NULL ) hS1F90->Delete();


  for (int region=0;region<2;region++)
    for (int signal=0;signal<3;signal++)
      {
	if ( hTBA_Z[signal][region]!=NULL )  hTBA_Z[signal][region]->Delete();
	if ( hTBA[signal][region]!=NULL   )  hTBA[signal][region]->Delete();  
      }
  
  if ( htEcho[0] !=NULL ) htEcho[0]->Delete();
  if ( htEcho[1] !=NULL ) htEcho[1]->Delete();
  if ( hEcho !=NULL ) hEcho->Delete();
  if ( hEcho1 !=NULL ) hEcho1->Delete();
  if ( hEcho2 !=NULL ) hEcho2->Delete();  
 


  if ( tS1S2_param[0]!=NULL ) tS1S2_param[0]->Delete();
  if ( tS1S2_param[1]!=NULL ) tS1S2_param[1]->Delete();
  
  for (int i=0;i<2;i++) 
    {
      if ( tS2_XY[i]!=NULL ) tS2_XY[i]->Delete();
      if ( tS1_XY[i]!=NULL ) tS1_XY[i]->Delete();
    }

  for (int ipar=0;ipar<4;ipar++)
    for (int region=0;region<2;region++)
      {
	for (int idrift=0;idrift<2;idrift++)
	  if ( S2Shape::hPar_PDF[idrift][region][ipar]!=NULL) S2Shape::hPar_PDF[idrift][region][ipar]->Delete();
	if ( S2Shape::hPar_Drift[region][ipar]!=NULL ) S2Shape::hPar_Drift[region][ipar]->Delete();
	if ( S2Shape::tPar_Drift[region][ipar]!=NULL ) S2Shape::tPar_Drift[region][ipar]->Delete();
      }
  if ( S2Shape::ht0_td!=NULL ) S2Shape::ht0_td->Delete();
  if ( S2Shape::htdrift!=NULL) S2Shape::htdrift->Delete();

  for (int idrift=0;idrift<2;idrift++)
    for (int ipar=0;ipar<3;ipar++)
      for (int ii=0;ii<2;ii++)
	if ( S2Shape::tPar_XY[idrift][ipar][ii]!= NULL) S2Shape::tPar_XY[idrift][ipar][ii]->Delete();

}


void SetHisto(int run)
{
  ResetHisto();

  for (int region=0;region<2;region++)
    for (int type=0;type<3;type++)
      {
	hdrift[type][region]=new TH1F(Form("hdrift_%d_%d_%d",run,type,region),"",100,0,tdriftMax);
	hdrift[type][region]->Sumw2();
      }

  for (int region=0;region<2;region++)
    {
      hS1[region] =new TH1F(Form("hS1_%d_%d",run,region ),"",50, 0.   ,s1max    );
      hS2[region] =new TH1F(Form("hS2_%d_%d",run,region ),"",500, 300.,35000.    );
      hS1Z[region] =new TH2F(Form("hS1Z_%d_%d",run,region) ,"",30.,0.,tdriftMax, 75 ,  0., s1max );
      hS2Z[region]=new TH2F(Form("hS2Z_%d_%d",run,region) ,"",30.,0.,tdriftMax ,100.,  1., 100.  );

      hS1S2[region]=new TH2F(Form("hS1S2_%d_%d",run,region),""  ,100,1.5    ,3.    ,200,-1.,2.);
      hS1S2p[region]=new TH2F(Form("hS1S2p_%d_%d",run,region),"",100,1.5    ,3.    ,200,-1.,2.);
      hS1[region]->Sumw2();
      hS2[region]->Sumw2();
      hS1S2[region]->Sumw2();
      hS1S2p[region]->Sumw2();
    }

  hS1F90=new TH2F(Form("hS1F90_%d",run),"",50,0.,s1max,100,0.,1.);
  hS1F90->Sumw2();

  for (int region=0;region<2;region++)  
    for (int signal=0;signal<3;signal++)
      {
	hTBA_Z[signal][region]=new TH2F(Form("hTBA_Z_%d_%d_%d",run,signal,region),"",50, 0., tdriftMax, 100, -1., 1. );
	hTBA_Z[signal][region]->Sumw2();
      
	hTBA[signal][region]=new TH1F(Form("hTBA_%d_%d_%d",run,signal,region),"", 100, -1., 1.);
	hTBA[signal][region]->Sumw2();
      
      }


  for (int i=0;i<2;i++) 
    {
      htEcho[i]=new TH1F(Form("htEcho_%d_%d",run,i),"",50,tdriftEcho[0],tdriftEcho[1]);         // dt
      htEcho[i]->Sumw2();
    }
  hEcho =new TH1F(Form("hEcho_%d",run),"",100,0.,s2maxEcho);          // chg
  hEcho1=new TH2F(Form("hEcho1_%d",run),"",50,0.,s2maxEcho,25,0.,1.); // f100 vs chg
  hEcho2=new TH2F(Form("hEcho2_%d",run),"",25,0.,35000.,25,0.,s2maxEcho);// chg vs s2
  hEcho->Sumw2(); hEcho1->Sumw2(); hEcho2->Sumw2();
  

  for (int region=0;region<2;region++)
    for (int type=0;type<2;type++)
      {
	tS1S2[region][type]=new TGraph();
	tS1S2[region][type]->SetName(Form("tS1S2_%d_%d_%d",run,region,type));
      }

  for (int ii=0;ii<2;ii++)
    {
      tS2_XY[ii]=new TGraphErrors();
      tS2_XY[ii]->SetName(Form("tS2_XY_%d",ii));
      tS1_XY[ii]=new TGraphErrors();
      tS1_XY[ii]->SetName(Form("tS1_XY_%d",ii));
    }

  for (int type=0;type<2;type++)
    {
      tS1S2_param[type]=new TGraphErrors();
      tS1S2_param[type]->SetName(Form("tS1S2_param_%d",type) );
    } 

}

//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------

void ReadHisto(int run, int FitS2Opt=1 )
{
 
  ResetHisto();
  S2Shape::FitS2Opt=FitS2Opt;

  //-----
  TFile *f1 = new TFile( Form("run_%d.histo.root",run)  ,  "READ" );  

  //-----
  for (int type=0;type<3;type++)
    for (int region=0;region<2;region++)
      { 
	std::string name=Form("hdrift_%d_%d_%d",run,type,region);  
	hdrift[type][region]=( f1->Get(name.c_str())==NULL ? new TH1F() : (TH1F*) f1->Get(name.c_str()) );
	hdrift[type][region]->SetDirectory(0);
      }

  //-----
  for (int region=0;region<2;region++)
    {   
      std::string name=Form("hS1_%d_%d",run,region);  
      hS1[region]=( f1->Get(name.c_str())==NULL ? new TH1F() : (TH1F*) f1->Get(name.c_str()) );
      hS1[region]->SetDirectory(0);

      name=Form("hS2_%d_%d",run,region);  
      hS2[region]=( f1->Get(name.c_str())==NULL ? new TH1F() : (TH1F*) f1->Get(name.c_str()) );
      hS2[region]->SetDirectory(0);

      name=Form("hS1S2_%d_%d",run,region);  
      hS1S2[region]=( f1->Get(name.c_str())==NULL ? new TH2F() : (TH2F*) f1->Get(name.c_str()) );
      hS1S2[region]->SetDirectory(0);

      name=Form("hS1S2p_%d_%d",run,region);  
      hS1S2p[region]=( f1->Get(name.c_str())==NULL ? new TH2F() : (TH2F*) f1->Get(name.c_str()) );
      hS1S2p[region]->SetDirectory(0);

      name=Form("hS1Z_%d_%d",run,region);  
      hS1Z[region]=( f1->Get(name.c_str())==NULL ? new TH2F() : (TH2F*) f1->Get(name.c_str()) );
      hS1Z[region]->SetDirectory(0);

      name=Form("hS2Z_%d_%d",run,region);  
      hS2Z[region]=( f1->Get(name.c_str())==NULL ? new TH2F() : (TH2F*) f1->Get(name.c_str()) );
      hS2Z[region]->SetDirectory(0);
    }
  
  for (int region=0;region<2;region++)
    for (int type=0;type<2;type++)
      {
	std::string name=Form("tS1S2_%d_%d_%d",run,region,type);
	tS1S2[region][type]= ( f1->Get(name.c_str())==NULL ? new TGraph() : (TGraph*) f1->Get(name.c_str()) );
      }

  {
    std::string name=Form("hS1F90_%d",run);  
    hS1F90=( f1->Get(name.c_str())==NULL ? new TH2F() : (TH2F*) f1->Get(name.c_str()) );
    hS1F90->SetDirectory(0);    
  }
  //----
  for (int region=0;region<2;region++)
    for (int signal=0;signal<3;signal++)
      {
	std::string name=Form("hTBA_Z_%d_%d_%d",run,signal,region);
	hTBA_Z[signal][region]=( f1->Get(name.c_str())==NULL ? new TH2F() : (TH2F*) f1->Get(name.c_str()) );
	hTBA_Z[signal][region]->SetDirectory(0);    
	
	name=Form("hTBA_%d_%d_%d",run,signal,region);
	hTBA[signal][region]=( f1->Get(name.c_str())==NULL ? new TH1F() : (TH1F*) f1->Get(name.c_str()) );
	hTBA[signal][region]->SetDirectory(0);    
		
    }


  //----
  for (int i=0;i<2;i++)
    {
      std::string name=Form("htEcho_%d_%d",run,i);
      htEcho[i]=( f1->Get(name.c_str())==NULL ? new TH1F() : (TH1F*) f1->Get(name.c_str()) );
      htEcho[i]->SetDirectory(0);    
    }
  {
    std::string name=Form("hEcho_%d",run);
    hEcho=( f1->Get(name.c_str())==NULL ? new TH1F() : (TH1F*) f1->Get(name.c_str()) );
    hEcho->SetDirectory(0);    

    name=Form("hEcho1_%d",run);
    hEcho1=( f1->Get(name.c_str())==NULL ? new TH2F() : (TH2F*) f1->Get(name.c_str()) );
    hEcho1->SetDirectory(0);    

    name=Form("hEcho2_%d",run);
    hEcho2=( f1->Get(name.c_str())==NULL ? new TH2F() : (TH2F*) f1->Get(name.c_str()) );
    hEcho2->SetDirectory(0);    
  }


  //-----
  for (int type=0;type<2;type++)
    {
      std::string name=Form("tS1S2_param_%d",type);
      tS1S2_param[type]=( f1->Get(name.c_str())==NULL ? new TGraphErrors() : (TGraphErrors*) f1->Get(name.c_str()) );     
    }
  
  for (int i=0;i<2;i++) 
    {
      std::string name=Form("tS2_XY_%d",i);
      tS2_XY[i]=( f1->Get(name.c_str())==NULL ? new TGraphErrors() : (TGraphErrors*) f1->Get(name.c_str()) );
      name=Form("tS1_XY_%d",i);
      tS1_XY[i]=( f1->Get(name.c_str())==NULL ? new TGraphErrors() : (TGraphErrors*) f1->Get(name.c_str()) );
    }

  f1->Close();


  //-----------------
  // Info file
  //-----------------
  std::ifstream infofile(Form("info.%d",run),std::ofstream::in);
  infofile>> tdrift0 >> tdrift1;
  infofile>> eLifetime ;
  printf("tdrift0=%5.2lf tdrift1=%5.2lf \n",tdrift0,tdrift1);
  printf("eLifetime=%5.2lf \n",eLifetime);
  for (int type=0;type<2;type++)
    {
      infofile>> parS2S1[type][0]>> parS2S1[type][1];
      printf(" S2 vs S1%s par=%5.2lf %5.2lf | Ratio at S1=1000=%5.2lf \n",(type==0?"":"p"),parS2S1[type][0],parS2S1[type][1],pow(10.,parS2S1[type][0]));      
    }
  infofile.close();
  


  //-----------------
  // (tdrift0,tdrift1)
  //-----------------
  tdrift1+=0.5; // Due to t0

  TGraph *t0=tS1S2[1][1];	
  std::string name=t0->GetName();
  TF1 *f=new TF1("f",FuncDrift0,0.,10.,4);
  f->SetParameter(0,1.); f->SetParameter(1,1.); f->SetParameter(2,2.); f->SetParameter(3,0.);
  TProfile *hp=new TProfile(Form("Prof_%s_tdrif0",name.c_str()),"",200,0.,tdrift1,-0.5,2.);
 
  for (int ii=0;ii<t0->GetN();ii++)
    { double x,y; t0->GetPoint(ii,x,y); if (y>-0.5&&y<1.5) hp->Fill(x,y); }
  hp->Fit(f,"QN","",0.,10.);
  tdrift0=f->GetParameter(2);
  printf("tdrift0=%5.2lf tdrift1=%5.2lf \n",tdrift0,tdrift1);

  //-----------------
  // S2 shape fits
  //-----------------

  if ( tTauS1!=NULL ) tTauS1->Delete();
  tTauS1=new TGraph();
  tTauS1->SetName(Form("tTauS1_%d",run));

  //-- Create histos
  S2Shape::ht0_td=new TH2F(Form("ht0_%d",run),"",22, 20.,65.,30.,0.,3.);
  S2Shape::ht0_td->Sumw2();
  S2Shape::htdrift=new TH1F(Form("hdriftFits_%d",run),"",20,35.,65.); 
  S2Shape::htdrift->Sumw2();

  for (int ipar=0;ipar<4;ipar++)
    for (int region=0;region<S2Shape::nRegions;region++)
      {
	for (int idrift=0;idrift<2;idrift++)
	  {
	    S2Shape::hPar_PDF[idrift][region][ipar]=new TH1F(Form("hpdf_%d_%d_%d_%d",idrift,ipar,region,run),"",90,
							     (ipar==3&&idrift==0?S2Shape::minR:S2Shape::ymin[ipar]),(ipar==3&&idrift==0?S2Shape::maxR:S2Shape::ymax[ipar]));
	    S2Shape::hPar_PDF[idrift][region][ipar]->Sumw2();
	  }
	S2Shape::hPar_Drift[region][ipar]=new TProfile(Form("hParprof_%d_%d_%d",ipar,region,run),"",30,S2Shape::xmin,S2Shape::xmax);
	S2Shape::hPar_Drift[region][ipar]->Sumw2();

	S2Shape::tPar_Drift[region][ipar]=new TGraph(); S2Shape::tPar_Drift[region][ipar]->SetName(Form("tPar_%d_%d_%d",ipar,region,run));
      }


  if ( FitS2Opt==0 ) return;


  //--- Load array
  std::ifstream s2file(Form("s2shape.fits.%d.%d",run,abs(FitS2Opt)),std::ifstream::in);
  if ( !s2file.is_open() ) return;
  printf("Reading %s\n",Form("s2shape.fits.%d.%d",run,abs(FitS2Opt)) );

  double cnt_s2[2][2][2]; // [Low/High S2][Low/High tdrift][ok/not ok]
  for (int i=0;i<2;i++) for (int j=0;j<2;j++) for (int k=0;k<2;k++) cnt_s2[i][j][k]=0.;
  double meanChi2[2][2]={ {0.,0.}, { 0.,0.} }; // [Low/High S2][All/Low/High tdrift] 


  std::vector<S2Shape::FitResult> fs2fit;
  fs2fit.clear();
  while (true)
    {
      S2Shape::FitResult s2fit;
      s2file >> s2fit.type >> s2fit.EvtNumber >> s2fit.td >> s2fit.s2>> s2fit.xy >>s2fit.region>> s2fit.rep
	     >> s2fit.tau1 >>s2fit.tau2 >> s2fit.p >> s2fit.T >> s2fit.sigma >> s2fit.norm >> s2fit.t0 >> s2fit.start
	     >> s2fit.chi2 >> s2fit.ndf >> s2fit.status >> s2fit.covstatus 
	     >> s2fit.etau1 >>s2fit.etau2 >> s2fit.ep >> s2fit.eT >> s2fit.esigma >> s2fit.enorm >> s2fit.et0 >> s2fit.estart
	     >> s2fit.s1>> s2fit.tauS1 >>s2fit.etauS1;

      if ( s2file.eof() ) break;	      

      //---
      int rS2=( s2fit.s2>5000 && s2fit.s2<30000. ? (s2fit.s2>5000 && s2fit.s2<10000 ? 0 : 1) : UNDEF );
      int rtd=( s2fit.td> 0.  && s2fit.td<60.    ? (s2fit.td>  0. && s2fit.td<10.   ? 0 : 1) : UNDEF );
      if ( rS2>=0 && rtd>=0 ) cnt_s2[rS2][rtd][0]+=1.;

      //---
      bool ok=( s2fit.covstatus==3 && s2fit.status==0&& s2fit.ndf>0 && s2fit.chi2/s2fit.ndf<S2Shape::chi2Max &&
		s2fit.eT    < S2Shape::eT_max && 
		s2fit.ep    < S2Shape::ep_max &&
		s2fit.esigma< S2Shape::esigma_max &&
		s2fit.et0   < S2Shape::et0_max && s2fit.p>S2Shape::p_min);
      if ( !ok ) continue;

      //---
      if ( rS2>=0 && rtd>=0 ) 
	{ cnt_s2[rS2][rtd][1]+=1.; meanChi2[rS2][rtd]+=(s2fit.chi2/s2fit.ndf);}
      

      if ( s2fit.s2< S2Shape::s2min ) continue; 

      fs2fit.push_back(s2fit);
      
      s2fit.td-=tdrift0;
    }
  s2file.close();

  printf("---------------------------------\n");
  printf(" S2 rec. efficiencies   \n");
  printf("   5000<S2<10000   |  td<10 : %5.2lf (%d/%d) <chi2>=%5.2lf |   td>10 : %5.2lf (%d/%d)  <chi2>=%5.2lf \n",
	 cnt_s2[0][0][1]/cnt_s2[0][0][0],  int(cnt_s2[0][0][1]),int(cnt_s2[0][0][0]),meanChi2[0][0]/cnt_s2[0][0][1],
	 cnt_s2[0][1][1]/cnt_s2[0][1][0],  int(cnt_s2[0][1][1]),int(cnt_s2[0][1][0]),meanChi2[0][1]/cnt_s2[0][1][1]);
  printf("  10000<S2<30000   |  td<10 : %5.2lf (%d/%d) <chi2>=%5.2lf |   td>10 : %5.2lf (%d/%d)  <chi2>=%5.2lf \n",
	 cnt_s2[1][0][1]/cnt_s2[1][0][0],  int(cnt_s2[1][0][1]),int(cnt_s2[1][0][0]),meanChi2[1][0]/cnt_s2[1][0][1],
	 cnt_s2[1][1][1]/cnt_s2[1][1][0],  int(cnt_s2[1][1][1]),int(cnt_s2[1][1][0]),meanChi2[1][1]/cnt_s2[1][1][1]);
  printf("---------------------------------\n");

  for (int iloop=0;iloop<2;iloop++) // 1st Loop (tdrift dependence)   2nd Loop ( PDF )
    for (int ipar=0;ipar<4;ipar++)
      for (int region=0;region<S2Shape::nRegions;region++)
	{
	  
	  TH1F *h0High=S2Shape::hPar_PDF[0][region][ipar];
	  TH1F *h0Low =S2Shape::hPar_PDF[1][region][ipar];
	  TProfile *h1=S2Shape::hPar_Drift[region][ipar];
	  TGraph *t0=S2Shape::tPar_Drift[region][ipar];
	  TF1 *f=new TF1("fDiffusion","sqrt(2.*x*[0]+[1]*[1])/[2]",0.,60.);
	
	  for (int i=0;i<(int) fs2fit.size();i++)
	    {
	      S2Shape::FitResult *x=&(fs2fit.at(i));
	      double val;
	      if ( ipar==0 ) val=x->tau2;
	      else if ( ipar==1 ) val=x->p;
	      else if ( ipar==2 ) val=x->T;
	      else if ( ipar==3 ) val=x->sigma;

	      bool ok=(( x->region==region || region==2 ) && val > S2Shape::ymin[ipar] && val<S2Shape::ymax[ipar] );		       
	      if ( !ok ) continue;
	      	      
	      
	      if ( iloop==1 ) 
		{
		  if ( ipar==3 && x->td>S2Shape::tdrift_min ) { f->SetParameters(S2Shape::parSigma[1]); val-=f->Eval(x->td); } // Use Inner to calculate residuals
		  if ( x->td>S2Shape::tdrift_min ) h0High->Fill(val);
		  else                             h0Low->Fill(val);
		}
	      else 
		{		  
		  h1->Fill(x->td,val);
		  t0->SetPoint(t0->GetN(),x->td,val);	
		  if ( ipar==0 && region==2 )
		    {
		      S2Shape::ht0_td->Fill(x->td,x->t0); 
		      S2Shape::htdrift->Fill(x->t0+x->td); 
		    }	
		}	      
	    }//events
	  
	  if ( iloop==0 && ipar==3 ) // Fit tdrift dependence
	    {
	      f->SetParameter(0,4.12e-4);
	      f->SetParameter(1,0.05);
	      f->FixParameter(2,50./(tdrift1-tdrift0));
	      h1->Fit(f,"QN","",S2Shape::tdrift_min,60.);
	      for (int kk=0;kk<3;kk++) S2Shape::parSigma[region][kk]=f->GetParameter(kk);
	      printf("Sigma Fit (region=%d): DL=%e +-%e Sigma0=%e +-%e  vd=%e \n",region,S2Shape::parSigma[region][0],f->GetParError(0),S2Shape::parSigma[region][1],f->GetParError(1),
		     S2Shape::parSigma[region][2]);
	    }
	  
	}//region/ipar

  //------- Mean of distributions
  for (int idrift=0;idrift<2;idrift++)
    for (int region=0;region<S2Shape::nRegions;region++)
      for (int ipar=0;ipar<4;ipar++)
	{
	  TH1F *h1=S2Shape::hPar_PDF[idrift][region][ipar];
	  h1->Scale(1./h1->Integral());
	  
	  double min=h1->GetMean()-2.*h1->GetRMS();
	  double max=h1->GetMean()+2.*h1->GetRMS();
	  TF1 *fgaus= new TF1("fgaus",FuncAsymGaus,min,max,4);
	  fgaus->SetParameter(0,h1->GetMaximum());
	  fgaus->SetParameter(1,h1->GetMean());
	  fgaus->SetParameter(2,h1->GetRMS());	  
	  fgaus->FixParameter(3,UNDEF);

	  h1->Fit(fgaus,"QN","",min,max);
	  if ( ipar==1 ) fgaus->SetParameter(3,fgaus->GetParameter(2));
	  min=fgaus->GetParameter(1)-2.*fgaus->GetParameter(2);
	  max=fgaus->GetParameter(1)+2.*fgaus->GetParameter(2);
	  h1->Fit(fgaus,"QN","",min,max);

	  S2Shape::norm[idrift][ipar][region]=fgaus->GetParameter(0);
	  S2Shape::mean[idrift][ipar][region]=fgaus->GetParameter(1);
	  S2Shape::rms1[idrift][ipar][region]=fabs(fgaus->GetParameter(2));
	  S2Shape::rms2[idrift][ipar][region]=(fgaus->GetParameter(3)==UNDEF?UNDEF:fabs(fgaus->GetParameter(3)));
	}

  int idrift=0;

  printf("--------------------------\n");
  printf(" Mean and RMS idrift=%d \n",idrift);
  printf("--------------------------\n");
  for (int ipar=0;ipar<4;ipar++)
    {
      double val,ff;
      if (ipar==0 ) { val=0.; ff=1.e3; }
      if (ipar==1 ) { val=0.; ff=1.;}
      if (ipar==2 ) { val=0.; ff=1.e3; }
      if (ipar==3 ) { val=0.; ff=1.e3; }
      printf("%s %s :",S2Shape::nameV[ipar].c_str(),(ipar==1?"":"[ns]"));
      if ( ipar==1 )  printf(" Mean %4.3lf/%4.3lf  RMS %4.3lf/%4.3lf \n", (S2Shape::mean[idrift][ipar][0]-val)*ff,(S2Shape::mean[idrift][ipar][1]-val)*ff,
			     S2Shape::rms1[idrift][ipar][0]*ff,S2Shape::rms1[idrift][ipar][1]*ff);
      else            printf(" Mean %3.0lf/%3.0lf  RMS %3.0lf/%3.0lf \n", (S2Shape::mean[idrift][ipar][0]-val)*ff,(S2Shape::mean[idrift][ipar][1]-val)*ff,
			     S2Shape::rms1[idrift][ipar][0]*ff,S2Shape::rms1[idrift][ipar][1]*ff);
    }
  


  //--- XY dependence
  for (int idrift=0;idrift<2;idrift++)
    for (int ipar=0;ipar<4;ipar++)
      {
	for (int ii=0;ii<2;ii++)
	  S2Shape::tPar_XY[idrift][ipar][ii]=new TGraphErrors();
	
	TGraphErrors *t0=S2Shape::tPar_XY[idrift][ipar][0];
	TGraphErrors *t1=S2Shape::tPar_XY[idrift][ipar][1];
	TF1 *f=new TF1("fDiffusion","sqrt(2.*x*[0]+[1]*[1])/[2]",0.,60.);
	
	for (int xy=0;xy<24;xy++)
	  {
	    std::vector<double> fy;
	    fy.clear();
	    for (int i=0;i<(int) fs2fit.size();i++)
	      {
		S2Shape::FitResult *x=&(fs2fit.at(i));
		double val;
		if ( ipar==0 ) val=x->tau2;
		else if ( ipar==1 ) val=x->p;
		else if ( ipar==2 ) val=x->T;
		else if ( ipar==3 ) val=x->sigma;
		
		bool Drift=( ( x->td>S2Shape::tdrift_min && idrift==0 ) || ( x->td<S2Shape::tdrift_min && idrift==1 )); 
		bool ok=( x->xy==xy  && Drift && (ipar==3 || ( val > S2Shape::ymin[ipar] && val<S2Shape::ymax[ipar] ) ) );			  
		if ( !ok ) continue;	      
		
		if ( ipar==3 && idrift==0 ) { f->SetParameters(S2Shape::parSigma[1]); val-=f->Eval(x->td); } // Use Inner to calculate residuals

		//Constrain data around +-3 sigma of the mean
		double mu=S2Shape::mean[idrift][ipar][x->region];
		double sigma=S2Shape::rms1[idrift][ipar][x->region];
		if ( val<mu-3.*sigma || val>mu+3.*sigma ) continue;

		fy.push_back( val );
	      }
	    double N=double(fy.size());
	    double mean=TMath::Mean((unsigned short) fy.size(),&fy[0]);
	    double rms =TMath::RMS((unsigned short) fy.size(),&fy[0]);
	    int ipoint=t0->GetN();
	    t0->SetPoint(ipoint,double(xy),mean);
	    t0->SetPointError(ipoint,0.,rms/sqrt(N));
	    
	    t1->SetPoint(ipoint,double(xy),rms);
	    t1->SetPointError(ipoint,0.,rms/sqrt(2.*N));
	    
	  }//xy
      }//ipar/idrift

  //--
}

//-------------------------------------------------------------------------------------------
void WriteHisto(int run)
{
  //--------------
  std::ofstream infofile(Form("info.%d",run),std::ofstream::out);
  infofile<< tdrift0 <<"  "<< tdrift1<<std::endl;
  infofile<< eLifetime <<std::endl;
  for (int type=0;type<2;type++)
    infofile<< parS2S1[type][0] << "  "<< parS2S1[type][1] <<std::endl;
  
  infofile.close();

  //--------------
  TFile *f1 = new TFile( Form("run_%d.histo.root",run)  ,  "RECREATE" );  
  
  for (int type=0;type<2;type++)
    tS1S2_param[type]->Write();
 
  for (int i=0;i<2;i++) 
    {
      tS1_XY[i]->Write();
      tS2_XY[i]->Write();
    }	

  for (int region=0;region<2;region++)
    {
      hS1[region]->Write();
      hS2[region]->Write();
      hS1S2[region]->Write();
      hS1S2p[region]->Write();
      hS1Z[region]->Write();
      hS2Z[region]->Write();
      for (int type=0;type<3;type++)
	hdrift[type][region]->Write();
    }

  for (int region=0;region<2;region++)
    for (int type=0;type<2;type++)
      tS1S2[region][type]->Write();
    

  for (int region=0;region<2;region++)
    for (int signal=0;signal<3;signal++)
      {
	hTBA_Z[signal][region]->Write();
	hTBA[signal][region]->Write();
      }


  hS1F90->Write();
  htEcho[0]->Write(); htEcho[1]->Write();
  hEcho->Write();
  hEcho1->Write();
  hEcho2->Write();

 
  f1->Close();



}



//---------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//    S2 Model
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------


namespace S2Model {

  static const double td=35./*us*/, sigma0=0.05/*mm*/,  DL=4.12e-4 /*mm^2/us*/, vd=0.93 /*mm/us*/;
  double Sigma=sqrt( sigma0*sigma0 + 2.*DL *td ) / vd *1.e3;// ns
 
  static const double tau1=11.,tau2=3.43e3;
  static const double eps2=20.; /*PE/e*/

  double fprompt=0.08,Transit=1.e3;/*ns*/

  double yIdealAux(double t, double tau, double T) {
    if(t < 0.0) {
      return 0.0;
    }
    else if(t <= T) {
      return (1.0 - TMath::Exp(-t/tau))/T;
    }
    else {
      return (TMath::Exp(T/tau) - 1.0) * TMath::Exp(-t/tau)/T;
    }
  }
  
  double idealModel(double* x, double* par) {
    double t    = x[0];
    double tau1 = par[0];
    double tau2 = par[1];
    double p    = par[2];
    double T    = par[3];
    return p * yIdealAux(t, tau1, T) + (1.0 - p) * yIdealAux(t, tau2, T);
  }
  
  double y2(double t, double tau, double sigma) {

    double A=t/(TMath::Sqrt2()*sigma);
    
    double B=-t/tau+TMath::Sq(sigma)/(2.0*TMath::Sq(tau));
    double C=(TMath::Sq(sigma)-t*tau)/(TMath::Sqrt2()*sigma*tau);
    double logD=B+log(TMath::Erfc(C));

    return TMath::Erf(A) - exp(logD);
  }
  
  double y1(double t, double tau, double T, double sigma) {
    return (y2(t, tau, sigma) - y2(t-T, tau, sigma))/(2.0*T);
  }
  

  double yfit(double*x, double* par) {
    
    double tau1_i    = par[0];
    double tau2_i  = par[1];
    double p_i     = par[2];
    double T_i     = par[3];
    double sigma_i = par[4];
    double A       = par[5];
    double t0      = par[6];
    
    double t       = x[0];

    return   A*(p_i * y1(t-t0, tau1_i, T_i, sigma_i) + (1.0 - p_i) * y1(t-t0, tau2_i, T_i, sigma_i));
  }
  
  double yfit_sa(double *x,double *par)
  {
    double sa=x[0];
    double t[1]={(sa-par[7])*FADCWidth/1.e3};
    
    double val=yfit(t,par);
    return val;
  }
}



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//         LAr/GAr properties
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

namespace LArProperties {
  
  static const double rhoLAr=1400.; //[kg/m3]
  static const double g=9.81;        //[m/s2]
  static const double dH=0.3;        // [m]  Depth of Gas Pocket
  static const double pLAr_ReD=1000.e2;                // Pa 
  static const double pGAr_ReD=pLAr_ReD+  rhoLAr*g*dH; // Pa 

  static const double GridToAnode=1.; // cm
  

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
    else         return tGArDriftSpeed_data->Eval(Td);
  }

  double GetGArDriftSpeed_kV(double E /*kV/cm*/ ) // [cm/us]
  {
    double Td=GetGArTd(E);
    return GetGArDriftSpeed(Td,false);
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



}


double GetELField(double GridToAnode_i /*cm*/, double h_g /*cm*/, double VAnode/*kV*/ ) // kV/cm
{
  const double h_e=GridToAnode_i-h_g;
  const double ratioEps=1.5;
  double EX=VAnode/(ratioEps*h_g+h_e);

  double EL=VAnode/(h_g+h_e/ratioEps);
  return EL;

}


void FitGasPocketAndGrid(TGraphErrors *transit, TGraphErrors *yield, double Vanode, bool UseS2  )
{
  // Important: Field Penetration from below the grid to the region above the grid has to be small. 
  //           Run 966: it seems that field above the grid is reduced by the First Ring voltage
  //           Run 978: at its naive value, field penetration is down not up. Using Run 978

  bool debug=false;

  //---Set up
  std::vector<double> ft ,fy,fy_e,   fxcm,fycm,   fhg, fga;
  std::vector<bool>   finner;
  if ( debug ) printf(" Setting Up \n-------\n");
  for (int i=0;i<transit->GetN();i++)
    {
      double x1,x2,t_i,y_i,ey_i;
      transit->GetPoint(i,x1,t_i);
      yield->GetPoint(  i,x2,y_i);
      ey_i=yield->GetErrorY(i);
      if ( x1!=x2 ) continue;
      int XY=int(x1);

      bool IsInner=false;
      for (int ixy=0;ixy<8;ixy++)
	if ( XYCentral[ixy]==XY ) IsInner=true;
     
      double xcm,ycm;
      GetXY(XY,xcm,ycm);

      double phi=atan2(ycm,xcm)*TMath::RadToDeg();
      double r  =sqrt(xcm*xcm+ycm*ycm);
      if ( debug ) printf(" XY=%d  Inner=%d  r=%5.2lf phi=%5.2lf  y_i=%5.2lf t_i=%5.2lf \n",XY,IsInner,r,phi,y_i,t_i);
      ft.push_back(t_i);
      fy.push_back(y_i);
      fy_e.push_back(ey_i);
      fxcm.push_back(xcm);
      fycm.push_back(ycm);
      finner.push_back(IsInner);
      fhg.push_back(0.);
      fga.push_back(LArProperties::GridToAnode);
    }

  int npixel=(int) ft.size();
  TGraphErrors *tYield[2][2];
  
  //---Gas Pocket with outer pixels

  if ( debug ) printf("\n\n Outer Pixels (h_g) \n-------\n");
  double mean_yield[2]={0.,0.};
 
  for (int pixel=0;pixel<npixel;pixel++)
    {
      if ( finner[pixel] && UseS2 ) continue;
      double T=ft[pixel];
      double h_g=UNDEF;
      
     
      for (int i=0;i<LArProperties::GridToAnode*300;i++)
	{
	  double h_g_i=double(i)/300.;
	  double EL=GetELField(LArProperties::GridToAnode, h_g_i, Vanode );
	  double Td=LArProperties::GetGArTd(EL);
	  double v=LArProperties::GetGArDriftSpeed_kV(EL);

	  double T_pred=h_g_i/v;
	  if ( T_pred>T ) 
	    { 
	      h_g=h_g_i;
	      double y_pred=LArProperties::GetLuminiscenceYield_kV(EL);
	      
	      if (debug )printf(" pixel=%d  T=%5.2lf T_pred=%5.2lf | v=%5.2lf [cm/us]  EL=%5.2lf [kV/cm] Td=%5.2lf|  h_g=%5.2lf | Yield=%5.2lf \n",pixel,T,T_pred,v,EL,Td,h_g,y_pred);
	      
	      fhg[pixel]=h_g;
	      mean_yield[0]+=fy[pixel];
	      mean_yield[1]+=y_pred;
	      break;
	    }
	}
    }
  //-Normalize predicted and measure yields to outer pixels
  mean_yield[0]/=(UseS2?16.:24.);
  mean_yield[1]/=(UseS2?16.:24.);
  if ( debug ) printf("\nMean Yield Outer=%5.2lf (data) %5.2lf (pred)\n\n", mean_yield[0],mean_yield[1]);
  for (int pixel=0;pixel<npixel;pixel++)
    {
      fy[pixel]/=mean_yield[0];
      fy_e[pixel]/=mean_yield[0];
    }

  //--Calculate grid to anode distance for inner pixels (Assumed that grid is below gas pocket)
  if ( debug ) printf("\n\n Inner Pixels (h_g/grid to anode) \n-------\n");
  if ( UseS2 )
    for (int pixel=0;pixel<npixel;pixel++)
      {
	//if ( !finner[pixel] ) continue;

	double T=ft[pixel];
	double yield=fy[pixel];
	double y_pred,T_pred;
	double min=1.e10,h_gFit=UNDEF,GridToAnodeFit=UNDEF;

	for (int i=0;i<LArProperties::GridToAnode*150;i++)
	  {
	    double GridToAnode_i=LArProperties::GridToAnode*1.5-double(i)/100.;
	    double h_g=UNDEF;
	    
	    for (int j=LArProperties::GridToAnode*300;j>=i;j--)
	      {
		double h_g_i=LArProperties::GridToAnode*1.5-double(j)/300.;
		
		double EL=GetELField(GridToAnode_i, h_g_i, Vanode );
		double Td=LArProperties::GetGArTd(EL);
		double v=LArProperties::GetGArDriftSpeed_kV(EL);
		
		y_pred=LArProperties::GetLuminiscenceYield_kV(EL)/mean_yield[1];	      
		T_pred=h_g_i/v;
		if (  T_pred>T ) 
		  { 
		    h_g=h_g_i;
		    break;
		  }
	      }// Gas pocket varying height
	    if ( fabs(y_pred-yield)<min )
	      {
		min=fabs(y_pred-yield);
		GridToAnodeFit=GridToAnode_i;
		h_gFit=h_g;
	      }
	  }// GridToAnode distance
      
	//----------
	if ( debug ) printf(" pixel=%d hg=%5.2lf ga=%5.2lf \n",pixel,h_gFit,GridToAnodeFit);
	fhg[pixel]=h_gFit;
	fga[pixel]=GridToAnodeFit;
	
      }//Pixel
  
  //---- Print out
  

  if ( debug ) printf("----------\n Summary \n---------\n");
  for (int inner=0;inner<2;inner++)
    {
      tYield[0][inner]=new TGraphErrors();
      tYield[1][inner]=new TGraphErrors();

      printf("\n");
      for (int pixel=0;pixel<npixel;pixel++)  
	{	
	  if ( finner[pixel]!=inner ) continue;
	  
	  double hg=fhg[pixel];
	  double ga=fga[pixel];
	  
	  double EL=GetELField(ga, hg, Vanode );
	  double Td=LArProperties::GetGArTd(EL);
	  double v=LArProperties::GetGArDriftSpeed_kV(EL);
	  
	  double y_pred=LArProperties::GetLuminiscenceYield_kV(EL)/mean_yield[1];	      
	  double t_pred=hg/v;


	  double phi=atan2(fycm[pixel],fxcm[pixel])*TMath::RadToDeg();
	  int ipoint=tYield[0][inner]->GetN();
	  tYield[0][inner]->SetPoint(ipoint,phi,fy[pixel]);
	  tYield[0][inner]->SetPointError(ipoint,0.,fy_e[pixel]);

	  tYield[1][inner]->SetPoint(ipoint,phi,y_pred);

	  if ( debug ) printf(" pixel=%d  (Fit) | GridToAnode=%5.2lf hg=%5.2lf |  yield=%5.2lf (%5.2lf)  T=%5.2lf (%5.2lf) | EL=%5.2lf v=%5.2lf Td=%5.2lf | Nph=%5.2lf (Opt. Eff=%5.2lf) \n",
			      pixel,ga,hg,fy[pixel],y_pred,ft[pixel],t_pred,EL,v,Td,fy[pixel]*mean_yield[1],30./fy[pixel]/mean_yield[1]);
	}
    }


  //----------
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,1.0*wh);
  cc[canvas]->Divide(3,1);
  cc[canvas]->Draw(); 

  for (int ipad=0;ipad<3;ipad++)
    {
      cc[canvas]->cd(ipad+1);

      if ( ipad<2 )
	{
	  int nx=4,ny=6;
	  TH2F *h0=new TH2F(Form("hXY_%d",ipad),"",nx,-2.5,2.5,ny,-2.5,2.5);
	  for (int pixel=0;pixel<npixel;pixel++)  
	    h0->Fill( fxcm[pixel],fycm[pixel], (ipad==0?fhg[pixel]*10.:fga[pixel]*10.));

	  h0->GetXaxis()->SetTitle("X[cm]");
	  h0->GetYaxis()->SetTitle("Y[cm]");
	  h0->DrawCopy("COLZ");
	  l.SetTextColor(kBlack);
	  l.DrawLatex(0.2,0.85,(ipad==0?"Gas Pocket height [mm]":"Grid to Anode distance [mm]"));
	}
      else
	{
	  double x0=0.20,y0=0.75;
	  TLegend *leg=new TLegend(x0,y0,x0+0.2,y0+0.15);  

	  double min=(ipad==2?0.5:0.1),max=(ipad==2?1.2:0.3);
	  TH2F *hFrame=new TH2F(Form("hEps2_XY_%d",ipad),"",1000,-180.,180.,1000,0.8, 1.45 );
	  hFrame->GetXaxis()->SetTitle("#Phi [deg]");
	  hFrame->GetYaxis()->SetTitle("#epsilon_{2}/<#epsilon_{2}>");
	  hFrame->DrawCopy();

	  int color[2]={kRed,kBlue};
	  for (int ipred=0;ipred<2;ipred++)
	    for (int inner=0;inner<2;inner++)
	      {
		TGraphErrors *t0=tYield[ipred][inner];
		t0->SetMarkerColor(color[inner]);
		t0->SetLineColor(color[inner]);
		t0->SetLineStyle(ipred+1);
		t0->SetMarkerStyle(ipred==0?20:24);
		t0->DrawClone("P");
		l.SetTextColor(color[inner]);
		if ( ipred==0 ) 
		  l.DrawLatex(0.65,0.85-double(inner)*0.05,(inner==0?"Outer pixels":"Inner pixels"));	      
		if ( inner==0 )
		  leg->AddEntry(t0,(ipred==0?"Data":"Model"),"P");
	      }
	  leg->DrawClone();
	}
    }

  canvas++;
  

}
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------



//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------


void Plot()
{
  SetStyle();

  TF1 *fdrift = new TF1("fdrif", FuncDrift, 25., 65.,3);
  fdrift->SetNpx(10000);
  double parD[3]={50.,10.,1.};

  double muKr[2];
  //-------------------------------------
  for (int iPlot=0;iPlot<3;iPlot++)
    {
      cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,1.2*wh);
      cc[canvas]->Divide(2,1);
      cc[canvas]->Draw();  
     
      //tdrift distribution
      if ( iPlot==0 )
	{
	  l.SetTextSize(0.025);
	  for (int region=0;region<2;region++)
	    {
	      cc[canvas]->cd(region+1);
	      double max=hdrift[0][region]->GetMaximum()*1.2;
	      for (int type=0;type<3;type++)
		{
		  int color[3]={kBlack,kBlue,kRed};
		  TH1F *h0=hdrift[type][region];
		  h0->GetXaxis()->SetRange(0,h0->GetXaxis()->FindBin(65.));
		  h0->GetXaxis()->SetTitle("t_{drift} [#mu s]");
		  h0->GetYaxis()->SetTitle("Number of Events");
		  h0->SetLineColor(color[type]); h0->SetMarkerColor(color[type]); h0->SetMarkerStyle(20);
		  h0->SetMaximum(max*1.2);
		  h0->DrawCopy(type==0?"E1":"HISTSAME"); 	
		  l.SetTextColor(color[type]);
		  l.DrawLatex(0.2,0.85-double(type)*0.05,(type==0?"All":(type==1?"n_{cluster}=2":"n_{cluster}=3 (echo)")));
		  if ( type == 0 )
		    {
		      fdrift->SetParameters(parD);
		      h0->Fit(fdrift,"QL0N","",25.,tdriftMax);
		      fdrift->DrawCopy("SAME");
		      l.DrawLatex(0.6,0.75,Form("t_{drift}^{max}=%5.2lf ",fdrift->GetParameter(0)));
		    }
		}
	      l.SetTextColor(kBlack);
	      l.DrawLatex(0.7,0.85,(region==0?"Outer":"Inner"));
	      l.DrawLatex(0.6,0.8,Form("%3.0lf<S1<%3.0lf",s1min_kr,s1max_kr));
	    }
	}
      // s1 (s1,f90) distributions
      if ( iPlot==1 )
	{
	  l.SetTextSize(0.035);

	  cc[canvas]->cd(1);
	  double max=TMath::Max(hS1[0]->GetMaximum(),hS1[1]->GetMaximum())*1.2;
	  for (int region=0;region<2;region++)
	    {
	      int color[2]={kBlue,kRed};
	      TH1F *h0=hS1[region];
	      h0->GetXaxis()->SetTitle("S1");
	      h0->GetYaxis()->SetTitle("Number of Events");
	      h0->SetLineColor(color[region]); h0->SetMarkerColor(color[region]); h0->SetMarkerStyle(20);
	      h0->SetMaximum(max*1.2);
	      h0->DrawCopy(region==0?"E1":"E1SAME"); 	
	      l.SetTextColor(color[region]);
	      l.DrawLatex(0.2,0.85-double(region)*0.05,(region==0?"Outer":"Inner"));	      

	      
	      double min=s1min_kr-100.,max=s1max_kr+100.;
	      TF1 *ff=new TF1("ff",FuncGausWithBG,min,max,5);
	      double mu=h0->GetXaxis()->GetBinCenter(h0->GetMaximumBin());
	      double parFit[5]={ h0->GetMaximum(),     mu,    sqrt(mu),  h0->GetBinContent(h0->GetXaxis()->FindBin(s1max_kr+50.)),0. };
	      ff->SetParameters(parFit);
	      TFitResultPtr rp=h0->Fit(ff,"QSN","",min,max);

	      TFitResult *r=rp.Get();
	      if ( (r->Status()==0 && r->CovMatrixStatus()==3 ) ) 
		{
		  ff->SetLineColor(color[region]); ff->DrawCopy("SAME");
		  l.DrawLatex(0.6,0.85-double(region)*0.05,Form("#mu=%5.1lf #sigma=%5.1lf ",ff->GetParameter(1),ff->GetParameter(2)));
		}	   
	      muKr[region]=ff->GetParameter(1);
	    }
	  l.SetTextColor(kBlack);
	  l.DrawLatex(0.2,0.75,"Rep=1");

	  cc[canvas]->cd(2);
	  TH2F *h0=hS1F90;
	  h0->GetXaxis()->SetTitle("S1");
	  h0->GetYaxis()->SetTitle("f100");
	  h0->DrawCopy("COLZ");
	  l.DrawLatex(0.5,0.85,Form("N=%d",int(h0->GetEntries())));
	  l.DrawLatex(0.6,0.80,"Rep=1");
	  gPad->SetLogz();
	}
      //  tdrift from tFits
      if ( iPlot==2 )
	{
	  cc[canvas]->cd(1);

	  TH2F *h0=S2Shape::ht0_td;
	  h0->GetXaxis()->SetTitle("t_{drift} [#mu s]");
	  h0->GetYaxis()->SetTitle("t_{0}[#mu s]");
	  h0->DrawCopy("COLZ");
	  l.SetTextColor(kBlack); l.DrawLatex(0.2,0.85,Form("s2>%5.0lf (Inner/Outer)",S2Shape::s2min));
	  gPad->SetLogz();

	  cc[canvas]->cd(2);

	  TH1F *h1=S2Shape::htdrift;
	  h1->GetXaxis()->SetTitle("t_{drift}+t_{0} [#mu s]");
	  h1->GetYaxis()->SetTitle("Number of events");
	  h1->DrawCopy("E1");
	  fdrift->SetParameters(parD);
	  h1->Fit(fdrift,"QLN0","",35.,tdriftMax);
	  fdrift->DrawCopy("SAME");
	  l.DrawLatex(0.2,0.25,Form("t_{drift}^{max}=%5.2lf ",fdrift->GetParameter(0)));	  

	}
      canvas++;
    }//iPlot



  //-------------------------------------
  //     Echo
  //-------------------------------------
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",1.5*ww,1.5*wh);
  cc[canvas]->Divide(2,2);
  cc[canvas]->Draw();  

  for (int ipad=0;ipad<4;ipad++)
    {
      cc[canvas]->cd(ipad+1);
      if ( ipad<2 )
	{
	  TH1F *h0=(ipad==0?htEcho[0]:hEcho);
	  h0->GetXaxis()->SetTitle(ipad==0?"t_{S3}-t_{S2} [#mu s]":"S3");
	  h0->GetYaxis()->SetTitle("Number of events");
	  h0->SetMarkerStyle(20); h0->SetMarkerColor(ipad==0?kBlue:kBlack); h0->SetLineColor(ipad==0?kBlue:kBlack);
	  h0->DrawCopy(ipad==0?"HIST":"E1");

	  if ( ipad==0 )
	    {
	      htEcho[1]->SetLineColor(kRed);
	      htEcho[1]->DrawCopy("HISTSAME");
	      l.SetTextColor(kBlue); l.DrawLatex(0.2,0.85,"t_{1st} (S2)");
	      l.SetTextColor(kRed); l.DrawLatex(0.2,0.80,"t_{50} (S2)");
	      l.SetTextColor(kBlack);
	    }
	  if ( ipad==1 )
	    {
	      double min=15.,max=s2maxEcho;
	      TF1 *feps2 = new TF1("feps2", FuncEps2, min, max,7);
	      double eps2=h0->GetXaxis()->GetBinCenter(h0->GetMaximumBin());
	      double par[7]={ eps2 , h0->GetMaximum(), 1., 1., 1. ,1. ,1.};	      
	      feps2->SetParameters(par);
	      feps2->SetParLimits(0,eps2-5.,100);
	      h0->Fit(feps2,"QN","",min,max);
	      feps2->DrawCopy("SAME");
	      l.DrawLatex(0.5,0.85,Form("#epsilon_{2}=%5.2lf +-%5.2lf",feps2->GetParameter(0),feps2->GetParError(0)));
	    }
	}
      else
	{
	  TH2F *h0=(ipad==2?hEcho1:hEcho2);
	  h0->GetXaxis()->SetTitle(ipad==2?"S3 [PE]":"S2 [PE]");
	  h0->GetYaxis()->SetTitle(ipad==2?"f100":"S3 [PE]");
	  h0->DrawCopy("COLZ");
	  gPad->SetLogz();
	}

    }
  canvas++;


  //-------------------------------------
  //     S1 Z corrections
  //-------------------------------------
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.0*ww,1.2*wh);
  cc[canvas]->Divide(2,2);
  cc[canvas]->Draw();  
  
  double tmin=2.,tmax=tdrift1-2.;
  for (int ipad=0;ipad<4;ipad++)
    {
      bool isS1=(ipad<2);
      int region=(ipad%2==0?1:0);

      cc[canvas]->cd(ipad+1);

      double min=(isS1?s1min_kr:UNDEF);
      double max=(isS1?s1max_kr:UNDEF);
      
      TH2F *h0=(isS1?hS1Z[region]:hS2Z[region]);
      TGraphErrors *t0=new TGraphErrors();
      if ( isS1 ) GetFittedProfile(h0,t0,min-50.,max+50.,true);
      else        GetProfile(h0,t0);
      if ( !isS1 ) { max=t0->Eval(0.)*1.2; min=t0->Eval(tdrift1)*0.2;}

      h0->GetXaxis()->SetRange(h0->GetXaxis()->FindBin(tmin),h0->GetXaxis()->FindBin(tmax));
      h0->GetYaxis()->SetRange(h0->GetYaxis()->FindBin(min),h0->GetYaxis()->FindBin(max));
      h0->GetXaxis()->SetTitle("t_{drift}");
      h0->GetYaxis()->SetTitle((isS1?"S1":"S2/S1"));
      h0->DrawCopy("COLZ");
      
      t0->DrawClone("P");
    
      if ( isS1 ) 
	{
	  TF1 *f=new TF1("f","[0]+(x-15.)*[1]",15.,tdrift1);
	  f->SetParameter(0,muKr[region]);
	  f->SetParameter(1,0.);
	  if ( IsAm ) t0->Fit(f,"QN","",15.,tdrift1);
	  double mu=f->GetParameter(0);
	  f->DrawCopy("SAME");
	  f->SetLineStyle(2); 
	  f->SetParameter(0,mu*0.98); f->DrawCopy("SAME");
	  f->SetParameter(0,mu*1.02); f->DrawCopy("SAME");
	}
      else
	{
	  TF1 *f=new TF1("f","[0]*exp(-(x-15.)/[1])",15.,tdrift1);
	  f->SetParameter(0,t0->Eval(15.));
	  f->FixParameter(1,eLifetime);
	  t0->Fit(f,"QN","",15.,tdrift1);
	  f->DrawCopy("SAME");
	  l.DrawLatex(0.5,0.85,Form("#tau=%5.1lf #mu s",eLifetime));
	}
      l.DrawLatex(0.2,0.85,"Rep=1 (S2)");
      l.DrawLatex(0.2,0.80,(region==0?"Outer":"Inner"));	      
    }

  canvas++;

  //-------------------------------------
  //     S2 XY corrections
  //-------------------------------------
  TGraphErrors *tDS50R=new TGraphErrors();
  GetDS50(2,tDS50R);

  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.*ww,2.0*wh);
  cc[canvas]->Divide(2,2);
  cc[canvas]->Draw();  

  bool UsePrompt=false;
  for (int ipad=0;ipad<4;ipad++)
    {
      cc[canvas]->cd(ipad+1);

      //--------
      if ( ipad<2 )
	{
	  int type=(UsePrompt?1:0);
	  int region=(ipad==0?1:0);
	  TH2F *h0=(UsePrompt?hS1S2p[region]:hS1S2[region]);
	  h0->GetXaxis()->SetTitle(UsePrompt?"log_{10} S1 (prompt)":"log_{10} S1");
	  h0->GetYaxis()->SetTitle("log_{10} S2 (zcorr)/S1");
	  //h0->GetYaxis()->SetRange(h0->GetYaxis()->FindBin(0.9),h0->GetYaxis()->FindBin(1.8));
	  
	  h0->DrawCopy("COLZ");
	  gPad->SetLogz();
	  
	  l.DrawLatex(0.2,0.85,"Rep=1");
	  l.DrawLatex(0.2,0.80,(region==0?"Outer":"Inner"));	      

	  if ( region==1)
	    {
	      l.SetTextColor(kBlue);
	      l.DrawLatex(0.2,0.25,Form("Ratio at S1=1000 : %3.2lf",pow(10.,parS2S1[type][0])));
	      l.SetTextColor(kBlack);
	    }

	  const int nQ=3;
	  double xq[3]={0.1,0.5,0.9};
	  TGraph *tg[nQ];
	  for (int iq=0;iq<3;iq++) 
	    { 
	      tg[iq]=new TGraph(); 
	      GetQuantiles(h0,tg[iq],xq[iq]);
	      tg[iq]->SetLineWidth(2.);
	      tg[iq]->SetLineStyle(iq==1?1:2);
	      tg[iq]->DrawClone("L");
	    }
	  
	  TF1 *f=new TF1("f","[0]+(x-3.)*[1]",log10(s1min_S2),log10(s1max_S2));    
	  f->SetParameters(parS2S1[type]);
	  f->SetLineColor(kBlue);f->DrawCopy("SAME");

	  tDS50R->SetLineWidth(2);
	  tDS50R->SetLineStyle(1); tDS50R->SetLineColor(kRed);
	  tDS50R->DrawClone("L");
	  
	}
      //--------.
      else
	{
	  double min=(ipad==2?0.5:0.1),max=(ipad==2?1.5:0.5);
	  //double min=(ipad==2?0.7:0.15),max=(ipad==2?1.2:0.3);
	  TH2F *hFrame=new TH2F(Form("hFrameXY_%d",ipad),"",1000,(PlotXY?-180.:-0.5),(PlotXY?180.:25),1000,min,max);
	  hFrame->GetXaxis()->SetTitle(PlotXY?"#Phi [deg]":"s2 max chan");
	  hFrame->GetYaxis()->SetTitle((ipad==2?"Mean  Ratio/Pred":"RMS Ratio/Pred"));
	  hFrame->DrawCopy();
	  
	  TF1 *f=new TF1("f","1.",(PlotXY?-180.:-0.5),(PlotXY?180.:25));
	  f->DrawCopy("SAME");

	  TGraphErrors *t0=tS2_XY[ipad-2];
	  TGraphErrors *t1[2];
	  if (PlotXY ) { t1[0]=new TGraphErrors(); t1[1]=new TGraphErrors(); GetMapXY(t0 , t1[0],t1[1]);}
	  else         { t1[0]=t0; t1[1]=t0; }

	  for (int ir=0;ir<(PlotXY?2:1);ir++)
	    {
	      TGraphErrors *tp=(PlotXY?t1[ir]:t0);
	      if ( PlotXY ) { tp->SetLineColor(ir==0?kBlue:kRed); tp->SetMarkerColor(ir==0?kBlue:kRed); l.SetTextColor(ir==0?kBlue:kRed); l.DrawLatex(0.5,0.85-double(ir)*0.05,(ir==0?"Outer":"Inner")); }
	      tp->DrawClone("P");
	    }
	  if ( !PlotXY )
	    for (int i=0;i<8;i++)
	      {
		int ch=XYCentral[i];
		TBox box(double(ch)-0.5,min,double(ch)+0.5,max);
		box.SetFillStyle(3001);
		box.SetFillColor(kGreen+2);
		box.DrawClone();
	      }
	}//
    }
  canvas++;


  //-------------------------------------
  //     S1 XY corrections
  //-------------------------------------
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.*ww,1.2*wh);
  cc[canvas]->Divide(2,1);
  cc[canvas]->Draw();  
  for (int ipad=0;ipad<2;ipad++)
    {
      cc[canvas]->cd(ipad+1);
      
      double min=(ipad==0?s1min_kr:0.0),max=(ipad==0?s1max_kr:0.2);
      TH2F *hFrame=new TH2F(Form("hFrameS1XY_%d",ipad),"",1000,(PlotXY?-180.:-0.5),(PlotXY?180.:25),1000,min,max);
      hFrame->GetXaxis()->SetTitle(PlotXY?"#Phi [deg]":"s2 max chan");
      hFrame->GetYaxis()->SetTitle((ipad==0?"<S1>":"RMS S1/<S1>"));
      hFrame->DrawCopy();
	  
      TF1 *f=new TF1("f","1.",(PlotXY?-180.:-0.5),(PlotXY?180.:25));
      f->DrawCopy("SAME");

      TGraphErrors *t0=tS1_XY[ipad];
      TGraphErrors *t1[2];
      if (PlotXY ) { t1[0]=new TGraphErrors(); t1[1]=new TGraphErrors(); GetMapXY(t0 , t1[0],t1[1]);}
      else         { t1[0]=t0; t1[1]=t0; }
      
      for (int ir=0;ir<(PlotXY?2:1);ir++)
	{
	  TGraphErrors *tp=(PlotXY?t1[ir]:t0);
	  if ( PlotXY ) { tp->SetLineColor(ir==0?kBlue:kRed); tp->SetMarkerColor(ir==0?kBlue:kRed); l.SetTextColor(ir==0?kBlue:kRed); l.DrawLatex(0.5,0.85-double(ir)*0.05,(ir==0?"Outer":"Inner")); }
	  tp->DrawClone("P");
	}
      if ( !PlotXY )
	for (int i=0;i<8;i++)
	  {
	    int ch=XYCentral[i];
	    TBox box(double(ch)-0.5,min,double(ch)+0.5,max);
	    box.SetFillStyle(3001);
	    box.SetFillColor(kGreen+2);
	    box.DrawClone();
	  }
      l.SetTextColor(kBlack);
      l.DrawLatex(0.6,0.8,Form("%3.0lf<S1<%3.0lf   Rep(S1)",s1min_kr,s1max_kr));      
    }
  canvas++;

  //-------------------------------------
  //    Extraction region
  //-------------------------------------

  cc[canvas]=new TCanvas(Form("c%d",canvas),"",1.5*ww,1.2*wh);
  cc[canvas]->Draw();  
  {
      int ipad=0;
      int type=0;

      TH2F *h0=new TH2F(Form("hlogResidual_%d",ipad),"",1000,(type==1?0.:15.),(type==1?5.:tdrift1),1000,(type==1?-1.:-0.5),(type==1?1.5:0.5) );
      h0->GetXaxis()->SetTitle("t_{drift} #mu s");
      h0->GetYaxis()->SetTitle("log_{10} S2/S1 - Pred");
      h0->DrawCopy();
      int color[2]={kBlue,kRed};
      l.SetTextColor(kBlack);
      if ( type == 1 ) l.DrawLatex(0.2,0.25,Form("Using S1_{prompt}  %3.0lf-%3.0lf PE",s1min,s1max));
      else             l.DrawLatex(0.2,0.25,Form("Using S1  %3.0lf-%3.0lf PE  Rep_{S1}=1",s1min_S2,s1max_S2));
	

      for (int region=0;region<2;region++)
	{
	  std::string name=tS1S2[region][type]->GetName();
	  TProfile *hp=new TProfile(Form("Prof_%s_%d",name.c_str(),ipad),"",(type==1?200:100),0.,tdrift1,-0.5,2.);

	  TGraph *t0=tS1S2[region][type];	
	  if ( t0->GetN()==0 ) continue;
	  t0->SetMarkerSize(0.5);
	  t0->SetMarkerColor(color[region]);
	  t0->DrawClone("P");
	  for (int ii=0;ii<t0->GetN();ii++)
	    { double x,y; t0->GetPoint(ii,x,y); if (y>-0.5&&y<1.5) hp->Fill(x,y); }
	    
	  hp->SetMarkerColor(kBlack);
	  hp->SetLineColor(kBlack); hp->SetMarkerStyle(region==1?20:30);
	  hp->SetMarkerSize(1.3);
	  hp->DrawClone("E1SAME");
	  l.SetTextColor(color[region]);
	  l.DrawLatex(0.65,0.85-double(region)*0.05,(region==0?"Outer":"Inner"));	      
	  if ( ipad==1 )
	    {
	      TF1 *f=new TF1("f",FuncDrift0,0.,10.,4);
	      f->SetParameter(0,1.); f->SetParameter(1,1.); f->SetParameter(2,2.); f->SetParameter(3,0.);
	      hp->Fit(f,"QN","",0.,10.);
	      f->SetLineColor(kBlack);
	      f->SetLineStyle(2-region);
	      f->DrawCopy("SAME");
	      l.SetTextColor(color[region]);
	      l.DrawLatex(0.65,0.85-double(region)*0.05,Form("%s t_{drift0}=%3.1lf #mus",(region==0?"Outer":"Inner"),f->GetParameter(2) ));	      
	    }
	  l.SetTextColor(kBlack);

	}
    }
  canvas++;

  //-------------------------------------
  //     S2
  //-------------------------------------

  cc[canvas]=new TCanvas(Form("c%d",canvas),"",1.5*ww,1.*wh);
  cc[canvas]->Draw();  
  
  {
    double max=TMath::Max(hS2[0]->GetMaximum(),hS2[1]->GetMaximum())*1.2;
    for (int region=0;region<2;region++)
      {
	int color[2]={kBlue,kRed};
	TH1F *h0=hS2[region];
	h0->GetXaxis()->SetTitle("S2 (zcorr)");
	h0->GetYaxis()->SetTitle("Number of Events");
	h0->SetLineColor(color[region]); h0->SetMarkerColor(color[region]); h0->SetMarkerStyle(20);
	h0->SetMaximum(max*1.3);
	h0->SetMinimum(0.5);
	h0->DrawCopy(region==0?"E1":"E1SAME"); 	
	l.SetTextColor(color[region]);
	l.DrawLatex(0.6,0.85-double(region)*0.05,Form("%s mean/rms=%3.0lf %3.2lf",(region==0?"Outer":"Inner"),h0->GetMean(),h0->GetRMS()/h0->GetMean() ) );	      
      }
    l.SetTextColor(kBlack);
    l.DrawLatex(0.6,0.55,"Rep=1");	 
    l.DrawLatex(0.6,0.50,Form("%3.0lf<S1<%3.0lf",s1min_kr,s1max_kr));
    gPad->SetLogy();
  } 
  canvas++;


  //-------------------------------------
  //     TBA
  //-------------------------------------

  //---- TBA distribution
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",1.5*ww,1.*wh);
  cc[canvas]->Draw();  
  
  double x0=0.20,y0=0.65;
  TLegend *leg=new TLegend(x0,y0,x0+0.25,y0+0.25);  
  
  std::string TBA=Form("(S^{top}-S^{bot})/(S^{top}+S^{bot})");
  for (int isS2=0;isS2<2;isS2++)
    for (int region=0;region<2;region++)
      {
	TH1F *h0=hTBA[isS2][region];
	if ( h0->GetEntries()==0 ) continue;
	h0->GetXaxis()->SetTitle(TBA.c_str());
	h0->GetYaxis()->SetTitle("A.U.");
	h0->Scale(1./h0->Integral());
	h0->SetMaximum(0.25);
	h0->SetLineColor(isS2==0?kRed:kBlue); h0->SetMarkerColor(isS2==0?kRed:kBlue);
	h0->SetLineStyle(2-region);
	h0->DrawCopy(isS2==0&&region==0?"HIST":"HISTSAME");
	leg->AddEntry(h0,Form("S%d (%s) %3.2lf",isS2+1,(region==1?"Inner":"Outer"),h0->GetMean()),"L");
      }
  leg->DrawClone();
  l.DrawLatex(0.6,0.85,Form("%3.0lf<s1<%3.0lf     s2>%3.0lf",s1min_tba,s1max_tba,s2min_tba));
  l.DrawLatex(0.55,0.82,Form("t_{drift}=   %3.0lf-%3.0lf         All",tdrift_min_tba,tdrift_max_tba));
  canvas++;
	      
  //-----TBA vs tdrift

  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2*ww,1.*wh);
  cc[canvas]->Divide(2,1);
  cc[canvas]->Draw();  

  for (int ipad=0;ipad<2;ipad++)
    {
      cc[canvas]->cd(ipad+1);
      std::string TBA=Form("(S%d^{top}-S%d^{bot})/(S%d^{top}+S%d^{bot})",ipad+1,ipad+1,ipad+1,ipad+1);      
  
      int region=1;

      TH2F *h0=hTBA_Z[ipad][region];
      h0->GetYaxis()->SetTitle(TBA.c_str());
      h0->GetXaxis()->SetTitle("t_{drift} #mu s");
      h0->GetXaxis()->SetRange(h0->GetXaxis()->FindBin(tdrift0),h0->GetXaxis()->FindBin(tdrift1));
      h0->GetYaxis()->SetRange(h0->GetYaxis()->FindBin(-0.5),h0->GetYaxis()->FindBin(0.5));
      h0->DrawCopy("COLZ");
      gPad->SetLogz();

      if ( ipad==0 )
	{	 
	  TF1 *f=new TF1("f","0.",tdrift0,tdrift1);
	  f->DrawCopy("SAME");
	}

      l.SetTextColor(kBlack);
      if ( ipad==0 ) l.DrawLatex(0.2,0.15,Form("%3.0lf<s1<%3.0lf ",s1min_tba,s1max_tba));
      else           l.DrawLatex(0.2,0.15,Form("s2>%3.0lf",s2min_tba));

      for (int region=0;region<2;region++)
	{
	  int color_i=(region==0?kRed:kBlue);
	  TGraphErrors *t0=new TGraphErrors();
	  GetProfile(hTBA_Z[ipad][region],t0);
	  t0->SetMarkerStyle(20);
	  t0->SetMarkerColor(color_i); t0->SetLineColor(color_i);	  
	  t0->DrawClone("P");	
	  l.SetTextColor(color_i);
	  
	  TF1 *f1=new TF1("f","(x-[0])*[1]",tdrift0,tdrift1);
	  f1->SetParameter(0,(tdrift1-tdrift0)/2.);
	  f1->SetParameter(1,-(tdrift1-tdrift0)/2.*0.2);
	  t0->Fit(f1,"QN","",tdrift0,tdrift1);
	  f1->SetLineColor(color_i);
	  f1->DrawCopy("SAME");

	  double z0=(f1->GetParameter(0)-tdrift0)/(tdrift1-tdrift0),ez0=f1->GetParError(0)/(tdrift1-tdrift0);
	  if ( ipad==0 )
	    l.DrawLatex(0.3,0.3-double(region)*0.05,Form("%s z0=%3.2lf +-%4.3lf [#mus] %3.2lf/%3.2lf",(region==0?"Outer":"Inner"),z0,ez0,f1->Eval(tdrift0),f1->Eval(tdrift1) ) );

	}
    }
  canvas++;

  //-------------------------------------
  //     Mean/RMS f100
  //-------------------------------------
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,1.2*wh);
  cc[canvas]->Divide(3,1);
  cc[canvas]->Draw();  

  cc[canvas]->cd(1);

  TH2F *h0=hS1F90;
  h0->GetXaxis()->SetTitle("S1");
  h0->GetYaxis()->SetTitle("f100");
  h0->DrawCopy("COLZ");
  l.DrawLatex(0.5,0.85,Form("N=%d",int(h0->GetEntries())));
  l.DrawLatex(0.6,0.80,"Rep=1");

  gPad->SetLogz();

    
  TGraphErrors *tRed[2],*tDS50[2];
  for (int i=0;i<2;i++)
    { tRed[i]=new TGraphErrors(); tDS50[i]=new TGraphErrors(); GetDS50(i,tDS50[i]); }
  
  int nbinsX=h0->GetNbinsX();
  for (int index=0;index<nbinsX;index++)
    {
      TString name=Form("ProjY_%d",index);
      TH1D *hY=h0->ProjectionY(name, index, index );
      double x=h0->GetXaxis()->GetBinCenter(index);
      if (hY->Integral()<10 ) continue;
      
      double S1=x;
      if ( S1<100. || S1>1000. ) continue;
      hY->GetXaxis()->SetRange(hY->GetXaxis()->FindBin(0.1),hY->GetXaxis()->FindBin(0.6));
      double N=hY->GetEntries();
      int ipoint=tRed[0]->GetN();
     
      double mean=hY->GetMean(), rms=hY->GetRMS();
      tRed[0]->SetPoint(ipoint,x,mean);
      tRed[0]->SetPointError(ipoint,0.,rms/sqrt(N));

      tRed[1]->SetPoint(ipoint,x,rms);
      tRed[1]->SetPointError(ipoint,0.,rms/sqrt(2.*N));
    }
  
  for (int ipad=1;ipad<3;ipad++)
    {
      cc[canvas]->cd(ipad+1);
      for (int i=0;i<2;i++)
	{
	  int color_i=(i==0?kRed:kBlue);
	  TGraphErrors *t0=(i==0?tRed[ipad-1]:tDS50[ipad-1]);
	  t0->GetYaxis()->SetTitle((ipad==1?"<f100>":"RMS f100"));
	  t0->GetXaxis()->SetTitle("S1");
	  t0->SetLineColor(color_i);  t0->SetMarkerColor(color_i); 
	  t0->DrawClone(i==0?"AP":"L");
	  l.SetTextColor(color_i);
	  l.DrawLatex(0.5,0.85-double(i)*0.05,i==0?"ReD":"DS-50");
	}
    }
  
  canvas++;


}





void PlotS2Shape()
{
  if ( S2Shape::FitS2Opt==0 ) return;

  std::string nameRegion[3]={"Outer","Inner","All"};
  //---------------------------------------------------------------------------
  //  tdrift dependence and PDF
  //---------------------------------------------------------------------------

  for (int iPlot=0;iPlot<2;iPlot++)
    {

      int nPad_x=2;
      int nPad_y=(S2Shape::FitS2Opt==1?2:1);
      double wx=(S2Shape::FitS2Opt==2?2.0:1.5)*ww;
      double wy=(S2Shape::FitS2Opt==2?1.2:1.5)*wh;

      cc[canvas]=new TCanvas(Form("c%d",canvas),"",wx,wy);
      if ( S2Shape::FitS2Opt<3 ) cc[canvas]->Divide(nPad_x,nPad_y);
      cc[canvas]->Draw(); 
     
      for (int ipar=(S2Shape::FitS2Opt==1?0:(S2Shape::FitS2Opt==2?2:3));ipar<4;ipar++)
	{
	  if ( S2Shape::FitS2Opt<3 )
	    {
	      if (S2Shape::FitS2Opt==2) cc[canvas]->cd(ipar-1);
	      else                      cc[canvas]->cd(ipar+1);
	    }

	  int idrift=UNDEF;
	  if ( iPlot>0 ) idrift=(iPlot==1?0:1);
	  
	  int color[3]={kBlue,kRed,kGreen+2};
	  
	  TH2F *h0;
	  if (iPlot==0) h0 =new TH2F(Form("hFit_tdrift_%d",ipar),"",1000, S2Shape::xmin      ,tdrift1,      1000, S2Shape::ymin[ipar],S2Shape::ymax[ipar]);
	  else          h0 =new TH2F(Form("hFit_%d_%d",ipar,idrift),""       ,1000,(ipar==3&&idrift==0?S2Shape::minR:S2Shape::ymin[ipar]),(ipar==3&&idrift==0?S2Shape::maxR:S2Shape::ymax[ipar]),1000,  0.  ,
				     S2Shape::hPar_PDF[idrift][1][ipar]->GetMaximum()*1.2);
	  h0->GetYaxis()->SetTitle((iPlot==0?S2Shape::nameV[ipar].c_str():"Number of Events"));
	  h0->GetXaxis()->SetTitle((iPlot==0?"t_{drift} [#mu s]": (ipar==3&&idrift==0?"#sigma-Pred [#mu s]":S2Shape::nameV[ipar].c_str())  ));
	  h0->DrawCopy();

	  int N[3]={ S2Shape::tPar_Drift[0][ipar]->GetN(), S2Shape::tPar_Drift[1][ipar]->GetN(), S2Shape::tPar_Drift[2][ipar]->GetN() };
	  if ( iPlot==0 ) { l.SetTextColor(kBlack); l.DrawLatex(0.2,0.90,Form("s2>%5.0lf N=%d %d (%d)",S2Shape::s2min,int(N[0]),int(N[1]),int(N[2])));}

	  for (int region=0;region<S2Shape::nRegions-1;region++)
	    {
	      // tdrift dependence 
	      if ( iPlot==0 )
		{
		  TProfile *h0=S2Shape::hPar_Drift[region][ipar];
		  TGraph *t0  =S2Shape::tPar_Drift[region][ipar];
		  if ( t0->GetN()==0 ) continue;
		  t0->SetMarkerSize(0.4);
		  t0->SetMarkerColor(color[region]);
		  t0->DrawClone("P");

		  if ( region==1 )
		    {
		      h0->SetMarkerColor(kBlack); h0->SetLineColor(kBlack); h0->SetMarkerStyle(20);
		      h0->DrawCopy("SAME");

		  
		      TF1 *f=new TF1("fDiffusion",(ipar==3?"sqrt(2.*x*[0]+[1]*[1])/[2]":"[0]+x*[1]"),S2Shape::tdrift_min,60.);
		      if ( ipar==3 ) f->SetParameters(S2Shape::parSigma[region]);
		      else           h0->Fit(f,"QN","",0.,60.);
		      f->SetLineColor(kBlack);
		      f->SetRange(0,60.); f->DrawCopy("SAME");		      
		      if ( ipar==3  ) 
			{
			  f->FixParameter(2,50./(tdrift1-tdrift0));
			  f->FixParameter(0,4.35e-4);
			  h0->Fit(f,"QN","",S2Shape::tdrift_min,60.);
			  f->SetLineColor(kRed);
			  f->DrawCopy("SAME");
			}
		    }

		}
	      else
		{
		  TH1F *h1=S2Shape::hPar_PDF[idrift][region][ipar];
		  h1->SetLineColor(color[region]); h1->SetMarkerColor(color[region]); h1->SetMarkerStyle(30); h1->SetMarkerSize(0.5);
		  h1->DrawCopy("E1SAME");

 		  TF1 *fgaus=new TF1("fgaus",FuncAsymGaus,(ipar==3&&idrift==0?S2Shape::minR:S2Shape::ymin[ipar]),(ipar==3&&idrift==0?S2Shape::maxR:S2Shape::ymax[ipar]),4);
		  fgaus->SetParameter(0,S2Shape::norm[idrift][ipar][region]);
		  fgaus->SetParameter(1,S2Shape::mean[idrift][ipar][region]);
		  fgaus->SetParameter(2,S2Shape::rms1[idrift][ipar][region]);
		  fgaus->SetParameter(3,S2Shape::rms2[idrift][ipar][region]);
		  fgaus->SetLineColor(color[region]);
		  fgaus->DrawCopy("SAME");
		}

	      l.SetTextColor(color[region]);
	      l.SetTextSize(0.035);
	      if ( iPlot>0 )
		{
		  l.DrawLatex(0.2,0.85-0.05*double(region),Form("%s #mu/#sigma=%5.3lf/%5.3lf",nameRegion[region].c_str(),S2Shape::mean[idrift][ipar][region],fabs(S2Shape::rms1[idrift][ipar][region]) ));
		  l.SetTextColor(kBlack); l.DrawLatex(0.2,0.90,Form("s2>%5.0lf t_{drift}%s%2.0lf ",S2Shape::s2min,(idrift==0?">":"<"),S2Shape::tdrift_min));
		}
	      else 
		{
		  l.DrawLatex(0.2,0.85-0.05*double(region),Form("%s",nameRegion[region].c_str()));
		}
	    }//region
	}//ipar
      canvas++;
    }//iPlot

  if ( S2Shape::FitS2Opt==3 ) return;
  //---------------------------------------------------------------------------
  //  XY dependence
  //---------------------------------------------------------------------------
  for (int idrift=0;idrift<1;idrift++)
    for (int iMeanRMS=0;iMeanRMS<2;iMeanRMS++)
      {
	cc[canvas]=new TCanvas(Form("c%d",canvas),"",(S2Shape::FitS2Opt==2?2.0:1.5)*ww,(S2Shape::FitS2Opt==2?1.2:1.5)*wh);
	cc[canvas]->Divide(2,S2Shape::FitS2Opt==2?1:2);
	cc[canvas]->Draw();  
	
	for (int ipar=(S2Shape::FitS2Opt==2?2:0);ipar<4;ipar++)
	  {
	    if (S2Shape::FitS2Opt==2) cc[canvas]->cd(ipar-1);
	    else                      cc[canvas]->cd(ipar+1);
	    
	    
	    double min=(iMeanRMS==0?S2Shape::mean[idrift][ipar][2]-1.5*fabs(S2Shape::rms1[idrift][ipar][2]) : fabs(S2Shape::rms1[idrift][ipar][2]*0.65) );
	    double max=(iMeanRMS==0?S2Shape::mean[idrift][ipar][2]+1.5*fabs(S2Shape::rms1[idrift][ipar][2]) : fabs(S2Shape::rms1[idrift][ipar][2]*1.55) );
	    
	    TH2F *hFrame=new TH2F(Form("hFrame_%d_%d_%d",iMeanRMS,ipar,idrift),"",1000,(PlotXY?-180.:-0.5),(PlotXY?180.:25),1000,min,max);
	    hFrame->GetYaxis()->SetTitle(Form("%s (%s)",(iMeanRMS==0?"Mean":"RMS"),(ipar==3&&idrift==0?"#sigma-Pred [#mu s]":S2Shape::nameV[ipar].c_str())));					     
	    hFrame->GetXaxis()->SetTitle(PlotXY?"#Phi [deg]":"s2 max chan");
	    hFrame->DrawCopy();
	    
	    TGraphErrors *t0=S2Shape::tPar_XY[idrift][ipar][iMeanRMS];	 
	    TGraphErrors *t1[2];
	    if (PlotXY ) { t1[0]=new TGraphErrors(); t1[1]=new TGraphErrors(); GetMapXY(t0 , t1[0],t1[1]);}
	    else         { t1[0]=t0; t1[1]=t0; }


	    for (int ir=0;ir<(PlotXY?2:1);ir++)
	      {
		int color_i=(ir==0?kBlue:kRed);
		TGraphErrors *tp=(PlotXY?t1[ir]:t0);
		if ( PlotXY ) { tp->SetLineColor(color_i); tp->SetMarkerColor(color_i);}
		tp->DrawClone("P");
		l.SetTextColor(color_i); l.DrawLatex(0.2,0.85-0.05*double(ir),Form("%s",nameRegion[ir].c_str()));
		
		TF1 *f=new TF1("f","[0]",-180.,180.);
		tp->Fit(f,"QN","");
		f->SetLineColor(kBlack);//color_i);
		//if ( iMeanRMS==0 ) f->DrawCopy("SAME");
	      }

	    if ( !PlotXY )
	      for (int i=0;i<8;i++)
		{
		  int ch=XYCentral[i];
		  TBox box(double(ch)-0.5,min,double(ch)+0.5,max);
		  box.SetFillStyle(3001);
		  box.SetFillColor(kGreen+2);
		  box.DrawClone();
		}
	    
	    l.SetTextColor(kBlack); l.DrawLatex(0.2,0.90,Form("s2>%5.0lf t_{drift}%s%2.0lf ",S2Shape::s2min,(idrift==0?">":"<"),S2Shape::tdrift_min));
	  }
	canvas++;
      }
}


void SaveCanvasesOnPDF(Int_t run)
{
  //Count how many good canvases!
  Int_t good=0;
  for (size_t l=0;l<MaxCanvas;l++)
    if (cc[l])
      good++;
    else
      break;

  std::cout << "Found " << good << " canvases to save" << std::endl;

  TString pdffilename;
  pdffilename.Form("run_%d_out.pdf",run); 
  for (size_t l=0;l<good;l++)
    {
      if (!l)
	cc[l]->Print(pdffilename+"(","pdf");
      else if (l == (good-1))
	cc[l]->Print(pdffilename+")","pdf");
      else
	cc[l]->Print(pdffilename,"pdf");    
    }
  return;
}
