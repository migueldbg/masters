#include "RDPulseFitter.hh"
#include "../Level1/RDconfig.h"

#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cmath>

#include <TF1.h>
#include <TFitResult.h>

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TStopwatch.h>

using namespace std;

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

RDPulseFitter::RDPulseFitter() :
    fUseBaker(false),fDebug(false)
{;}

//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

void RDPulseFitter::Clear() 
{
    fpulsefit.clear();
}

bool RDPulseFitter::IsDS50Binning()
{ return UseDS50Binning;}

void RDPulseFitter::ReBin( std::vector<double> *wf,  int start, int stop , TGraphErrors *tpeak, double rmsBaseline, 
                           std::vector<double> &tLow, std::vector<double> &tHigh, std::vector<double> &tMean,  std::vector<double> &fChg )
{
    tLow.clear(); tHigh.clear(); tMean.clear(); fChg.clear();

    //------------------
    // Rebin data
    //------------------

    int firstPE=UNDEF;
    double chg=0.;
    for (int ii=start; ii<stop;ii++)
    {
        chg+=(-wf->at(ii));
        if ( chg>0.5 && firstPE==UNDEF ) firstPE=ii;
    }
    if ( firstPE==UNDEF ) return;

    const double minW=50;
    const double minChg=(UseDS50Binning?chg/nBinsDS50:5.);


    for (int ii=firstPE; ii<stop;ii++)
    {
        double t_i=double(ii);
        double chg_i=0.,mtime=0.,rtime=0.;
        for ( int jj=ii;jj<stop;jj++)
        {
            double chg_j=-wf->at(jj);
            double t_j  =double(jj);
            chg_i+=  chg_j;
            mtime+= (chg_j*t_j);
            rtime+= (chg_j*t_j*t_j);

            if (  t_j-t_i >minW && chg_i>minChg )
            {
                mtime/=chg_i;
                rtime/=chg_i;
                double x=mtime;
                double rms=sqrt(rtime-mtime*mtime)/sqrt(chg_i-1);
                double ex=rms;

                double dt=(t_j-t_i+1);

                double y=chg_i/dt;
                double ep=fpoisson*sqrt(chg_i/fpoisson)/dt;
                double ey=sqrt(ep*ep+rmsBaseline*rmsBaseline);

                int ipoint=tpeak->GetN();
                tpeak->SetPoint(ipoint,x,-y);
                tpeak->SetPointError(ipoint,ex,ey);

                tLow.push_back(t_i); tHigh.push_back(t_j);  tMean.push_back(mtime); fChg.push_back(chg_i);

                ii=jj;
                break;
            }
        }
    }//ii
}

void RDPulseFitter::ReBin( std::vector<double> *wf,  int start, int stop , TGraphErrors *tpeak, double rmsBaseline )
{
    std::vector<double> tLow,tHigh,tMean,Chg;
    ReBin(wf,start,stop,tpeak,rmsBaseline,tLow,tHigh,tMean,Chg);
}

void RDPulseFitter::FitS2( std::vector<double> *wf,  int start, int stop, double rmsBaseline ) // Cumulative trace
{
    if ( fDebug ) printf("Fitting S2, %d-%d \n",start,stop);
    int n_samp=wf->size();
    if ( start<0 ) start=0;
    if ( stop>=n_samp ) stop=n_samp;

    TStopwatch w;
    w.Start();

    PulseFitInfo_t tfit;

    tfit.start=start;
    tfit.end=(stop-start>10000?start+10000:stop);
    tfit.status=-1; tfit.covstatus=-1; tfit.ndf=0; tfit.chi2=-1;
    for (int ipar=0;ipar<nParPulseFitter;ipar++) { tfit.par[ipar]=0.; tfit.epar[ipar]=0.;}

    tfit.type=1;
    tfit.sipm=-1;

    double chg=0.;
    for (int ii=tfit.start;ii<tfit.end;ii++)
        chg+=wf->at(ii);

    //------------------
    // ReBin data
    //-----------------
    TGraphErrors *tpeak=new TGraphErrors();
    std::vector<double> tLow,tHigh,tMean,fChg;
    ReBin(wf,tfit.start,tfit.end, tpeak, rmsBaseline,tLow,tHigh,tMean,fChg);

    if ( tpeak->GetN()==0 || tLow.size()==0 ) { delete tpeak; return; }

    if ( fDebug )
        for (int i=0;i<tLow.size();i++)
            printf(" t=%d-%d  <t>=%5.2lf Chg=%5.2lf \n",int(tLow[i]),int(tHigh[i]),tMean[i],fChg[i]);    

    //------------------
    // Function to minimize
    //-----------------
    auto chi2Function = [&](const double *par)
    {
        Double_t f=0;

        for (size_t i=0;i<tLow.size();i++)
        {
            if ( !fUseBaker  )
            {
                double t[1]={ tMean[i] };
                double dt=tHigh[i]-tLow[i]+1;
                double rho=fChg[i]/dt;

                double rho_exp=fabs(yfit_sa(t,par));
                double erho=fpoisson*sqrt(rho*dt/fpoisson)/dt;

                f+= pow( (rho-rho_exp)/erho,2);// Minimization
            }
            else
            {
                double expected=0.;
                //const int nSteps=10;
                const int nSteps=2;
                double step=(tHigh[i]-tLow[i])/double(nSteps);
                for ( int ii=0;ii<nSteps;ii++)
                {
                    double x[2]={ tLow[i]+double(ii)*step, x[0]+step };
                    double y[2]={  fabs(yfit_sa(&(x[0]),par)),  fabs(yfit_sa(&(x[1]),par)) }; // Numerically should be always>0
                    //if ( y[0]==0. || y[1]==0. ) printf("Warning: numerical problem in chi2\n");
                    expected+= (y[0]+y[1])/2.;
                }
                expected*=(step/fpoisson);
                double observed=fChg[i]/fpoisson;
                double logExpected=(expected<=std::numeric_limits<double>::min()?log(std::numeric_limits<double>::min()):log(expected));
                f+=(  expected- observed + observed*( log(observed)-logExpected ) ); //Poisson
            }
        }
        return (fUseBaker?2.:1.)*f;
    };

    //------------------
    // Minimizer
    //------------------
    //       minName                   algoName
    //   Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
    //   Minuit2                     Fumili2
    //   Fumili
    //   GSLMultiMin                  ConjugateFR, ConjugatePR, BFGS,
    //                                BFGS2, SteepestDescent
    //   GSLMultiFit
    //   GSLSimAn
    //   Genetic

    ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
    minimum->SetMaxFunctionCalls(1000); // for Minuit/Minuit2
    minimum->SetMaxIterations(10000);  // for GSL
    minimum->SetTolerance(0.001);
    minimum->SetPrintLevel(fDebug?1:0);
    minimum->SetErrorDef(1.);

    //-----
    // Variable setting and minimization
    //-----

    const int npar=11;
    ROOT::Math::Functor fcn(chi2Function,npar);
    minimum->SetFunction(fcn);

    double sigmaL=0.1;
    double t0    = TMath::Max(0.10,0.10+2.00*sigmaL);
    double chgS2 = 2.*chg/1.e3;
    double tau1  =11.e-3;
    double tau2  =3.43;
    double pgar  =0.1;
    double T     =1.2;

    double tau_spe=575.e-3, p_spe=0.06, sigma_spe=7.5e-3;
    /*
    // Use the extended ser file mean parameters for SiPM
    auto tau_sipm   = RDconfig::GetInstance()->GetSER_Tau();
    auto p_sipm     = RDconfig::GetInstance()->GetSER_p();
    auto sigma_sipm = RDconfig::GetInstance()->GetSER_Sigma();
    
    double tau_spe=0., p_spe=0., sigma_spe=0.;
    for (size_t i=0; i<tau_sipm.size(); i++)
    {
        tau_spe   += tau_sipm.at(i);
        p_spe     += p_sipm.at(i);
        sigma_spe += sigma_sipm.at(i);
    }
    
    tau_spe = tau_spe/(1.0*tau_sipm.size())*1.e-3;
    p_spe = p_spe/(1.0*p_sipm.size());
    sigma_spe = sigma_spe/(1.0*sigma_sipm.size())*1.e-3;
        
    //printf("tau_spe=%E p_spe=%E sigma_spe=%E\n", tau_spe, p_spe, sigma_spe);
    */
    double init[npar]=           {        tau1    ,    tau2  ,    pgar ,    T ,   sigmaL,   chgS2   ,  t0  ,  double(start)     , tau_spe  , p_spe   ,  sigma_spe    };
    double low[npar] =           {      UNDEF     ,    UNDEF  ,     0. ,   0.1,    0.01 ,   UNDEF   ,  0.  ,      UNDEF         ,   UNDEF  ,  UNDEF  ,  UNDEF        };
    double high[npar]=           {      UNDEF     ,    UNDEF  ,     1. ,    5.,      5. ,   UNDEF   ,  5.  ,      UNDEF         ,   UNDEF  ,  UNDEF  ,  UNDEF        };
    double step[npar]=           {        0.      ,    0.1    ,   0.01 ,   0.1,    0.01 , 0.1*fabs(chgS2) , 0.01 ,       0.           ,      0.  ,     0.  ,      0.       };
    const std::string name[npar]={      "tau1"    ,   "tau2"  ,  "pgar",   "T",  "sigma", "chgS2"   , "t0" ,   "startS2"        , "tau_spe", "p_spe" , "sigma_spe"   };

    for (int ipar=0;ipar<npar;ipar++)
    {
        bool IsFixed  =( ipar==0  ||  ipar>=7 );
        bool IsLimited=( !IsFixed && low[ipar]!=UNDEF );

        if ( IsFixed )           minimum->SetFixedVariable(ipar,name[ipar] ,init[ipar]);
        else if ( IsLimited )    minimum->SetLimitedVariable(ipar,name[ipar] ,init[ipar], step[ipar], low[ipar], high[ipar] );
        else                     minimum->SetVariable(ipar,name[ipar] ,init[ipar],step[ipar] );
    }


    minimum->FixVariable(2);
    minimum->Minimize();
    minimum->ReleaseVariable(2);
    minimum->Minimize();
    
    const double chi2=(std::isnan(minimum->MinValue())? UNDEF: minimum->MinValue() );
    const double ndf=double(tLow.size());
    const int status=minimum->Status();
    const int covstatus=minimum->CovMatrixStatus();

    //----------------------------
    //   Save results
    //----------------------------

    const double *par = minimum->X(), *epar=minimum->Errors();

    tfit.status=status;
    tfit.covstatus=covstatus;
    tfit.chi2=chi2;
    tfit.ndf=ndf;

    for (int ipar=0;ipar<11;ipar++)
    { tfit.par[ipar]=par[ipar]; tfit.epar[ipar]=epar[ipar]; }

    fpulsefit.push_back(tfit);

    if (fDebug)
    {
        printf("--------------------------\n");
        printf("Chi2=%5.2lf ndf=%d | Status=%d %d  | Edm=%5.2e \n",chi2,int(ndf),status,covstatus,minimum->Edm());
        printf(" tau1=%5.2lf (%5.2lf)  tau2=%5.2lf (%5.2lf)  p=%5.2lf (%5.2lf)   T=%5.2lf  (%5.2lf)   Sigma=%5.2lf (%5.2lf)    norm=%5.2lf (%5.2lf)  t0=%5.2lf (%5.2lf)  start=%5.0lf (%5.0lf) \n",
               par[0],epar[0],	   par[1],epar[1],	   par[2],epar[2],	   par[3],epar[3],
                par[4],epar[4],	   par[5],epar[5],	   par[6],epar[6],	   par[7],epar[7]);
        printf(" tau_spe=%5.2lf (%5.2lf)    p_spe=%5.2lf  (%5.2lf)  sigma_spe=%5.0lf (%5.0lf)   \n",
               par[8],epar[8],	   par[9],epar[9],	   par[10],epar[10] );

        printf("\n");
        printf("Niter= %d   Time=%lf [sec]  CPU time=%lf [sec] \n",minimum->NIterations(), w.RealTime(),  w.CpuTime());
        printf("----------------------------------------------\n");
    }

    delete minimum;
    delete tpeak;

    //------
    w.Stop();
}

void RDPulseFitter::FitSPE( std::vector<double> *wf, int start, int stop,  int isipm, double rms )
{
    PulseFitInfo_t tfit;
    tfit.start=0; tfit.end=0;
    tfit.status=-1; tfit.covstatus=-1; tfit.ndf=0; tfit.chi2=-1;
    for (int ipar=0;ipar<nParPulseFitter;ipar++) { tfit.par[ipar]=0.; tfit.epar[ipar]=0.;}

    tfit.type=0;
    tfit.sipm=isipm;

    int n_samp=wf->size();
    if ( start<0 ) start=0;
    if ( stop>=n_samp ) stop=n_samp;

    int window=30;
    double thr=rms*sqrt(double(window))*10.;

    bool HasPE=false;
    // Find One PE
    for (int ii=start;ii<stop;ii++)
    {
        double sum=0.;
        for (int sa=ii;sa<ii+window;sa++)
            sum+=-wf->at(sa);
        
        if ( fDebug )
            printf("Partial sum=%+.5E, thr=%.5E for ii=%d\n", sum, thr, ii);
        
        if ( sum > thr )
        {
            tfit.start=ii-30;
            tfit.end=ii+1500;
            HasPE=true;
            break;
        }
    }

    if ( fDebug )
        printf("%s\n", (HasPE)?"YES PE":"NO PE");
    if ( !HasPE ) return;

    //-------------------------------
    TGraphErrors *tg=new TGraphErrors();
    double tmax=UNDEF, ymax=0.;
    for (int sa=tfit.start;sa<tfit.end;sa++)
    {
        int ipoint=tg->GetN();
        double val=wf->at(sa);
        if ( val>ymax ) { ymax=val; tmax=double(sa); }
        tg->SetPoint(ipoint,double(sa),val);
        tg->SetPointError(ipoint,0.,rms*1.25);
    }

    //--------------------------------
    double chg=0;
    for (int sa=tfit.start;sa<tfit.end;sa++)
        chg+=wf->at(sa);

    double par[7]= {  double(tfit.start+60), 575., 0.1, 10., chg*FADCWidth, chg*FADCWidth , double(tfit.end+tfit.start) };
    double epar[7]={   0., 0., 0., 0., 0., 0., 0. };
    //TF1 *fmodel = new TF1("fmodel", ImpulseFunc, double(tfit.start), double(tfit.end), 7);
    TF1 *fmodel = new TF1("fmodel", this, &RDPulseFitter::ImpulseFunc, double(tfit.start), double(tfit.end), 7);
    fmodel->SetNpx(10000);

    fmodel->SetParLimits(0, double(tfit.start), double(tfit.end));
    fmodel->SetParLimits(1, 1.0, 2.e3);
    fmodel->SetParLimits(2, 0. , 1.  );
    fmodel->SetParLimits(3, 1.0, 50. );
    fmodel->SetParameters(par);
    fmodel->FixParameter(5,par[5]);
    fmodel->FixParameter(6,par[6]);

    TFitResultPtr rp=tg->Fit("fmodel", "QNS", "", double(tfit.start), double(tfit.end));
    TFitResult*    r=rp.Get();

    for (int ipar=0;ipar<7;ipar++)
    {
        par[ipar]=fmodel->GetParameter(ipar);
        epar[ipar]=fmodel->GetParError(ipar);
    }

    //----------------------------
    //   Save results
    //----------------------------

    tfit.status=r->Status();
    tfit.covstatus=r->CovMatrixStatus();
    if ( fDebug >1 ) r->Print("V");
    tfit.chi2=fmodel->GetChisquare();
    tfit.ndf=double(fmodel->GetNDF());

    //---Fit parameters
    for (int ipar=0;ipar<7;ipar++)
    { tfit.par[ipar]=par[ipar]; tfit.epar[ipar]=epar[ipar]; }

    //-- Other usefull variables
    tfit.par[7]=fmodel->Eval(tfit.par[0]);            // Value of model at t0
    tfit.par[8]=tg->Eval(tfit.par[0]);                // Value of data  at t0
    tfit.par[9]=tg->Eval(tmax);                       // Maximum value of data in ROI
    tfit.par[10]=tmax;                                 // Position of maximum value of data in ROI

    fpulsefit.push_back(tfit);
    delete tg;
    delete fmodel;

    if ( fDebug )
    {
        printf(" sipm=%d  chg=%5.2lf   rms=%5.2lf sa=%d-%d t0=%5.2lf (%5.2lf ) tau=%5.2lf (%5.2lf) p=%5.2lf (%5.2lf) sigma=%5.2lf (%5.2lf) Q=%5.2lf (%5.2lf) status=%d covstatus=%d chi2=%5.4lf ndf=%d \n", tfit.sipm,chg,rms, tfit.start,tfit.end, tfit.par[0], tfit.epar[0], tfit.par[1], tfit.epar[1], tfit.par[2],tfit.epar[2], tfit.par[3],tfit.epar[3], tfit.par[4]/FADCWidth, tfit.epar[4]/FADCWidth, tfit.status, tfit.covstatus, tfit.chi2,tfit.ndf);

        printf("sipm=%d  Max(Data)= %5.4lf at %5.3lf [sa]  , Max(Model)=%5.4lf (%5.4lf) at %5.3lf [sa] \n",tfit.sipm,tfit.par[8],tfit.par[9],  tfit.par[6],tfit.par[7], tfit.par[0]);
    }
}

double RDPulseFitter::exp_erfc(double t, double tau, double sigma)
{
    return Exp(-t/tau + Sq(sigma)/(2.0*Sq(tau)) + Log(Erfc((Sq(sigma)-t*tau)/(Sqrt2()*sigma*tau))));
}

double RDPulseFitter::y2(double t, double tau, double sigma, double tau_spe, double p_spe, double sigma_spe)
{
    double sigma_total = Hypot(sigma, sigma_spe), p = p_spe;
    
    if(tau_spe == 0.0 || p_spe == 1.0) {
        return Erf(t/(Sqrt2()*sigma_total)) - exp_erfc(t, tau, sigma_total);
    }
    return Erf(t/(Sqrt2()*sigma_total)) - ((tau-p*tau_spe)*exp_erfc(t, tau, sigma_total) -
                                           (1-p)*tau_spe*exp_erfc(t, tau_spe, sigma_total))/(tau-tau_spe);
}

double RDPulseFitter::y1(double t, double tau, double T, double sigma, double tau_spe, double p_spe, double sigma_spe)
{
    return (y2(t, tau, sigma, tau_spe, p_spe, sigma_spe) - y2(t-T, tau, sigma, tau_spe, p_spe, sigma_spe))/(2.0*T);
}

double RDPulseFitter::s2_model(double*x, double* par)
{
    double t         = x[0];
    double tau1      = par[0];
    double tau2      = par[1];
    double p         = par[2];
    double T         = par[3];
    double sigma     = par[4];
    double tau_spe   = par[5];
    double p_spe     = par[6];
    double sigma_spe = par[7];
    
    return p*y1(t, tau1, T, sigma, tau_spe, p_spe, sigma_spe) + (1-p)*y1(t, tau2, T, sigma, tau_spe, p_spe, sigma_spe);
}

double RDPulseFitter::fit_model(double*x, const double* par)
{
    double t         = x[0];
    double tau1      = par[0];
    double tau2      = par[1];
    double p         = par[2];
    double T         = par[3];
    double sigma     = par[4];
    double A         = par[5];
    double t0        = par[6];
    double start     = par[7];
    double tau_spe   = par[8];
    double p_spe     = par[9];
    double sigma_spe = par[10];
    
    return A*(p*y1(t-t0, tau1, T, sigma, tau_spe, p_spe, sigma_spe) + (1-p)*y1(t-t0, tau2, T, sigma, tau_spe, p_spe, sigma_spe));
}

double RDPulseFitter::yfit_sa(double *x,const double *par)
{
    double sa=x[0];
    double t[1]={(sa-par[7])*FADCWidth/1.e3};
    
    double val=fit_model(t,par);
    return val;
}

double RDPulseFitter::ImpulseFunc(double *v, double *par) // The SPE pulse
{
    const double t0=par[0],tau=par[1],fslow=1.-par[2],sigma=par[3],Q=par[4],Q0=par[5],dsa=par[6];
    double x=(v[0]-t0)*FADCWidth;

    const double shift=(Q0-Q)/dsa;
    const double MPI=TMath::Pi();

    double alpha = 1./tau;
    double sigmasq = sigma*sigma;
    double sigmasq_alpha = sigmasq*alpha;
    double slow = 0.5 * alpha * exp(-alpha * (x - 0.5*sigmasq_alpha)) * (1 + TMath::Erf((x - sigmasq_alpha)/(sqrt(2)*sigma))) ;

    double fast = 1/(sqrt(2*MPI)*sigma)*exp(-x*x/(2*sigmasq));

    return Q*(fast*(1.-fslow)+slow*fslow)+shift;   
}
