// Clusterization algorithm
#ifndef RDPULSEFITTER_H
#define RDPULSEFITTER_H

#include <vector>

#include <TMath.h>
#include <TGraphErrors.h> 

#ifndef UNDEF
#define UNDEF -100
#endif

using namespace TMath;

static const double FADCWidth=2.;
static const bool UseDS50Binning=false;
static const double nBinsDS50=100;
static const double fpoisson=0.15;
static const int nParPulseFitter=11;

//-------------------------------------------------------------

typedef struct PulseFitInfo {
  double par[nParPulseFitter],epar[nParPulseFitter];
  int type;
  int status,covstatus,ndf;
  double chi2;
  int sipm, start,end;
} PulseFitInfo_t;

//--------------------------------------------------------------

class RDPulseFitter {
private:

  std::vector<PulseFitInfo_t> fpulsefit;

public:
  RDPulseFitter();
  ~RDPulseFitter(){;};
  void Clear();

  void SetUseBaker(bool val){fUseBaker = val;};
  bool IsUsingBaker(){return fUseBaker;};

  bool IsDS50Binning();

  void SetDebug(bool val){fDebug = val;};
  bool GetDebug(){return fDebug;};

  void ReBin( std::vector<double> *wf,  int start, int stop , TGraphErrors *tpeak, double rmsBaseline );
  void FitS2( std::vector<double> *wfc, int start, int stop , double rms);
  void FitSPE(std::vector<double> *wfc, int start, int stop,  int isipm, double rms);
  
  double yfit_sa(double *x,const double *par);
  double ImpulseFunc(double *v, double *par);

  int    GetnFits() const { return fpulsefit.size(); }  
  int    fit_status(int i) const   { return fpulsefit.at(i).status;}
  int    fit_covstatus(int i) const   { return fpulsefit.at(i).covstatus;}
  int    fit_ndf(int i) const   { return fpulsefit.at(i).ndf;}
  int    fit_type(int i) const   { return fpulsefit.at(i).type;}
  int    fit_start(int i) const   { return fpulsefit.at(i).start;}
  int    fit_end(int i) const   { return fpulsefit.at(i).end;}
  int    fit_sipm(int i) const   { return fpulsefit.at(i).sipm;}
  double fit_chi2(int i) const   { return fpulsefit.at(i).chi2;}

  int    GetnPar(int i)       const  { return fpulsefit.size(); }
  double fit_par(int i,int j) const  
  { 
    if  (j>=0 && j<nParPulseFitter )  return fpulsefit.at(i).par[j];
    else                              return -1;
  }
  double fit_epar(int i,int j) const 
  { 
    if  (j>=0 && j< nParPulseFitter  ) return fpulsefit.at(i).epar[j];
    else                               return -1;
  }

private:
  bool fUseBaker;
  bool fDebug;
  void ReBin( std::vector<double> *wf, int start, int stop, TGraphErrors* tpeak, double rmsBaseline, 
              std::vector<double> &tLow, std::vector<double> &tHigh, std::vector<double> &tMean, std::vector<double> &fChg);
  double exp_erfc(double t, double tau, double sigma);
  double y2(double t, double tau, double sigma, double tau_spe, double p_spe, double sigma_spe);
  double y1(double t, double tau, double T, double sigma, double tau_spe, double p_spe, double sigma_spe);
  double s2_model(double*x, double* par);
  double fit_model(double*x, const double* par);
};

#endif
