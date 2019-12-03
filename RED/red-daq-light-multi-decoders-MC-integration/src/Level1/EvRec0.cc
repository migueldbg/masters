#include "EvRec0.hh"


EvRec0::EvRec0() :
  evheader(0),clusters(0),fits(0)
{
  //cout << "EvRec0::EvRec0()" << endl;  
}

//=============================================================================

void EvRec0::Init()
{
  //cout << "EvRec0::Init" << endl;
 
  evheader = 0;
  event_number = 0;
  charge_total = 0;
  f90_tot = 0;
  charge_total_lsci = 0.; 
  lsci_psd_tot = 0.;
  baseline_mean_all = 0;
  baseline_rms_all = 0;
  for (size_t t=0; t<baseline_mean.size(); t++) {
    baseline_mean[t] = 0.;
  }
  for (size_t t=0; t<baseline_rms.size(); t++) {
    baseline_rms[t] = 0.;
  }
  for (size_t t=0; t<f90.size(); t++) {
    f90[t] = 0.;
  }
  for (size_t t=0; t<start_time.size(); t++) {
    start_time[t] = 0.;
  }
  for (size_t t=0; t<xmin.size(); t++) {
    xmin[t] = 0.;
  }
  for (int t=0; t<xmax.size(); t++) {
    xmax[t] = 0.;
  }
  for (size_t t=0; t<ymin.size(); t++) {
    ymin[t] = 0.;
  };
  for (int t=0; t<ymax.size(); t++) {
    ymax[t] = 0.;
  }
  for (int t=0; t<charge.size(); t++) {
    charge[t] = 0.;
  }
  /*
  for (int t=0; t<number_of_peaks.size(); t++) {
    number_of_peaks[t] = 0.;
  };
  */
};

//=============================================================================

void EvRec0::Clear()
{    
  if (evheader)
    delete  evheader;
  evheader = 0;

  for (size_t i=0;i<clusters.size();i++)
    delete clusters.at(i);
  clusters.clear();
  number_of_clusters = 0;

  for (size_t i=0;i<fits.size();i++)
    delete fits.at(i);
  fits.clear();
  nfits = 0;


  event_number = 0;
  charge_total = 0;
  f90_tot = 0;
  charge_total_lsci = 0; 
  lsci_psd_tot = 0;
  baseline_mean_all = 0;
  baseline_rms_all = 0;

  baseline_mean.clear();
  baseline_rms.clear();
  f90.clear();
  start_time.clear();
  xmin.clear();
  xmax.clear();
  ymin.clear();
  ymax.clear();
  charge.clear();
  //number_of_peaks.clear(); 
}

//=============================================================================
void EvRec0::AddCluster(RDCluster* aCluster)
{
  clusters.push_back(aCluster);
  number_of_clusters = (int) clusters.size();
}

//=============================================================================
void EvRec0::AddFit(RDPulseFit* aFit)
{
  fits.push_back(aFit);
  nfits = (int) fits.size();
}

//============================================================================

void EvRec0::Fill(EvHeader * evh, int evn, double charget, double f90_t, vector <double> basel_mean,
        vector <double> basel_rms, vector <double> charge, vector <double> af90,
        vector <double> star_t, vector <double> axmin, vector <double> axmax,
        vector <double> aymin, vector <double> aymax, vector <int> number_p) {

        evheader = evh;
        event_number = evn;
        charge_total = charget;
        f90_tot = f90_t;
        baseline_mean = basel_mean;
        baseline_rms = basel_rms;
        f90 = af90;
        start_time = star_t;
        xmin = axmin;
        xmax = axmax;
        ymin = aymin;
        xmax = aymax;
        //number_of_peaks = number_p;

    
};

void EvRec0::Dump() {
    
    for (size_t t=0; t<baseline_mean.size(); t++) {
        
        cout << "baseline_mean[" << t << "]=" << baseline_mean[t] << endl;
    }
    for (int t=0; t<baseline_rms.size(); t++) {
        
        cout << "baseline_rms[" << t << "]=" << baseline_rms[t] << endl;
    }
    for (int t=0; t<f90.size(); t++) {
        
        cout << "f90[" << t << "]=" << f90[t] << endl;
    }
    for (int t=0; t<start_time.size(); t++) {
        
        cout << "start_time[" << t << "]=" << start_time[t] << endl;
    }
    for (int t=0; t<xmin.size(); t++) {
        
        cout << "xmin[" << t << "]=" << xmin[t] << endl;
    }    
    for (int t=0; t<xmax.size(); t++) {
        
        cout << "xmax[" << t << "]=" << xmax[t] << endl;
    };   
    for (int t=0; t<ymin.size(); t++) {
        
        cout << "ymin[" << t << "]=" << ymin[t] << endl;
    }   
    for (int t=0; t!=ymax.size(); t++) {
        
        cout << "ymax[" << t << "]=" << ymax[t] << endl;
    };   
    /*
    for (int t=0; t!=_number_of_peaks.size(); t++) {
        
        cout << "number_of_peaks[" << t << "]=" << _number_of_peaks[t] << endl;
    };
    */

};

ClassImp ( EvRec0 )



