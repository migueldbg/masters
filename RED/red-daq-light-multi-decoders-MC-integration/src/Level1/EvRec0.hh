#ifndef EvRec0_H
#define EvRec0_H
#include <stdlib.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>


#include<list>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<stdlib.h>
#include<stdio.h>

#include "EvHeader.hh"
#include "../Modules/RDCluster.hh"
#include "../Modules/RDPulseFit.hh"

#include "TROOT.h"
#include "TObject.h"

using namespace std;

//class EvHeader;

/*
  This is the persistent Rec1 Event class
  it should be the event that is dispached to the EvRec0 observers.
*/

class EvRec0 : public TObject {


private:

  EvHeader* evheader;
  int event_number;
  double charge_total, f90_tot;
  double charge_total_lsci, lsci_psd_tot;
  double baseline_mean_all;
  double baseline_rms_all;
  vector <double> baseline_mean;
  vector <double> baseline_rms;
  vector <double> charge;
  vector <double> f90;
  vector <double> start_time;
  vector <double> xmin;
  vector <double> xmax;
  vector <double> ymin;
  vector <double> ymax;
  int number_of_clusters;
  vector<RDCluster*> clusters;

  int nfits;
  vector<RDPulseFit*> fits;

public:

  EvRec0();
  ~EvRec0(){delete evheader;};

  // Copy constructor: default only with C++11
  //EvRec0(const EvRec0&);

  // set methods
  void SetEvHeader(EvHeader *e)         {evheader =  e;};
  void SetEvNumber(int n)               {event_number = n;};
  void SetChargeTot(double qt)          {charge_total = qt;};
  void SetChargeTotLSci(double qt)      {charge_total_lsci = qt;};
  void SetLSciPSDTot(double val)        {lsci_psd_tot = val;};
  void SetF90Tot(double ft)             {f90_tot = ft;};
  //void SetNPeaks(vector <int> np)       {_number_of_peaks = np;};
  void SetBaseMean(vector <double> bm)  {baseline_mean = bm;};
  void SetBaseRMS(vector <double> rms)  {baseline_rms = rms;};
  void SetCharge(vector <double> q)     {charge = q;};
  void SetF90(vector<double> f)         {f90 = f;};
  void SetStartTime(vector <double> st) {start_time = st;};
  void SetXmin(vector <double> val)     {xmin = val;};
  void SetXmax(vector <double> val)     {xmax = val;};
  void SetYmin(vector <double> val)     {ymin = val;};
  void SetYmax(vector <double> val)     {ymax = val;};
  void SetBaselineMeanAll(double val)   {baseline_mean_all = val;};
  void SetBaselineRMSAll(double val)    {baseline_rms_all = val;};
  void AddCluster(RDCluster* aCluster);
  void AddFit(RDPulseFit* aFit);

  // get methods
  EvHeader* GetEvHeader()         {return evheader;};
  int GetEvNumber()               {return event_number;};
  double GetChargeTot()           {return charge_total;};
  double GetF90Tot()              {return f90_tot;};
  double GetChargeTotLSci()       {return charge_total_lsci;};
  double GetLSciPSDTot()          {return lsci_psd_tot;};
  double GetBaselineMeanAll()     {return baseline_mean_all;};
  double GetBaselineRMSAll()      {return baseline_rms_all;};
  vector <double> GetBaseMean()   {return baseline_mean;};
  vector <double> GetBaseRMS()    {return baseline_rms;};
  vector <double> GetCharge()     {return charge;};
  vector <double> GetF90()        {return f90;};
  vector <double> GetStartTime()  {return start_time;};
  vector <double> GetXmin()       {return xmin;};
  vector <double> GetXmax()       {return xmax;};
  vector <double> GetYmin()       {return ymin;};
  vector <double> GetYmax()       {return ymax;};

  int GetNClusters()     {return number_of_clusters;};
  vector<RDCluster*> GetClusters(){return clusters;};
  RDCluster* GetCluster(int id)   {return clusters.at(id);};

  int GetNfits()     {return nfits;};
  vector<RDPulseFit*> GetFits(){return fits;};
  RDPulseFit* GetFits(int id)   {return fits.at(id);};


  void Fill(EvHeader * evh, int evn, double charget, double f90_t, vector <double> basel_mean,
	    vector <double> basel_rms, vector <double> charge, vector <double> f90,
	    vector <double> star_t, vector <double> xmin, vector <double> xmax,
	    vector <double> ymin, vector <double> ymax, vector <int> number_p);
  void Init();
  void Clear();
  void Dump();

  ClassDef( EvRec0, 6 );

};

#endif
