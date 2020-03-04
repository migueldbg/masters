#ifndef _RDcluster_
#define _RDcluster_

#include "TObject.h"

#include <vector>

using namespace std;

class RDCluster : public TObject
{ 
 public:
  RDCluster();
  ~RDCluster();
  RDCluster(const RDCluster&); //copy constructor
  const RDCluster& operator=(const RDCluster&); //assignment
  int operator==(const RDCluster&) const; //equality
  

  int start_time;
  int start;
  int stop;
  int rep; 
  double charge;
  double f90;
  double mean_time;
  double rms_time;
  double min_x;
  double min_y;
  double max_x;
  double max_y;
  double bar_x; //barycenter
  double bar_y;
  double pos_x; //s1/s2max
  double pos_y;
  double cdf_time;
  double fixed_time;
  double tot_charge_top;
  double tot_charge_bottom;
  vector<double> charge_top;
  vector<double> charge_bottom;

  ClassDef(RDCluster, 7);
  
};

#endif
