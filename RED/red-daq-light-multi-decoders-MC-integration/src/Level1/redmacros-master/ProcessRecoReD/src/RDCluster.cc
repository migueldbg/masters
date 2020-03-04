#include "RDCluster.hh"

//============================================================
RDCluster::RDCluster()
{
  start_time = 0;
  start = 0;
  stop = 0;
  rep = 0; 
  charge = 0.;
  f90 = 0.;
  mean_time = 0.;
  rms_time = 0.;
  min_x = 0.;
  min_y = 0.;
  max_x = 0.;
  max_y = 0.;
  cdf_time = 0.;
  fixed_time = 0.;  
  tot_charge_top = 0.;
  tot_charge_bottom = 0.;
  bar_x = 0;
  bar_y = 0;
  pos_x = 0;
  pos_y = 0;
  charge_top.clear();
  charge_bottom.clear();
}

//============================================================
RDCluster::~RDCluster()
{;}

//============================================================
RDCluster::RDCluster(const RDCluster& right)
{
  start_time = right.start_time;
  start = right.start;
  stop = right.stop;
  rep = right.rep; 
  charge = right.charge;
  f90 = right.f90;
  mean_time = right.mean_time;
  rms_time = right.rms_time;
  min_x = right.min_x;
  min_y = right.min_y;
  max_x = right.max_x;
  max_y = right.max_y;
  cdf_time = right.cdf_time; 
  fixed_time = right.fixed_time; 
  tot_charge_top = right.tot_charge_top;
  tot_charge_bottom = right.tot_charge_bottom;
  bar_x = right.bar_x;
  bar_y = right.bar_y;
  pos_x = right.pos_x;
  pos_y = right.pos_y;
  charge_top = right.charge_top;
  charge_bottom = right.charge_bottom;
}

//============================================================
const RDCluster& RDCluster::operator=(const RDCluster& right)
{
  start_time = right.start_time;
  start = right.start;
  stop = right.stop;
  rep = right.rep; 
  charge = right.charge;
  f90 = right.f90;
  mean_time = right.mean_time;
  rms_time = right.rms_time;
  min_x = right.min_x;
  min_y = right.min_y;
  max_x = right.max_x;
  max_y = right.max_y;
  cdf_time = right.cdf_time; 
  fixed_time = right.fixed_time; 
  tot_charge_top = right.tot_charge_top;
  tot_charge_bottom = right.tot_charge_bottom;
  bar_x = right.bar_x;
  bar_y = right.bar_y;
  pos_x = right.pos_x;
  pos_y = right.pos_y;
  charge_top = right.charge_top;
  charge_bottom = right.charge_bottom;
  return *this;
}

//============================================================
int RDCluster::operator==(const RDCluster& right) const
{
  return (start_time != right.start_time || 
	  start != right.start || 
	  stop != right.stop ||
	  rep != right.rep || 
	  charge != right.charge ||
	  f90 != right.f90 ||
	  mean_time != right.mean_time ||
	  rms_time != right.rms_time ||
	  min_x != right.min_x || 
          cdf_time != right.cdf_time ||
          fixed_time != right.fixed_time)
    ? 0 : 1;
}

ClassImp ( RDCluster)
