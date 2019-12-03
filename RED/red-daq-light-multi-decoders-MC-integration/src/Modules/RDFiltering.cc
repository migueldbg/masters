// Filtering Module
#include "RDFiltering.hh"
#include <iostream>
using namespace std;

//==============================================================
// Constructor
RDFiltering::RDFiltering(const vector<double> *wf, int nsamp)
{
  nsamp_ = nsamp;
  wf_ = new vector<double>(wf->size());
  for (size_t i=0;i<wf->size();i++)
    wf_->at(i) = wf->at(i);
  ma_ = new vector<double>(nsamp_);
  sc_ = new vector<double>(nsamp_);
  hp_ = new vector<double>(nsamp_);
  sn_ = new vector<double>(nsamp_);
}

//==============================================================
// Destructor

RDFiltering::~RDFiltering()
{ 
  delete ma_;
  delete sc_;
  delete hp_;
  delete sn_;
  delete wf_;
}

//==============================================================
// Moving average
vector<double>* RDFiltering::mavg(int step) 
{
  vector<double> ma_temp(nsamp_ + step);
   
  int n0 = (step/2);

  for (int i = 0; i < nsamp_ + step; i++) 
    {
      if (i - n0 >= 0 && i - n0 < nsamp_) 
	ma_temp[i] = wf_->at(i - n0);
      else if 
	(i - n0 < 0) ma_temp[i] = wf_->at(0); 
      else 
	ma_temp[i] = wf_->at(nsamp_ - 1);
   }
    
  ma_->at(0) = 0;
  for (int k = 0; k < step; k++) 
    ma_->at(0) += ma_temp[k]; 
  ma_->at(0) /= step;
  for (int i = 1; i < nsamp_; i++) 
    ma_->at(i) = ma_->at(i - 1) + 
      (ma_temp[i + step] - ma_temp[i])/step; 
  
  return ma_;    
}

//==============================================================
// Cumulative Sum
vector<double>* RDFiltering::cumsum() 
{
  sc_->at(0) = wf_->at(0);
  for (int i = 1; i < nsamp_; i++ ) 
    {
      sc_->at(i) = sc_->at(i - 1) + wf_->at(i);
    }
   
   return sc_;
}

//==============================================================
// High pass filter
vector<double>* RDFiltering::highpass(int dt, double RC) 
{   
  hp_->at(0) = wf_->at(0);

  double p = RC/(RC + dt);
   
   for (int i = 1; i < nsamp_; i++) 
     hp_->at(i) = p*(hp_->at(i - 1) + wf_->at(i) - wf_->at(i - 1));
   return hp_;    
}

//==============================================================
// SNIP backgrond detection 
vector<double>* RDFiltering::snip(int iter, int av) 
{   
  // moving avg with padding 
  vector<double> sn_temp(nsamp_+av);
   
   int n0 = (av/2);

   for (int i = 0; i < nsamp_ + av; i++) 
     {
       if (i - n0 >= 0 && i - n0 < nsamp_) 
	 sn_temp[i] = wf_->at(i - n0);
      else if (i - n0 < 0) 
	sn_temp[i] = wf_->at(0); 
      else 
	sn_temp[i] = wf_->at(nsamp_ - 1);
   }
    
   sn_->at(0) = 0;
   for (int k = 0; k < av; k++) 
     sn_->at(0) += sn_temp[k]; 
   sn_->at(0) /= av;
   for (int i = 1; i < nsamp_; i++) 
     sn_->at(i) = sn_->at(i - 1) + (sn_temp[i + av] - sn_temp[i])/av; 

   //snip
   vector<double> y0(nsamp_);

   for (int k = 0; k < iter; k = k + av) {     
      for (int i = 1; i < nsamp_; i++) {
	double m1 = sn_->at(i); 
	double m2; 
	if (i - k < 0) 
	  m2 = 0.5*(sn_->at(0) + sn_->at(i + k));
	else if 
	  (i + k > nsamp_ - 1) m2 = 0.5*(sn_->at(i - k) + sn_->at(nsamp_ - 1));
	else 
	  m2 = 0.5*(sn_->at(i - k) + sn_->at(i + k)); 	
         y0[i] = max(m1, m2);         
      }
   }      

   for (int i = 1; i < nsamp_; i++) 
     sn_->at(i) = y0[i];

   return sn_;    
}


