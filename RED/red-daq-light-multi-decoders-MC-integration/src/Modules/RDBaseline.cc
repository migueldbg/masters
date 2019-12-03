#include "RDBaseline.hh"

RDBaseline::RDBaseline() : 
  q_(0),m_(0),b_(0),s_(0)
{;}

//=============================================================================
void RDBaseline::DoIt(const vector<double> *wf, int pre_min, int pre_max, 
		   int pos_min, int pos_max) 
{ // interpolation baseline
  b_ = 0; 
  s_ = 0; 
  int dn (0);
  
  for (int i = pre_min; i < pre_max; i++) 
    {
      double val =  wf->at(i);
      b_ += val;
      s_ += val*val;
      dn++;
    }
  
  int n = dn; 
  double y0 = b_;
  
  for (int i = pos_min; i < pos_max; i++) 
    {
      double val =  wf->at(i);
      b_ += val;
      s_ += val*val;
      dn++;
    }
  
  double y1 = b_ - y0;
  
  b_ = b_/dn;
  
  s_ = s_/dn;
  s_ = sqrt(s_ - b_*b_);
  
  y0 = y0/n;
  y1 = y1/(dn - n);
  
  double x0 = 0.5*(pre_max + pre_min); 
  double x1 = 0.5*(pos_max + pos_min);
  
  m_ = (y1 - y0)/(x1 - x0);
  q_ = y0 - m_*x0;    
}

//=============================================================================

void RDBaseline::DoIt(const vector<double> *wf, int pre_min, int pre_max) 
{ 
  // std baseline
  b_ = 0; 
  s_ = 0; 
  int n (0);
    
  for (int i = pre_min; i < pre_max; i++) 
    {
      double val = wf->at(i);
      b_ += val;
      s_ += val*val;
      n++;
  }

  b_ = b_/n;
  
  s_ = s_/n;
  s_ = sqrt(s_ - b_*b_);
  
  m_ = 0;
  q_ = b_;  
}

//=============================================================================
vector<double> RDBaseline::GetBaseline(const vector<double> *wf, int pre_min, 
           int pre_max,int pos_min, int pos_max)
{
  //First (re)calculate parameters
  DoIt(wf,pre_min,pre_max,pos_min,pos_max);

  //Check how many samples. If fBaseline do nothing: vector components 
  //are simply overwritten
  if (wf->size() != fBaseline.size())
    {
      fBaseline.clear();
      fBaseline.resize(wf->size(),0);
    }    

  //Fill it. At the moment, using linear model
  for (size_t i=0;i<wf->size();i++)
    fBaseline.at(i) = i*m_ + q_;

  return fBaseline;
}


