#include "RDIntegral.hh"

using namespace std;

RDIntegral::RDIntegral() : 
  e_(0.),f90_(0.),mt_(0.),rmst_(0.)
{;}

//======================================================================

void RDIntegral::DoIt(const vector<double> *wf, 
		     double m, double q, 
		     int min, int max, int max90) 
{
   e_ = 0;
   f90_ = 0;
   mt_ = 0;
   rmst_ = 0;
   double mt2 = 0;

   for (int i = min; i < max; i++) {
     double ee = m*i + q - wf->at(i);
      e_ += ee;
      mt_ += (i - min)*ee;
      mt2 += (i - min)*(i - min)*ee;
   
      if (i < max90) f90_ += m*i + q - wf->at(i);      
   }

   if(e_) { 
      mt_ =  mt_/e_;
      mt2 = mt2/e_;
      rmst_= mt2 - mt_*mt_;
   }
   else {
      mt_ = 0;
      rmst_ = 0;
   }
   

   if(e_) f90_ = f90_/e_;
   else f90_ = -100;   

}

