#include "RDFind_min.hh"
#include <iostream> 

RDFind_min::RDFind_min() :
  x_(0.),y_(0.),x1(0.),y1(0.)
{;}

//====================================================================

void RDFind_min::DoIt(const std::vector<double> *wf, int min, int max, double b, double cdf, double xf) 
{     

  x_ = 0;
  y_ = wf->at(0);   
  x1 = 0;
  y1 = wf->at(0); 

  if (min < 0 || max > (wf->size()-1))
    {
      std::cout << "find_min::find_min() invalid range" << std::endl;
      std::cout <<"min = " << min << "; max = " << max << std::endl; 
     return;
    }
  
  
  // find min x and y 
  for (int i = min; i < max; i++) {
    if (wf->at(i) < y_) {
      x_ = i;  
      y_ = wf->at(i);
    }
    if (wf->at(i) > y1) {
      x1 = i;
      y1 = wf->at(i);
    }
  }
  
  // find start time with CFD
  t_ = min;
  for(int i = min; i < max; i++) {
    if (wf->at(i) < y_ - (y_ - b)*cdf) {
      t_ = i;
      if (wf->at(i-1) > y_ - (y_ - b)*cdf) 
	{ //interpolate t_ b/w the two samples to improve resol.
	  t_ = i-1 +
            ((y_ - (y_ - b)*cdf) - wf->at(i-1)) / (wf->at(i)-wf->at(i-1));
	}
      break;
    }
  }

  // find start time with fixed threshold
  ft_ = min;   
  for(int i = min; i < max; i++) {
       if (wf->at(i) - b < xf) {
          ft_ = i; 
          break;
       }
  }  

}


