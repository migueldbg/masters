// Charge and f90 reconstruction
#ifndef RDINTEGRAL_H
#define RDINTEGRAL_H

#include <vector>

class RDIntegral {
  double e_, f90_;
  double mt_, rmst_;
 public:
  RDIntegral();
  ~RDIntegral(){;};

  void DoIt(const std::vector<double> *wf, double m, 
	   double q, int min, int max, int max90);
  double charge() {return e_;}
  double f90() {return f90_;}  
  double mean_time() {return mt_;}      
  double rms_time() {return rmst_;}      
};

#endif
