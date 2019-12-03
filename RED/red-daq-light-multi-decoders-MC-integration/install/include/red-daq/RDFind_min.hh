// Minimum finder module
#ifndef RDFIND_MIN_H
#define RDFIND_MIN_H

#include <vector>

class RDFind_min {
 private:
  double x_, y_, t_, ft_;
  double x1,y1;
 public:
  RDFind_min();
  ~ RDFind_min(){;};

  void DoIt(const std::vector<double> *wf, 
	    int min, int max, double b, double cdf, double xf);
  double x_min() {return x_;}
  double y_min() {return y_;}
  double s_time() {return t_;}
  double f_time() {return ft_;}     
  double x_max() {return x1;};
  double y_max() {return y1;};    
};

#endif
