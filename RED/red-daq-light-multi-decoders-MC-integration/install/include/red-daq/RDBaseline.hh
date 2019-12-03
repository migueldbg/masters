// Baseline module
#ifndef RDBASELINE_H
#define RDBASELINE_H

#include <math.h> 
#include <vector>

using namespace std;

class RDBaseline {
 private:
  double q_, m_, b_, s_;
 public:
  RDBaseline();
  ~RDBaseline(){;};

  void DoIt(const vector<double> *wf, int pre_min, int pre_max, 
	    int pos_min, int pos_max);
  std::vector<double> GetBaseline(const vector<double> *wf, int pre_min, 
           int pre_max,int pos_min, int pos_max);

  void DoIt(const vector<double> *wf, int pre_min, int pre_max);

  double mean() {return b_;}
  double rms() {return s_;}        
  double q() {return q_;}        
  double m() {return m_;}        

  std::vector<double> fBaseline;
};


#endif
