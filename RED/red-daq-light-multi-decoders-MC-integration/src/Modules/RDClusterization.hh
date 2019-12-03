// Clusterization algorithm
#ifndef RDCLUSTERIZATION_H
#define RDCLUSTERIZATION_H

#include <vector>

class RDClusterization {
private:
  std::vector<int> stime_;
  std::vector<int> str_;
  std::vector<int> stp_;
  std::vector<int> rp_; 
  
public:
  RDClusterization(const std::vector<double> *wf, int start, int n_samp, 
		     double thr, int sigma,  int pre_cl, 
		     int post_cl);  
  const std::vector<int>& start_time() const {return stime_;}
  const std::vector<int>& start() const {return str_;}
  const std::vector<int>& stop() const {return stp_;}
  const std::vector<int>& rep() const {return rp_;} 
};
#endif
