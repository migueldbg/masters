// Filtering Module
#ifndef RDFILTERING_H
#define RDFILTERING_H

#include <algorithm>
#include <vector>

class RDFiltering {  
  int nsamp_;
  std::vector<double> *wf_;
  std::vector<double> *ma_;
  std::vector<double> *sc_;
  std::vector<double> *hp_;
  std::vector<double> *sn_;

public:
  RDFiltering(const std::vector<double> *wf, int nsamp);
  ~RDFiltering(); 

  //Filtering methods
  std::vector<double> *mavg(int step);
  std::vector<double> *cumsum();
  std::vector<double> *highpass(int dt, double RC);
  std::vector<double> *snip(int iter, int av);  
  
};


#endif
