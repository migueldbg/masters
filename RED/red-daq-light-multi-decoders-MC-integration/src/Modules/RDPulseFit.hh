#ifndef _RDPulseFit_
#define _RDPulseFit_

#include "TObject.h"

#include <vector>

using namespace std;

class RDPulseFit : public TObject
{ 
 public:
  RDPulseFit();
  ~RDPulseFit();
  RDPulseFit(const RDPulseFit&); //copy constructor
  const RDPulseFit& operator=(const RDPulseFit&); //assignment
  int operator==(const RDPulseFit&) const; //equality
 
  //data members
  int status,covstatus,ndf;
  double chi2;
  std::vector<double> par,epar;
  int start,end; 
  int sipm;      // [-1] All [isipm]
  int type;      // [0] SPE fit [1] S2 fit

  ClassDef(RDPulseFit, 10);
  
};

#endif
