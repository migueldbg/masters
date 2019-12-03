#ifndef _Analyzer_hh_
#define _Analyzer_hh_

#include <map>
#include <vector>
#include <string>
#include "EvRec0.hh"

class EvRaw0;
class RDIntegral;
class RDBaseline;
class RDFind_min;
class RDPulseFitter;

class Analyzer
{  
public:
  Analyzer(bool islaser, const std::string& customfile="NULL");   
  ~Analyzer();
  
  //returns the event number (-1 if failed)
  EvRec0* const ProcessEvent(const EvRaw0*); 
  //
  EvRec0* const GetRecEvent() const;

  void SetVerbosity(int val) {fVerbosity=val;};

private:
  EvRec0* fRecEvent;
  
  std::map<std::string, double> cfg; 
  int fVerbosity;
  
  std::string fCustomFile;

  RDIntegral* fIntegral;
  RDBaseline* fBaseline;
  RDFind_min* fFindMin;
  RDPulseFitter* fPFitter;

  bool fFitS2;

};

#endif
