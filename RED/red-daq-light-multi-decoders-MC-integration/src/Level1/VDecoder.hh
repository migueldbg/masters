#ifndef _VDecoder_hh_
#define _VDecoder_hh_

#include <string>
#include <iostream>
#include "EvRaw0.hh"
#include "RunHeader.hh"

class VDecoder
{  
public:
  VDecoder();   
  ~VDecoder();

  //file management
  virtual bool OpenFiles(std::string,int) = 0;
  virtual bool CloseFiles(bool) = 0;
  //Main processing method
  virtual EvRaw0* const ReadEvent() = 0; 
  //header info, event
  virtual bool ReadEventHeader();
  virtual bool ReadRunHeader();
  virtual bool ReadRunTrailer();

  //
  EvRaw0* const GetRawEvent() const;
  //
  RunHeader* GetRunHeader();

  void SetVerbosity(int vl){fVerbosity = vl;};
  bool IsOver(){return fIsOver;};
  int GetEventNumber(){return fEventNumber;};
  void SetAlignTimes(bool val){fAlignTimes = val;};
 
  void SetEventsToSkip(int ntoskip){fEventsToSkip = ntoskip;};
  void SetCustomNSamples(int ns);
  void ReleaseCustomNSamples(){fCustomNSamples = false;};

protected:
  EvRaw0* fRawEvent;
  RunHeader* fRunHeader;
  
  int fVerbosity;
  std::string fDataDir;

  bool fIsOver;
  int fEventNumber;
  int fRunNumber;
  bool fAlignTimes; 
  int fEventsToSkip;
  bool fEventSkipped;

  bool fCustomNSamples;
  int fNSamples;

};

#endif
