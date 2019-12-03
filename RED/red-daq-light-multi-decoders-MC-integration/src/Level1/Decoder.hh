#ifndef _Decoder_hh_
#define _Decoder_hh_

#include <string>
#include <iostream>
#include "EvRaw0.hh"
#include "RunHeader.hh"
#include "VDecoder.hh"

#include "ADCBoard.hh"

class Decoder : public VDecoder
{  
public:
  Decoder(std::string datadir);   
  ~Decoder();

  //file management
  // Implement virtual interface.
  bool OpenFiles(std::string,int);
  bool CloseFiles(bool);
  EvRaw0* const ReadEvent(); 

  //Inherited methods
  /*
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
  */

  size_t GetNumberOfBoards(){return fBoards->size();};

private:
  std::vector<ADCBoard*> *fBoards;
  std::string fDataDir;
  std::string fListFiles;
  /* //Inherited members 
      EvRaw0* fRawEvent;
      RunHeader* fRunHeader;
      
      int fVerbosity;
      
      bool fIsOver;
      int fEventNumber;
      int fRunNumber;
      bool fAlignTimes; 
      int fEventsToSkip;
      bool fEventSkipped;
      
      bool fCustomNSamples;
      int fNSamples;
  */


  //Time stamps
  vector<UInt_t> TTg; // Time tag of last good event
  vector<UInt_t> TT;  // Time tag of current event
  vector<Int_t> dT;   // dT wrt last good event
  vector<Int_t> dTB0; // dT wrt board 0 in current event

  std::vector<ADCBoard*> * const GetBoards() const
  {return fBoards;};
};

#endif
