#ifndef ADCBoard_H
#define ADCBoard_H

#include<vector>
#include<string>
#include<iostream>
#include<fstream>

#include "Rtypes.h"

#include "DBService.hh"

#include "ADCEvent.hh"
#include "ADCFile.hh"

class ADCBoard
{

 public:

  ADCBoard(int);
  ~ADCBoard();

  void Reset();

  int GetBoardId() { return fBoardId; }

  ADCFile* File(int i);
  ADCFile* CurrentFile(){return File(fCurrentFile);};
  ADCFile* AddFile(std::string);
  int GetNFiles() { return fFiles.size(); }
  std::string GetFileName(int);
  std::string GetCurrentFileName(){return GetFileName(fCurrentFile);};  

  ADCEvent* Event() { return fADCEvent; }
  ADCEvent* NextEvent();

  void SetVerbose(Int_t v) { fVerbose = v; }
  void ForceNumberOfSamplesTo(int nsamples);
  bool IsNSamplesForced(){return fCustomNSamples;};
  int GetForcedNSamples(){return fNSamples;};

 private:

  int ReadFileHead();
  int ReadNextEvent();
  int UnpackEvent(UInt_t);
  int UnpackFileHead();
  int UnpackFileTail();

 private:

  int fBoardId;

  Int_t fVerbose;

  //DBService* fDB;

  std::vector<ADCFile> fFiles;

  UInt_t fCurrentFile;
  std::ifstream fFileHandle;
  void* fBuffer;

  ADCEvent* fADCEvent;
  int fNSamples;
  bool fCustomNSamples;
  

};
#endif
