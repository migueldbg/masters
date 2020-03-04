#ifndef EvHeader_H
#define EvHeader_H
#include <stdlib.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>


#include<list>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<stdlib.h>
#include<stdio.h>


#include "TROOT.h"
#include "TObject.h"

using namespace std;


/* 
  This is the persistent Raw1 Event class
  it should be the event that is dispached to the EvRaw0 observers.
*/

// class EvRaw0 : public TObject {
class EvHeader : public TObject {
    
protected:
  Int_t   event_num;     // event number
  string date;
  string run_type;
  string rec_ver;
  string input_file;
  Int_t run_number;
  //Int_t n_samples;
  Int_t n_channels;
  Int_t peak_buffer;
  UInt_t time;
  vector<UInt_t> boardtimes;
  
  
public:
  
  EvHeader();
  ~EvHeader(){};
  
  // set method
  void SetEventNum(int en)          { event_num = en; };
  void SetDate(string d)            { date = d; };
  void SetRunType(string rt)        { run_type = rt; };
  void SetRecVer(string rv)         { rec_ver = rv; };
  void SetInputFile(string ifi)     { input_file = ifi; };
  void SetRunNumber(Int_t num)      { run_number = num; };
  //void SetNSamples(Int_t num)       { n_samples = num; };
  void SetNChannels(Int_t num)      { n_channels = num; };
  void SetPeakBuffer(Int_t num)     {peak_buffer = num;};
  void SetTime(UInt_t num) {time = num;};
  void SetNumberOfBoards(size_t nboards);
  void SetBoardTime(int id,size_t time);

  // get method
  Int_t GetEventNum()       { return event_num; };
  string GetDate()          { return date; };
  string GetRunType()       { return run_type; };
  string GetRecVer()        { return rec_ver; };
  string GetInputFile()     { return input_file; };
  Int_t GetRunNumber()      { return run_number; };
  //Int_t GetNSamples()       { return n_samples; };
  Int_t GetNChannels()      { return n_channels; };
  Int_t GetPeakBuffer()     { return peak_buffer;};
  UInt_t GetTime()          { return time;};
  UInt_t GetBoardTime(int id)          { return boardtimes.at(id);};
  
  void Fill(string date, string run_type, string rec_ver, string input_file, int run_num, int n_sam,
	    int n_ch, int p_buff, int ev_num);
    
    
    
  // Dump method
  
  void Dump();
  
  ClassDef( EvHeader, 3 )
    
};

#endif
