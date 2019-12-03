#ifndef EvRaw0_HH
#define EvRaw0_HH
#include <stdlib.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <map>

#include<list>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<stdlib.h>
#include<stdio.h>

#include "EvHeader.hh"


#include "TROOT.h"
#include "TObject.h"
#include "TH1D.h"

using namespace std;

//class EvRaw0;

/* 
  This is the persistent Raw0 Event class
  it should be the event that is dispached to the EvRaw0 observers.
*/

class EvRaw0 : public TObject {


private:

  EvHeader * evheader; 
  map<int,vector<double>* > wf;
       
public:
  EvRaw0();
  ~EvRaw0();
  //Copy constructor: automatic only with C++11
  //EvRaw0(const EvRaw0&);
  
  // set methods
  void SetEvHeader(EvHeader *e)         {evheader = e;};
  void SetWF(vector<double> *wform, int channel);
  
  // get methods
  EvHeader* const GetEvHeader() const   {return evheader;};
  int GetChannelNumber() const {return (int) wf.size();};
  int GetSampleNumber(int channumber=0) const;
  map<int,vector<double>* > const GetWFs() const {return wf;};
  vector<double>* const GetWF(int channumber) const;

  //string GetName() {return Subject::GetName();};
  //bool Notify();

  //void Fill(EvHeader * evh, int n_ch, int n_sam, vector <int> wform);
  void Init();
  void Clear();
  void Dump();
  void FillHeader(int,string,string,string,string,size_t);
  void AddBoardTime(int boardid, size_t time);
  void SetNumberOfBoards(size_t nboards);
  
  TH1D* GimmeHist(int channumber);
  void Draw(int channumber, Option_t* option = "") { GimmeHist(channumber)->Draw(option); };

  ClassDef( EvRaw0, 2 );

};

#endif
