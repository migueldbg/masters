#include "EvRaw0.hh"
#include <iostream>
#include <map>

using namespace std;

//constructor
EvRaw0::EvRaw0()
{
  evheader = new EvHeader;
  //cout << "EvRaw0::EvRaw0()" << endl; 
}

//destructor
EvRaw0::~EvRaw0()
{
  if (evheader)
    delete evheader;
}

/*
bool EvRaw0::Notify() 
{   
  bool myb;
  return myb;
};
*/

void EvRaw0::Init(){};


void EvRaw0::Clear()
{
  //Clean up Header
  if (evheader)
    delete evheader;
  evheader = new EvHeader;

  //Clear wfs
  
  for (map<int,vector<double>* >::iterator it=wf.begin(); it!=wf.end(); ++it)
    {
      delete it->second;
    }
  
  wf.clear();
 
};

void EvRaw0::SetWF(vector<double> *wform, int channel)  
{
  //Already exists: replace it
  if ( wf.count(channel))
    {
      delete wf.find(channel)->second;
      wf.find(channel)->second = wform;
    }      
  else //add it
    wf.insert(std::pair<int,vector<double>* >(channel,wform) );
  return;
}

int EvRaw0::GetSampleNumber(int channel) const
{  
  return (wf.count(channel)) ? wf.find(channel)->second->size() : -1;
}

vector<double>* const  EvRaw0::GetWF(int channel) const
{
  return (wf.count(channel)) ? wf.find(channel)->second : 0;
}


void EvRaw0::Dump() 
{  
  for (map<int,vector<double>* >::iterator it=wf.begin(); it!=wf.end(); ++it)
    {
      vector<double>* theWF = it->second;
      for (int t=0; t<theWF->size(); t++) 
	{	  
	  cout <<" wf" << it->first << "[" << t << "]=" << theWF->at(t) 
	       << endl;
	
      }
    }
}

void EvRaw0::SetNumberOfBoards(size_t nboards)
{
  evheader->SetNumberOfBoards(nboards);
}


void EvRaw0::AddBoardTime(int boardid, 
			  size_t time)
{
  evheader->SetBoardTime(boardid,time);
}

void EvRaw0::FillHeader(int eventnumber,
			string date,
			string runtype,
			string recver,
			string inputfile,
			size_t time)
{
  evheader->SetEventNum(eventnumber);
  evheader->SetDate(date);
  evheader->SetRunType(runtype);
  evheader->SetRecVer(recver);
  evheader->SetInputFile(inputfile);
  evheader->SetTime(time);
}



TH1D* EvRaw0::GimmeHist(int channumber)
{
   vector<double>* wf = GetWF(channumber);	
   TH1D* hh = new TH1D("histo","hh",wf->size(),0,wf->size()-1);
   for (size_t i=0;i<wf->size();i++)
     hh->SetBinContent(i+1,wf->at(i));   
  return hh;
}

ClassImp ( EvRaw0 )



