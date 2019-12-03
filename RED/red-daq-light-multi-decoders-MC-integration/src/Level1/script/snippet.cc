#include <iostream>
#include <cstring>
#include <map>
#include <vector>

#include "red-daq/Decoder.hh"
#include "red-daq/EvRaw0.hh"
#include "red-daq/RDconfig.h"

int main()
{
  //input information: either runlist or run number (for database)
  std::string listfile = "run_xxx.lst";
  int run = 0;
  
  //Instantiate decoder
  Decoder* theDecoder = new Decoder("");
  theDecoder->SetVerbosity(0); 
  theDecoder->SetAlignTimes(true);

  if (!(theDecoder->OpenFiles(listfile,run)))
    {
      fprintf (stderr,"Unable to open file. Exiting.\n");
      exit(1);
    }

  // Loop over all events in files  
  for(int iloop=0; !(theDecoder->IsOver());iloop++)
    { 
      //Read next event
      EvRaw0* theEvent = theDecoder->ReadEvent();      
      if (theDecoder->IsOver())
	continue;
      
      //infos here!
      int evNumber = theEvent->GetEvHeader()->GetEventNum();   
      //all waveforms here
      map<int,vector<double>* > wfs = theEvent->GetWFs();

      //Loop over waveforms
      for (map<int,vector<double>* >::iterator it=wfs.begin(); 
	   it!=wfs.end(); ++it)
	{
	  // channel number
	  int k = it->first;
	  //waveform
	  vector<double>  *wf = it->second;    
	  
	  cout << "Wf for channel: " << k << " has " << 
	    wf->size() << " samples " << endl;
	}
    }
  theDecoder->CloseFiles(false);
  

  return 0;
}
