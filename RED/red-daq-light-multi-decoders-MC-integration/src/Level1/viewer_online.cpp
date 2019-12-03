#include <iostream>
#include <cstring>
#include <map>
#include <vector>

#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TTree.h> 
#include <TF1.h>
#include <TApplication.h>
#include <TTimer.h>
#include <TSystem.h>

#include "getopt.h"

#include "Decoder.hh"
#include "RDconfig.h"
#include "EvRaw0.hh"
#include "RDconfig.h"
#include "RDBaseline.hh"
#include "RDFind_min.hh"
#include "RDIntegral.hh"
#include "RDFiltering.hh"
#include "RDClusterization.hh"


using namespace std;

void viewer(int run, string listfile,int customnsamples=-1) {  
   
   TApplication* app = new TApplication("App",0,0);
   int ev = 0; 
   std::string datadir = "";
   Decoder* theDecoder = new Decoder(datadir);
   if (customnsamples>0)
     theDecoder->SetCustomNSamples(customnsamples);
   if (!(theDecoder->OpenFiles(listfile,run)))
    {
      fprintf (stderr,"Unable to open file. Exiting.\n");
      exit(1);
    }

   double scale = 0.75;
   size_t nBoards = theDecoder->GetNumberOfBoards();
   std::cout << "Found " << nBoards << " boards in this run" << std::endl;
   std::vector<TCanvas*> c_top;
   for (size_t ic=0;ic<nBoards; ic++)
     {
       c_top.push_back(new TCanvas(Form("Display%d",ic), "Display", 1000*scale, 1000*scale));
       c_top.at(ic)->Divide(4,4); //16 per canvas
       c_top.at(ic)->SetBatch(kFALSE);
     }
     
   TTimer *timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE); 
   
   std::map<int,TGraph*> g;
   bool firsttime=true;
  
   Int_t evNumber = 0;
 
   //cout << "Sono qui: "<< ev << " " << sch << endl;
   while (ev > -1) {     
     // plot waveform with basic reconstruction
   
      EvRaw0* evRaw = theDecoder->ReadEvent();    
      if (theDecoder->IsOver()) //File is over: rewind!
	{
	  if (evNumber < ev)
	    {
	      cout << "The event #" << ev << " does not exist" << endl;
	      exit(1);
	    }
	  theDecoder->CloseFiles(false);
	  theDecoder->OpenFiles(listfile,run);	  
	  continue;
	}
      evNumber = evRaw->GetEvHeader()->GetEventNum();       
      if (evNumber < ev)
	continue;
      else if (evNumber > ev) //rewind
	{
	  delete theDecoder;
	  theDecoder = new Decoder(datadir);
	  //theDecoder->CloseFiles(false);
	  theDecoder->OpenFiles(listfile,run);
	  continue;
	}
      if (!(ev%10))
	cout << "Display event #" << ev << endl;
      
      //all waveforms here
      map<int,vector<double>* > wfs = evRaw->GetWFs();
      size_t n_samp = evRaw->GetSampleNumber(wfs.begin()->first);     

      Int_t n_channels = evRaw->GetChannelNumber();
     
      for (map<int,vector<double>* >::iterator it=wfs.begin(); 
	   it!=wfs.end(); ++it)
	{
	  int k = it->first;
	  vector<double>  *wf;
	  if (wfs.count(k))
	    wf = wfs.find(k)->second;
	  else 
	    continue;
	  //if g missing, create
	  if (!(g.count(k)))
	    {
	      g.insert(std::pair<int,TGraph*>(k,new TGraph));
	      g.find(k)->second->SetLineColor(kBlue); 
	      g.find(k)->second->SetLineWidth(2);
	    }
	  
	  for (int i = 0; i < n_samp; i++) 
	    g.find(k)->second->SetPoint(i, i, wf->at(i));
	  g.find(k)->second->SetTitle(Form("Event %d - Ch %d;samples;ADC", ev, k));    
	  
	  Int_t index1 = k/16;
	  Int_t index2 = k%16;
	 
	  c_top.at(index1)->cd(1 + index2); 
	  g.find(k)->second->Draw("AL");   

	}

      // draw canvas
      for (size_t ij=0;ij<c_top.size();ij++)
	{
	  c_top.at(ij)->Modified();
	  c_top.at(ij)->Update();
	}

      timer->TurnOn(); // for intective modification of the canvas
      timer->Reset();
             
      //sleep(5);
      timer->TurnOff();
      
      sleep(0.7);
      ev++;
   }

   //Cleanup
   for (std::map<int,TGraph*>::iterator it=g.begin(); it!=g.end(); ++it)
     {
       delete it->second;
     }
  
   delete theDecoder;
   delete app;
   return;
}


int main(int argc, char *argv[]) {

  int c;
  std::string listfile = "";
  int customnsamples = -1;
  int runnr = 0; 
   static struct option long_options[] =
        {            
          {"run",required_argument,0,'r'},
          {"list",required_argument,0,'l'},
	  {"nsamples",required_argument,0,'a'},
          {"help",no_argument,0,'h'},	  
          {0, 0, 0, 0}
        };

  // getopt_long stores the option index here. 
  int option_index = 0;

   // Parse options
  while ((c = getopt_long (argc, argv, "r:l:a:h",long_options,
			   &option_index)) != -1) {
    switch (c)
      {	
      case 'r':
	if (runnr!=0) {
          fprintf (stderr, "Error while processing option '-r'. Multiple runs specified.\n");
          exit(1);
	}
        if ( sscanf(optarg,"%d",&runnr) != 1 ) {
          fprintf (stderr, "Error while processing option '-r'. Wrong parameter '%s'.\n", optarg);
          exit(1);
        }
        if (runnr<0) {
          fprintf (stderr, "Error while processing option '-r'. Run number set to %d (must be >=0).\n", runnr);
          exit(1);
        }
        fprintf(stdout,"Merging files from run %d\n",runnr);
        break;      
      case 'l':
        listfile = optarg;
        fprintf(stdout,"Data will be read from files listed in '%s'\n",listfile.c_str());
        break;     
      case 'a':
	customnsamples = atoi(optarg);
	if (customnsamples<=0)
	  {
	    fprintf (stderr, "Error while processing option '--nsamples'. Wrong parameter '%s'.\n", optarg);
	    exit(1);
	  }	
        fprintf(stdout,"Set custom number of samples %d\n",customnsamples);
        break;
      case 'h':
        fprintf(stdout,"\nviewer ([-r run_number]|[-l list_file]) [-h]\n\n");
        fprintf(stdout,"  -r: define run to process\n");
        fprintf(stdout,"  -l: define file with list of data files to process\n");
        fprintf(stdout,"      n.b. either -r or -l must be specified\n");
        fprintf(stdout,"  -h: show this help message and exit\n\n");
        exit(0);
      case '?':
	fprintf (stderr,"Unknown option character `\\x%x'.\n",optopt);
	exit(1);
      default:
        abort();
      }
  }

  // Verify that some input was specified
  if ( listfile.compare("")==0 && runnr==0 ) {
    fprintf (stderr,"No run number and no list file were specified. Exiting.");
    exit(1);
  }
  if ( listfile.compare("")!=0 && runnr!=0 ) {
    fprintf (stderr,"Both run number and list file were specified. Exiting.");
    exit(1);
  }


  //Main call
  viewer(runnr,listfile,customnsamples);
  

  return 0; 
}







