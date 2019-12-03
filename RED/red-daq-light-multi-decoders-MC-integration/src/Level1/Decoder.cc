#include "Decoder.hh"
#include "DBService.hh"
#include "ADCEventTags.hh"
#include "inttypes.h"
#include "unistd.h"
#include <cstdio>
#include <ctime>
#include <fstream>
#include <cstring>
#include <vector>


using namespace std;

//=============================================================================
Decoder::Decoder(string datadir) : VDecoder(),
				   fDataDir(datadir)
{
  fListFiles = "";
  fBoards = new std::vector<ADCBoard*>;
}

//=============================================================================

Decoder::~Decoder() 
{
  if (fBoards)
    {
      for (size_t i=0;i<fBoards->size();i++)
	delete fBoards->at(i);
      delete fBoards;
    } 
}

//=============================================================================
bool Decoder::OpenFiles(string listfile,int runnr)
{
  ADCBoard* board;
  fRunNumber = runnr;
  fListFiles = listfile;

  // If no list is specified, get board and file lists from DB
  if (fListFiles.compare("")==0) {

    // Get handle to DB
    DBService* db = DBService::GetInstance();

    // Get from DB list of board ids used in current run
    std::vector<int> boardList;
    int rc = db->GetBoardList(boardList,fRunNumber);
    if (rc != DBSERVICE_OK) {
      printf("ERROR retrieving from DB list of files for run %d. Aborting\n",fRunNumber);
      return false;
    }
    if (boardList.size() == 0) {
      printf("ERROR no boards associated to run %d found in DB. Aborting\n",
	     fRunNumber);
      return false;
    }

    // Create vector of boards
    for(unsigned int b=0; b<boardList.size(); b++) 
      {
	printf("Board %d\tReading data from board id %d\n",b,boardList[b]);
	board = new ADCBoard(boardList[b]);
	board->SetVerbose(fVerbosity);
	if (fCustomNSamples)
	  board->ForceNumberOfSamplesTo(fNSamples);

	// Get list of files created for each board during run
	vector<string> fileList;
	rc = db->GetFileList(fileList,fRunNumber,boardList[b]);
	if (rc != DBSERVICE_OK) {
	  printf("ERROR retrieving from DB list of files for run %d board id %d. Aborting\n",
		 runnr,boardList[b]);
	  return false;
	}
	// Add files to board
	for(unsigned int f=0; f<fileList.size(); f++) {
	  string filePath = fDataDir+"/"+fileList[f];
	  printf("\tFile %d - Reading from file %s\n",f,filePath.c_str());
	  board->AddFile(filePath);
	}

	// Add board to boards vector
	fBoards->push_back(board);      
      }
  } 
  else 
    {
      // Get list of boards and files from file
      ifstream list;
      string line;
      int bid;
      char bfile[1024];
      list.open(fListFiles.c_str());
      if (!list.good())
	{
	  cout << "Problem in opening the file " << fListFiles << endl;
	  return false;
	}
      while(!list.eof())
	{
	  getline(list,line);
	  if (line.compare("")!=0) 
	    {
	      if ( sscanf(line.c_str(),"%d %s",&bid,bfile) != 2 ) {
		printf("ERROR while parsing list '%s'. Aborting\n",
		       fListFiles.c_str());
		return false;
	      }

	      // Find board with correct board id or add it to list if new
	      std::vector<ADCBoard*>::iterator it;
	      for (it = fBoards->begin(); it != fBoards->end(); ++it) {
		board = *it;
		if (board->GetBoardId() == bid) break;
	      }
	      if (it == fBoards->end()) {
		printf("Board id %d\n",bid);
		board = new ADCBoard(bid);
		board->SetVerbose(fVerbosity);
		if (fCustomNSamples)
		  board->ForceNumberOfSamplesTo(fNSamples);
		fBoards->push_back(board);	      
	      }

	      // Add file to board
	      std::string filePath = fDataDir+"/"+bfile;
	      printf("\tBoard %d - File %s\n",bid,filePath.c_str());
	      board->AddFile(filePath);
	    }
	}
      list.close();      
    }

  // Show list of known boards/files
  printf("Reading %d board(s)\n",(int)fBoards->size());
  for (std::vector<ADCBoard*>::iterator it = fBoards->begin(); 
       it != fBoards->end(); ++it) 
    {
      board = *it;
      printf("Board %d Files %d\n",board->GetBoardId(),board->GetNFiles());
      for(int f=0; f<board->GetNFiles(); f++) 	
	printf("File %d %s\n",f,board->GetFileName(f).c_str());     	
    }

  //Create and fill RunHeader
  if (fRunHeader)    
    delete fRunHeader;
  fRunHeader = new RunHeader();
  fRunHeader->SetRunNumber(fRunNumber);
  
  time_t starttime = (time_t) fBoards->at(0)->File(0)->GetStartTime();
  //cout << "Start time: " << starttime << endl;
  std::tm * ptm = std::localtime(&starttime);
  char buffer[32];
  // Format: 15.06.2009 20:20:00
  std::strftime(buffer, 32, "%d.%m.%Y %H:%M:%S", ptm);  
  fRunHeader->SetDate(string(buffer));
  //cout << " Time: " << buffer << endl;

  //Initialize vectors for timing
  TTg.clear();
  TT.clear();
  dT.clear();
  dTB0.clear();
  TTg.resize(fBoards->size(),0);
  TT.resize(fBoards->size(),0);
  dT.resize(fBoards->size(),0);
  dTB0.resize(fBoards->size(),0);

  fEventSkipped = false;

  return true;
}

//=============================================================================
bool Decoder::CloseFiles(bool updatedb)
{
  // If input was from a real run, finalize output and update DB if required
  if (fListFiles.compare("")==0) 
    {
      printf("Run %d closed after writing %d events\n",fRunNumber,fEventNumber);
      if (updatedb == 0) {
	printf("Option -u was not specified: DB not updated.\n");
      }
      else 
	{
	  // Get handle to DB
	  DBService* db = DBService::GetInstance();

	  // Update DB with number of events for this run
	  int n_events;
	  int rc = db->GetRunEvents(n_events,fRunNumber);
	  if (rc != DBSERVICE_OK) {
	    printf("ERROR retrieving from DB old number of events for run %d. Aborting\n",fRunNumber);
	    return false;
	  }
	  if (n_events == fEventNumber) {
	    printf("Run %d Events %d. DB is up-to-date: no action.\n",
		   fRunNumber,fEventNumber);
	  } else {
	    if (n_events == 0) {
	      printf("Run %d Events %d. Writing number of events to DB.\n",
		     fRunNumber,fEventNumber);
	    } else {
	      printf("Run %d Events %d. WARNING - DB returns %d events: updating DB.\n",
		     fRunNumber,fEventNumber,n_events);
	    }
	    int rc = db->UpdateRunEvents(n_events,fRunNumber);
	    if (rc != DBSERVICE_OK) {
	      printf("ERROR updating DB for run %d. Aborting\n",fRunNumber);
	      return false;
	    }
	  }
	  
	}
    } 
  else 
    {  
      printf("A total of %d merged events were written to the output file\n",fEventNumber);
      sleep(2);      
    }

  //close files here
  return true;
}

//=============================================================================
EvRaw0* const  Decoder::ReadEvent() 
{ 
  fRawEvent->Clear();

  // Load next event for all boards
  UInt_t nEOR = 0;

  std::vector<int> eventToSkip(fBoards->size(),0);
  //Check if there are any events to skip
  if (fEventsToSkip>0 && !fEventSkipped)
    for (size_t b=1;b<eventToSkip.size();b++)
      {
	eventToSkip[b] = fEventsToSkip;  
	cout << "Skip " << eventToSkip[b] << " events on board " << b << endl;
      }

  if (!fEventSkipped && fEventsToSkip) {   
      // Skip first events for selected boards
      for(unsigned int b=0; b<fBoards->size(); b++) {
        while (eventToSkip[b] > 0) {
          if ( fBoards->at(b)->NextEvent() == 0 ) {
            nEOR++;
            printf("*** Board %2d - End of Run.\n",
		   fBoards->at(b)->GetBoardId());
            break;
          }
          eventToSkip[b]-- ;
        }
      }
      fEventSkipped = true;
    }

  for(unsigned int b=0; b<fBoards->size(); b++) {
    if ( fBoards->at(b)->NextEvent() == 0 ) {
      nEOR++;
      printf("*** Board %2d - End of Run.\n",fBoards->at(b)->GetBoardId());
    }
  }
  if (nEOR != 0) 
    {
      if (nEOR == fBoards->size()) 
	printf("All boards reached end of run.\n");
      else 
	printf("WARNING: only %d board(s) reached end of run.\n",nEOR);  
      fIsOver = true;
      return fRawEvent;
    }
  
  nEOR = 0;
  UInt_t in_time = 0;
  while ( (! in_time) && (nEOR==0) ) 
    {    
      // Get timing information
      UInt_t bmax = 0;
      for(unsigned int b=0; b<fBoards->size(); b++) {
	TT[b] = 
	  fBoards->at(b)->Event()->GetEventTimeTag() & 0x7FFFFFFF;
	/*
	  cout << "Reading timestamp of event " << 
	  fBoards->at(b)->Event()->GetEventCounter() <<
	  " from board: " << b << "--> " << 
	  fBoards->at(b)->Event()->GetEventTimeTag() << endl;
	*/
	if(!fEventNumber) 
	  {
	    //TT0[b] = TT[b];
	    TTg[b] = TT[b]; // Event 0 is good by definition
	  }
	dT[b] = TT[b]-TTg[b];
	if (dT[b]<0) dT[b] += (1<<31);
	if (dT[b]>dT[bmax]) bmax = b;
	dTB0[b] = TT[b]-TT[0];
	
	if (fVerbosity>=2) {
	  if (b==0) {
	    printf("- Board %2d n. %2d NEv %8u Tabs %u Dt %u \n",
		   fBoards->at(b)->GetBoardId(),b,
		   fBoards->at(b)->Event()->GetEventCounter(),TT[b],dT[b]);
	    /*
	      printf("- Board %2d NEv %8u Tabs %f (0x%08x) Dt %f (0x%08x)\n",b,
	      fBoards->at(b)->Event()->GetEventCounter(),
	      TT[b]*8.5E-9,TT[b],dT[b]*8.5E-9,dT[b]);
	    */
	  } else {
	    printf("- Board %2d n. %2d NEv %8u Tabs %u Dt %u DtB0 %u \n",
		   fBoards->at(b)->GetBoardId(),b,
		   fBoards->at(b)->Event()->GetEventCounter(),TT[b],dT[b],dTB0[b]);

	    /*
	    printf("- Board %2d NEv %8u Tabs %f (0x%08x) Dt %f (0x%08x) DtB0 %f (0x%08x)\n",b,
		   fBoards->at(b)->Event()->GetEventCounter(),
		   TT[b]*8.5E-9,TT[b],dT[b]*8.5E-9,dT[b],dTB0[b]*8.5E-9,dTB0[b]);
	    */
	  }
	}
      }
      
      // Verify if all events are in time (0x100 clock cycles tolerance).
      // Tolerance is now 2 us (0x100= 256 clocks)
      in_time = 1;
      for(unsigned int b=0; b<fBoards->size(); b++) {
	if ( (dT[bmax]-dT[b]) >= 0x100 && fAlignTimes) {
	  in_time = 0;
	  printf("*** Board %2d - Skipping event %8u\n",
		 fBoards->at(b)->GetBoardId(),
		 fBoards->at(b)->Event()->GetEventCounter());
	  if ( fBoards->at(b)->NextEvent() == 0 ) {
	    printf("*** Board %2d - End of Run.\n",fBoards->at(b)->GetBoardId());
	    nEOR++;
	  }
	}
      }
    }

    // If one or more files reached EOR, stop processing events
    if (nEOR != 0) {
      if (nEOR == fBoards->size()) {
	printf("All boards reached end of run.\n");
      } else {
	printf("WARNING: only %d board(s) reached end of run.\n",nEOR);
      }
      fIsOver = true;
      return fRawEvent;
    }

    // Tag current set of times ad "good"
    //Information kept in the header
    fRawEvent->SetNumberOfBoards(fBoards->size());

    // Tag current set of times ad "good"
    for (unsigned int b=0; b<fBoards->size(); b++) 
      {
	TTg[b] = TT[b];
	fRawEvent->AddBoardTime(b,TTg[b]);
      }

    int n_channels=0;

    //Readout channels
    for (size_t b=0; b<fBoards->size(); b++) 
      {
	for(int i=0; i<16; i++)
	  {    
	    int chanID = b*16+i;
	    if(fBoards->at(b)->Event()->GetAcceptedChannelMask() & (0x1<<i))
	      {
		int n_samp = (fBoards->at(b)->IsNSamplesForced()) ? 
		  fBoards->at(b)->GetForcedNSamples() : //overridden
		  fBoards->at(b)->Event()->GetAcceptedChannelMask() >> 16 ; //read from the evheader
                if (n_samp > ADCEVENT_NSAMPLES)
                   { 
		     cout << "The waveform should have " << n_samp << " samples, but " << endl;
                     cout << "the maximum dimension allowed in ADCEventTags (ADCEVENT_NSAMPLES) is " 
                         << ADCEVENT_NSAMPLES << endl;
                     cout << "** check it **" << endl;
                     fIsOver = true;   
                     return fRawEvent;
                   }

		n_channels++;     

		//The memory will be released by EvRaw0
		vector<double> *theWF = new vector<double>(n_samp);
		for (size_t is=0;is<theWF->size();is++)		  
		  theWF->at(is) = fBoards->at(b)->Event()->GetADCChannelSample(i,is);
		//the fRawEvent interface will take care of releasing the memory
		fRawEvent->SetWF(theWF,chanID);
			
		if (fVerbosity)
		  {
		    cout << "Channel found: " << n_samp << " " <<
		      fBoards->at(b)->Event()->GetAcceptedChannelMask() << " " << 
		      (0x1<<i) << " " << 
		      chanID << endl;
		  }
	      }
	}
      }
    //fill header
    //ADCFile* theFile = fBoards->at(0)->GetCurrentFile();
    fRawEvent->FillHeader(fEventNumber, //eventnumber
			  "NULL", //date
			  "NULL", //run type
			  "NULL", //Rev Version
			  fBoards->at(0)->GetCurrentFileName(), //input file
			  TT.at(0)); //time stamp
  
			      
    
    //Done
    fEventNumber++;
    return fRawEvent;
}

