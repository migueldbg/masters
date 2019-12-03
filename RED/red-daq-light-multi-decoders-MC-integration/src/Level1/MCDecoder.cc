#include "MCDecoder.hh"
#include "MCTraces.hh"

#include <inttypes.h>
#include <unistd.h>

#include <cstdio>
#include <ctime>
#include <fstream>
#include <cstring>
#include <vector>

using namespace std;

//=============================================================================
MCDecoder::MCDecoder(std::string configfile, int evts) 
    : VDecoder(), fMaxEvents(evts)
{
    fTraces = new MCTraces(configfile);
}

//=============================================================================

MCDecoder::~MCDecoder() 
{
    delete fTraces;
}


//=============================================================================
bool MCDecoder::OpenFiles(string listfile, int runnr)
{
    //Open the file and do any other operation required at the beginning
    //of the file
    
    //Create and fill RunHeader
    if (fRunHeader)
        delete fRunHeader;
    fRunHeader = new RunHeader();
    fRunHeader->SetRunNumber(fRunNumber);

    return true;
}

//=============================================================================
bool MCDecoder::CloseFiles(bool updatedb)
{
    //Do any operation required at the closeout of file

    //close files here
    return true;
}

//=============================================================================
EvRaw0* const  MCDecoder::ReadEvent() 
{ 
    //Clear and re-fill the fRawEvent and fTraces objects
    fRawEvent->Clear();
    fTraces->Reset();    
    
    // Get PE times
    fTraces->GetEvent();
    
    // Add PE times to trace
    for (int icl=0; icl<fTraces->GetNClusters(); icl++)
        fTraces->AddToTrace(icl);
    
    // Add fluctuations and load into EvRaw0
    for (int isipm=0; isipm<28; isipm++) 
    {
        fTraces->AddFluctuations(isipm);
        auto trace = fTraces->GetTrace(isipm);
        int  nsa   = trace.size();
        auto theWF = new vector<double>(nsa);
        for (int sa=0; sa<nsa; sa++)
            theWF->at(sa) = -trace.at(sa);
        fRawEvent->SetWF(theWF, isipm);
    }

    //fill header
    fRawEvent->FillHeader(fEventNumber, //eventnumber
                          "NULL", //date
                          "NULL", //run type
                          "NULL", //Rev Version
                          "NULL", //input file
                          0); //time stamp

    //Done
    fEventNumber++;
    
    if (fEventNumber >= fMaxEvents)
        fIsOver = true;
    
    return fRawEvent;
}

//=============================================================================

bool MCDecoder::ReadEventHeader()
{
    //Dummy! Remove it if not necessary (dummy implementation inherited)
    return true;
}

//=============================================================================

bool MCDecoder::ReadRunHeader()
{
    //Dummy! Remove it if not necessary (dummy implementation inherited)
    return true;
}

//=============================================================================

bool MCDecoder::ReadRunTrailer()
{
    //Dummy! Remove it if not necessary (dummy implementation inherited)
    return true;
}
