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
VDecoder::VDecoder() :
  fRunHeader(0),fVerbosity(0),
  fIsOver(false),fEventNumber(0),fRunNumber(0),fAlignTimes(false),
  fEventsToSkip(0),fEventSkipped(false),fCustomNSamples(false),
  fNSamples(-1)
{
  fRawEvent = new EvRaw0();  
}

//=============================================================================

VDecoder::~VDecoder() 
{  
  if (fRawEvent)
    delete fRawEvent; 
  if (fRunHeader)
    delete fRunHeader;
}

//=============================================================================
bool VDecoder::OpenFiles(string /*listfile*/,int /*runnr*/)
{
  //Dummy implementation
  return true;
}

//=============================================================================
bool VDecoder::CloseFiles(bool /*updatedb*/)
{
  //Dummy implementation
  return true;
}

//=============================================================================
bool VDecoder::ReadRunHeader()
{
  //dummy implementation
  return true;
}

//=============================================================================
bool VDecoder::ReadRunTrailer()
{
  //dummy implementation
  return true;
}

//=============================================================================
EvRaw0* const VDecoder::GetRawEvent() const
{
  return fRawEvent;
}
//=============================================================================
RunHeader* VDecoder::GetRunHeader() 
{
  return fRunHeader;
}

//=============================================================================
bool VDecoder::ReadEventHeader()
{
  //dummy implementation
  return true;
}

//=============================================================================
void VDecoder::SetCustomNSamples(int ns)
{
  fCustomNSamples = true;
  fNSamples = ns;
}
