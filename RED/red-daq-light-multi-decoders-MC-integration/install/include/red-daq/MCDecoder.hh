#ifndef _MCDecoder_hh_
#define _MCDecoder_hh_

#include <string>
#include <iostream>

#include "EvRaw0.hh"
#include "RunHeader.hh"
#include "VDecoder.hh"

class MCTraces;

class MCDecoder : public VDecoder
{  
public:
    MCDecoder(std::string configfile, int evts=100000);
    ~MCDecoder();

    //file management
    // Implement virtual interface.
    bool OpenFiles(std::string,int);
    bool CloseFiles(bool);
    EvRaw0* const ReadEvent();
    bool ReadEventHeader();
    bool ReadRunHeader();
    bool ReadRunTrailer();

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

private:

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
    MCTraces* fTraces;
    int       fMaxEvents;
};

#endif
