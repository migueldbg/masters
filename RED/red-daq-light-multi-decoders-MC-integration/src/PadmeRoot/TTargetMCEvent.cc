// --------------------------------------------------------------
// History:
//
// Created by Emanuele Leonardi (emanuele.leonardi@roma1.infn.it) 2016-02-02
//
// --------------------------------------------------------------

#include "TTargetMCEvent.hh"

#include "TTargetMCHit.hh"

ClassImp(TTargetMCEvent)

TTargetMCEvent::TTargetMCEvent() : TMCVEvent(TTargetMCHit::Class())
{;}

TTargetMCEvent::~TTargetMCEvent()
{;}
