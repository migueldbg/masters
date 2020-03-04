#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TRint.h"
#include "TProof.h"
#include "TProofMgr.h"
#include "TProofLog.h"
//#include "TROOT.h"
#include "TDSet.h"

using namespace std;

#include "src/RedAna.h"

TRint* theApp;
void Process(TChain* chain, TString fHistOutName);
void ProcessRecoReD(TString fInNameEvent, TString outDirname ="");


void ProcessRecoReD(TString fInNameEvent, TString outDirname) {
  gROOT->SetBatch(kTRUE);

  TChain* chain = new TChain("reco");
  chain->Add(fInNameEvent.Data());

  TString fHistOutName(fInNameEvent);
  if (outDirname.CompareTo("")!=0) fHistOutName.Replace(0, fHistOutName.Last('/'), outDirname); // replace directory name to output only if out_dir is not "".
  fHistOutName.Replace(fHistOutName.Length()-5, 5, "_Hist.root");

  Process(chain, fHistOutName);

  delete chain;

}

void ProcessRecoReD() {
  gROOT->SetBatch(kTRUE);

  TChain* chain = new TChain("reco");
  chain->Add("data/run_1501.root");
  chain->Add("data/run_1502.root");
  chain->Add("data/run_1503.root");
  chain->Add("data/run_1504.root");
  chain->Add("data/run_1506.root");
  chain->Add("data/run_1508.root");
  chain->Add("data/run_1509.root");
  chain->Add("data/run_1510.root");
  chain->Add("data/run_1511.root");
  chain->Add("data/run_1512.root");
  chain->Add("data/run_1513.root");
  chain->Add("data/run_1514.root");
//  chain->Add("data/run_1515.root");
//  chain->Add("data/run_1516.root");

  Process(chain, "AllRun_Hist.root");

  delete chain;

}

void Process(TChain* chain, TString fHistOutName){
  std::cout << "outhistfile: " << fHistOutName.Data() << std::endl;

  TString home(gSystem->Getenv("PWD"));
  TString option = fHistOutName;

  RedAna *selector = new RedAna();
  selector->SetOutputFile(fHistOutName.Data());
  chain->Process(selector, option.Data());
//  chain->Process(selector, option.Data(), 1000000, 0);//nEvents,0);

  std::cout << "outhistfile: " << fHistOutName.Data() << std::endl;

}

#ifndef __CINT__
int main(int argc, char **argv) {
	theApp = new TRint("App", &argc, argv, NULL, 0);
	theApp->Connect("KeyPressed(Int_t)","TSystem",gSystem,"ExitLoop()");
	if ( theApp->Argc() == 1 ) {
		ProcessRecoReD();
	} else if ( theApp->Argc() == 2 ) {
		std::cout << "\n==========> ProcessRecoReD <=============" << std::endl;
		ProcessRecoReD(theApp->Argv(1));

	} else if ( theApp->Argc() == 3 ) {
		std::cout << "\n==========> ProcessRecoReD <=============" << std::endl;
		ProcessRecoReD(theApp->Argv(1), theApp->Argv(2));
	} else {
		std::cout << "Usage:" <<theApp->Argc()<< std::endl;
        std::cout << "./ProcessRecoReD fInNameEvent.root" << std::endl;
        std::cout << "./ProcessRecoReD fInNameEvent.root outDir" << std::endl;
		return 0;
	}

	std::cout << "==> Application finished." << std::endl;
	return 0;
}
#endif /* __CINT __ */
