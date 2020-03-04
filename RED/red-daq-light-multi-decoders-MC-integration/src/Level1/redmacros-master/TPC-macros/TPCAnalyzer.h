//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Aug  2 16:39:32 2019 by ROOT version 6.04/14
// from TTree reco/Reco events
// found on file: reco/naples/camp_VIII/run_1137.root
//////////////////////////////////////////////////////////

#ifndef TPCAnalyzer_h
#define TPCAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH2.h>

// Header file for the classes stored in the TTree if any.
#include "EvRec0.hh"
#include "TObject.h"

class TPCAnalyzer : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
 //EvRec0          *recoevent;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Int_t           event_number;
   Double_t        charge_total;
   Double_t        f90_tot;
   Double_t        charge_total_lsci;
   Double_t        lsci_psd_tot;
   Double_t        baseline_mean_all;
   Double_t        baseline_rms_all;
   vector<double>  baseline_mean;
   vector<double>  baseline_rms;
   vector<double>  charge;
   vector<double>  f90;
   vector<double>  start_time;
   vector<double>  xmin;
   vector<double>  xmax;
   vector<double>  ymin;
   vector<double>  ymax;
   Int_t           number_of_clusters;
   vector<RDCluster*> clusters;

   // List of branches
   TBranch        *b_recoevent_fUniqueID;   //!
   TBranch        *b_recoevent_fBits;   //!
   TBranch        *b_recoevent_event_number;   //!
   TBranch        *b_recoevent_charge_total;   //!
   TBranch        *b_recoevent_f90_tot;   //!
   TBranch        *b_recoevent_charge_total_lsci;   //!
   TBranch        *b_recoevent_lsci_psd_tot;   //!
   TBranch        *b_recoevent_baseline_mean_all;   //!
   TBranch        *b_recoevent_baseline_rms_all;   //!
   TBranch        *b_recoevent_baseline_mean;   //!
   TBranch        *b_recoevent_baseline_rms;   //!
   TBranch        *b_recoevent_charge;   //!
   TBranch        *b_recoevent_f90;   //!
   TBranch        *b_recoevent_start_time;   //!
   TBranch        *b_recoevent_xmin;   //!
   TBranch        *b_recoevent_xmax;   //!
   TBranch        *b_recoevent_ymin;   //!
   TBranch        *b_recoevent_ymax;   //!
   TBranch        *b_recoevent_number_of_clusters;   //!
   TBranch        *b_recoevent_clusters;   //!

   //Histos
   TH1F * hS1;
   TH1F * hTBA;
   TH2F * h2S1vsTBA;
   TH1F * hS1cor;
   TH1F * hS2;
   TH1F * hDriftT;
   TH1F * hS3Time;
   TH1F * hS3;
   TH2F * h2S1vsT;
   TH2F * h2S2vsT;


   Int_t isDoublePhase=0;
   Int_t isDoublePhaseDAQ=0;
   Double_t LYFit[2];
   Double_t FanoFit[2];
   Double_t resFit[2];

   //EndHistos

   TPCAnalyzer(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~TPCAnalyzer() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(TPCAnalyzer,0);
};

#endif

#ifdef TPCAnalyzer_cxx
void TPCAnalyzer::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_recoevent_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_recoevent_fBits);
   fChain->SetBranchAddress("event_number", &event_number, &b_recoevent_event_number);
   fChain->SetBranchAddress("charge_total", &charge_total, &b_recoevent_charge_total);
   fChain->SetBranchAddress("f90_tot", &f90_tot, &b_recoevent_f90_tot);
   fChain->SetBranchAddress("charge_total_lsci", &charge_total_lsci, &b_recoevent_charge_total_lsci);
   fChain->SetBranchAddress("lsci_psd_tot", &lsci_psd_tot, &b_recoevent_lsci_psd_tot);
   fChain->SetBranchAddress("baseline_mean_all", &baseline_mean_all, &b_recoevent_baseline_mean_all);
   fChain->SetBranchAddress("baseline_rms_all", &baseline_rms_all, &b_recoevent_baseline_rms_all);
   fChain->SetBranchAddress("baseline_mean", &baseline_mean, &b_recoevent_baseline_mean);
   fChain->SetBranchAddress("baseline_rms", &baseline_rms, &b_recoevent_baseline_rms);
   fChain->SetBranchAddress("charge", &charge, &b_recoevent_charge);
   fChain->SetBranchAddress("f90", &f90, &b_recoevent_f90);
   fChain->SetBranchAddress("start_time", &start_time, &b_recoevent_start_time);
   fChain->SetBranchAddress("xmin", &xmin, &b_recoevent_xmin);
   fChain->SetBranchAddress("xmax", &xmax, &b_recoevent_xmax);
   fChain->SetBranchAddress("ymin", &ymin, &b_recoevent_ymin);
   fChain->SetBranchAddress("ymax", &ymax, &b_recoevent_ymax);
   fChain->SetBranchAddress("number_of_clusters", &number_of_clusters, &b_recoevent_number_of_clusters);
   fChain->SetBranchAddress("clusters", &clusters, &b_recoevent_clusters);
}

Bool_t TPCAnalyzer::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef TPCAnalyzer_cxx
