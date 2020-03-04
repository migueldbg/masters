//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 28 11:02:23 2020 by ROOT version 6.12/06
// from TTree reco/Reco events
// found on file: run_1397.root
//////////////////////////////////////////////////////////

#ifndef RedAna_h
#define RedAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph2D.h>

// Headers needed by this particular selector
#include "EvHeader.hh"
#include "RDCluster.hh"
#include "RDPulseFit.hh"

#ifndef UNDEF
#define UNDEF -100
#endif

struct RunInfo {
	RunInfo() { clear(); }

	Double_t tdrift_max;

	//For S1 correction
	vector<Double_t> para_s1_tba_corr, para_s1_tdrift_corr;

	//For S2 correction
	Double_t electron_lifetime;

	inline void clear(){
		tdrift_max=UNDEF;
		para_s1_tba_corr.clear();
		para_s1_tdrift_corr.clear();
		electron_lifetime=UNDEF;
	}
};


class RedAna : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTreeReader     fReaderForGetEntries; // To get Entries. Avoid clush due to a bug.
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
//   TTreeReaderValue<EvRec0> evReco = {fReader, "recoevent"};
   TTreeReaderValue<unsigned int> fUniqueID = {fReader, "fUniqueID"};
   TTreeReaderValue<unsigned int> fBits = {fReader, "fBits"};
   TTreeReaderArray<EvHeader> evheader = {fReader, "evheader"};
   TTreeReaderValue<Int_t> event_number = {fReader, "event_number"};
   TTreeReaderValue<Double_t> charge_total = {fReader, "charge_total"};
   TTreeReaderValue<Double_t> f90_tot = {fReader, "f90_tot"};
   TTreeReaderValue<Double_t> charge_total_lsci = {fReader, "charge_total_lsci"};
   TTreeReaderValue<Double_t> lsci_psd_tot = {fReader, "lsci_psd_tot"};
   TTreeReaderValue<Double_t> baseline_mean_all = {fReader, "baseline_mean_all"};
   TTreeReaderValue<Double_t> baseline_rms_all = {fReader, "baseline_rms_all"};
   TTreeReaderArray<double> baseline_mean = {fReader, "baseline_mean"};
   TTreeReaderArray<double> baseline_rms = {fReader, "baseline_rms"};
   TTreeReaderArray<double> charge = {fReader, "charge"};
   TTreeReaderArray<double> f90 = {fReader, "f90"};
   TTreeReaderArray<double> start_time = {fReader, "start_time"};
//   TTreeReaderArray<double> start_time_th = {fReader, "start_time_th"};
   TTreeReaderArray<double> xmin = {fReader, "xmin"};
   TTreeReaderArray<double> xmax = {fReader, "xmax"};
   TTreeReaderArray<double> ymin = {fReader, "ymin"};
   TTreeReaderArray<double> ymax = {fReader, "ymax"};
   TTreeReaderValue<Int_t> number_of_clusters = {fReader, "number_of_clusters"};
   TTreeReaderArray<RDCluster> clusters = {fReader, "clusters"};
   TTreeReaderArray<RDPulseFit> fits = {fReader, "fits"};
   TTreeReaderValue<Int_t> nfits = {fReader, "nfits"};
 //  TTreeReaderArray<int> fits = {fReader, "fits"};


   RedAna(TTree * /*tree*/ =0) : fRunInfo() { }
   virtual ~RedAna() { }
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


   //Custom functions
   Bool_t          isAmpeak();
   Double_t        getS1corr(Double_t t_drift);
   Double_t        getS1corr(Double_t t_drift, Int_t channel);
   void            getTopMaxChan_and_Frac(RDCluster clstr, Int_t &max_chan, Double_t &frac);
   Double_t        getS13Dcorr(Double_t bar_x, Double_t bar_y, Double_t t_drift, Double_t topCharge, Double_t botCharge);
   Double_t        getS2Zcorr(Double_t t_drift);
   Double_t        getS2MaxS2Chancorr(Int_t max_s2_chan);
   Int_t		   GetMaxChannel(double x, double y);
   void            ReadRunConfig();
   void            SetRunNumber(Int_t runnumber){}
   Bool_t          IsInnerChannel(Int_t ch);
   void            BookHistograms();
   void            SetOutputFile(TString name){fOutName = name;}

private:
   RunInfo fRunInfo;
   TString fOutName;

   //Function for corrections
   TF1 *fpol;

   //Histograms
   //For SiTel
   TH2D *hBanana, *hCBanana;
   TH1D *hdT_Si_TPC;

   //For LSci
   TH1D **hdT_Si_LSci, **hdT_TPC_LSci;
   TH2D **hLSci_PSD_Charge, **hdT_Si_LSci_TPC, **hdT_Si_LSci_LSciPSD, **hF90_S1_3Coinci_LSci, **hS2_S1_3Coinci_LSci;

   //For TPC
   TH1D *hEventConter;
//   TH2D *hTotalS1corr_TDrift, *hTotalF90_TDrift, *hRun_dT_Trigger;
   TH1D *hTDrift, *hF90, *hS1, *hS1zcorr, *hS2, *hS2_zcorr, *hS2_xycorr, *hS2_xyzcorr, *hNclusters;
   TH1D *hTBAS1, *hTBAS1All, *hTBAS1wS2s, *hTBAS2, *hS2_efficiency, *hS2MaxChan;
   TH2D *hF90_S1, *hF90_S1_SiTelCoinci, *hS1_tdrift, *hS1zcorr_tdrift, *hS2_tdrift, *hS2ovS1_tdrift, *hS2_Am_tdrift, *hS2_Am_s2maxchan;
   TH2D *hS1Topcorr_tdrift, *hS13Dcorr_tdrift, *hS2Topcorr_tdrift;
   TH2D *hNclusters_evtnum, *htdrift_evtnum, *hS2RMS_evtnum, *hS2TR_evtnum, *hS2ovS1_evtnum;
   TH2D *hYBary_XBary, *hYBary_XBary_corr, *hXBary_tdrift, *hYBary_tdrift;
   TH2D *hTBAS1_tdrift, *hS1_TBAS1, *hS1zcorr_TBAS1, *hS2_TBAS2, *hS2_S1, *hS2_S1_SiTelCoinci;
   TH2D *hdT_Si_TPC_TPCF90, *hStartTime_ChNum, *hS2_frac_chnum;

   TH3D *hYBary_XBary_s2maxchan, *hS1_tdrift_s2maxchan;

   TH2D *hChTopS1Frac;
   TMultiGraph *gMYBary_XBary_corr;
   TGraph *gYBary_XBary_corr[22];
   TH3D *hS1_barXY[10], *hS1Top_barXY[10], *hS1Bottom_barXY[10];

   ClassDef(RedAna,0);

   long int AtEntry;
   long int NEntries;
	void TopChCorr(vector<double> top, vector<double> bot, double &xx, double &yy, double &chargeTopcorr);
	double xch[22];
	double ych[22];
	int chS1FracFill[22];

	int runtype;
	TGraph2D* gS1map[10][2];

};

#endif

#ifdef RedAna_cxx
void RedAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

	tree->LoadTree(-1);
   fReader.SetTree(tree);
   fReaderForGetEntries.SetTree(tree);
   NEntries = fReader.GetEntries(true);
}

Bool_t RedAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef RedAna_cxx
