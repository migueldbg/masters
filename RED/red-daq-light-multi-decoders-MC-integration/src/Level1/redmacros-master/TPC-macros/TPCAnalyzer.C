#define TPCAnalyzer_cxx
// The class definition in TPCAnalyzer.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("TPCAnalyzer.C")
// root> T->Process("TPCAnalyzer.C","some options")
// root> T->Process("TPCAnalyzer.C+")
//

#include "TPCAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TStyle.h"

#include <iostream>
#include <vector>

#define TIME 0.002

using namespace std ;
using namespace TMath ;

Double_t f90cut=0.2;
Int_t    nClus = 1;

TH1F *ham ;
double fam241(double *x, double *p) {
  double LY  = p[0];
  double sigma0 = p[1] ;
  double fano   = p[2] ;
  double width = ham->GetBinWidth(3);
  double tot = 0 ;
  for(int i=1;i<ham->GetXaxis()->GetNbins();++i) {
    double ene = ham->GetBinCenter(i);
    if(ene > 60) continue; 
    double val   = ham->GetBinContent(i);
    double npe   = LY*ene ;
    double sigma = sqrt(pow(sigma0,2)+npe*fano);
    //cout << npe << " " << sigma << " " << fano << endl;
    tot += val*Gaus(x[0],npe,sigma,1)/width;    
  }
  
  return p[3]*tot ;
  
}


void TPCAnalyzer::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();


   //Create Histograms
   hS1    = new TH1F("hS1","S1 Raw Charge",400,0.,1000.);
   hTBA   = new TH1F("hTBA","Tob Bottom Asymmetry",200,-0.5,0.5);
   h2S1vsTBA = new TH2F ("h2S1vsTBA","S1 vs TBA",200,-0.5,0.5,400,0,1000.);
   // applying cuts
   hS1cor    = new TH1F("hS1cor","S1 Corrected",400,0.,1000.);
   // more for double phase operation
   hS2    = new TH1F("hS2","S2 Raw Charge",400,0.,30000.);
   hDriftT = new TH1F("hDriftT","Drift Time [us]",180,0,90);
   hS3Time = new TH1F("hS3Time","S2-S2 Time [us]",90,0,90);
   hS3    = new TH1F("hS3","S3 Raw Charge",200,0.,400.);
   h2S1vsT = new TH2F ("h2S1vsT","S1 vs Drift Time",180,0,90,400,0,1000.);
   h2S2vsT = new TH2F ("h2S2vsT","S2 vs Drift Time",180,0,90,400,0,30000.);
}

void TPCAnalyzer::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t TPCAnalyzer::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either TPCAnalyzer::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

  //RunFlag
  // RunFlag==1  Am241 runs in SP
  // RunFlag==2  Am241 runs in DP

  fChain->GetEntry(entry);
  
  // Define Cuts (depending on RunFlag)
  
  if (isDoublePhase==0) {
    if (number_of_clusters!=1) return kTRUE;
    if (clusters[0]->f90<=f90cut || clusters[0]->f90>0.6) return kTRUE;
    if (clusters[0]->rep==0) return kTRUE;
  } else if (isDoublePhase==1) {
    if (number_of_clusters<2) return kTRUE;
    if (clusters[0]->f90<=f90cut || clusters[0]->f90>0.6) return kTRUE;
    if (clusters[0]->rep==0) return kTRUE;
    if (clusters[1]->f90>f90cut) return kTRUE;
  }

  Double_t S1time=clusters[0]->cdf_time * TIME;
  if (isDoublePhase==0 && isDoublePhaseDAQ==0) {
    if (S1time <5. || S1time >6.) return kTRUE;
  } else if (isDoublePhase==1 && isDoublePhaseDAQ==0) {
    if (S1time <9. || S1time >10.) return kTRUE;
  } else if (isDoublePhase==0 && isDoublePhaseDAQ==1) {
    if (S1time <9. || S1time >10.) return kTRUE;
  } else if (isDoublePhase==1 && isDoublePhaseDAQ==2) { //long daq window
    if (S1time <59. || S1time >60.) return kTRUE;
  }
  Double_t S2time=-1;
  Double_t S3time=-1;
  Double_t driftTime=-1;
  Double_t S3S2Time=-1;
  Double_t S3S1Time=-1;
  if (number_of_clusters>1) {
    S2time=clusters[1]->fixed_time * TIME;
    driftTime=S2time-S1time;
  }
  if (number_of_clusters>2) {
    S3time=clusters[1]->fixed_time * TIME;
    S3S1Time=S3time-S1time;
    S3S2Time=(clusters[2]->cdf_time-clusters[1]->cdf_time) * TIME;
  }

  Double_t TBA = (clusters[0]->tot_charge_top-clusters[0]->tot_charge_bottom)/
    (clusters[0]->tot_charge_top+clusters[0]->tot_charge_bottom);
  Double_t S1 = clusters[0]->charge;

  //  cout <<TBA<<" "<<S1<<" "<<clusters[0]->tot_charge_top<<" "<<clusters[0]->tot_charge_bottom<<endl;


  //Fill S1 for quenching studies only for nclus==2 in double phase case

  if (isDoublePhase==0) {
    hS1->Fill(S1);
    hTBA->Fill(TBA);
    h2S1vsTBA->Fill(TBA,S1);
    //S1 cor for now only a cut on TBA
    if (TBA>-0.4 && TBA<-0.05) hS1cor->Fill(S1);

  } else if (isDoublePhase==1) {

    if (number_of_clusters==2) {
      hS1->Fill(S1);
      hTBA->Fill(TBA);
      h2S1vsTBA->Fill(TBA,S1);
      //S1 cor for now only a cut on TBA
      if (TBA>-0.4 && TBA<-0.05) hS1cor->Fill(S1);
    }

  }

  
  // return if no double phase
  if (!isDoublePhase) return kTRUE; 
  
  Double_t S2=clusters[1]->charge;
  
  hS2->Fill(S2);
  hDriftT->Fill(driftTime);
  
  h2S1vsT->Fill(driftTime,S1);
  h2S2vsT->Fill(driftTime,S2);

  if (number_of_clusters>2 && clusters[2]->charge<400) {
    Double_t S3=clusters[2]->charge;
    hS3->Fill(S3);
    //    cout <<S3S2Time<<endl;
    hS3Time->Fill(S3S2Time);
  }

    

  return kTRUE;
}

void TPCAnalyzer::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void TPCAnalyzer::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.


  //int myfitter(TH1D * htop, float min=480., float max=800.) {
  gStyle->SetOptFit(1111);
  hS1->Draw();
  
  TFile *fmc = TFile::Open("am241.root");
  ham = (TH1F*) fmc->Get("ham");
  TF1 * fun = new TF1("fun",fam241,0,1000,4);
  fun->SetParameters(hS1->GetMaximumBin()*hS1->GetBinWidth(1)/59.6,5,3.,1e4);
  fun->SetLineColor(64);
  //  fun->Draw();
  fun->SetParNames("LY","#sigma_0","Fano","A");
  fun->FixParameter(1,0);
  fun->FixParameter(2,3);
  //  fun->SetParLimits(2,0.,10.);
  //  TH1D *hd = htop ;
  hS1cor->GetXaxis()->SetTitle("[pe]");
  //  float min = 480;
  //  float max = 800;
  
  float min= (float) hS1cor->GetMaximumBin()*hS1cor->GetBinWidth(1) * 0.8; //lower fit range peak-20%
  float max= (float) hS1cor->GetMaximumBin()*hS1cor->GetBinWidth(1) * 1.15; //upper fit range peak+10%

  hS1cor->Fit("fun","R","",min, max);
  fun->ReleaseParameter(2);
  hS1cor->Fit("fun","R","",min, max);
  hS1cor->Fit("fun","R","",min, max);

  LYFit[0] = fun->GetParameter(0);   LYFit[1] = fun->GetParError(0);
  FanoFit[0] = fun->GetParameter(2);   FanoFit[1] = fun->GetParError(2);
  resFit[0] = sqrt(FanoFit[0]/(LYFit[0]*59.6));  
  resFit[1] = 0.5 * sqrt(pow(FanoFit[1]/FanoFit[0],2) + pow(LYFit[1]/LYFit[0],2))*resFit[0];

  delete(fun);
}


