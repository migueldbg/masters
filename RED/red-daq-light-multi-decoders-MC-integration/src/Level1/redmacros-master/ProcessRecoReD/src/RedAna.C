#define RedAna_cxx
// The class definition in RedAna.h has been generated automatically
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
// root> T->Process("RedAna.C")
// root> T->Process("RedAna.C","some options")
// root> T->Process("RedAna.C+")
//


#include "RedAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>
#include <TF1.h>
#include <iostream>
using namespace std;

const Int_t NLSciChan = 6;
//#define BF_RUN1628
#ifdef BF_RUN1628
//Original Configuration up to Run1628
const TString LSciName[NLSciChan]            = {"LSci0", "LSci1", "LSci3_7", "LSci4_6", "LSci5", "LSci8"};
const TString LSciAngle[NLSciChan]           = {  "4.3",    "90",      "40",      "20",     "0",    "90"};
const Int_t LSci_ChID[NLSciChan]             = {      0,       1,         3,         4,       5,       8};
const Double_t LSci_PSD_threshold[NLSciChan] = {   0.16,    0.14,      0.14,      0.12,    0.15,    0.15};
double LSciOffset[NLSciChan]                 = {      6,       4,         6,         4,       0,       2}; //LSci time offset in ns //the ch0 have larger seperation in fast and slow
//double LSciOffset[NLSciChan]                 = {      0,     2.4,       4.4,       3.8,       0,       0}; //From Marco, LSci time offset in ns //the ch0 have larger seperation in fast and slow

#else
//Channel Configuration after Run1629
const TString LSciName[NLSciChan]            = {"LSci0", "LSci1", "LSci3_7", "LSci4_6", "LSci5", "LSci8"};
const TString LSciAngle[NLSciChan]           = {  "4.3",    "90",      "40",      "20",     "0",    "90"};
const Int_t LSci_ChID[NLSciChan]             = {      0,       1,         2,         3,       4,       5};
const Double_t LSci_PSD_threshold[NLSciChan] = {   0.16,    0.14,      0.14,      0.12,    0.15,    0.15};
double LSciOffset[NLSciChan]                 = {      6,       4,         6,         4,       0,       2}; //LSci time offset in ns //the ch0 have larger seperation in fast and slow
#endif

const Int_t ch_SiDeltaE(30), ch_SiE(31);

void RedAna::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   AtEntry = -1;
}

void RedAna::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   ReadRunConfig();
   BookHistograms();
    double xchinit[22] = {1, 2, 3, 2,   0, 1, 1, 3, 0, 1, 2,   3, 2, 0, 1, 2, 3, 0, 1, 3, 2, 0};
	double ychinit[22] = {3, 5, 4, 4, 4.5, 4, 5, 5, 0, 0, 0, 0.5, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3};
	for(int i=0; i<22; i++){
	  xch[i] = (xchinit[i]+0.5)*5./4.;
	  ych[i] = (ychinit[i]+0.5)*5./6.;
	  chS1FracFill[i] = 0;
	}

}

Bool_t RedAna::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

	Int_t nentries = fReaderForGetEntries.GetEntries(true);
	if (entry % (nentries / 10) == 0) printf("Processing Entry number %ld [%ld %% of %d]\n", (long int) entry, (long int) (entry / (nentries / 100)), nentries);
	if (gDebug) std::cout << Form("============================Event %d Starts!!=============================",(int) entry + 1) << std::endl;

	fReader.SetLocalEntry(entry);

	if(*event_number < AtEntry){
		std::cout<<"get the end of one data file."<<std::endl;
	}
	AtEntry = *event_number;

	for (Int_t ich =0; ich<32; ich++) hStartTime_ChNum->Fill(ich, start_time[ich]);
	//This is the info of the Si detectors
	double Si_deltaE, Si_E, Si_T;
	if(runtype==1){
		Si_deltaE = baseline_mean[ch_SiDeltaE] - ymin[ch_SiDeltaE];
		Si_E = baseline_mean[ch_SiE] - ymin[ch_SiE];
		Si_T = start_time[ch_SiDeltaE];
		hBanana->Fill(Si_E, Si_deltaE);

		Bool_t isBeLow = false;
		if (Si_deltaE >700 && Si_deltaE<900 && Si_E>2200 && Si_E<2600) isBeLow = true;
	//	if (!isBeLow) return kTRUE;
	}

	hNclusters->Fill(*number_of_clusters);
	hNclusters_evtnum->Fill(entry, *number_of_clusters);

	Bool_t isTPCCoinci(false);
	Double_t f90cut(0.1);
	if (*number_of_clusters==0) return kTRUE; // no clusters
	if (clusters[0].f90<f90cut) return kTRUE; // select only S1 for the first pulse

	Double_t dT_Si_TPC;
	if(runtype==1){
		dT_Si_TPC = (clusters[0].cdf_time - Si_T)*2.; // in ns  cdf_time gives better correlation with SiTel
		hdT_Si_TPC->Fill(dT_Si_TPC);
		if(TMath::Abs(dT_Si_TPC+25.)<10.) isTPCCoinci = true;
		if (isTPCCoinci) hCBanana->Fill(Si_E, Si_deltaE);
	}

//	if(clusters[0].f90<0.4) return kTRUE; // neutron selection in TPC
//	if(clusters[0].f90>0.4) return kTRUE; // gamma selection in TPC

	Int_t LSci_ch(0);
	vector<Bool_t> isLSciCoinci(NLSciChan);
	if(runtype==1){
		for(Int_t ilsci=0; ilsci<NLSciChan; ilsci++){
			isLSciCoinci.push_back(false);
			LSci_ch = LSci_ChID[ilsci];
			Double_t charge_LSci = charge[LSci_ch];//
			Double_t amp_LSci = -(ymin[LSci_ch]-baseline_mean[LSci_ch]); // minus because wf is negative

			hLSci_PSD_Charge[ilsci]->Fill(charge_LSci, f90[LSci_ch]);

			if(charge_LSci < 2000.) continue; // LSci cuts
			Double_t dT_Si_LSci = (start_time[LSci_ch] - Si_T)*2. - LSciOffset[ilsci]; // in ns
			Double_t dT_TPC_LSci = (start_time[LSci_ch] - clusters[0].cdf_time)*2. - LSciOffset[ilsci]; // in ns

			hdT_Si_LSci_LSciPSD[ilsci]->Fill(dT_Si_LSci, f90[LSci_ch]);
	//		if(f90[LSci_ch] < LSci_PSD_threshold[ilsci]) continue; // LSci cuts

			hdT_Si_LSci[ilsci]->Fill(dT_Si_LSci);
			hdT_TPC_LSci[ilsci]->Fill(dT_TPC_LSci);
			hdT_Si_LSci_TPC[ilsci]->Fill(dT_Si_TPC, dT_Si_LSci);

			if(TMath::Abs(dT_Si_LSci-25.)<10.) isLSciCoinci[ilsci] = true;
		}
	}

//////////////////  TPC Part  /////////////////////////
//	if(clusters[0].rep!=1) return kTRUE; // pile up cuts
	Double_t S1 = clusters[0].charge;
	Double_t TBAS1 = (clusters[0].tot_charge_top - clusters[0].tot_charge_bottom) / S1;
	hTBAS1All->Fill(TBAS1);

	if (*number_of_clusters > 1) hTBAS1wS2s->Fill(TBAS1);
	if (*number_of_clusters != 2) return kTRUE; // only event with 2 pulses
//	if (*number_of_clusters!=1) return kTRUE; // only one pulse event
//	if(clusters[1].rep!=1) return kTRUE; // pile up cuts
//	if(clusters[0].rep!=1 || clusters[1].rep!=1) return kTRUE; // pile up cuts
//	if(clusters[0].rep==1 && clusters[1].rep==1) return kTRUE;

	hF90->Fill(clusters[0].f90);
	hTBAS1All->Fill(TBAS1);

//	   if (clusters[0].f90<f90cut ||  clusters[1].f90>f90cut) return kTRUE;
	if (*number_of_clusters == 2 && clusters[1].f90>f90cut) return kTRUE; // S2 only for the second pulse.

	Double_t S2 = (*number_of_clusters == 2)? clusters[1].charge: -1.;

	Double_t s2_max_frac(0.);
	Int_t s2_max_chan(-1);
	if (*number_of_clusters == 2) getTopMaxChan_and_Frac(clusters[1], s2_max_chan, s2_max_frac);
	hS2MaxChan->Fill(s2_max_chan);


	Double_t t_drift = (*number_of_clusters == 2)? (clusters[1].cdf_time - clusters[0].cdf_time) / 500.: -1.;
	Double_t S1_zcorr = getS1corr(t_drift, s2_max_chan) * S1;
	Double_t S2_zcorr = getS2Zcorr(t_drift) * S2;
	Double_t S2XYcorrFactor = getS2MaxS2Chancorr(s2_max_chan);
	Double_t S2_xycorr = S2XYcorrFactor * S2;
	Double_t S2_xyzcorr = S2XYcorrFactor * S2_zcorr;

	Double_t S1_topcorr, S2_topcorr, S1_3Dcorr;
	Double_t S1_bar_x, S1_bar_y, S2_bar_x, S2_bar_y;
	TopChCorr(clusters[0].charge_top, clusters[0].charge_bottom, S1_bar_x, S1_bar_y, S1_topcorr);
	TopChCorr(clusters[1].charge_top, clusters[1].charge_bottom, S2_bar_x, S2_bar_y, S2_topcorr);
	//S1_topcorr *= S1_zcorr/S1;

	hS1->Fill(S1);
	hS1zcorr->Fill(S1_zcorr);
//	hS1zcorr->Fill(S1_3Dcorr);
	hS2->Fill(S2);
	hF90_S1->Fill(S1, clusters[0].f90);
	hS2_S1->Fill(S1, S2);
	hS2_zcorr->Fill(S2_zcorr);
	hS2_xycorr->Fill(S2_xycorr);
	hS2_xyzcorr->Fill(S2_xyzcorr);

	if(runtype==1){
		if (isTPCCoinci) {
			hF90_S1_SiTelCoinci->Fill(S1, clusters[0].f90);
			hS2_S1_SiTelCoinci->Fill(S1, S2);
		}
		for(Int_t ilsci=0; ilsci<NLSciChan; ilsci++) {
			if (!isTPCCoinci || !isLSciCoinci[ilsci]) continue;
	//		hF90_S1_3Coinci_LSci[ilsci]->Fill(S1, clusters[0].f90);
	//		hS2_S1_3Coinci_LSci[ilsci]->Fill(S1, S2);
			hF90_S1_3Coinci_LSci[ilsci]->Fill(S1_zcorr, clusters[0].f90);
			hS2_S1_3Coinci_LSci[ilsci]->Fill(S1_zcorr, S2_zcorr);
	//		cout<<"S1: "<<S1<<" S1 corr: "<<S1_zcorr<<endl;
	//		cout<<"S2: "<<S2<<" S2 corr: "<<S2_zcorr<<endl;
		}
		hdT_Si_TPC_TPCF90->Fill(dT_Si_TPC, clusters[0].f90);
	}

	hTDrift->Fill(t_drift);

	htdrift_evtnum->Fill(entry, t_drift);
	hS2ovS1_evtnum->Fill(entry, S2 / S1);

	Double_t TBAS2 = (*number_of_clusters == 2)? (clusters[1].tot_charge_top - clusters[1].tot_charge_bottom) / S2: -1.;
	hTBAS1->Fill(TBAS1);
	hTBAS2->Fill(TBAS2);
	hS1_TBAS1->Fill(TBAS1, S1);
	hS2_TBAS2->Fill(TBAS2, S2);
	hS1zcorr_TBAS1->Fill(TBAS1, S1_zcorr);
	hTBAS1_tdrift->Fill(t_drift, TBAS1);

	hS1_tdrift->Fill(t_drift, S1);
	hS2_tdrift->Fill(t_drift, S2);
	hS2Topcorr_tdrift->Fill(t_drift, S2_topcorr);
	hS1zcorr_tdrift->Fill(t_drift, S1_zcorr);
	hS1Topcorr_tdrift->Fill(t_drift, S1_topcorr);
	if (isAmpeak()) {
		hS2_Am_tdrift->Fill(t_drift, S2);
		hS2_Am_s2maxchan->Fill(s2_max_chan, S2);//S2_zcorr);
	}
//	if (isAmpeak() && !IsInnerChannel(s2_max_chan)) hS2_Am_tdrift->Fill(t_drift, S2);
	hS2ovS1_tdrift->Fill(t_drift, S2 / S1);

	if (*number_of_clusters == 2 && clusters[1].f90<0.1 &&clusters[0].f90>0.1){ // fine corrections.
		double x, y, S2_topcorr;
		//if(clusters[1].charge>100){
		//	cout<<"bar:"<<clusters[1].bar_x<<" "<<clusters[1].bar_y<<endl;
			x = clusters[1].bar_x;
			y = clusters[1].bar_y;
			TopChCorr(clusters[1].charge_top,clusters[1].charge_bottom, x, y, S2_topcorr);
		//}

		hYBary_XBary_corr->Fill(x,y);//clusters[1].bar_x, clusters[1].bar_y);
		if(clusters[0].f90>0.1 && clusters[1].f90<0.1 && clusters[0].charge>200 && clusters[0].charge<800){
			int maxChannel = GetMaxChannel(clusters[1].pos_x, clusters[1].pos_y);
			TGraph *gTemp = (TGraph*)(gMYBary_XBary_corr->GetListOfGraphs()->At(maxChannel));
			gTemp->SetPoint(gTemp->GetN(), x, y);//clusters[1].bar_x, clusters[1].bar_y);
			//if(chS1FracFill[maxChannel]<80*(maxChannel==4||maxChannel==11? 2 : 1) && t_drift>35){
				for(int ich=0; ich<22; ich++)
					hChTopS1Frac->Fill(xch[ich], ych[ich], clusters[0].charge_top[ich]);
				chS1FracFill[maxChannel]++;
			//}
		}
		hXBary_tdrift->Fill(t_drift, clusters[1].bar_x);
		hYBary_tdrift->Fill(t_drift, clusters[1].bar_y);
		hYBary_XBary->Fill(clusters[1].bar_x, clusters[1].bar_y);
		hYBary_XBary_s2maxchan->Fill(s2_max_chan, clusters[1].bar_x, clusters[1].bar_y);
//		if(!IsInnerChannel(s2_max_chan))
			hS1_tdrift_s2maxchan->Fill(s2_max_chan, t_drift, S1);

		//3D map corrected S1
		S1_3Dcorr = getS13Dcorr(x,y,t_drift,S1_topcorr-clusters[0].tot_charge_bottom, clusters[0].tot_charge_bottom);
		hS13Dcorr_tdrift->Fill(t_drift, S1_3Dcorr);

		Double_t S2RMS = TMath::Sqrt(clusters[1].rms_time) / 500.;
		hS2RMS_evtnum->Fill(entry, S2RMS);
		if(*nfits>0){
			Double_t S2TR = fits[0].par[3] ;
			hS2TR_evtnum->Fill(entry, S2TR);
		}
		if(runtype==0){
			if(clusters[0].charge>400&&clusters[0].charge<900){
				int iz = t_drift/7.;
				if(iz>=0&&iz<10){
					hS1_barXY[iz]->Fill(x,y,S1_topcorr);
					hS1Top_barXY[iz]->Fill(x,y,S1_topcorr-clusters[0].tot_charge_bottom);
					hS1Bottom_barXY[iz]->Fill(x,y,clusters[0].tot_charge_bottom);
				}
			}
		}
		if(runtype==1){//make final plots


		}
	}


	for(Int_t ich=0; ich<clusters[1].charge_top.size(); ich++){
		(ich==4||ich==11)? hS2_frac_chnum->Fill(ich, 0.5*clusters[1].charge_top[ich]/clusters[1].tot_charge_top): // for summed ch, divide by 2.
				           hS2_frac_chnum->Fill(ich, clusters[1].charge_top[ich]/clusters[1].tot_charge_top);

	}

	return kTRUE;
}

void RedAna::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

	TList *list = GetOutputList();

	hS1_tdrift->FitSlicesY(0, 0, -1, 20);
	TH1D *hmean = (TH1D*)gDirectory->Get(Form("%s_1",hS1_tdrift->GetName()));
	list->Add(hmean);

	hS2_tdrift->FitSlicesY(0, 0, -1, 20);
	hmean = (TH1D*)gDirectory->Get(Form("%s_1",hS2_tdrift->GetName()));
	list->Add(hmean);

	hS2ovS1_tdrift->FitSlicesY(0, 0, -1, 20);
	hmean = (TH1D*)gDirectory->Get(Form("%s_1",hS2_tdrift->GetName()));
	list->Add(hmean);

	hS2_efficiency->Divide(hTBAS1wS2s, hTBAS1All);

}

void RedAna::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

	  Info("RedAna::Terminate()","closining ...");

	  Int_t nhists = 0;
	  TList *list = GetOutputList();

	  TIter next(list);
	  TFile *hOutF = new TFile(fOutName.Data(), "RECREATE");
	  TH1 *h;
	  TObject* obj;
	  while ((obj = (TObject*) next()))
	  {
	      if(!obj) continue;
	      if (obj->InheritsFrom(TH1::Class()))
	      {
	          h = (TH1*) obj;
			  if(std::strcmp(h->GetName(),"hChTopS1Frac")==0){
				h->Scale(1./h->Integral()*24);
			  }
	          h->Write();
	      }else if (obj->InheritsFrom(TMultiGraph::Class()))
	      {
	          TMultiGraph* gM = (TMultiGraph*) obj;
	          gM->Write();
	      }
	  }
	  cout<<" S1Frac fill numbers: ";
	  for(int ich=0; ich<22; ich++){
		cout<<" "<<chS1FracFill[ich];
	  }
	  cout<<endl;

	  hOutF->Save();
	  hOutF->Close();
}


Bool_t RedAna::isAmpeak() {
	   Double_t S1 = clusters[0].charge;
	   return S1>400. && S1<800.;
}


//                      Ch: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21
const Int_t chan_group[] = {0, 2, 2, 0, 4, 0, 5, 2, 3, 3,  3,  4,  1,  3,  0,  1,  2,  3,  0,  1,  0,  2}; // group number for channels
const Int_t ngroups = 6;

Double_t par[ngroups][4]={{0.954998, 0.208693, -0.153149, -0.168458},
		{0.920263, 0.53024, -0.933853, 0.416495},
		{1.02041, -0.117352, 0.516841, -0.540358},
		{1.0067, -0.222406, 1.03052, -1.00699},
		{1.01371, -0.0177725, 0.462407, -0.551185},
		{0.992104, 0.174704, -0.128643, -0.0957555}};

TF1 *fpol_norm[ngroups];


Double_t RedAna::getS1corr(Double_t t_drift) { // This should be depends on t_drift_max. Need to be improve once field configuration is fixed.
//	Double_t t_drift_max(76.); // from run1432
	Double_t refpoint = fpol->Eval(0.5*fRunInfo.tdrift_max);
	return refpoint/fpol->Eval(t_drift);
}

Double_t RedAna::getS1corr(Double_t t_drift, Int_t channel) {
	return 1./fpol_norm[chan_group[channel]]->Eval(t_drift/fRunInfo.tdrift_max);
}

Double_t RedAna::getS13Dcorr(Double_t bar_x, Double_t bar_y, Double_t t_drift, Double_t S1_top, Double_t S1_bot){
	double x = bar_x<1.6?1.6:bar_x>3.4?3.4:bar_x;
	double y = bar_y<1.6?1.6:bar_y>3.4?3.4:bar_y;
	double iz = t_drift/7.;
	iz = iz>10?10:iz<0?0:iz;
	double corrtop, corrbot;
	if(iz<=0.5){
		corrtop = gS1map[0][0]->Interpolate(x,y);
		corrbot = gS1map[0][1]->Interpolate(x,y);
	}else if(iz<9.5){
		int iz1 = (int)(iz-0.5);
		int iz2 = iz1+1;
		double corrtop1 = gS1map[iz1][0]->Interpolate(x,y);
		double corrtop2 = gS1map[iz2][0]->Interpolate(x,y);
		double corrbot1 = gS1map[iz1][1]->Interpolate(x,y);
		double corrbot2 = gS1map[iz2][1]->Interpolate(x,y);
		corrtop = corrtop1*(iz2-iz+0.5) + corrtop2*(iz-0.5-iz1);
		corrbot = corrbot1*(iz2-iz+0.5) + corrbot2*(iz-0.5-iz1);
	}else{
		corrtop = gS1map[9][0]->Interpolate(x,y);
		corrbot = gS1map[9][1]->Interpolate(x,y);
	}
	return corrtop*S1_top+corrbot*S1_bot;
}


Double_t RedAna::getS2Zcorr(Double_t t_drift) {
	return TMath::Exp(t_drift/fRunInfo.electron_lifetime);
}
//With electron lifetime correction applied before obtaining xy correction
//1, 0.936023, 0.872466, 0.991714, 0.883636, 1.03647, 0.962989, 0.838443, 0.84197, 1.01422, 0.90065, 0.850821, 0.940773, 0.860827, 1.04908, 0.950591, 0.848482, 0.877947, 0.989487, 0.897336, 0.945046, 0.909548,

//Without electron lifetime correction applied before obtaining xy correction
//1, 0.931492, 0.868636, 0.986842, 0.879336, 1.03835, 0.958531, 0.833352, 0.838891, 1.00982, 0.896783, 0.846638, 0.936747, 0.857295, 1.04543, 0.94881, 0.844424, 0.874508, 0.988424, 0.892424, 0.940913, 0.906368,
//                      Ch: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21
const Double_t maxS2ChanCorrFac[] = {1, 0.931492, 0.868636, 0.986842, 0.879336, 1.03835, 0.958531, 0.833352, 0.838891, 1.00982, 0.896783, 0.846638, 0.936747, 0.857295, 1.04543, 0.94881, 0.844424, 0.874508, 0.988424, 0.892424, 0.940913, 0.906368};

Double_t RedAna::getS2MaxS2Chancorr(Int_t max_s2_chan) {
	return maxS2ChanCorrFac[max_s2_chan];
}


void RedAna::getTopMaxChan_and_Frac(RDCluster clstr, Int_t &max_chan, Double_t &frac) {
	Int_t nch_top = clstr.charge_top.size();
	max_chan=-1; frac=-1;
	Double_t max_charge(0.), tot_charge(0.);
	for(Int_t ich=0; ich<nch_top; ich++){
		tot_charge+=clstr.charge_top[ich];
		Double_t weighted_charge = (ich==4||ich==11)? 0.5*clstr.charge_top[ich]:clstr.charge_top[ich]; // for summed channel, divided by 2.
		if(weighted_charge> max_charge){
			max_charge=weighted_charge;
			max_chan = ich;
		}
	}
	frac = max_charge/tot_charge;
	return ;
}

Bool_t RedAna::IsInnerChannel(Int_t ch){
	return ch==0 || ch==3 || ch==5 || ch==12 || ch==14 || ch==15 || ch==18 || ch==20;
}

Int_t RedAna::GetMaxChannel(double x, double y){
//					  0  1  2  3    4  5  6  7  8  9 10   11 12 13 14 15 16 17 18 19 20 21
//    double x[22] = {1, 2, 3, 2,   0, 1, 1, 3, 0, 1, 2,   3, 2, 0, 1, 2, 3, 0, 1, 3, 2, 0};
//	  double y[22] = {3, 5, 4, 4, 4.5, 4, 5, 5, 0, 0, 0, 0.5, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3};
	int maxchannel;
	x = x*4./5.-0.5;
	y = y*6./5.-0.5;
	maxchannel = x > 2.9 ?	y > 4.5 ? 7 :
							y > 3.5 ? 2 :
							y > 2.5 ? 19:
							y > 1.5 ? 16: 11 :
				 x > 1.9 ?  y > 4.5 ? 1 :
							y > 3.5 ? 3 :
							y > 2.5 ? 20:
							y > 1.5 ? 15:
							y > 0.5 ? 12: 10 :
				 x > 0.9 ?  y > 4.5 ? 6 :
							y > 3.5 ? 5 :
							y > 2.5 ? 0 :
							y > 1.5 ? 18:
							y > 0.5 ? 14: 9  :
						    y > 4.4 ? 4 :
							y > 2.5 ? 21:
							y > 1.5 ? 17:
							y > 0.5 ? 13: 8 ;
	return maxchannel;
}


void RedAna::ReadRunConfig() {
//	fReader.SetEntry(0); // get runnumber from input tree.

	if (fOutName.Contains("1432")){
		fRunInfo.tdrift_max = 76.; // us
//		fRunInfo.para_s1_tba_corr.push_back(606.9);
//		fRunInfo.para_s1_tba_corr.push_back(3.879);
//		fRunInfo.para_s1_tba_corr.push_back(-0.0669);
		fRunInfo.para_s1_tdrift_corr.push_back(643.7);
		fRunInfo.para_s1_tdrift_corr.push_back(1.875);
		fRunInfo.para_s1_tdrift_corr.push_back(-0.02751);
		fRunInfo.electron_lifetime = 272.4;
	} else if (fOutName.Contains("1428")){
		fRunInfo.tdrift_max = 68.; // us
		fRunInfo.para_s1_tba_corr.push_back(606.9);
		fRunInfo.para_s1_tba_corr.push_back(3.879);
		fRunInfo.para_s1_tba_corr.push_back(-0.0669);
		fRunInfo.para_s1_tdrift_corr.push_back(606.9);
		fRunInfo.para_s1_tdrift_corr.push_back(3.879);
		fRunInfo.para_s1_tdrift_corr.push_back(-0.0669);
		fRunInfo.electron_lifetime = 335.8;
	} else if (fOutName.Contains("1525")){
		fRunInfo.tdrift_max = 68.; // us
		fRunInfo.para_s1_tba_corr.push_back(606.9);
		fRunInfo.para_s1_tba_corr.push_back(3.879);
		fRunInfo.para_s1_tba_corr.push_back(-0.0669);
		fRunInfo.para_s1_tdrift_corr.push_back(595.2); //From Run1525
		fRunInfo.para_s1_tdrift_corr.push_back(3.146); //From Run1525
		fRunInfo.para_s1_tdrift_corr.push_back(-0.05389); //From Run1525
		fRunInfo.electron_lifetime = 313.1; //From Run1525

	} else {
		fRunInfo.tdrift_max = 68.; // us
		fRunInfo.para_s1_tba_corr.push_back(606.9);
		fRunInfo.para_s1_tba_corr.push_back(3.879);
		fRunInfo.para_s1_tba_corr.push_back(-0.0669);
		fRunInfo.para_s1_tdrift_corr.push_back(606.9);
		fRunInfo.para_s1_tdrift_corr.push_back(3.879);
		fRunInfo.para_s1_tdrift_corr.push_back(-0.0669);
		fRunInfo.electron_lifetime = 1095;
	}
	// get run type
	if(fOutName.Contains("scan")){
		// calibration
		runtype=0;
	}else{
		// beam
		runtype=1;
	}
	//load S1 correction maps.
	for(int iz=0; iz<10; iz++){
		gS1map[iz][0] = (TGraph2D*)(TFile::Open("S1Map_fieldScan.root")->Get(Form("S1CorrMapTop_z%d",iz)))->Clone();
		gS1map[iz][1] = (TGraph2D*)(TFile::Open("S1Map_fieldScan.root")->Get(Form("S1CorrMapBot_z%d",iz)))->Clone();
		cout<< gS1map[iz][1]->GetN()<<endl;
	}

	cout<<"runtype is "<<(runtype==1 ? "beam":"calibration")<<endl;

	fpol = new TF1("fpol", "pol2(0)", 0, 100);
	fpol->SetParameters(fRunInfo.para_s1_tdrift_corr[0], fRunInfo.para_s1_tdrift_corr[1], fRunInfo.para_s1_tdrift_corr[2]); // from run1432

	for (Int_t i=0; i<ngroups; i++) {
		fpol_norm[i] = new TF1(Form("fpol_norm_%d", i), "pol3(0)", 0, 1);
		fpol_norm[i]->SetParameters(par[i]);
	}


	cout<<"t_drift_max: "<<fRunInfo.tdrift_max<<endl;
}



enum EvntCut_t {TOTAL, ACCEPTED, NCHANNEL, BASELINE, BASIC, LIVETIME, NPULSE, TRIGERWINDOW, SMALLF90, VALIDS2};
const Int_t nEventBins = 20;
const TString eventlabels[nEventBins] = {"Total", "Accepted", "Number of Channels", "Baseline Cut", "After Basic Cuts", "Short livetime", "# of pulse Cut", "Trigger Time Cut", "Small F90", "Valid S2"};

void RedAna::BookHistograms()//runtype: 0 calibration, 1 beam
{
    Info("SladRead_multihits_Selector::BookHistograms()","booking histogram...");

    TList  *list = GetOutputList();
	Int_t nentries = 210184;//fReaderForGetEntries.GetEntries(true);



    Int_t ndTTPCBins(100), ndTLSciBins(150), ndTTPCLSci(100);
//    Float_t mindTTPC(-1.05), maxdTTPC(-.8),  mindTLSci(-3.5), maxdTLSci(-3);
    Float_t mindTTPC(-50), maxdTTPC(0),  mindTLSci(-50), maxdTLSci(100); // ns
    Float_t mindTTPCLSci(0), maxdTTPCLSci(100); // ns

    Int_t nBinX(100), nBinY(100), nVarianceXY(280);
    Double_t xmin(-20.), xmax(20.), ymin(-20.), ymax(20.), varmin(80.), varmax(220.);

	Int_t nS1Bins(400), nS2Bins(100), nS2ovS1Bins(100), nF90Bins(240), nDriftT(50);//nDriftT((int)(fDriftTimeMax*1.2));
	Float_t minS1(0.), maxS1(1.5e3);
	Float_t minS2(0.), maxS2(25e+3);
//	Float_t minS2(0.), maxS2(5e+3);
	Float_t minS2ovS1(0.), maxS2ovS1(50.), minF90(0.), maxF90(1.2);
    Double_t minDriftT(0.), maxDriftT((int)(80.*1.2));
    Int_t nBaryBins(100), nTBABins(100), nS2RMSBins(200);
    Float_t minBary(1.5), maxBary(3.5), minTBA(-0.6), maxTBA(0.2), minS2RMS(1), maxS2RMS(4);

    Int_t nChanSiPMs(22);
	Float_t minChan(-0.5), maxChan(21.5);


    Int_t nEvtNumBins(nentries/1000+1);
    Float_t minEvtNum(0), maxEvtNum((nentries/1000+1)*1000);//70000);
 cout<<"nentries: "<<nentries<<" nentries/1000+1="<<nentries/1000+1<<" (nentries/1000+1)*1000="<<(nentries/1000+1)*1000<<endl;

    TString hname, htitle;
	if(runtype==1){
		hBanana = new TH2D("hBanana", "Full banana; E_{SiTel} [A.U.]; #DeltaE_{SiTel} [A.U.]", 400,1500,4000,400,400,1500);
		list->Add(hBanana);

		hCBanana = new TH2D("hCBanana", "Banana with TPC Coincidence; E_{SiTel} [A.U.]; #DeltaE_{SiTel} [A.U.]", 400,1500,4000,400,400,1500);
		list->Add(hCBanana);
    hdT_Si_TPC = new TH1D("hdT_Si_TPC","DeltaT between SiTel and TPC; dT (TPC-SiTel) [ns]", ndTTPCBins, mindTTPC, maxdTTPC);
    list->Add(hdT_Si_TPC);

    hdT_Si_LSci = new TH1D*[NLSciChan];
    hdT_TPC_LSci = new TH1D*[NLSciChan];
    hdT_Si_LSci_TPC = new TH2D*[NLSciChan];
    hLSci_PSD_Charge = new TH2D*[NLSciChan];
    hdT_Si_LSci_LSciPSD = new TH2D*[NLSciChan];
    hF90_S1_3Coinci_LSci = new TH2D*[NLSciChan];
    hS2_S1_3Coinci_LSci = new TH2D*[NLSciChan];
    for(Int_t ilsci=0; ilsci<NLSciChan; ilsci++){
        hname = Form("hdT_Si_LSci_%d", ilsci);
        htitle = Form("DeltaT between SiTel and %s at %s#circ; dT (LSci-SiTel) [ns]", LSciName[ilsci].Data(), LSciAngle[ilsci].Data());
        hdT_Si_LSci[ilsci] = new TH1D(hname, htitle, ndTLSciBins, mindTLSci, maxdTLSci);
        list->Add(hdT_Si_LSci[ilsci]);

        hname = Form("hdT_TPC_LSci_%d", ilsci);
        htitle = Form("DeltaT between SiTel and %s at %s#circ; dT (LSci-TPC) [ns]", LSciName[ilsci].Data(), LSciAngle[ilsci].Data());
        hdT_TPC_LSci[ilsci] = new TH1D(hname, htitle, ndTTPCLSci, mindTTPCLSci, maxdTTPCLSci);
        list->Add(hdT_TPC_LSci[ilsci]);

        hname = Form("hdT_Si_LSci_TPC_%d", ilsci);
        htitle = Form("DeltaT between SiTel, TPC, and %s at %s#circ; dT (TPC-SiTel) [ns]; dT (LSci-SiTel) [ns]", LSciName[ilsci].Data(), LSciAngle[ilsci].Data());
        hdT_Si_LSci_TPC[ilsci] = new TH2D(hname, htitle, ndTTPCBins, mindTTPC, maxdTTPC, ndTLSciBins, mindTLSci, maxdTLSci);
        list->Add(hdT_Si_LSci_TPC[ilsci]);

        hname = Form("hLSci_PSD_Charge_%d", ilsci);
        htitle = Form("PSD vs Charge (%s at %s#circ); Charge [A.U.]; PSD", LSciName[ilsci].Data(), LSciAngle[ilsci].Data());
        hLSci_PSD_Charge[ilsci] = new TH2D(hname, htitle, 100, 0, 5.e+5,100,0,0.5);
        list->Add(hLSci_PSD_Charge[ilsci]);

        hname = Form("hdT_Si_LSci_LSciPSD_%d", ilsci);
        htitle = Form("LSci PSD vs DeltaT between SiTel and %s at %s#circ; dT (LSci-SiTel) [ns]; PSD_{Lsci}", LSciName[ilsci].Data(), LSciAngle[ilsci].Data());
        hdT_Si_LSci_LSciPSD[ilsci] = new TH2D(hname, htitle, ndTLSciBins, mindTLSci, maxdTLSci, 100, 0, 0.5);
        list->Add(hdT_Si_LSci_LSciPSD[ilsci]);

        hname = Form("hF90_S1_3Coinci_LSci_%d", ilsci);
        htitle = Form("TPC F90 vs S1 with triple coincidence (%s at %s#circ); S1 [PE]; f90;", LSciName[ilsci].Data(), LSciAngle[ilsci].Data());
        hF90_S1_3Coinci_LSci[ilsci] = new TH2D(hname, htitle, nS1Bins, minS1, maxS1, nF90Bins, minF90, maxF90);
        list->Add(hF90_S1_3Coinci_LSci[ilsci]);

        hname = Form("hS2_S1_3Coinci_LSci_%d", ilsci);
        htitle = Form("S2 vs S1 with triple coincidence (%s at %s#circ); S1 [PE]; S2 [PE];", LSciName[ilsci].Data(), LSciAngle[ilsci].Data());
        hS2_S1_3Coinci_LSci[ilsci] = new TH2D(hname, htitle, nS1Bins, minS1, maxS1, nS2Bins, minS2, maxS2);
        list->Add(hS2_S1_3Coinci_LSci[ilsci]);

    }
	}


    hEventConter = new TH1D("hEventConter", "Event Counter", nEventBins, -0.5, nEventBins-0.5);
//    for (Int_t i=0; i<nEventBins; i++) hEventConter->GetXaxis()->SetBinLabel(i+1, eventlabels[i].Data());
    for (Int_t i=0; i<nEventBins; i++) hEventConter->GetXaxis()->SetBinLabel(i+1, Form("CX%d",i));
    hEventConter->SetFillColor(kCyan+1);
    list->Add(hEventConter);

    hTDrift = new TH1D("hTDrift", "Drift time; t_drift [#mus]", nDriftT, minDriftT, maxDriftT);
    list->Add(hTDrift);

    hF90 = new TH1D("hF90", "F90; f90", nF90Bins, minF90, maxF90);
    list->Add(hF90);

    hS1 = new TH1D("hS1", "Charge of the first pulse; clusters[0].charge [PE]", nS1Bins, minS1, maxS1);
    list->Add(hS1);

    hS1zcorr = new TH1D("hS1zcorr", "Charge of the first pulse (z corrected); S1zcorr [PE]", nS1Bins, minS1, maxS1);
    list->Add(hS1zcorr);

    hS2 = new TH1D("hS2", "Charge of the second pulse; clusters[1].charge [PE]", nS2Bins, minS2, maxS2);
    list->Add(hS2);

    hS2_zcorr = new TH1D("hS2_zcorr", "Charge of the second pulse (electron lifetime corrected); S2_zcorr [PE]", nS2Bins, minS2, maxS2);
    list->Add(hS2_zcorr);

    hS2_xycorr = new TH1D("hS2_xycorr", "S2 (x-y position dependency corrected); S2_xycorr [PE]", nS2Bins, minS2, maxS2);
    list->Add(hS2_xycorr);

    hS2_xyzcorr = new TH1D("hS2_xyzcorr", "S2 (x-y-z position dependency corrected); S2_xyzcorr [PE]", nS2Bins, minS2, maxS2);
    list->Add(hS2_xyzcorr);

    hNclusters = new TH1D("hNclusters", "Number of clusters in each event; number_of_clusters", 11, -0.5, 10.5);
    list->Add(hNclusters);

    hTBAS1 = new TH1D("hTBAS1", "TBA of S1; (top-bot)/S1", nTBABins, minTBA, maxTBA);
    list->Add(hTBAS1);

    hTBAS1wS2s = new TH1D("hTBAS1wS2s", "TBA of S1 from more than two pulses events; (top-bot)/S1", nTBABins, minTBA, maxTBA);
    list->Add(hTBAS1wS2s);

    hTBAS1All = new TH1D("hTBAS1All", "TBA of S1 for all events; (top-bot)/S1", nTBABins, minTBA, maxTBA);
    list->Add(hTBAS1All);

    hTBAS2 = new TH1D("hTBAS2", "TBA of S2; (top-bot)/S2", nTBABins, minTBA, maxTBA);
    list->Add(hTBAS2);

    hS2_efficiency = new TH1D("hS2_efficiency", "# of events w more than 2 pulses / # of nonzero pulse events  ; (top-bot)/S1", nTBABins, minTBA, maxTBA);
    list->Add(hS2_efficiency);

    hS2MaxChan = new TH1D("hS2MaxChan", "Max Channel for S2 (only top channels); s2_max_chan; Events", 32, -0.5, 31.5);
    list->Add(hS2MaxChan);

    hF90_S1 = new TH2D("hF90_S1", "F90 vs S1; S1 [PE]; f90", nS1Bins, minS1, maxS1, nF90Bins, minF90, maxF90);
    list->Add(hF90_S1);

    hF90_S1_SiTelCoinci = new TH2D("hF90_S1_SiTelCoinci", "TPC F90 vs S1 with SiTel coincidence; S1 [PE]; f90", nS1Bins, minS1, maxS1, nF90Bins, minF90, maxF90);
    list->Add(hF90_S1_SiTelCoinci);

    hS1_tdrift = new TH2D("hS1_tdrift", "S1 vs t_drift; t_drift [#mus]; S1 [PE]", nDriftT, minDriftT, maxDriftT, nS1Bins, minS1, maxS1);
    list->Add(hS1_tdrift);

    hS1zcorr_tdrift = new TH2D("hS1zcorr_tdrift", "S1 (z corrected) vs t_drift; t_drift [#mus]; S1zcorr [PE]", nDriftT, minDriftT, maxDriftT, nS1Bins, minS1, maxS1);
    list->Add(hS1zcorr_tdrift);
    hS1Topcorr_tdrift = new TH2D("hS1Topcorr_tdrift", "S1 (top tile corrected) vs t_drift; t_drift [#mus]; S1topcorr [PE]", nDriftT, minDriftT, maxDriftT, nS1Bins, minS1, maxS1);
    list->Add(hS1Topcorr_tdrift);
    hS13Dcorr_tdrift = new TH2D("hS13Dcorr_tdrift", "S1 (3D corrected) vs t_drift; t_drift [#mus]; S13Dcorr [PE]", nDriftT, minDriftT, maxDriftT, nS1Bins, minS1, maxS1);
    list->Add(hS13Dcorr_tdrift);

    hS2_tdrift = new TH2D("hS2_tdrift", "S2 vs t_drift; t_drift [#mus]; S2 [PE]", nDriftT, minDriftT, maxDriftT, nS2Bins, minS2, maxS2);
    list->Add(hS2_tdrift);
    hS2Topcorr_tdrift = new TH2D("hS2Topcorr_tdrift", "S2 vs t_drift; t_drift [#mus]; S2 [PE]", nDriftT, minDriftT, maxDriftT, nS2Bins, minS2, maxS2);
    list->Add(hS2Topcorr_tdrift);

    hS2_Am_tdrift = new TH2D("hS2_Am_tdrift", "S2 from Am peak vs t_drift; t_drift [#mus]; S2 [PE]", nDriftT, minDriftT, maxDriftT, nS2Bins, minS2, maxS2);
    list->Add(hS2_Am_tdrift);

    hS2_Am_s2maxchan = new TH2D("hS2_Am_s2maxchan", "S2 from Am peak vs max S2 channel; S2 max chan; S2 [PE]", nChanSiPMs, minChan, maxChan, nS2Bins, minS2, maxS2);
    list->Add(hS2_Am_s2maxchan);

    hS2_S1 = new TH2D("hS2_S1", "S2 vs S1; S1 [PE]; S2 [PE]", nS1Bins, minS1, maxS1, nS2Bins, minS2, maxS2);
    list->Add(hS2_S1);

    hS2_S1_SiTelCoinci = new TH2D("hS2_S1_SiTelCoinci", "S2 vs S1 with SiTel coincidence; S1 [PE]; S2 [PE]", nS1Bins, minS1, maxS1, nS2Bins, minS2, maxS2);
    list->Add(hS2_S1_SiTelCoinci);

    hS2ovS1_tdrift = new TH2D("hS2ovS1_tdrift", "S2/S1 vs t_drift; t_drift [#mus]; S2/S1", nDriftT, minDriftT, maxDriftT, nS2ovS1Bins, minS2ovS1, maxS2ovS1);
    list->Add(hS2ovS1_tdrift);

    hNclusters_evtnum = new TH2D("hNclusters_evtnum", "Number of clusters in each event vs event number; event_number; number_of_clusters", nEvtNumBins, minEvtNum, maxEvtNum, 11, -0.5, 10.5);
    list->Add(hNclusters_evtnum);

    htdrift_evtnum = new TH2D("htdrift_evtnum", "t_drift vs event number; event_number; t_drift", nEvtNumBins, minEvtNum, maxEvtNum, nDriftT, minDriftT, maxDriftT);
    list->Add(htdrift_evtnum);

    hS2RMS_evtnum = new TH2D("hS2RMS_evtnum", "S2 Waveform RMS vs event number; event_number; S2 rms_time [#mus]", nEvtNumBins, minEvtNum, maxEvtNum, nS2RMSBins, minS2RMS, maxS2RMS);
    list->Add(hS2RMS_evtnum);

    hS2TR_evtnum = new TH2D("hS2TR_evtnum", "S2 Waveform rise time vs event number; event_number; S2 rise_time [#mus]", nEvtNumBins, minEvtNum, maxEvtNum, nS2RMSBins, minS2RMS, maxS2RMS);
    list->Add(hS2TR_evtnum);

    hS2ovS1_evtnum = new TH2D("hS2ovS1_evtnum", "S2/S1 vs event number; event_number; S2/S1", nEvtNumBins, minEvtNum, maxEvtNum, nS2ovS1Bins, minS2ovS1, maxS2ovS1);
    list->Add(hS2ovS1_evtnum);

    hYBary_XBary = new TH2D("hYBary_XBary", "Barycenter Y vs Barycenter X; bar_x; bar_y", nBaryBins, minBary, maxBary, nBaryBins, minBary, maxBary);
    list->Add(hYBary_XBary);
    hYBary_XBary_corr = new TH2D("hYBary_XBary_corr", "Barycenter Y vs Barycenter X CorrS1; bar_x; bar_y", nBaryBins, minBary, maxBary, nBaryBins, minBary, maxBary);
    list->Add(hYBary_XBary_corr);
	gMYBary_XBary_corr = new TMultiGraph("gXY_by_ch","XY_by_ch");
	for(int i=0;i<22; i++){
		gYBary_XBary_corr[i] = new TGraph();
		gYBary_XBary_corr[i]->SetNameTitle(Form("gxy_ch%d",i),Form("gxy_ch%d",i));
		gMYBary_XBary_corr->Add(gYBary_XBary_corr[i],"p");
	}
	list->Add(gMYBary_XBary_corr);

    hXBary_tdrift = new TH2D("hXBary_tdrift", "Barycenter X vs t_drift; t_drift; bar_x", nDriftT, minDriftT, maxDriftT, nBaryBins, minBary, maxBary);
    list->Add(hXBary_tdrift);

    hYBary_tdrift = new TH2D("hYBary_tdrift", "Barycenter Y vs t_drift; t_drift; bar_y", nDriftT, minDriftT, maxDriftT, nBaryBins, minBary, maxBary);
    list->Add(hYBary_tdrift);

    hTBAS1_tdrift = new TH2D("hTBAS1_tdrift", "TBA of S1 vs t_drift; t_drift [#mus];  (top-bot)/S1", nDriftT, minDriftT, maxDriftT, nTBABins, minTBA, maxTBA);
    list->Add(hTBAS1_tdrift);

    hS1_TBAS1 = new TH2D("hS1_TBAS1", "S1 vs TBA of S1; (top-bot)/S1; S1 [PE]", nTBABins, minTBA, maxTBA, nS1Bins, minS1, maxS1);
    list->Add(hS1_TBAS1);

    hS1zcorr_TBAS1 = new TH2D("hS1zcorr_TBAS1", "S1 (z corrected) vs TBA of S1; (top-bot)/S1; S1zcorr [PE]", nTBABins, minTBA, maxTBA, nS1Bins, minS1, maxS1);
    list->Add(hS1zcorr_TBAS1);

    hS2_TBAS2 = new TH2D("hS2_TBAS2", "S2 vs TBA of S2; (top-bot)/S2; S2 [PE]", nTBABins, minTBA, maxTBA, nS2Bins, minS2, maxS2);
    list->Add(hS2_TBAS2);

    hdT_Si_TPC_TPCF90 = new TH2D("hdT_Si_TPC_TPCF90","F90_{TPC} vs DeltaT between SiTel and TPC; dT (TPC-SiTel) [ns]; F90_{TPC}", ndTTPCBins, mindTTPC, maxdTTPC, 400, 0, 1);
    list->Add(hdT_Si_TPC_TPCF90);

    hStartTime_ChNum = new TH2D("hStartTime_ChNum","Pulse start_time vs Global channel number; channel; start_time [sample]", 32, -0.5, 31.5, 200, 4400, 5400);
//    hStartTime_ChNum = new TH2D("hStartTime_ChNum","Pulse start_time vs Global channel number; channel; start_time [sample]", 32, -0.5, 31.5, 200, 4800, 5400);
    list->Add(hStartTime_ChNum);

    hYBary_XBary_s2maxchan = new TH3D("hYBary_XBary_s2maxchan", "Barycenter Y vs Barycenter X vs max S2 channel; S2 max chan; bar_x; bar_y; ", nChanSiPMs, minChan, maxChan, nBaryBins, minBary, maxBary, nBaryBins, minBary, maxBary);
    list->Add(hYBary_XBary_s2maxchan);

    hS2_frac_chnum = new TH2D("hS2_frac_chnum", "S2 fraction vs channel; ch number; S2_ch/S2_tot; ", nChanSiPMs, minChan, maxChan, 100, 0, 0.5);
    list->Add(hS2_frac_chnum);

    hS1_tdrift_s2maxchan = new TH3D("hS1_tdrift_s2maxchan", "S1 vs t_drift vs max S2 channel; S2 max chan; t_drift [#mus]; S1 [PE]", nChanSiPMs, minChan, maxChan, nDriftT, minDriftT, maxDriftT, nS1Bins, minS1, maxS1);
    list->Add(hS1_tdrift_s2maxchan);

	//channel efficiency correction map generator
	hChTopS1Frac = new TH2D("hChTopS1Frac","Channel light faction from uniform S1 in the bottom half of TPC; x [cm]; y[cm]", 4,0,5,6,0,5);
	list->Add(hChTopS1Frac);

	if(runtype==0){
		//hist for Am241 S1 3D correction map. 
		for(int iz=0; iz<10; iz++){
			hS1_barXY[iz] = new TH3D(Form("hS1_barXY_z%d",iz),Form("S1 distribution across XY in slice z[%0.1lf,%0.1lf] [us]",iz*0.7,iz*0.7+1),15,1.5,3.5,15,1.5,3.5,100,400,800);
			hS1Top_barXY[iz] = new TH3D(Form("hS1Top_barXY_z%d",iz),Form("S1 top distribution across XY in slice z[%0.1lf,%0.1lf] [us]",iz*0.7,iz*0.7+1),15,1.5,3.5,15,1.5,3.5,100,100,500);
			hS1Bottom_barXY[iz] = new TH3D(Form("hS1Bottom_barXY_z%d",iz),Form("S1 bottom distribution across XY in slice z[%0.1lf,%0.1lf] [us]",iz*0.7,iz*0.7+1),15,1.5,3.5,15,1.5,3.5,150,200,800);
			list->Add(hS1_barXY[iz]);
			list->Add(hS1Top_barXY[iz]);
			list->Add(hS1Bottom_barXY[iz]);
		}
	}

}

void RedAna::TopChCorr(vector<double> top, vector<double> bot, double &xx, double &yy, double &S1corr){
	if(top.size()!=22||bot.size()!=2){cout<<"channel not matching! exit!"<<endl; return;}
	//double x[20] = {1, 0, 2, 3, 2, 1, 0, 1, 2,  3, 2, 0, 1, 2,   0, 3, 1, 3, 2, 1};
	//double y[20] = {5, 3, 5, 4, 4, 4, 2, 0, 0,0.5, 1, 1, 1, 2, 4.5, 2, 2, 3, 3, 3};
    //double xch[22] = {1, 2, 3, 2,   0, 1, 1, 3, 0, 1, 2,   3, 2, 0, 1, 2, 3, 0, 1, 3, 2, 0};
	//double ych[22] = {3, 5, 4, 4, 4.5, 4, 5, 5, 0, 0, 0, 0.5, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3};
	//  xch[i] = (xch[i]+0.5)*5./4.;
	//  ych[i] = (ych[i]+0.5)*5./6.;

	//corr from run1609
	//					D1		D2		D3		D4		D5E3	E2		E4		E5		A1		A2		A3		A4B4	A5		B1		B2		B3		B5		C1		C2		C3		C4		C5
	//double corr[22] = {0.8841, 1.02, 1.0166, 0.8942, 1.3159, 0.8754, 1.0519, 1.1080, 1.1766, 1.0671, 1.1347, 1.4512, 0.9404, 1.0182, 0.9022, 0.9034, 1.0385, 1.0363, 0.9054, 1.0336, 0.9116, 1.0194};
	// corr from bkg, uniform S1 80events/max ch. 200-800peS1 >35us. 
	//double corr[22] =       {1./1.18694, 1./0.98288, 1./0.971182, 1./1.17896, 2./1.83117,
    //                         1./1.17237, 1./1.01723, 1./0.827953, 1./0.801992, 1./0.95377,
    //                         1./0.886194,2./1.66258, 1./1.07892 , 1./0.982774, 1./1.12815,
    //                         1./1.17234, 1./0.953743,1./0.966676, 1./1.16131 , 1./0.958728,
    //                         1./1.12457, 1./0.999582};

	// corr from bkg, uniform S1 no max ch selection. 20-800peS1 >35us.
	double corr[22] =       {1./1.19552, 1./0.97908, 1./0.97473 , 1./1.17607 , 2./1.82954,
                             1./1.18200, 1./1.00148, 1./0.823792, 1./0.811283, 1./0.943718,
                             1./0.883644,2./1.64743, 1./1.07917 , 1./0.98764 , 1./1.12608,
                             1./1.17155, 1./0.959034,1./0.978697, 1./1.16663 , 1./0.960965,
                             1./1.13294, 1./0.989004};


	double meancorr = 0;
	for(int ich=0; ich<22; ich++) meancorr+=corr[ich];
	meancorr = 22./meancorr;
	for(int ich=0; ich<22; ich++) corr[ich]*=meancorr;

	double tot_top = 0;
	for(int i=0; i<22; i++){
	  tot_top+=top[i];
	}
	double x_bar=0, y_bar=0, tot = 0;
	for(int id=0; id<22; id++){
		x_bar += xch[id]*corr[id]*top[id];
		y_bar += ych[id]*corr[id]*top[id];
		tot += corr[id]*top[id];
	}
	xx = x_bar/tot;
	yy = y_bar/tot;
	S1corr = tot + bot[0]+bot[1];

}
