#include "sidutility.cc"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

/* ************************************************************************************************************************* *
 * File: TimeOfFlight_SiTelLSci.C (ROOT macro).                                                                               *
 *                                                                                                                           *
 * Author: Miguel Del Ben Galdiano (miguel.galdiano@gmail.com).                                                              *
 * Date of creation : June 28 2020.                                                                                          *
 *                                                                                                                           *
 * Summary of File:                                                                                                          *
 *                                                                                                                           *
 *    This macro should contain all functions that relate to the time of flight between the SiTel and LSci. The purpouse of  *
 *    each function will be properly explained in a comment located right above the function. Some concepts are used by all  *
 *    functions, and as such will be explained here.                                                                         *
 *                                                                                                                           *
 *    > TIME OF FLIGHT PARAMETER >                                                                                           *
 *      Definition: 2*(start_time[] - 0.5*(start_time[30] + start_time[31] - 7.45)).                                         *
 *                                                                                                                           *
 *      The parameter start_time[] is used instead of xmin[] because the later is quantized in 2ns steps, while the former   *
 *      is continuous. This allows for a better histogram. Also, instead of using only the timing of the thinner silicon     *
 *      detector (start_time[30]), we used a modified avarage of the thinner and larger silicon detectors. The value         *
 *      subtracted from start_time[31] is obtained from constructing a histogram of the difference between the two           *
 *      start_time[] values and fitting it with a gaussian function and taking its mean. This subtraction is made to match   *
 *      the timing of the larger silicon detector, that happens later, with that of the smaller silicon detector.            *
 *                                                                                                                           *
 *    > LSci CHANNEL MAPPING >                                                                                               *
 *      LSci:       0    1    3 & 7    4 & 6    5    8                                                                       *
 *      Channel:    0    1      3        4      5    8                                                                       *
 * ************************************************************************************************************************* */

using namespace std;

TCut DefineLSciPSDCut(Int_t chanID, Double_t psdMin, Double_t psdMax);
TCut DefineLSciChargeCut(Int_t chanID, Double_t chargeMin, Double_t chargeMax);


/* THStack* CutAnalysis( Int_t run, Int_t chanID, TString stackTitle = "SiTel-LSci ToF", TTree* reco = NULL, bool draw = true )
 *
 *  Summary of Function:
 *
 *    For a given TTree and a given channel number, this function calculates the time of flight between
 *    the SiTel detector and the given channel. Taking the appropriate channel numbers, the function
 *    calculates the ToF between the SiTel and a LSci. Three histograms are generated, each with
 *    progressively stricter cuts. The first has no cuts applied, the second considers TPC signal quality
 *    cuts and the third also selects the low energy berilium blob. The function then returns a THStack
 *    object containing all three histograms.
 *
 *  Parameters   :  run        >> the run containing the desired data.
 *                  chanID     >> then number indicating which channel to consider
 *                  stackTitle >> the title of the stack generated
 *                  reco       >> the TTree containing the desired data. If null, the code finds the TTree acording
 *                                to the run.
 *                  draw       >> a bool that tells the code wether the resulting THStack is to be drawn or not.
 *
 *  Return Value :  THStack* tofStack
 */
THStack* CutAnalysis( Int_t run, Int_t chanID, TString stackTitle = "SiTel-LSci ToF", TTree* reco = NULL, bool draw = true ){

  if (reco == NULL){
    TStyle* sidStyle = SetSidStyle();   sidStyle -> cd();

    TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
    TFile* file = CheckFile(file_name);
    reco;  file -> GetObject("reco", reco);
  }

  Int_t cfg = 2;
  Double_t s1Min  = 50;      Double_t s1Max  = 1000;
  Double_t f90Min = 0.1;     Double_t f90Max = 1.0;
  Double_t tofMin = -700;    Double_t tofMax = 100 ;

  Double_t tofBinSize = 0.5;
  Double_t tofBinNumber = (tofMax - tofMin)/tofBinSize;


  TCut   tpcCut   = DefineCuts(cfg, f90Min, f90Max, s1Min, s1Max);
  TCutG* lowBeCut = LowBeGraphCut(run, "LowBeCut");

  TH1F*       tofHist;
  const char* histNames[3] = {"tpc+lowBeCut"       , "tpcCut", "noCut" };
  TCut        cuts[3]      = { tpcCut && "LowBeCut",  tpcCut , "" };
  EColor      colors[3]    = { kRed                ,  kBlue  , kBlack};

  TString stackName = Form("tofStack_%d", chanID);

  THStack* tofStack = new THStack(stackName, stackTitle);


  for (Int_t i = 0; i < 3; i++){

    tofHist = new TH1F(histNames[i], "Time of Flight (SiTel-LSci)", tofBinNumber, tofMin, tofMax);
    reco -> Project(histNames[i], Form("2*(start_time[%d] - 0.5*(start_time[30] + start_time[31] - 7.45))", chanID), cuts[i]);

    tofHist -> SetLineColor(colors[i]);
    tofHist -> SetLineWidth(2);

    tofStack -> Add(tofHist);

  }

  if (draw){
    TCanvas* canvas = new TCanvas("canvas", "canvas", 400*gdRatio, 400);
    tofStack -> Draw("nostack");
  }

  return tofStack;
}

/* void CutAnalysisAllChan( Int_t run )
 *
 *  Summary of Function:
 *
 *    This function constructs the SiTel-LSci ToF histogram for each LSci and returns a canvas with all six
 *    plots. Important to note: this function assumes a specific channel mapping (Feb2020_corr_channelmapping).
 *    and the chanIDs and chanNames variables are defined with that in mind.
 *
 *  Parameters   :  run >> the run containing the desired data.
 *
 *  Return Value :  void
 */
void CutAnalysisAllChan( Int_t run ){

  TStyle* sidStyle = SetSidStyle();   sidStyle -> cd();

  TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);

  const Int_t chanAmount = 6;
  Int_t chanIDs[chanAmount] = {0, 1, 3, 4, 5, 8};
  TString chanNames[chanAmount] = {"LSci 0", "LSci 1", "LSci 3 & 7", "LSci 4 & 6", "LSci 5", "LSci 8"};

  THStack* tofStacks[chanAmount];

  Double_t height = 250;    Double_t width = gdRatio*height;
  TCanvas* canvas = new TCanvas("canvas", "ToF (SiTel-LSci)", 2*width, 3*height);
  canvas -> Divide(2,3);

  for (Int_t i = 0; i < chanAmount; i++){

    cout << chanIDs[i] << endl;
    tofStacks[i] = CutAnalysis(run, chanIDs[i], chanNames[i], reco, false);

    canvas -> cd(i+1);
    tofStacks[i] -> Draw("NOSTACK");
    tofStacks[i] -> GetXaxis() -> SetTitle("ToF (SiTel-LSci) [ns]");
    tofStacks[i] -> GetYaxis() -> SetTitle("Counts/0.5 ns");
  }

  canvas -> Modified();

  canvas -> cd(1);

  TLegend* legend = new TLegend(0.2, 0.2);
  TString  legendNames[3] = {"No Cut", "TPC Cut", "TPC + Low Be Cut"};

  TList* histList = tofStacks[0] -> GetHists();
  TH1F*  hist[histList -> GetEntries()];

  for (Int_t i = 0; i < histList -> GetEntries(); i++){
    hist[i] = (TH1F*) histList -> At(histList -> GetEntries() - (i+1));
    legend -> AddEntry(hist[i], legendNames[i], "l");
  }

  legend -> Draw();

}



void SiTelLSciToFNeutron( Int_t run, Int_t chanID ){

  TStyle* sidStyle = SetSidStyle();   sidStyle -> cd();

  TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);


  Double_t psdMin      = 0.0;   Double_t psdMax      = 1.0;
  Double_t lsChargeMin = 0;     Double_t lsChargeMax = 200e3;
  Double_t tofMin      = -200;  Double_t tofMax      = 800 ;

  Double_t tofBinSize = 0.5;
  Double_t tofBinNumber = (tofMax - tofMin)/tofBinSize;


  TH1F* siTelLSciTof[3];
  const char* histNames[3] = {"noCut", "neutronCut", "gammaCut"};
  TCut cuts[3] = {"", Form("lsci_psd_tot[%d] >= 0.15 && lsci_psd_tot[%d] <= 1.0 && charge[%d] > 400", chanID, chanID, chanID),
                      Form("lsci_psd_tot[%d] >= 0.0 && lsci_psd_tot[%d] <= 0.15 && charge[%d] > 400", chanID, chanID, chanID)};
  EColor colors[3] = {kBlack, kBlue, kRed};

  for (Int_t i = 0; i < 3; i++){

    siTelLSciTof[i] = new TH1F(histNames[i], "Time of Flight (SiTel-LSci)", tofBinNumber, tofMin, tofMax);
    reco -> Project(histNames[i], Form("2*(0.5*(start_time[30] + start_time[31] - 7.45) - start_time[%d])", chanID), cuts[i]);

    siTelLSciTof[i] -> SetLineColor(colors[i]);
    siTelLSciTof[i] -> SetLineWidth(2);

    if (i == 0) {
      siTelLSciTof[i] -> GetXaxis() -> SetTitle("ToF (SiTel-LSci) [ns]");
      siTelLSciTof[i] -> GetXaxis() -> SetTitle("Counts/0.5 ns");
      siTelLSciTof[i] -> Draw("HIST");
    } else {
      siTelLSciTof[i] -> Draw("HIST SAME");
    }
  }

  TLegend*legend = new TLegend(0.85,0.74,0.95,0.95);
  legend -> SetTextFont(102);
  legend -> SetTextSize(0.03);
  legend -> AddEntry(siTelLSciTof[0], "No Cuts", "l");
  legend -> AddEntry(siTelLSciTof[1], "NR Selection", "l");
  legend -> AddEntry(siTelLSciTof[2], "ER Selection", "l");
  legend -> Draw();

}


void LSciPSDvChargeStudy( Int_t run, Int_t chanID ){

  TStyle* sidStyle = SetSidStyle();   sidStyle -> cd();

  TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);


  Double_t psdMin     = 0;   Double_t psdMax     = 1.0;
  Double_t chargeMin  = 0;   Double_t chargeMax  = 400e3;
  Double_t chargeMinZ = 0;   Double_t chargeMaxZ = 2e3;   //Z means 'zoomed'

  Double_t psdBinSize     = 2E-3;   Int_t psdBinNumber     = (psdMax - psdMin)/psdBinSize;
  Double_t chargeBinSize  = 4E3;    Int_t chargeBinNumber  = (chargeMax - chargeMin)/chargeBinSize;
  Double_t chargeBinSizeZ = 10;     Int_t chargeBinNumberZ = (chargeMaxZ - chargeMinZ)/chargeBinSizeZ;

  TCut chargeCut  = DefineLSciChargeCut(chanID, chargeMin, chargeMax);
  TCut chargeCutZ = DefineLSciChargeCut(chanID, chargeMinZ, chargeMaxZ);
  TCut psdCut = DefineLSciPSDCut(chanID, psdMin, psdMax);


  TH2F* PSDvCharge  = new TH2F("psd_charge" , "PSD vs Charge (LSci 0)", chargeBinNumber, chargeMin, chargeMax,
                                                                        psdBinNumber, psdMin, psdMax);
  TH2F* PSDvChargeZ = new TH2F("psd_chargeZ", "PSD vs Charge (LSci 0)", chargeBinNumberZ, chargeMinZ, chargeMaxZ,
                                                                        psdBinNumber, psdMin, psdMax);

  PSDvCharge  -> GetXaxis() -> SetTitle("Charge [PE]");    PSDvCharge  -> GetYaxis() -> SetTitle("PSD");
  PSDvChargeZ -> GetXaxis() -> SetTitle("Charge [PE]");    PSDvChargeZ -> GetYaxis() -> SetTitle("PSD");

  reco -> Project("psd_charge" , Form("f90[%d]:charge[%d]", chanID, chanID), chargeCut  && psdCut);
  reco -> Project("psd_chargeZ", Form("f90[%d]:charge[%d]", chanID, chanID), chargeCutZ && psdCut);


  Double_t height = 500; Double_t width = gdRatio * height;
  TCanvas* canvas = new TCanvas("canvas", "canvas", 2*width, height);
  canvas -> Divide(2,1);

  canvas -> cd(1);  PSDvCharge  -> Draw("HIST COLZ");
  canvas -> cd(2);  PSDvChargeZ -> Draw("HIST COLZ");
}

TH2* PSDvCharge( Int_t run, Int_t chanID, TString histTitle = "PSD v Charge", TTree* reco = NULL, bool draw = true ){

  if (reco == NULL){
    TStyle* sidStyle = SetSidStyle();   sidStyle -> cd();
    sidStyle -> SetOptStat(0);

    TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
    TFile* file = CheckFile(file_name);
    reco;  file -> GetObject("reco", reco);
  }


  Double_t psdMin     = 0.;   Double_t psdMax     = 0.5;
  Double_t chargeMin  = 0.;   Double_t chargeMax  = 400e3;

  Double_t psdBinSize     = 4E-3;   Int_t psdBinNumber     = (psdMax - psdMin)/psdBinSize;
  Double_t chargeBinSize  = 2E3;    Int_t chargeBinNumber  = (chargeMax - chargeMin)/chargeBinSize;

  TCut chargeCut  = DefineLSciChargeCut(chanID, chargeMin, chargeMax);
  TCut psdCut     = DefineLSciPSDCut(chanID, psdMin, psdMax);


  TH2F* PSDvCharge  = new TH2F(Form("psd_charge_%d", chanID), histTitle, chargeBinNumber, chargeMin, chargeMax,
                                                                         psdBinNumber, psdMin, psdMax);

  PSDvCharge -> GetXaxis() -> SetTitle("Charge [PE]");
  PSDvCharge -> GetYaxis() -> SetTitle("PSD");

  reco -> Project(Form("psd_charge_%d", chanID), Form("f90[%d]:charge[%d]", chanID, chanID), chargeCut  && psdCut);


  if (draw){
    TCanvas* canvas = new TCanvas("canvas", "canvas", gdRatio*500, 500);
    canvas -> SetLogz();
    PSDvCharge  -> Draw("HIST COLZ");
  }

  return PSDvCharge;
}

void PSDvChargeAllChan( Int_t run ){
  
}


TCut DefineLSciPSDCut(Int_t chanID, Double_t psdMin, Double_t psdMax){

  TCut psdMinCut = Form("f90[%d] >= %f", chanID, psdMin);
  TCut psdMaxCut = Form("f90[%d] <= %f", chanID, psdMax);
  TCut psdCut = psdMinCut && psdMaxCut;

  return psdCut;
}

TCut DefineLSciChargeCut(Int_t chanID, Double_t chargeMin, Double_t chargeMax){

  TCut chargeMinCut = Form("charge[%d] >= %f", chanID, chargeMin);
  TCut chargeMaxCut = Form("charge[%d] <= %f", chanID, chargeMax);
  TCut chargeCut = chargeMinCut && chargeMaxCut;

  return chargeCut;
}
