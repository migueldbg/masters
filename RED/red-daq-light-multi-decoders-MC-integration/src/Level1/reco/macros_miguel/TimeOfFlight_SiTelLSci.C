#include "sidutility.cc"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

/* ************************************************************************************************************************* *
 * File: TimeOfFlight_SiTelTPC.C (ROOT macro).                                                                               *
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

/* THStack* TimeOfFlight( Int_t run, TTree* reco, Int_t chanID, TString stackTitle = "stack" )
 *
 *  Summary of Function:
 *
 *    For a given TTree and a given channel number, this function calculates the time of flight between
 *    the SiTel detector and the given channel. Taking the appropriate channel numbers, the function
 *    calculates the ToF between the SiTel and a LSci. Three histograms are generated, each with
 *    progressively stricter cuts. The first has no cuts applie, the second considers tpc signal quality
 *    cuts and the third also selects the low energy berilium blob. The function then returns a THStack
 *    object containing all three histograms.
 *
 *  Parameters   :  run        >> the run containing the desired data
 *                  reco       >> the TTree containing the desired data
 *                  chanID     >> then number indicating which channel to consider
 *                  stackTitle >> the title of the stack generated
 *
 *  Return Value : THStack* tofStack
 */
THStack* ToFCutAnalysis( Int_t run, TTree* reco, Int_t chanID, TString stackTitle = "stack" ){
  /*
  TStyle* sidStyle = SetSidStyle();   sidStyle -> cd();

  TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);
  */

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

  return tofStack;
}


void AllLSciToFCutAnalysis(Int_t run){

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
    tofStacks[i] = ToFCutAnalysis(run, reco, chanIDs[i], chanNames[i]);

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


void LSciPSD( Int_t run, Int_t chanID ){

  TStyle* sidStyle = SetSidStyle();   sidStyle -> cd();

  TString file_name = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* file = CheckFile(file_name);
  TTree* reco;  file -> GetObject("reco", reco);


  Double_t psdMin       = 0;   Double_t psdMax       = 1.0;
  //Double_t psdMinZ      = 0;   Double_t psdMaxZ      = 1.0;
  Double_t lsChargeMin  = 0;   Double_t lsChargeMax  = 200e3;
  Double_t lsChargeMinZ = 0;   Double_t lsChargeMaxZ = 2e3;

  Double_t psdBinSize       = 0.002;  Int_t psdBinNumber       = (psdMax - psdMin)/psdBinSize;
  Double_t lsChargeBinSize  = 200;    Int_t lsChargeBinNumber  = (lsChargeMax - lsChargeMin)/lsChargeBinSize;
  Double_t lsChargeBinSizeZ = 10;     Int_t lsChargeBinNumberZ = (lsChargeMaxZ - lsChargeMinZ)/lsChargeBinSizeZ;

  TCut lsChargeCut  = Form("charge[%d] >= %f && charge[%d] <= %f", chanID, lsChargeMin,  chanID, lsChargeMax);
  TCut lsChargeCutZ = Form("charge[%d] >= %f && charge[%d] <= %f", chanID, lsChargeMinZ, chanID, lsChargeMaxZ);
  TCut psdCut = Form("f90[%d] >= %f && f90[%d] <= %f", chanID, psdMin, chanID, psdMax);

  TH2F* PSDvCharge = new TH2F("PSDvCharge", "PSD vs Charge (LSci 0)",
                              lsChargeBinNumber, lsChargeMin, lsChargeMax, psdBinNumber, psdMin, psdMax);
  TH2F* PSDvChargeZ = new TH2F("PSDvChargeZ", "PSD vs Charge (LSci 0)",
                               lsChargeBinNumberZ, lsChargeMinZ, lsChargeMaxZ, psdBinNumber, psdMin, psdMax);

  PSDvCharge  -> GetXaxis() -> SetTitle("Charge");    PSDvCharge  -> GetYaxis() -> SetTitle("PSD");
  PSDvChargeZ -> GetXaxis() -> SetTitle("Charge");    PSDvChargeZ -> GetYaxis() -> SetTitle("PSD");

  reco -> Project("PSDvCharge",  Form("f90[%d]:charge[%d]", chanID, chanID), lsChargeCut  && psdCut);
  reco -> Project("PSDvChargeZ", Form("f90[%d]:charge[%d]", chanID, chanID), lsChargeCutZ && psdCut);


  Double_t height = 500; Double_t width = gdRatio * height;
  TCanvas* canvas = new TCanvas("canvas", "canvas", 2*width, height);
  canvas -> Divide(2,1);

  canvas -> cd(1);  PSDvCharge  -> Draw("HIST COLZ");
  canvas -> cd(2);  PSDvChargeZ -> Draw("HIST COLZ");
}
