#include "sidutility.cc"

#include "TFile.h"
#include "TH1.h"
#include "TStyle.h"
#include "TTree.h"

using namespace std;
using namespace RooFit;



TCut SiTelTPCToFCut(Double_t tofMin, Double_t tofMax);

// GOAL: Learning the basics of RooFit and how to apply it.

void Gaussian(){

  RooRealVar x("x", "x", -10, 10);
  RooRealVar mean("mean", "Mean of Gaussian", 0, -10, 10);
  RooRealVar sigma("sigma", "Width of Gaussian", 3, -10, 10);

  /* Declaration of the variables that will be used in the model.
   *
   * x: the independet variable, allowed to float within the range (-10 to 10). Its initial value is set to the center of the range.
   * mean: the mean of the gaussian, also allowed to float, but with initial value explicitely defined to be '0'.
   * sigma: the sigma of the gaussian, also allowed to float, but with initial value explicitely defined to be '3'
   */

  RooGaussian gauss("gauss", "Gaussian (x, mean, sigma)", x, mean, sigma);

  // Initialization of the gaussian PDF using the previously defined variable.

  RooPlot* frame = x.frame();
  // Creation of a empty frame with the 'x' variable along the x-axis.

  gauss.plotOn(frame);
  // Tells the gaussian model to plot itself in the frame created previously.

  frame -> Draw();
  // Draws the frame, with the gaussian model included, in a TCanvas object.


  sigma = 2;
  // This changes the initial value of the sigma parameter.

  gauss.plotOn(frame, LineColor(kRed));
  frame -> Draw();
  // Tells the model (with new sigma value) to draw itself in the frame. Also has a new parameter, changing the color of the curve.
  // Note: the frame still contains the previous plot of the model with sigma=3, even after the model is changed.
}

void Histogram( Int_t run ){

  // <code-fold> PRELIMINARY STEPS FOR HISTOGRAM GENERATION FROM DATA
  TStyle* sidStyle = SetSidStyle();
  sidStyle -> cd();
  // Defining the visual style of the output.

  TString runFileName = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* runFile = CheckFile(runFileName);
  TTree* reco;  runFile -> GetObject("reco", reco);
  // Getting the TTree object containing the unbinned data.

  Double_t f90Min = 0.4;    Double_t f90Max = 0.6;
  Double_t s1Min  = 200.0;   Double_t s1Max  = 220.0;
  Double_t tpcToFMin  = 22;   Double_t tpcToFMax  = 27;
  // Defining parameter values that will be used to make the appropriate event selection.

  Double_t f90BinSize = 0.005;   Double_t f90BinCount = (f90Max - f90Min)/f90BinSize;
  // Generated histogram bin number and bin size.

  TCut tpcCut    = DefineCuts(2, f90Min, f90Max, s1Min, s1Max);
  TCut tpcToFCut = SiTelTPCToFCut(tpcToFMin, tpcToFMax);
  // Selection criteria applied for event selection.

  // </code-fold>

  TH1F* f90Histogram = new TH1F("f90Hist", "f_{90} Histogram", f90BinCount, f90Min, f90Max);
  reco -> Project("f90Hist", "clusters[0].f90", tpcCut);
  // Getting a histogram generated from unbinned data located in the TTree object 'reco'.

  RooRealVar f90("f90", "f_{90}", f90Min, f90Max);
  RooDataHist data("data", "Dataset with f_{90}", f90, f90Histogram);
  // Declaration of the variable and binned data containing the histogram.

  RooPlot* frame = f90.frame();
  data.plotOn(frame);
  frame -> Draw();

  RooRealVar mean("mean", "Mean of Gaussian", 0.55, f90Min, f90Max);
  RooRealVar sigmaRighQ("sigmaRight", "Right Width of Gaussian", 0.05, 0, 1);
  RooRealVar sigmaLeft("sigmaLeft", "Left Width of Gaussian", 0.05, 0, 1);

  RooGaussian asymGauss("gauss", "Assymetric Gaussian (x, mean, sigma)", f90, mean, sigmaLeft);
  // Defining the model function.

  asymGauss.fitTo(data);
  // Fits the biined data (histogram) to the model using a binned log-likelihood method and MINUIT.

  mean.Print();
  //sigmaRight.Print();
  sigmaLeft.Print();
  // Prints the obtained values from the fit for the model's parameters.

  asymGauss.plotOn(frame);
  frame -> Draw();
}

void Analysis( Int_t run, Double_t s1Bin = 10, Double_t s1Min = 0., Double_t s1Max = 1000. ){

  // Defining the visual style of the output.
  TStyle* sidStyle = SetSidStyle();
  sidStyle -> cd();

  // Getting the TTree object containing the unbinned data.
  TString runFileName = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* runFile = CheckFile(runFileName);
  TTree* reco;  runFile -> GetObject("reco", reco);

  // Parameters used for event selection.
  Double_t f90Min    = 0.4;    Double_t f90Max    = 0.6;
  Double_t tpcToFMin = 22.;    Double_t tpcToFMax = 27.;


  // Selection criteria definition.
  TCut tpcToFCut = SiTelTPCToFCut(tpcToFMin, tpcToFMax);

  // Quantities necessary to generate multiple fits, each for a different S1 region.
  const Int_t histNumber = (s1Max - s1Min)/s1Bin;
  Int_t s1Low;
  Int_t s1Upp;

  // Declaration of the variable that will hold the generated histograms.
  Double_t f90BinSize  = 0.005;
  Double_t f90BinCount = (f90Max - f90Min)/f90BinSize;

  TH1F* f90Histogram[histNumber];

  RooRealVar f90("f90", "f_{90}", f90Min, f90Max);
  RooDataHist* data[histNumber];

  RooRealVar mean("mean", "Mean of Gaussian", 0.55, f90Min, f90Max);
  RooRealVar sigmaRighQ("sigmaRight", "Right Width of Gaussian", 0.05, 0, 1);
  RooRealVar sigmaLeft("sigmaLeft", "Left Width of Gaussian", 0.05, 0, 1);

  RooGaussian asymGauss("gauss", "Assymetric Gaussian (x, mean, sigma)", f90, mean, sigmaLeft);

  RooPlot* f90Frame = f90.frame();
  RooPlot* drawFrame;

  Double_t height = 350; Double_t width = gdRatio*height;
  TCanvas* canvas1 = new TCanvas("c1", "f90 Distribution (1/2)", 3*width, 3*height);
  TCanvas* canvas2 = new TCanvas("c2", "f90 Distribution/ (2/2)", 3*width, 3*height);
  canvas1 -> Divide(3,3);
  canvas2 -> Divide(3,3);

  for (Int_t i = 0; i < histNumber; i++){

    s1Low = s1Min + i*s1Bin;
    s1Upp = s1Min + (i + 1)*s1Bin;

    TCut tpcCut = DefineCuts(1, f90Min, f90Max, s1Low, s1Upp);

    f90Histogram[i] = new TH1F(Form("f90Hist%d", i+1), "f_{90} Histogram", f90BinCount, f90Min, f90Max);
    reco -> Project(Form("f90Hist%d", i+1), "clusters[0].f90", tpcCut);

    data[i] = new RooDataHist(Form("data%d", i+1), "Dataset with f_{90}", f90, f90Histogram[i]);

    asymGauss.fitTo(*data[i]);
    drawFrame = f90Frame -> emptyClone("drawFrame");
    drawFrame -> SetTitle(Form("f90 (S1: %d - %d PE); f90; Events / 0.05", s1Low, s1Upp));
    if (i < 9){
      canvas1 -> cd(i+1);
    } else {
      canvas2 -> cd(i+1 - 9);
    }
    data[i] -> plotOn(drawFrame);
    asymGauss.plotOn(drawFrame);
    drawFrame -> Draw();
  }
}

TCut SiTelTPCToFCut(Double_t tofMin, Double_t tofMax){

  TCut tofMinCut = Form("2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time) >= %f", tofMin);
  TCut tofMaxCut = Form("2*(0.5*(start_time[30] + start_time[31] - 7.45) - clusters[0].cdf_time) <= %f", tofMax);
  TCut tofCut = tofMinCut && tofMaxCut;

  return tofCut;
}
