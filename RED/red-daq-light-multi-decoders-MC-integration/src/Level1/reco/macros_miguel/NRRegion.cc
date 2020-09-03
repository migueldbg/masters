#include "sidutility.cc"

#include "TFile.h"
#include "TH1.h"
#include "TStyle.h"
#include "TTree.h"

using namespace std;
using namespace RooFit;

/* ************************************************************************************************************************* *
 * File: NRRegion.cc (ROOT macro).                                                                                             *
 *                                                                                                                           *
 * Author: Miguel Del Ben Galdiano (miguel.galdiano@gmail.com).                                                              *
 * Date of creation : August 29 2020.                                                                                        *
 *                                                                                                                           *
 * Summary of File:                                                                                                          *
 *                                                                                                                           *
 *    This macro should contain all functions that pertain to determining the NR region of interest considering different    *
 *    parameters spaces. The NR region of interest should be better thought of as the NR 'candidate' region, as there is no  *
 *    Number guarantee that a event within a given NR region isn't simply a background event.                                  *                                                                       *
 *                                                                                                                           *
 * ************************************************************************************************************************* */

 //<code-fold> AUXILIARY FUNCTIONS DECLARATION.
 Int_t NumberOfHistograms(TDirectory* directory);
 Int_t NumberOfCanvas(Int_t plotNumber, Int_t padNumber);
 //</code-fold>

/* void GenerateF90Histograms( Int_t run, Double_t s1BinSize = 20, Double_t s1Min = 0, Double_t s1 = 1e4, bool write = false )
 *
 *  Summary of Function:
 *
 *    This function is intended to be used for determining the NR region in the F90 vs S1 parameter space.
 *
 *    This function generates a set of f90 histograms for a given run data set. The function considers the S1 range determined
 *    by the user and within that range, generates a number of f90 histograms. The number of the histgograms is defined by the
 *    S1 bin size, also determined by the user. For each S1 region, a f90 histogram is generated and then saved to a root file
 *    for later use.
 *
 *  Parameters   : run       >> the run containing the reconstructed data.
 *                 s1BinSize >> size of each s1 region for which a f90 histogram will be generated.
 *                 s1Min     >> minimum value of s1 to be considered.
 *                 s1Max     >> maximum value of s1 to be considered.
 *                 write     >> indicates whether or not the generated histograms will be saved to a root file.
 *
 *  Return Value : void.
 */
void GenerateF90Histograms( Int_t run, Double_t s1BinSize = 20, Double_t s1Min = 0, Double_t s1Max = 1e4, bool write = false ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> cd();

  // Getting the TTree object containing the unbinned data.
  TString runFileName = runsDirectoryPath + Form("/run_%d.root", run);
  TFile* runFile = CheckFile(runFileName);
  TTree* reco;  runFile -> GetObject("reco", reco);


  // Variables to determine the histogram's relevant parameters (number of bins, size of bin, minimum and maximum value).
  Double_t f90Min = 0.4;
  Double_t f90Max = 0.6;

  Double_t f90BinSize = 0.005;   Double_t f90BinNumber = (f90Max - f90Min)/f90BinSize;

  const Int_t s1BinNumber = (const Int_t)(s1Max - s1Min) / s1BinSize;

  // Declaring the variables that will hold the generated histograms.
  TH1F*   f90Hist[s1BinNumber];
  TString histName[s1BinNumber];
  TString histTitle[s1BinNumber];


  // Variables to be used when selecting the appropriate s1 bin for which to construct a f90 histogram.
  TCut s1Range;
  Double_t s1Low;
  Double_t s1Upp;


  // Generates the f90 histograms for different s1 regions.
  for (Int_t i = 0; i < s1BinNumber; i++){

    s1Low =  i      * s1BinSize + s1Min;
    s1Upp = (i + 1) * s1BinSize + s1Min;

    s1Range = DefineS1Range(s1Low, s1Upp);

    histName[i]  = Form("f90_histogram_%d", i+1);
    histTitle[i] = Form("f_{90} Distribution (S1 Interval: %d - %d PE)", (Int_t) s1Low, (Int_t) s1Upp);
    f90Hist[i] = new TH1F(histName[i], histTitle[i], f90BinNumber, f90Min, f90Max);

    reco -> Project(histName[i], "clusters[0].f90", s1Range);

  }

  // Saves the generated histograms to root file for later use.
  if (write){

    TString outputFileName = analysisDirectoryPath + Form("/fit_%d.root", run);
    TFile* outputFile = CheckFile(outputFileName);

    TDirectory* histDir = MakeDirectory("histograms", "histograms");
    histDir -> cd();
    TDirectory* f90Dir = MakeDirectory("f90", "f90");
    f90Dir -> Delete("*;*");

    for (Int_t i = 0; i < s1BinNumber; i++)  f90Dir -> WriteObject(f90Hist[i], histName[i], "OverWrite");

    outputFile -> Close();

  }
}

/* FitF90Histograms( Int_t run )
 *
 *  Summary of Function:
 *
 *  Parameters :
 *
 *  Return Value:
 */
void FitF90Histograms( Int_t run, bool draw = true ){

  TStyle* sidStyle = SetSidStyle();
  sidStyle -> cd();

  // Opening the file containing the histograms generated by 'GenerateF90Histograms'
  TString outputFileName = analysisDirectoryPath + Form("/fit_%d.root", run);
  TFile*  outputFile = CheckFile(outputFileName);

  TDirectory* histDir = MakeDirectory("histograms", "histograms");
  histDir -> cd();
  TDirectory* f90Dir = MakeDirectory("f90", "f90");


  // Declaring variables to hold the histograms.
  const Int_t histNumber = NumberOfHistograms(f90Dir);

  TH1F*   f90Hist[histNumber];
  TString histName[histNumber];

  Int_t binNumber;

  // Getting the histograms located in the 'fit' file.
  for (Int_t i = 0; i < histNumber; i++){

    histName[i] = Form("f90_histogram_%d", i+1);
    f90Hist[i] = (TH1F*) f90Dir -> Get(histName[i]);

  }

  // Getting the lower and upper boundaries of the f90 histogram, as they will be used to define the RooFit model.
  Double_t f90Min = f90Hist[0] -> GetXaxis() -> GetBinLowEdge(1);
  Double_t f90Max = f90Hist[0] -> GetXaxis() -> GetBinLowEdge(f90Hist[0] -> GetNbinsX() + 1);


  // <code-fold> SETTING UP THE MODEL WITH ROOFIT USED TO FIT THE F90 HISTOGRAMS.
  RooDataHist* f90Data[histNumber];

  RooRealVar f90("f90", "f_{90}", f90Min, f90Max);
  RooRealVar mean("mean", "Mean of Gaussian", 0.55, f90Min, f90Max);
  RooRealVar sigma("sigma", "Width of Gaussian", 0.05, 0, 1);

  RooGaussian* gauss[histNumber];
  //</code-fold>

  // Fitting the binned data (f90 histograms) to the model defined.
  for (Int_t i = 0; i < histNumber; i++){

    gauss[i]   = new RooGaussian(Form("f90_gauss_%d", i+1), "Gaussian (f_{90}, mean, sigma)", f90, mean, sigma);
    f90Data[i] = new RooDataHist(Form("f90_data_%d", i+1), "Dataset with f_{90}", f90, f90Hist[i]);

    gauss[i] -> fitTo(*f90Data[i]);
  }

  if (draw){

    Int_t padTrue = 9; // The desired number of pads for each canvas.
    const Int_t canvasNumber = NumberOfCanvas(histNumber, padTrue);
    TCanvas* canvas[canvasNumber];

    // <code-fold> SETTING UP THE NUMBER OF PADS FOR EACH CANVAS
    // Variables to divide each canvas in the appropriate number of pads.
    Int_t padNumber[canvasNumber];

    Int_t remainder;

    // Dividing each canvas with the appropriate number of pads.
    for (Int_t i = 0; i < canvasNumber; i++){

      remainder = histNumber - padTrue*i;

      if ( remainder >= padTrue)
        padNumber[i] = padTrue;
      else if (remainder < padTrue && remainder != 0){
        padNumber[i] = remainder;
      }

      // This statement determines how the canvas will be divided depending on the number of pads given.
      switch(padNumber[i]){
        case 1:
          break;
        case 2:
          canvas[i] -> Divide(2,1);
          break;
        case 3:
          canvas[i] -> Divide(3,1);
          break;
        case 4:
          canvas[i] -> Divide(2,2);
          break;
        case 5:
        case 6:
          canvas[i] -> Divide(2,3);
          break;
        case 7:
        case 8:
          canvas[i] -> Divide(2,4);
          break;
        case 9:
          canvas[i] -> Divide(3,3);
          break;
        case 10:
          canvas[i] -> Divide(2,5);
          break;
      }
    }
    // </code-fold>

    for ( Int_t i = 0; i < histNumber; i++ ){

      canvas[i/padTrue] -> cd(i+1 - (i/padTrue * padTrue));
    }
  }

}

void test(Int_t histNumber, Int_t padTrue){

  for (Int_t i = 0; i < histNumber; i++){
    cout << "i = " << i << "    c(i) = " << i/padTrue << "    d(i) = " << i + 1 - (i/padTrue * padTrue) << endl;
  }
}

// =================================== AUXILIARY FUNCTIONS USED IN THE MACRO =================================== //

/* Int_t NumberOfHistograms ( TDirectory* directory )
 *
 * Summary of Function:
 *
 *    The NumberofHistograms function counts the number of one dimenstional histograms in a given directory
 *    by interating over all objects located in it and checking which inherit from TH1. It then returns the
 *    number of 1D histograms in the directory.
 *
 * Parameters   : directory >> a pointer to the the directory where the histograms are located.
 *
 * Return value : Int_t histNumber.
 */
Int_t NumberOfHistograms( TDirectory* directory ){

  Int_t histNumber = 0;
  TKey* key;

  TIter next((TList*) directory -> GetListOfKeys());

  while ( (key = (TKey*) next()) ){
    TClass *objectClass = gROOT -> GetClass(key->GetClassName());

    if ( objectClass -> InheritsFrom("TH1") ) {
      histNumber++;
    }
  }

  return histNumber;
}

/* Int_t NumberOfCanvas( Int_t plotNumber, Int_t padNumber )
 *
 */
Int_t NumberOfCanvas( Int_t plotNumber, Int_t padNumber ){

  Int_t canvasNumber;

  if (plotNumber % padNumber == 0){
    canvasNumber = plotNumber/padNumber;
  } else if (plotNumber % padNumber != 0){
    canvasNumber = plotNumber/padNumber + 1;
  }

  return canvasNumber;
}
