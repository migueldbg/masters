#include <iostream>
#include <fstream>
  #include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TRint.h"
#include "TList.h"
#include "TCollection.h"
#include "TKey.h"
#include "TClass.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TColor.h"
#include "TLegend.h"
#include "TVirtualFitter.h"
#include <TGaxis.h>
#include <TLatex.h>

using namespace std;
#ifndef __CINT__
  int main(int argc, char **argv);
#endif

TRint* theApp;
TCanvas *c1;
TLegend* leg;
vector<int> colors;
vector<TString> labels;
void SetColors();
void LoopOverKeys(TList *sourcelist);
void DrawHistograms(TList *hlist);
Double_t getNormalization(const TH1& h);

void CompareRuns(const char** runs, Int_t nruns) {
  SetColors();
//  gStyle->SetPalette(kBird);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
//  gStyle->SetOptTitle(0);
  TGaxis::SetMaxDigits(3);

//  gROOT->SetBatch(kTRUE);

  Double_t btmMargin(0.11), topMargin(0.05), ritMargin(0.08), lftMargin(0.11);
  c1 = new TCanvas("c1","c1",1000, 750);
  c1->SetBottomMargin(btmMargin);
  c1->SetTopMargin(topMargin);
  c1->SetRightMargin(ritMargin);
  c1->SetLeftMargin(lftMargin);
  c1->cd();
  c1->SetGrid();
  c1->Draw();

  leg = new TLegend(0.7, 0.2, 0.9, 0.35);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetFillColor(10);
//  leg->SetHeader("Drift Time");

  TList *fileList = new TList();
  for (Int_t irun = 1; irun <= nruns; irun++){
	  TString filename = Form("data/run_%s_Hist.root", runs[irun]);
	  TFile *file = new TFile(filename.Data());
	  fileList->Add(file);
	  labels.push_back(runs[irun]);
	  cout<<"run: "<<runs[irun]<<endl;
  }

  LoopOverKeys(fileList);


  delete c1;
  delete leg;
}

void LoopOverKeys(TList *sourcelist){
  TFile *first_source = (TFile*)sourcelist->First();
  first_source->cd();
  TDirectory *current_sourcedir = gDirectory;
  cout << "path: " << current_sourcedir->GetPath() << endl;
  //gain time, do not add the objects in the list in memory
  Bool_t status = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  // loop over all keys in this directory
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key, *oldkey=0;
  while ( (key = (TKey*)nextkey())) {
    //keep only the highest cycle number for each key
    if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

    TString tmp(key->GetName());
    if (tmp.Contains("hRunTime") || tmp.Contains("hEventConter")) continue;
    if (tmp.Contains("hS2OverS2TotVsR_ch")) continue;
    if (tmp.Contains("pS2OverS2TotVsXY_ch")) continue;
    if (tmp.Contains("LSci")) continue;


    // read object from first source file
//    first_source->cd();
    TObject *obj = key->ReadObj();

    cout << "Obtain "<<key->GetName()<<endl;
    if ( obj->IsA()->InheritsFrom( "TH1" ) ) {
      if(obj->IsA()->InheritsFrom(TH2::Class()) || obj->IsA()->InheritsFrom(TH3::Class())) continue;
      // descendant of TH1 -> merge it
      TH1 *h1 = (TH1*)obj;

      TList *hlist = new TList();
      hlist->Add(h1);
      // loop over all source files and add the content of the
      // correspondant histogram to the one pointed to by "h1"
      TFile *nextsource = (TFile*)sourcelist->After( first_source );
      while ( nextsource ) {

        // make sure we are at the correct directory level by cd'ing to path
        nextsource->cd();
        TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
        if (key2) {
           TH1 *h2 = (TH1*)key2->ReadObj();
           hlist->Add(h2);
        }

        nextsource = (TFile*)sourcelist->After( nextsource );
      }

      DrawHistograms(hlist);
    } else {

      // object is of no type that we know or can handle
      cout << "Unknown object type, name: "
           << obj->GetName() << " title: " << obj->GetTitle() << endl;
    }


  } // while ( ( TKey *key = (TKey*)nextkey() ) )
  TH1::AddDirectory(status);
}

void DrawHistograms(TList *hlist){
#define NORM
  TIter next(hlist);
  TH1 *obj;
  Int_t counter(0);
  c1->cd();
  TString hname;
  while ((obj = (TH1*) next())) {
      hname = obj->GetName();
      obj->SetName(Form("%s_%d", hname.Data(), counter));
      obj->SetLineColor(colors.at(counter+1));
      obj->SetMarkerColor(colors.at(counter+1));
//      obj->SetLineStyle(counter+1);
      obj->SetLineWidth(1);
      if ( obj->IsA()->InheritsFrom( "TH2" ) ){
          c1->SetLogy(0);
//          (counter==0)? obj->DrawClone("colz") : obj->DrawClone("colz same");
          c1->Clear();
          c1->SetRightMargin(0.11);
          obj->DrawClone("colz");
          c1->Print(Form("%s_%s.pdf",hname.Data(),labels.at(counter).Data()));

      }else {
          c1->SetRightMargin(0.08);

#ifdef NORM
          obj->Scale( getNormalization(*obj) );
#else
          obj->SetYTitle("counts");
#endif
          TH1* htmp = (counter==0)? (TH1*)obj->DrawClone("HIST") : (TH1*)obj->DrawClone("same HIST");
          leg->AddEntry(htmp, labels.at(counter).Data(), "L");
          if(counter==0){
#ifdef NORM
  TLatex L;
  L.SetNDC();
  L.SetTextSize(0.03);
  Double_t xposi = 0.75;
  L.DrawLatex(xposi, 0.75, "Normalized");
#endif
          }
//          c1->SetLogy();
      }
      counter++;
//      theApp->Run();
  }
  leg->Draw();


  c1->Update();
  theApp->Run();
  c1->Print(Form("%s.pdf",hname.Data()));
  leg->Clear();
  c1->Clear();
}

//____________________________________________________________________________________________________
Double_t getNormalization(const TH1& h)
{
  /// Get normalization of histogram
  /// Normalization = 1/(Nevts * binWidth)

  const Double_t ntrack   = h.Integral() ;
  if( ntrack == 0 ) return 1.0 ;

  const Double_t binWidth = h.GetBinWidth(1) ;
  if( binWidth == 0.0 ) return 1.0 ;

  return 1.0/(ntrack*binWidth);
}


void SetColors(){
    colors.push_back(TColor::GetColor("#E6855E"));
    colors.push_back(TColor::GetColor("#44A5CB"));
    colors.push_back(TColor::GetColor("#D45D87"));
    colors.push_back(TColor::GetColor("#40BFB0"));
    colors.push_back(TColor::GetColor("#F9DB57"));
    colors.push_back(TColor::GetColor("#9D73BB"));
    colors.push_back(TColor::GetColor("#009F8C"));
    colors.push_back(TColor::GetColor("#8B90BE"));
    colors.push_back(TColor::GetColor("#F3C0AB"));

}

#ifndef __CINT__
int main(int argc, char **argv) {
	theApp = new TRint("App", &argc, argv, NULL, 0);
        theApp->Connect("KeyPressed(Int_t)","TSystem",gSystem,"ExitLoop()");
	if ( theApp->Argc() != 1 ) {
		std::cout << "==> Application start." << std::endl;
		CompareRuns((const char**)theApp->Argv(), theApp->Argc()-1);

	} else {
		std::cout << "Usage:" <<theApp->Argc()<< std::endl;
		std::cout << "./comp fInNameHiVariance fInNameLowVariance" << std::endl;
		return 0;
	}


	std::cout << "==> Application finished." << std::endl;
//	theApp->Run();

	return 0;
}
#endif /* __CINT __ */
