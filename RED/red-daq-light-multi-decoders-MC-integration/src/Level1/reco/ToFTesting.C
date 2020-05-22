#include "sidutility.cc"

#include <TCut.h>
#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <TTree.h>

/* *********************************************************************************************************************** *
 *                                                                                                                         *
 * File: ToFDefinitionTest.C                                                                                               *
 *                                                                                                                         *
 * Author: Miguel Del Ben Galdiano.                                                                                        *
 * Date of Creation : April 24 2020.                                                                                       *
 *                                                                                                                         *
 * Summary of File:                                                                                                        *
 *                                                                                                                         *
 *    The goal of this macro is to test the differences, if any, between different possible definitions of the time of     *
 *    flight parameter. This macro will test how different ToFs (such as that between SiTel and TPC or SiTel and LSci)     *
 *    are different. It will also test how the definition of the parameter based on the reconstructed values obtained      *
 *    from the reconstructed run files show different behaviours, if any. This code shall test the following definitions:  *
 *                                                                                                                         *
 *      > 2*(xmin[] - clusters[0].min_x).   <-- currently the one in use.                                                  *
 *      > 2*(start_time[] - clusters[0].cdf_time).                                                                         *
 *                                                                                                                         *
 *    This macro will also be used to study other aspects of the time of flight analysis. These will be based on studies   *
 *    made by members of the collaboration, such as to better understand their results and develop a more uniform          *
 *    definition of the parameter, such as to allow better comparisions.                                                   *
 *                                                                                                                         *
 * *********************************************************************************************************************** */


void ToFSiTelTPCDefinitionTest( int run ){

  TFile* file = CheckFile( Form("runs/run_%d.root", run) );
  TTree* reco = 0; file -> GetObject("reco", reco);

  Int_t cfg = 2;
  Double_t f90_min = 0.0;    Double_t f90_max = 1.0;
  Double_t s1_min  = 0.0;    Double_t s1_max = 10000;
  Double_t s2_min  = 0.0;    Double_t s2_max = 30000;
  Double_t tof_min = -100.;  Double_t tof_max = 100.;

  TH1F* tof_histograms[2] ;

  TString tof_definitions[2] = {"2*(xmin[30] - clusters[0].min_x)",
                                "2*(start_time[30] - clusters[0].cdf_time)"   };

  TCut general_cut = DefineCuts(cfg, f90_min, f90_max, s1_min, s1_max, s2_min, s2_max, 0, 0);
  TCut tof_cuts[2];

  THStack* tof_stack = new THStack("tof_stack", "Time of Flight (SiTel and TPC)");

  for ( Int_t i = 0; i < 2; i++ ){

    tof_cuts[i] = tof_definitions[i] + Form(" > %f", tof_min) + " && " + tof_definitions[i] + Form(" < %f", tof_max);

    reco -> Draw( tof_definitions[i], tof_cuts[i] && general_cut, "goff");
    tof_histograms[i] = (TH1F*) gDirectory -> Get("htemp");
    tof_histograms[i] = (TH1F*) NormalizeHistogram(tof_histograms[i]);
    tof_histograms[i] -> SetLineColor(2 + 1*i);

    tof_stack -> Add(tof_histograms[i]);
  }

  tof_stack -> Draw("hist nostack");
  tof_stack -> GetXaxis() -> SetTitle("ns");

  auto legend = new TLegend(0.710396,0.751579,0.888614,0.882105);
  legend -> AddEntry(tof_histograms[0], tof_definitions[0], "l");
  legend -> AddEntry(tof_histograms[1], tof_definitions[1], "l");
  legend -> Draw();
}

void ToFBinTest( int run ){

  TFile* file = CheckFile( Form("runs/run_%d.root", run) );
  TTree* reco = 0; file -> GetObject("reco", reco);

  Int_t cfg = 2;
  Double_t f90_min = 0.0;    Double_t f90_max = 1.0;
  Double_t s1_min  = 0.0;    Double_t s1_max = 10000;
  Double_t s2_min  = 0.0;    Double_t s2_max = 30000;
  Double_t tof_min = -100.;  Double_t tof_max = 100.;

  TH1F* tof_histograms[2][4];
  Int_t num_bins[4] = {100, 200, 400, 800};

  TString tof_definitions[2] = {"2*(xmin[30] - clusters[0].min_x)",
                                "2*(start_time[30] - clusters[0].cdf_time)"};

  TCut general_cut = DefineCuts(cfg, f90_min, f90_max, s1_min, s1_max, s2_min, s2_max, 0, 0);
  TCut tof_cuts[2];

  for (Int_t i = 0; i < 2; i++){
    tof_cuts[i] = tof_definitions[i] + Form(" > %f", tof_min) + " && " + tof_definitions[i] + Form(" < %f", tof_max);
    for(Int_t j = 0; j < 4; j++){
      tof_histograms[i][j] = new TH1F(Form("tof_hist%d%d", i, j), "Timing (SiTel and TPC); ns", num_bins[j], tof_min, tof_max);
      reco -> Project( Form("tof_hist%d%d", i, j), tof_definitions[i], tof_cuts[i] && general_cut );
    }
  }

  Double_t gold_ratio = 1.61803398875;
  Double_t height = 500;    Double_t width = gold_ratio * height;

  TCanvas* canvas1 = new TCanvas("canvas1", "Timing Bin Size", 2*width, height);
  canvas1 -> Divide(2,1);

  TCanvas* canvas2 = new TCanvas("canvas2", "Timing Bin Size", 2*width, height);
  canvas2 -> Divide(2,1);

  TCanvas* canvas3 = new TCanvas("canvas3", "Timing Bin Size", 2*width, height);
  canvas3 -> Divide(2,1);

  TCanvas* canvas4 = new TCanvas("canvas4", "Timing Bin Size", 2*width, height);
  canvas4 -> Divide(2,1);

  canvas1 -> cd(1); tof_histograms[0][0] -> Draw(); canvas1 -> cd(2); tof_histograms[1][0] -> Draw();
  canvas2 -> cd(1); tof_histograms[0][1] -> Draw(); canvas2 -> cd(2); tof_histograms[1][1] -> Draw();
  canvas3 -> cd(1); tof_histograms[0][2] -> Draw(); canvas3 -> cd(2); tof_histograms[1][2] -> Draw();
  canvas4 -> cd(1); tof_histograms[0][3] -> Draw(); canvas4 -> cd(2); tof_histograms[1][3] -> Draw();

}
