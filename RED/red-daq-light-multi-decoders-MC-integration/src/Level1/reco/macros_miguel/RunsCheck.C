#include "sidutility.cc"

#include <TCut.h>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <iterator>
#include <vector>

/* *********************************************************************************************************************** *
 *
 * File: RunsCheck.C
 *
 * Author: Miguel Del Ben Galdiano.
 * Date of Creation : April 18 2020.
 *
 * Summary of File:
 *
 *    The goal of this macro is to run a preliminary analysis to a set of runs (from 1489 up to 1521). It will also serve
 *    as a learning experience to figure out how to join TTrees with the same structure that come from different root
 *    files. The runs will be separated in two groups: those with drift field turned off (1489 to 1500) and those with
 *    drift field turned on (1501 to 1521). The first group will be further divided based on the DAQ configuration, while
 *    the second will be divide based on the CH2 target density. I believe these groups will show similar behaviour and
 *    the subdivisions will prove to be irrelevant, but nonetheless.
 *
 *    TABLE OF DIVISIONS:
 *
 *    ::                      FIELD OFF                      ::                   FIELD ON                     ::
 *    ::  SiMaster_TPCslave   Si_Standalone    Si_and_TPC    ::    401 ug/cm2     426 ug/cm2     423 ug/cm2    ::
 *    ::    1489 to 1494      1495 to 1496    1499 to 1500   ::   1501 to 1510   1511 to 1516   1519 to 1521   ::
 *
 *    DISCARDED RUNS:
 *      > 1497: Au (90 um/cm2) run.
 *      > 1498: Faulty channels, told not to use.
 *      > 1505: Au (90 um/cm2) run.
 *      > 1507: JUNK.
 *      > 1517: JUNK.
 *      > 1518: JUNK.
 *
 *    SPECIFIC RUN NOTES:
 *      > 1501: Rise field! Beam not stable.
 *
 * *********************************************************************************************************************** */

TDirectory* CH2_dirs[3];
TDirectory* s1_dirs[3];
TDirectory* s2_dirs[3];
TDirectory* dt_dirs[3];
TDirectory* f90_dirs[3];

const char* dir_name[3]    = {"CH2_401_ug_cm2", "CH2_426_ug_cm2", "CH2_423_ug_cm2"};
const char* subdir_name[4] = {"s1", "s2", "drift_time", "f90"};

Int_t    num_bins = 100;  Int_t cfg        = 2;
Double_t s1_min   = 0.;   Double_t s1_max  = 5000.;
Double_t s2_min   = 0.;   Double_t s2_max  = 10000.;
Double_t dt_min   = 0.;   Double_t dt_max  = 80;
Double_t f90_min  = 0.;   Double_t f90_max = 1.0;
Double_t tof_min  = 0.;   Double_t tof_max = 0.;
bool normalize = true;

void CheckRunFieldOnS1(){

  std::vector<int> runs_on {1501, 1502, 1503, 1504, 1506, 1508, 1509, 1510, 1511, 1512, 1513, 1514, 1515, 1516, 1519, 1520, 1521};

  TString file_name = Form( "analysis_%d%d.root", runs_on.at(0), runs_on.at(runs_on.size() - 1) );
  TFile* output_file = CheckFile(file_name);

  for (Int_t i  =  0; i < 3; i++){
    output_file -> cd();  CH2_dirs[i] = MakeDirectory(dir_name[i], dir_name[i]);        CH2_dirs[i] -> Delete("*;*");
    CH2_dirs[i] -> cd();  s1_dirs[i]  = MakeDirectory(subdir_name[0], subdir_name[0]);  s1_dirs[i] -> Delete("*;*");
  }

  TCut histogram_cuts = DefineCuts(cfg, f90_min, f90_max, s1_min, s1_max, s2_min, s2_max, tof_min, tof_max);

  TString run_name;  TString hist_name;  TString hist_title;

  for (Int_t i = 0; i < runs_on.size(); i++){

    run_name   = Form("runs/run_%d.root", runs_on.at(i));
    hist_name  = Form("s1_%d", runs_on.at(i));
    hist_title = Form("S1 Charge (Run: %d); S1 (PE)", runs_on.at(i));

    if ( i < 8 ) {
      WriteS1Histogram(s1_dirs[0], run_name, histogram_cuts, num_bins, s1_min, s1_max, normalize, hist_name, hist_title);
    } else if ( i < 14 ) {
      WriteS1Histogram(s1_dirs[1], run_name, histogram_cuts, num_bins, s1_min, s1_max, normalize, hist_name, hist_title);
    } else {
      WriteS1Histogram(s1_dirs[2], run_name, histogram_cuts, num_bins, s1_min, s1_max, normalize, hist_name, hist_title);
    }
  }
  delete output_file;
}

void CheckRunFieldOnS2(){

  std::vector<int> runs_on {1501, 1502, 1503, 1504, 1506, 1508, 1509, 1510, 1511, 1512, 1513, 1514, 1515, 1516, 1519, 1520, 1521};

  TString file_name = Form( "analysis_%d%d.root", runs_on.at(0), runs_on.at(runs_on.size() - 1) );
  TFile* output_file = CheckFile(file_name);

  for (Int_t i  =  0; i < 3; i++){
    output_file -> cd();  CH2_dirs[i] = MakeDirectory(dir_name[i], dir_name[i]);        CH2_dirs[i] -> Delete("*;*");
    CH2_dirs[i] -> cd();  s2_dirs[i]  = MakeDirectory(subdir_name[1], subdir_name[1]);  s2_dirs[i] -> Delete("*;*");
  }

  TCut histogram_cuts = DefineCuts(cfg, f90_min, f90_max, s1_min, s1_max, s2_min, s2_max, tof_min, tof_max);

  TString run_name;  TString hist_name;  TString hist_title;

  for (Int_t i = 0; i < runs_on.size(); i++){

    run_name   = Form("runs/run_%d.root", runs_on.at(i));
    hist_name  = Form("s2_%d", runs_on.at(i));
    hist_title = Form("S2 Charge (Run: %d); S2 (PE)", runs_on.at(i));

    if ( i < 8 ) {
      WriteS2Histogram(s2_dirs[0], run_name, histogram_cuts, num_bins, s2_min, s2_max, normalize, hist_name, hist_title);
    } else if ( i < 14 ) {
      WriteS2Histogram(s2_dirs[1], run_name, histogram_cuts, num_bins, s2_min, s2_max, normalize, hist_name, hist_title);
    } else {
      WriteS2Histogram(s2_dirs[2], run_name, histogram_cuts, num_bins, s2_min, s2_max, normalize, hist_name, hist_title);
    }
  }
  delete output_file;
}

void CheckRunFieldOnDrift(){

  std::vector<int> runs_on {1501, 1502, 1503, 1504, 1506, 1508, 1509, 1510, 1511, 1512, 1513, 1514, 1515, 1516, 1519, 1520, 1521};

  TString file_name = Form( "analysis_%d%d.root", runs_on.at(0), runs_on.at(runs_on.size() - 1) );
  TFile* output_file = CheckFile(file_name);

  for (Int_t i  =  0; i < 3; i++){
    output_file -> cd();  CH2_dirs[i] = MakeDirectory(dir_name[i], dir_name[i]);        CH2_dirs[i] -> Delete("*;*");
    CH2_dirs[i] -> cd();  dt_dirs[i]  = MakeDirectory(subdir_name[2], subdir_name[2]);  dt_dirs[i] -> Delete("*;*");
  }

  TCut histogram_cuts = DefineCuts(cfg, f90_min, f90_max, s1_min, s1_max, s2_min, s2_max, tof_min, tof_max);

  TString run_name;  TString hist_name;  TString hist_title;

  for (Int_t i = 0; i < runs_on.size(); i++){

    run_name   = Form("runs/run_%d.root", runs_on.at(i));
    hist_name  = Form("drift_time_%d", runs_on.at(i));
    hist_title = Form("Drift Time (Run: %d); Time (#mus)", runs_on.at(i));

    if ( i < 8 ) {
      WriteDriftTimeHistogram(dt_dirs[0], run_name, histogram_cuts, num_bins, dt_min, dt_max, normalize, hist_name, hist_title);
    } else if ( i < 14 ) {
      WriteDriftTimeHistogram(dt_dirs[1], run_name, histogram_cuts, num_bins, dt_min, dt_max, normalize, hist_name, hist_title);
    } else {
      WriteDriftTimeHistogram(dt_dirs[2], run_name, histogram_cuts, num_bins, dt_min, dt_max, normalize, hist_name, hist_title);
    }
  }
  delete output_file;
}

void CheckRunFieldOnF90(){

  std::vector<int> runs_on {1501, 1502, 1503, 1504, 1506, 1508, 1509, 1510, 1511, 1512, 1513, 1514, 1515, 1516, 1519, 1520, 1521};

  TString file_name = Form( "analysis_%d%d.root", runs_on.at(0), runs_on.at(runs_on.size() - 1) );
  TFile* output_file = CheckFile(file_name);

  for (Int_t i  =  0; i < 3; i++){
    output_file -> cd();  CH2_dirs[i] = MakeDirectory(dir_name[i], dir_name[i]);         CH2_dirs[i] -> Delete("*;*");
    CH2_dirs[i] -> cd();  f90_dirs[i]  = MakeDirectory(subdir_name[3], subdir_name[3]);  f90_dirs[i] -> Delete("*;*");
  }

  TCut histogram_cuts = DefineCuts(cfg, f90_min, f90_max, s1_min, s1_max, s2_min, s2_max, tof_min, tof_max);

  TString run_name;  TString hist_name;  TString hist_title;

  for (Int_t i = 0; i < runs_on.size(); i++){

    run_name   = Form("runs/run_%d.root", runs_on.at(i));
    hist_name  = Form("f90_%d", runs_on.at(i));
    hist_title = Form("f90 (Run: %d); f90", runs_on.at(i));

    if ( i < 8 ) {
      WriteF90Histogram(f90_dirs[0], run_name, histogram_cuts, num_bins, f90_min, f90_max, normalize, hist_name, hist_title);
    } else if ( i < 14 ) {
      WriteF90Histogram(f90_dirs[1], run_name, histogram_cuts, num_bins, f90_min, f90_max, normalize, hist_name, hist_title);
    } else {
      WriteF90Histogram(f90_dirs[2], run_name, histogram_cuts, num_bins, f90_min, f90_max, normalize, hist_name, hist_title);
    }
  }
  delete output_file;
}

void CheckRunFieldOn(bool s1 = true, bool s2 = true, bool drift_time = true, bool f90 = true){

  if ( s1 )         CheckRunFieldOnS1();
  if ( s2 )         CheckRunFieldOnS2();
  if ( drift_time ) CheckRunFieldOnDrift();
  if ( f90 )        CheckRunFieldOnF90();
}


void CheckRunFieldOffS1() {

  std::vector<int> runs_off { 1489, 1490, 1491, 1492, 1493, 1494, 1495, 1496, 1499, 1500 };

  TString file_name   = Form( "analysis_%d-%d.root", runs_off.at(0), runs_off.at(runs_off.size()-1) );
  TFile*  output_file = CheckFile(file_name);

  output_file -> cd();
  TDirectory* Si_TPCslave   = MakeDirectory("TPC_SiMaster_TPCSlave", "TPC_SiMaster_TPCSlave");
  TDirectory* Si_Standalone = MakeDirectory("Si_Standalone", "Si_Standalone");
  TDirectory* Si_and_TPC    = MakeDirectory("Si_and_TPC", "Si_and_TPC");


  TCut histogram_cuts = DefineCuts(cfg, f90_min, f90_max, s1_min, s1_max, 0., 0.);

  TString run_name;  TString hist_name;  TString hist_title;
  for (Int_t i = 0; i < runs_off.size(); i++){

    run_name  = Form("runs/run_%d.root", runs_off.at(i));
    hist_name = Form("s1_%d", runs_off.at(i));
    hist_title = Form("S1 Charge (Run: %d); S1 (PE)", runs_off.at(i));

    if ( i < 6 ) {
      WriteS1Histogram(Si_TPCslave,   run_name, histogram_cuts, num_bins, s1_min, s1_max, normalize, hist_name, hist_title);
    } else if ( i < 8 ) {
      WriteS1Histogram(Si_Standalone, run_name, histogram_cuts, num_bins, s1_min, s1_max, normalize, hist_name, hist_title);
    } else {
      WriteS1Histogram(Si_and_TPC,    run_name, histogram_cuts, num_bins, s1_min, s1_max, normalize, hist_name, hist_title);
    }
  }
  delete output_file;
}
