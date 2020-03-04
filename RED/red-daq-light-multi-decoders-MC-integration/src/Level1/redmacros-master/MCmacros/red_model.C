#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF2.h"
#include "TStyle.h"

//List of the parameters
double LEff       = 0.315 ;  // visible energy fraction 
double alpha      = 1.;      // Nex/Ni - fixed in this model 
double g1         = 0.2 ;    // detection probability of scintillation light ( in DS50, 0.16 gives a LY of 7.1 PE/keV @200 V/cm. Scales linearly. Expect ~0.2?) 
double g2         = 28.1 ;     // PE/e- ( in DS50 is 23. In ReD? ) 
double sigmaS1    = 1.4 ;    // resolution of S1. It may include a gaussian smearing? 
double sigmaS2    = 10. ;     // resolution of S2. It may include a gaussian smearing? 
double aparameter = 0.08;   // The amplitude of the "directional effect"
TRandom3 * ran ; 

///reference 

/*
double LEff       = 0.315 ;  // visible energy fraction 
double alpha      = 1.;      // Nex/Ni - fixed in this model 
double g1         = 0.2 ;    // detection probability of scintillation light ( in DS50, 0.16 gives a LY of 7.1 PE/keV @200 V/cm. Scales linearly. Expect ~0.2?) 
double g2         = 10 ;     // PE/e- ( in DS50 is 23. In ReD? ) 
double sigmaS1    = 2. ;    // resolution of S1. It may include a gaussian smearing? 
double sigmaS2    = 2. ;     // resolution of S2. It may include a gaussian smearing? 
double aparameter = 0.08;   // The amplitude of the "directional effect"


*/






//////////////////////////////////////////////////////////////////////////////////////////
// Recombination Probability  (fraction of e-/ion pairs that recombine). 
//
// Extracted from ARIS, extrapolated and approximated to the value at 70 keV
//
// The model for the effect of the track orientation is invented. Theta is the track azimuthal angle, 
// the amplitude of the affect (aparameter) is set above
//////////////////////////////////////////////////////////////////////////////////////////
double recop (double E, double theta) {
  double   rr = 0.9 *(1. + aparameter*abs(cos(theta))) ;
  return   rr>1?1.:rr  ;  
}









//////////////////////////////////////////////////////////////////////////////////////////
// Model for S1 and S2 
//
// Quenching is extracted from ARIS (0.3) and it is assumed as fraction of visible energy
// 
// The effective model is based on DarkSide-50 (DSLight3). Pay attention to the fluctutations and 
// to the correlations, since the original version is optimized for S1 only (The S2 is cross
// calibrated based on DS50 data and it is in good agreement with ARIS) 
// 
// 
//////////////////////////////////////////////////////////////////////////////////////////
void S1S2 (double E /*keV*/, double theta, double &S1, double &S2) {

  double Nq   = LEff * E / 19.5e-3 ;
  //Fluctuate here? 
  Nq          = ran->Gaus(Nq,  sqrt(Nq)) ;      

  double Ni   = Nq/(1 + alpha) ; 
  
  //compute recombination as function of angle. E is redundant for the moment
  double reco = recop (E, theta)  ; 
  
  //compute Nex (excitons) and Ni (ions)
  
  //A - Fluctuate? 
  double Nex  = Nq -  Ni + ran->Binomial( Ni, reco )  ; 
  //B - Do not fluctuate
  //double Nex  = Nq - ( 1 - reco ) * Ni  ; 
  Ni   = Nq - Nex  ; 

  //In DSLigh3 we fluctuate here. If fluctuations are above or below, may be skipped?  
  double Nexg =  ran->Binomial ( Nex,  g1 ) ;   
  double Nig  =  Ni*g2 ; //  ran->Gaus(Ni*g2 ,  sqrt(Ni*g2)) ;   

  //Detection probability and multilpication. Fluctuations? 
  //Nexg  *= g1  ; 
  //Nig   *= g2 ; 

  //Detector S1 and S2 resolutions
  S1 = ran->Gaus( Nexg, sigmaS1  * sqrt ( Nexg )   ) ;    
  S2 = ran->Gaus( Nig , sigmaS2  * sqrt ( Nig  )   ) ;   
    
}



















//This function prepares the data samples. 
//It reads a dstree, it produces 2D distributions of Ar40 recoil energy and angle wrt to vertical. 
//The preliminary selection is based on TOF, no selection based on Si detector in the scattering chamber 
//
//Some extra histograms are filled for validation

void prepare_MC_samples () {

   ran = new TRandom3 ; 

   //TFile *f = TFile::Open("/storage/gpfs_ds50/darkside/users/kussds/test/ReD_LNS_Neutrons_v2.root");
   //TFile *f = TFile::Open("/storage/gpfs_ds50/darkside/users/kussds/test/ReD_LNS_Neutrons_v2.root");
   //TTree * dstree = (TTree*)  f->Get("dstree");

   TChain * dstree = new TChain("dstree") ; 
   //for (int i=0;i<30;++i) dstree->Add(Form ("g4ds10/Linux-g++/Red/out_Red_%i.root", i) ) ; 
   //dstree->Add( "g4ds10/Linux-g++/Red/merged.root"); 
   dstree->Add( "merged.root"); 

//Declaration of leaves types
   Int_t           ev;
   Int_t           pdg;
   Float_t         ene0;
   Float_t         s1ene;
   Float_t         s2ene;
   Float_t         veto_visene;
   Float_t         mu_visene;
   Float_t         tpcene;
   Float_t         vetoene;
   Float_t         muene;
   Float_t         x;
   Float_t         y;
   Float_t         z;
   Float_t         radius;
   Float_t         px;
   Float_t         py;
   Float_t         pz;
   Float_t         bx;
   Float_t         by;
   Float_t         bz;
   Int_t           npeS1;
   Int_t           npeS2;
   Int_t           npe;
   Int_t           munpe;
   Int_t           vnpe;
   Int_t           nph;
   Int_t           ndaughters;
   Int_t           ndeposits;
   Int_t           nusers;
   Int_t           dau_id[0];
   Int_t           dau_pdg[0];
   Int_t           dau_pid[0];
   Int_t           dau_process[0];
   Double_t        dau_time[0];
   Float_t         dau_ene[0];
   Float_t         dau_x[0];
   Float_t         dau_y[0];
   Float_t         dau_z[0];
   Float_t         dau_r[0];
   Float_t         dau_px[0];
   Float_t         dau_py[0];
   Float_t         dau_pz[0];
   Int_t           dep_pdg[10000];
   Int_t           dep_mat[10000];
   Double_t        dep_time[10000];
   Float_t         dep_ene[10000];
   Float_t         dep_step[10000];
   Float_t         dep_x[10000];
   Float_t         dep_y[10000];
   Float_t         dep_z[10000];
   Float_t         dep_r[10000];
   Int_t           userint1[300];
   Int_t           userint2[300];
   Float_t         userfloat1[300];
   Float_t         userfloat2[300];
   Double_t        userdouble0[300];
   Double_t        pe_time[00];
   Int_t           pe_pmt[00];
   Double_t        vpe_time[0];
   Int_t           vpe_pmt[0];
   Double_t        mupe_time[0];
   Int_t           mupe_pmt[0];
   Int_t           ph_volume[0];
   Int_t           ph_pid[0];
   Float_t         ph_wl[0];
   Float_t         ph_x[0];
   Float_t         ph_y[0];
   Float_t         ph_z[0];
   Double_t        ph_time[0];
   Int_t           nTPC;
   Double_t        TPCtime;
   Float_t         TPCene;
   Float_t         TPCf90;
   Int_t           nAr40;
   Float_t         Ar40eneCalc;
   Float_t         Ar40ene;
   Float_t         Ar40x;
   Float_t         Ar40y;
   Float_t         Ar40z;
   Float_t         Ar40px;
   Float_t         Ar40py;
   Float_t         Ar40pz;
   Int_t           LSci_ndeps[9];
   Double_t        LSci_time[9];
   Float_t         LSci_x0[9];
   Float_t         LSci_y0[9];
   Float_t         LSci_z0[9];
   Float_t         LSci_x[9];
   Float_t         LSci_y[9];
   Float_t         LSci_z[9];
   Float_t         LSci_eneRaw[9];
   Float_t         LSci_ene[9];
   Float_t         LSci_eneToF[9];
   Int_t           nLSci_trg;

   // Set branch addresses.
   dstree->SetBranchAddress("ev",&ev);
   dstree->SetBranchAddress("pdg",&pdg);
   dstree->SetBranchAddress("ene0",&ene0);
   dstree->SetBranchAddress("s1ene",&s1ene);
   dstree->SetBranchAddress("s2ene",&s2ene);
   dstree->SetBranchAddress("veto_visene",&veto_visene);
   dstree->SetBranchAddress("mu_visene",&mu_visene);
   dstree->SetBranchAddress("tpcene",&tpcene);
   dstree->SetBranchAddress("vetoene",&vetoene);
   dstree->SetBranchAddress("muene",&muene);
   dstree->SetBranchAddress("x",&x);
   dstree->SetBranchAddress("y",&y);
   dstree->SetBranchAddress("z",&z);
   dstree->SetBranchAddress("r",&radius);
   dstree->SetBranchAddress("px",&px);
   dstree->SetBranchAddress("py",&py);
   dstree->SetBranchAddress("pz",&pz);
   dstree->SetBranchAddress("bx",&bx);
   dstree->SetBranchAddress("by",&by);
   dstree->SetBranchAddress("bz",&bz);
   dstree->SetBranchAddress("npeS1",&npeS1);
   dstree->SetBranchAddress("npeS2",&npeS2);
   dstree->SetBranchAddress("npe",&npe);
   dstree->SetBranchAddress("munpe",&munpe);
   dstree->SetBranchAddress("vnpe",&vnpe);
   dstree->SetBranchAddress("nph",&nph);
   dstree->SetBranchAddress("ndaughters",&ndaughters);
   dstree->SetBranchAddress("ndeposits",&ndeposits);
   dstree->SetBranchAddress("nusers",&nusers);
   dstree->SetBranchAddress("dau_id",&dau_id);
   dstree->SetBranchAddress("dau_pdg",&dau_pdg);
   dstree->SetBranchAddress("dau_pid",&dau_pid);
   dstree->SetBranchAddress("dau_process",&dau_process);
   dstree->SetBranchAddress("dau_time",&dau_time);
   dstree->SetBranchAddress("dau_ene",&dau_ene);
   dstree->SetBranchAddress("dau_x",&dau_x);
   dstree->SetBranchAddress("dau_y",&dau_y);
   dstree->SetBranchAddress("dau_z",&dau_z);
   dstree->SetBranchAddress("dau_r",&dau_r);
   dstree->SetBranchAddress("dau_px",&dau_px);
   dstree->SetBranchAddress("dau_py",&dau_py);
   dstree->SetBranchAddress("dau_pz",&dau_pz);
   dstree->SetBranchAddress("dep_pdg",dep_pdg);
   dstree->SetBranchAddress("dep_mat",dep_mat);
   dstree->SetBranchAddress("dep_time",dep_time);
   dstree->SetBranchAddress("dep_ene",dep_ene);
   dstree->SetBranchAddress("dep_step",dep_step);
   dstree->SetBranchAddress("dep_x",dep_x);
   dstree->SetBranchAddress("dep_y",dep_y);
   dstree->SetBranchAddress("dep_z",dep_z);
   dstree->SetBranchAddress("dep_r",dep_r);
   dstree->SetBranchAddress("userint1",userint1);
   dstree->SetBranchAddress("userint2",userint2);
   dstree->SetBranchAddress("userfloat1",userfloat1);
   dstree->SetBranchAddress("userfloat2",userfloat2);
   dstree->SetBranchAddress("userdouble0",userdouble0);
   dstree->SetBranchAddress("pe_time",&pe_time);
   dstree->SetBranchAddress("pe_pmt",&pe_pmt);
   dstree->SetBranchAddress("vpe_time",&vpe_time);
   dstree->SetBranchAddress("vpe_pmt",&vpe_pmt);
   dstree->SetBranchAddress("mupe_time",&mupe_time);
   dstree->SetBranchAddress("mupe_pmt",&mupe_pmt);
   dstree->SetBranchAddress("ph_volume",&ph_volume);
   dstree->SetBranchAddress("ph_pid",&ph_pid);
   dstree->SetBranchAddress("ph_wl",&ph_wl);
   dstree->SetBranchAddress("ph_x",&ph_x);
   dstree->SetBranchAddress("ph_y",&ph_y);
   dstree->SetBranchAddress("ph_z",&ph_z);
   dstree->SetBranchAddress("ph_time",&ph_time);
   dstree->SetBranchAddress("nTPC",&nTPC);
   dstree->SetBranchAddress("TPCtime",&TPCtime);
   dstree->SetBranchAddress("TPCene",&TPCene);
   dstree->SetBranchAddress("TPCf90",&TPCf90);
   dstree->SetBranchAddress("nAr40",&nAr40);
   dstree->SetBranchAddress("Ar40eneCalc",&Ar40eneCalc);
   dstree->SetBranchAddress("Ar40ene",&Ar40ene);
   dstree->SetBranchAddress("Ar40x",&Ar40x);
   dstree->SetBranchAddress("Ar40y",&Ar40y);
   dstree->SetBranchAddress("Ar40z",&Ar40z);
   dstree->SetBranchAddress("Ar40px",&Ar40px);
   dstree->SetBranchAddress("Ar40py",&Ar40py);
   dstree->SetBranchAddress("Ar40pz",&Ar40pz);
   dstree->SetBranchAddress("LSci_ndeps",LSci_ndeps);
   dstree->SetBranchAddress("LSci_time",LSci_time);
   dstree->SetBranchAddress("LSci_x0",LSci_x0);
   dstree->SetBranchAddress("LSci_y0",LSci_y0);
   dstree->SetBranchAddress("LSci_z0",LSci_z0);
   dstree->SetBranchAddress("LSci_x",LSci_x);
   dstree->SetBranchAddress("LSci_y",LSci_y);
   dstree->SetBranchAddress("LSci_z",LSci_z);
   dstree->SetBranchAddress("LSci_eneRaw",LSci_eneRaw);
   dstree->SetBranchAddress("LSci_ene",LSci_ene);
   dstree->SetBranchAddress("LSci_eneToF",LSci_eneToF);
   dstree->SetBranchAddress("nLSci_trg",&nLSci_trg);

//     This is the loop skeleton
//       To read only selected branches, Insert statements like:
//       
  dstree->SetBranchStatus("*",0);  // disable all branches
  dstree->SetBranchStatus("nAr40",1);        // activate branchname
  dstree->SetBranchStatus("Ar40ene",1);      
  //dstree->SetBranchStatus("Ar40px",1);     
  //dstree->SetBranchStatus("Ar40py",1);     
  //dstree->SetBranchStatus("Ar40pz",1);     
  dstree->SetBranchStatus("nusers",1);       
  dstree->SetBranchStatus("userdouble0",1);  

  dstree->SetBranchStatus("TPCtime",1);      
  dstree->SetBranchStatus("LSci_time",1);    
  dstree->SetBranchStatus("nLSci_trg", 1) ; 
  
   Long64_t nentries = dstree->GetEntries();


  TH2D * h2[9]  ; 
  TH2D * hEneLSC_theta[9]  ; 
  for (int i=0;i<9;++i)  { 
    h2[i] = new TH2D (Form ( "h2_%i",i) ,Form ( "h2_%i",i) ,  500, 0, 500, 500, 30, 80	) ;  
    h2[i]->SetLineColor(i+1) ; 
    h2[i]->SetMarkerStyle(7) ; 
    h2[i]->SetMarkerColor(i+1) ; 
    hEneLSC_theta[i] = new TH2D (Form ( "hEneLSC_theta_%i",i) ,Form ( "hEneLSC_theta_%i",i) ,  2000, 0, 200, 3140, 0, 3.14);  
  } 
  TH1D * hTheta = new TH1D ("htheta", "", 500, -1 , 4 ) ;  
  TH1D * hEne   = new TH1D ("hEne", "", 500, 0 , 500 ) ;  
  TH1D * hEneTOF= new TH1D ("hEneTOF", "", 500, 0 , 500 ) ;  
  hEneTOF->SetLineColor(2) ; 
  TH1D * hTPCt  = new TH1D ("hTPCt", "", 500, -100 , 400 ) ;  
  TH1D * hLSCt  = new TH1D ("hLSCt", "", 500, -100 , 400 ) ;  

  TH1D * hEneLSC[9], *hS1[9], *hS2[9] , *htheta[9]; 
  for (int i=0;i<9;++i)  { 
    hEneLSC[i] = new TH1D (Form ( "hEneLSc%i",i) ,Form ( "hEneLSc%i",i) ,  500, -100,400 ) ;  
    hS1[i]  = new TH1D (Form ( "hS1_%i",i) ,Form ( "hS1_%i",i) ,  100, 0, 500 ) ;  
    hS2[i]  = new TH1D (Form ( "hS2_%i",i) ,Form ( "hS2_%i",i) ,  100, 0, 15000 ) ;  
    htheta[i] = new TH1D (Form ( "htheta_%i",i) ,Form ( "htheta_%i",i) ,  515, -1, 4.15 ) ;  
    hEneLSC[i]->SetLineColor(i+1) ; 
    hS2[i]->SetLineColor(i+1) ; 
    hS1[i]->SetLineColor(i+1) ; 
    htheta[i]->SetLineColor(i+1) ; 
    htheta[i]->SetLineWidth(2) ; 
    hS2[i]->SetLineWidth(2) ; 
    hS1[i]->SetLineWidth(2) ; 
  } 

  for (Long64_t i=0; i<nentries;i++) {
     dstree->GetEntry(i);

     //if (nAr40 > 1 )  continue ; 
     double theta = acos ( userdouble0[3]  )    ;  

     if (theta == 0 )  continue ; 
     if (nLSci_trg>1)  continue ; 

     hTheta->Fill(theta) ;
     hEne->Fill(Ar40ene) ; 

     hTPCt->Fill(TPCtime) ; 

     int nch = -1 ; 
     if (TPCtime>35 && TPCtime<41) { 
       hEneTOF->Fill(Ar40ene) ; 
       for (int j=0;j<9;++j) {
	 if (LSci_time[j]>-1 ) hLSCt->Fill(LSci_time[j]) ; 
	 if (LSci_time[j]>20 && LSci_time[j]<26) { 
	   hEneLSC[j]->Fill(Ar40ene) ; 
	   htheta[j]->Fill(theta) ; 
	   nch = j ;   
	   //the following histo is what is going to be stored
	   hEneLSC_theta[j]->Fill(Ar40ene, theta) ;
	 }
       }
     }

     if (nch<0) continue ; 

     double S1, S2 ; 
     //preliminary. Do the S1 and S2 simulations on the fly. 
     //A dedicated function is below 
     S1S2 (Ar40ene, theta, S1, S2)  ; 
     h2[nch]->Fill(S1, 1.*S2/S1) ;   
     hS1[nch]->Fill(S1) ;   
     hS2[nch]->Fill(S2) ;   

  }

  TFile * fout = new TFile("spectra.root", "recreate") ; 
  for (int i=0;i<9;++i)  { 
    hEneLSC_theta[i]->Write() ; 
  }
  fout->Close() ; 
}















//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//This is the function that simulate S1 and S2, given a 2D distribution of 40Ar recoil energy vs angle
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
void simulate (double myStat=-1) {
  
  TFile * fin = TFile::Open("spectra.root") ; 
  TH2D * hEneLSC_theta[8], *hQY[8] ; 
  TH2D * h2[8]  ; 
  TH1D * hS1[8]  ; 
  TH1D * hTh[8]  ; 
  TH1D * hEne[8]  ; 
  TH1D * hS2[8]  ; 
  
  for (int i=0;i<8;++i)  { 
    hEneLSC_theta[i] = (TH2D*) fin->Get(Form("hEneLSC_theta_%i",i) ) ; 
    h2[i]  = new TH2D (Form ( "S2_vs_S1_%i",i) ,Form ( "S2vsS1_%i",i) ,  500, 0, 500, 300, 0, g2	) ;  
    h2[i]->GetYaxis()->SetTitle("S2/S1") ; 
    h2[i]->GetXaxis()->SetTitle("S1 [PE]") ; 
    hQY[i] = new TH2D (Form ( "hQY_%i",i) ,Form ( "inization yield for validation %i",i) ,  500, 0, 500, 400, 0, 4	) ;  
    hQY[i]->GetYaxis()->SetTitle("ionization yield [e-/keV]") ; 
    hQY[i]->GetXaxis()->SetTitle("Energy [keV]") ; 
    hS1[i] = new TH1D (Form ( "hS1_%i",i) ,Form ( "hS1_%i",i) ,  200, 0, 500	   ) ;  
    hS1[i]->GetXaxis()->SetTitle("S1 [PE]") ; 
    hEne[i]= new TH1D (Form ( "hEne_%i",i) ,Form ( "hEne_%i",i) ,  200, 0, 200	   ) ;  
    hEne[i]->GetXaxis()->SetTitle("Energy [keV_NR]") ; 
    hS2[i] = new TH1D (Form ( "hS2_%i",i) ,Form ( "hS2_%i",i) ,  200, 0, g2*200	   ) ;  
    hS2[i]->GetXaxis()->SetTitle("S2 [PE]") ; 
    hTh[i] = new TH1D (Form ( "hTh_%i",i) ,Form ( "hTh_%i",i) ,  100, 0, 3.14   ) ;  
    hTh[i]->GetXaxis()->SetTitle("Azimuthal angle [rad]") ; 
    
    h2[i]->SetLineColor(i+1) ; 
    h2[i]->SetMarkerStyle(7) ; 
    h2[i]->SetMarkerColor(i+1) ; 
    hS1[i]->SetLineColor(i+1) ; 
  } 
  
  
  //decide aboout the statistics. At the moment, it is limited by the MC statistics
  double MaxStat = 5 * hEneLSC_theta[1]->GetEntries() ; 
  if (myStat>0 && myStat<MaxStat) MaxStat = myStat ; 
      	
  ran = new TRandom3 ; 
  double S1, S2 ; 
  
  double Ar40ene, theta ; 
  
  //limit the plots to 5 LSci for the moment   
  for (int nHisto=1;nHisto<6;++nHisto) {
  
    for (int i=0;i<MaxStat;++i){
    
      //Get a pair of 40Ar energy and angle from the MC  
      hEneLSC_theta[nHisto]->GetRandom2(Ar40ene, theta) ; 
      hEne[nHisto]->Fill(Ar40ene) ; 
      hTh[nHisto]->Fill(theta) ; 
      
      //simulate S1 and S2 
      S1S2 (Ar40ene, theta, S1, S2)  ; 
            
      hS1[nHisto]->Fill(S1) ; 
      hS2[nHisto]->Fill(S2) ; 
      hQY[nHisto]->Fill(Ar40ene, S2/(Ar40ene*LEff)/g2) ; 
      h2 [nHisto]->Fill(S1,1.*S2/S1) ; 
    }  
    
    double meanS1 = hS1[nHisto]->GetMean();   
    double meanS2 = hS2[nHisto]->GetMean(); 
     
  }
  
  
  
  
  
  //results - 2D distributions S2/S1 vs S1   
  TCanvas * c2 = new TCanvas("c2", "S2 vs S1 for the first 5 LSCi", 800, 600)   ; 
  
  h2[1]->Draw("*") ;   
  h2[2]->Draw("same,*") ;   
  h2[3]->Draw("same,*") ;   
  h2[4]->Draw("same,*") ;   
  h2[5]->Draw("same,*") ; 

  



  //S1 for a few LSci 
  TCanvas * cS1 = new TCanvas("cS1", "S1", 800, 600)    ; 
  hS1[1]->Draw() ; 
  hS1[2]->Draw("same") ; 
  hS1[3]->Draw("same") ; 
  hS1[4]->Draw("same") ; 
  hS1[5]->Draw("same") ; 

  TF1 * ffit = new TF1 ("ffit", "[0]* ( [1]*[2]*exp(-x/[2]) + (1-[1]) * TMath::Gaus(x,[3],[4]) ) ", 55, 200) ; 
  ffit->SetParameters(1000, .3 ,30,200,22) ;
  ffit->SetParNames("A", "ratio" ,"tau","MEAN S1","sigma S1") ;


  //validation - the SCENE effect is 10% difference between ortogonal tracks. 
  //this TGraph contains the mean vaklues of the S1 peaks
  TCanvas * cS1bis = new TCanvas("cS1bis", "S1 vs theta", 800, 600)    ; 
  TGraphErrors * S1vsTheta = new TGraphErrors ; 
  hS1[1]->Fit(ffit)  ; 
  S1vsTheta->SetPoint(0, 90 , ffit->GetParameter(3) ) ; 
  S1vsTheta->SetPointError(0, 0,  ffit->GetParError(3) ) ;   
  double meanS1_1 =  ffit->GetParameter(3) ; 
  double sigmaS1_1 =  ffit->GetParameter(4) ; 
  hS1[2]->Fit(ffit)  ; 
  S1vsTheta->SetPoint(1, 60 , ffit->GetParameter(3) ) ; 
  S1vsTheta->SetPointError(1, 0,  ffit->GetParError(3) ) ;   
  hS1[3]->Fit(ffit)  ; 
  S1vsTheta->SetPoint(2, 40 , ffit->GetParameter(3) ) ; 
  S1vsTheta->SetPointError(2, 0,  ffit->GetParError(3) ) ;   
  hS1[4]->Fit(ffit)  ; 
  S1vsTheta->SetPoint(3, 20 , ffit->GetParameter(3) ) ; 
  S1vsTheta->SetPointError(3, 0,  ffit->GetParError(3) ) ;   
  hS1[5]->Fit(ffit)  ; 
  S1vsTheta->SetPoint(4, 0 , ffit->GetParameter(3) ) ; 
  S1vsTheta->SetPointError(4, 0,  ffit->GetParError(3) ) ;   
  double meanS1_5 = ffit->GetParameter(3) ; 
  double sigmaS1_5 = ffit->GetParameter(4) ; 
    
  S1vsTheta->Draw("a*") ; 
  
  //Additional plots 
  TCanvas * cE = new TCanvas  ; 
  hEne[1]->Draw() ; 
  TCanvas * cS2 = new TCanvas  ; 
  hS2[1]->Draw() ; 
  
  //validation - ionization yield, to be compared with what is extracted from DS50/ARIS 
  TCanvas * cQY = new TCanvas  ; 
  hQY[1]->Draw("colz") ; 
  
  //validation - ionization yield, to be compared with what is extracted from DS50/ARIS 
  TCanvas * cQY3 = new TCanvas  ; 
  hQY[3]->Draw("colz") ; 
    
  //validation - ionization yield, to be compared with what is extracted from DS50/ARIS 
  TCanvas * cQY2 = new TCanvas  ; 
  hQY[5]->Draw("colz") ; 








  //Fitting procedure to be standardized
  //results - 2D fit S2/S1 vs S1    and computation of "distance" 
  TCanvas * c2cc= new TCanvas("c2cc", "S2 vs S1 for the first LSC", 800, 600)   ; 
  TH2D * hbi2 = (TH2D*) h2[1]->Clone("hbi2") ; 
  hbi2->RebinX(4) ; 
  hbi2->RebinY(4) ; 
  
  double meanS2_1 = 1130/2. * (0.1 ) * g2; 

  TF2 *f2 = new TF2("f2","bigaus(0)",  meanS1_1*.7,   meanS1_1*1.4,   meanS2_1*.25/meanS1_1 ,    meanS2_1/meanS1_1*1.8  ); 

  f2->SetParameters(MaxStat,meanS1_1,sigmaS1_1,meanS2_1/meanS1_1,meanS2_1*0.1/meanS1_1, -0.6); 
  f2->SetParLimits(0,0.01*MaxStat,10*MaxStat) ; 
  f2->SetLineColor(4) ; 
  f2->SetParLimits(1,100,300) ; 
  f2->SetParLimits(2,0,300) ; 
  f2->SetParLimits(5,-1,1) ; 
  hbi2->Fit(f2,"RL"); 

  hbi2->Draw("*"); 
  f2->Draw("same") ;


  double meanS2_5 = 1130/2. *(1 -  (0.9+aparameter )) * g2; 
  TCanvas * c5cc= new TCanvas("c5cc", "S2 vs S1 for the fifth LSC", 800, 600)   ; 
  TF2 *f25 = new TF2("f25","bigaus",  meanS1_5*.7,   meanS1_5*1.4,   meanS2_5*.25/meanS1_5 ,    meanS2_5/meanS1_5*3.4); 
  TH2D * hbi5 = (TH2D*) h2[5]->Clone("hbi5") ; 
  hbi5->RebinX(4) ; 
  hbi5->RebinY(4) ; 
  f25->SetParameters(MaxStat,meanS1_5,sigmaS1_5,meanS2_5/meanS1_5,meanS2_5*0.1/meanS1_5, -.6);  
  f25->SetParLimits(0,0.01*MaxStat,100*MaxStat) ; 
  f25->SetParLimits(1,100,300) ; 
  f25->SetParLimits(2,0,300) ; 
  f25->SetParLimits(5,-1,1) ; 
  f25->SetLineColor(1) ; 
  hbi5->Fit(f25,"RL"); 
  hbi5->Draw("*"); 
  f25->Draw("same") ;  
  
  cout <<endl << endl << endl <<"DISTANCE" << endl ; 
  TF2 *f2final = new TF2("f2final","bigaus(0) * bigaus(6) ",0,400,0.001,meanS2_1/meanS1_1*3); 
  for (int i=0;i<6;++i)  f2final->SetParameter(i,  f25->GetParameter(i)); 
  for (int i=6;i<12;++i) f2final->SetParameter(i,  f2->GetParameter(i-6)); 
  f2final->SetParameter(0,  1); 
  f2final->SetParameter(6,  1); 

  TCanvas *cfinal = new TCanvas ; 
  f2final->Draw() ; 
  f2->Draw("same") ; 
  f25->Draw("same") ; 
  cout << "overlap:     " <<  f2final -> Integral(0,400,0.001,meanS2_1/meanS1_1*3) << endl ; 
  cout << "overlap*1e5: " << f2final -> Integral(0,400,0.001,meanS2_1/meanS1_1*3)/1e-5 << endl ; 
    
}




//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//Main
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
void red_model () {
  gStyle->SetOptStat(0) ; 
  gStyle->SetNumberContours (8) ; 
  prepare_MC_samples() ; 
  simulate() ; 
  
}
