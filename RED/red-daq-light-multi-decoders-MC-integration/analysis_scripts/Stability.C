#include "Utils.C"
#include "Math/SpecFunc.h"
#include "Math/Functor.h"
#include "Fit/Fitter.h"
#include "TF2.h"

#define UNDEF -100

//  E=28 MeV

//
// Banana runs with Low Threshold   (Already compared in calibSi.C )
//          Coll   Run
//           0    624 
//           1    676

//
// Banana runs with High Threshold (no access to the LiC elastic peak, only low energy LiP and the whole Be band )
//          Coll   Run
//           0    645/661 
//           1    672/680 and  715/717/719/720

// Goal       :  estimate beam variability and validate cuts on LowBe blob
// Procedure  :  use BeBand SiE, SiTelEnergy and SiMonitor for the two collimators and different runs

// Conclusions : 1) Beam Energy seems stable within 0.2 MeV according to Be Blob displacement.
//               2) Monitor LiC elastic peak seems to be jumping uncorrelated with the Be Blob displacement. Problem with Monitor Amplifier, or high rate of 
//                                       low energy BG that displaces baseline?
//               3) BeEx total rate is ~30-50% of LowBe. Using energy cuts  18.5-20.5 (Coll 0) and 18.-19.5 (Coll 1), contamination should be much smaller than 10%
void PlotEnergies()  
{
  SetStyle();

  int color[nChannels]={kGreen+2, kOrange, kRed, kBlue, kMagenta  };
  
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,1.*wh);
  cc[canvas]->Divide(3,1);
  cc[canvas]->Draw();  
  l.SetTextSize(0.03);

  for (int ipad=0;ipad<3;ipad++)
    {
      cc[canvas]->cd(ipad+1);
  
      for (int icase=0;icase<(ipad==1?1:2);icase++)    
	for (int ie=0;ie<2;ie++)
	  for (int SelForward=1;SelForward>=0;SelForward--)
	    {
	      int ichan;
	      if ( ipad==1 ) ichan=0;
	      else           ichan=(icase==0?2:4);
	      
	      double beamEne=(ie==0?27.8:28.2);
	      TGraph *tg=new TGraph();
	      if ( ipad<2 ) GetEnergy(beamEne,ichan,SelForward,tg,(ipad==0?7.:25.));
	      else          GetThetaNeutron(beamEne,ichan,SelForward,tg,(ipad==1?25:7.));
	      if ( tg->GetN()==0 ) continue;
	      
	      tg->SetMaximum((ipad<2?28.5:55.));
	      tg->SetMinimum((ipad<2?16.5: 0.));
 	      tg->SetLineColor(color[ichan]);
	      tg->SetLineStyle(ie+1);
	      tg->GetXaxis()->SetTitle("#theta_{Si Detector}");
	      tg->GetYaxis()->SetTitle("KE [MeV]");
	      tg->DrawClone((icase==0&& ie==0&&SelForward==1?"AL":"L"));
	      l.SetTextColor(color[ichan]);
	      l.DrawLatex((ipad==2?0.6:0.2),(ipad==2?0.35:0.55)-double(ichan)*0.03,ChannelName[ichan].c_str());
	    }
      l.SetTextColor(kBlack);
      l.DrawLatex(0.5,0.85,Form("E_{beam}=%5.2lf+-0.2 MeV",28.));
    }

  canvas++;
  
}


void Stability(int Set)  // [0]  Coll 0  (>=12 hour diff)  [1] Coll 1 (>=12 hour diff)  [2]  Coll 1 ( < 12 hour diff)
{
  SetStyle();
  std::string NameSet[3]={ "Coll. 0 (Runs >10 hour apart)", "Coll. 1 (Runs >10 hour apart)", "Coll. 1 (Runs <10 hour apart)"};
  int nF_in[3]={ 2, 3, 4 };
  int run_i[3][4]={ {645,661,UNDEF,UNDEF}, { 672, 680, 715, UNDEF }, { 715,717,719,720 } };
  
  int color[4]={kBlue, kRed, kGreen+2, kOrange };
  
  nF=nF_in[Set];
  for (int iF=0;iF<nF;iF++)
    {
      iF_plot.push_back(iF);
      frun_i[iF]=run_i[Set][iF]; ftime_i[iF]=10.;  IsTPC=false;

      TTree *t;
      std::string FileName=Form("reco/run_%d.root",run_i[Set][iF]);
      TFile *f=new TFile(FileName.c_str());
      f->GetObject("reco",t);
      Summary(t,iF);
      SetHistos(t,iF);
    }



  //----------------
  // Get references
  //-----------------
  
  const int nF_ref=2;
  TH1F   *hE_Ref[2][nF_ref];
  TH1F   *hSiE_Ref[2][nF_ref];
  TH2F   *hBananas_Ref[nF_ref];
  
  TFile *f1 = new TFile("bananas.root", "READ" );
  for (int iF=0;iF<nF_ref;iF++)  
    {
      std::string name=Form("h2D%d",iF);
      if ( f1->Get(name.c_str())==NULL )    hBananas_Ref[iF]=new TH2F();
      else                                  hBananas_Ref[iF]= (TH2F*)   f1->Get(name.c_str());
      hBananas_Ref[iF]->SetDirectory(0);

      for (int iBand=0;iBand<2;iBand++)
	{
	  std::string name=Form("h1D_%d%d",iF,iBand);
	  if ( f1->Get(name.c_str())==NULL )   {printf("NULL\n"); hE_Ref[iBand][iF]=new TH1F();}
	  else                                  hE_Ref[iBand][iF]=(TH1F*) f1->Get(name.c_str());
	  hE_Ref[iBand][iF]->SetDirectory(0);
	  hE_Ref[iBand][iF]->SetLineColor(kBlack);
	  hE_Ref[iBand][iF]->SetLineStyle(2);
	  hE_Ref[iBand][iF]->SetLineWidth(2);
	}
      for (int iBand=0;iBand<2;iBand++)
	{
	  std::string name=Form("h1D_SiE_%d%d",iF,iBand);
	  if ( f1->Get(name.c_str())==NULL )    hSiE_Ref[iBand][iF]=new TH1F();
	  else                                  hSiE_Ref[iBand][iF]=(TH1F*) f1->Get(name.c_str());
	  hSiE_Ref[iBand][iF]->SetDirectory(0);
	  hSiE_Ref[iBand][iF]->SetLineColor(kBlack);
	  hSiE_Ref[iBand][iF]->SetLineStyle(2);
	  hSiE_Ref[iBand][iF]->SetLineWidth(2);
	}
    }
  f1->Close();

  //----------------------------------------
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,1.*wh);
  cc[canvas]->Divide(3,1);
  cc[canvas]->Draw();  
      

  for (int ipad=0;ipad<3;ipad++)
    {
      cc[canvas]->cd(ipad+1);

      std::string xaxis[3]={ "SiE [MeV]", "SiTel [MeV] (E+dE)","SiMonit [MeV]" };

      double xmin[3]={ 11.5,18.,23. }, xmax[3]={ 20., 25., 24.5};
      TH2D *hFrame=new TH2D(Form("hX%d",ipad),"",1000,xmin[ipad],xmax[ipad] ,1000,0.,0.12);
      hFrame->GetXaxis()->SetTitle(xaxis[ipad].c_str());
      hFrame->GetYaxis()->SetTitle("A.U");
      hFrame->DrawCopy();

      for (int iF=0;iF<nF;iF++)
	{
	  int color_i=color[iF];
	  
	  TH1F *h0;
	  if (ipad==0 )      h0=hE[iF]; 
	  else if (ipad==1 ) h0=hE_C[iF]; 
	  else if (ipad==2 ) h0=hMon[iF]; 
		  		  
	  h0->SetLineColor(color_i);
	  h0->Scale(1./h0->Integral());
	  h0->SetMaximum(0.12);
	  h0->DrawCopy("HISTSAME");

	  // Plot Ref
	  if ( ipad<2 )
	    {
	      int iref=(Set==0?0:1);
	      if (ipad==0 ) hSiE_Ref[1][iref]->DrawCopy("HISTSAME");			
	      else          hE_Ref[1][iref]->DrawCopy("HISTSAME");			
			
	    }
	}//iF
      l.SetTextColor(kRed);
      l.SetTextSize(0.03);
      l.DrawLatex(0.2,0.85,NameSet[Set].c_str());     
      l.DrawLatex(0.2,0.80,"E=28 MeV");
      l.SetTextSize(0.035);
      l.SetTextColor(kBlack);
      l.DrawLatex(0.2,0.75,"Black dashed is ref");

    }//ipad
}
 
