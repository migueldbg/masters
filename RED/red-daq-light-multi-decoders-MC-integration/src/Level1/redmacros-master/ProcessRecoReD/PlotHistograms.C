#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TRint.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"
#include "TKey.h"
#include "TMath.h"
#include "TPDF.h"
#include <TGaxis.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TGraphErrors.h>
#include <Rtypes.h>
#include "TH2D.h"
#include "THStack.h"
#include "TPaveStats.h"



using namespace std;

#ifndef __CINT__
  int main(int argc, char **argv);
#endif

TRint* theApp;
TCanvas *c1;
void PlotdT_Si_LSci_LSciPSD (TFile *file);
void PlotS2Efficiency(TFile *file);
void PlotS2vsS1(TFile *file);
void PlotS1S2Projection(TFile *file);
void PlotF90vsS1(TFile *file);
void CompareS1_S1zcorr(TFile *file);
void CompareS2corrections(TFile *file);
void PlotHistograms(TString fInName, TString outDirname="");
TF1* FitS1Zdependency(TH2* h2); // Fit S1 vs t_drift histogram and return fit function
void FitElectronlifetime(TH2* h2);
void DrawS2Fraction(TFile *file);
void DrawChannelMap(TH3 *h3);
void Fit_s1_tdrift_with_s2maxchan_slice(TH3 *h3); // Function to fit each channel and plot by inner and outer channels
void Fit_s1_tdrift_with_group(TH3 *h3); // Function to get S1 z-correction by group
void Get_S2_MaxS2Chan_Correction(TH2* h2);
void DrawSlice3D(TH3 *h3);

Bool_t IsInnerChannel(Int_t ch){
	return ch==0 || ch==3 || ch==5 || ch==12 || ch==14 || ch==15 || ch==18 || ch==20;
}

void PlotHistograms(TString fInName, TString outDirname) {
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  TGaxis::SetMaxDigits(3);
  gStyle->SetTitleOffset(1.2, "YZ");
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatY(0.88);
  gStyle->SetStatH(0.075);

//  gStyle->SetPalette(kBird);
  gROOT->SetBatch(kTRUE);

  c1 = new TCanvas("c1","c1",1000, 750);
  c1->cd();
  c1->Draw();
  c1->SetGrid();

  cout<<"Input file: "<<fInName.Data()<<endl;
  TFile *file = new TFile(fInName.Data());
  if(!file) cout<<"file: "<<fInName.Data()<<" is not found."<<endl;

  if (outDirname.CompareTo("")!=0) {
	  gSystem->Exec(Form("mkdir %s", outDirname.Data()));
	  gSystem->cd(outDirname.Data());
  }

  Int_t nhists = 0;
  TIter next(file->GetListOfKeys());
//  file->ls();
  TH1 *h;
  TObject* obj;
  TKey* key;
  while ((key = (TKey*) next()))
  {
      obj = (TObject*)key->ReadObj();
      if (obj->InheritsFrom(TH1::Class())){
          h = (TH1*) obj;
          TString hname(h->GetName());
          if(h->InheritsFrom(TH2::Class())){
        	  if(h->GetEntries()==0) continue;
        	  Bool_t setlogz(false);
              if(hname.Contains("hStartTime_ChNum")) setlogz=true;
              if(hname.Contains("hS2_Am_s2maxchan")){
            	  setlogz=true;
                  Get_S2_MaxS2Chan_Correction((TH2*)h);
              }
//              Draw2DHist(h);
              h->Draw("colz");

              if(setlogz) c1->SetLogz(1);
              c1->Update();
//              theApp->Run();
              c1->Print(Form("%s.pdf", h->GetName()));
              if(setlogz) c1->SetLogz(0);

          } else if(h->InheritsFrom(TH3::Class())){
//              Draw3DHist(h);
//              h->Draw("iso");
//              c1->Update();
//              theApp->Run();
//              c1->Print(Form("%s.pdf", h->GetName()));
//              DrawSlice3D((TH3*)h);
        	  if(hname.Contains("hYBary_XBary_chnum")) DrawChannelMap((TH3*)h);
        	  if(hname.Contains("hS1_tdrift_s2maxchan")) {
        		  Fit_s1_tdrift_with_s2maxchan_slice((TH3*)h);
        		  Fit_s1_tdrift_with_group((TH3*)h);
        	  }
          } else {//TH1
        	  if(h->GetEntries()==0) continue;
              h->Draw();
              c1->Update();
//              theApp->Run();
              c1->Print(Form("%s.pdf", h->GetName()));
          }
          if(hname.Contains("hS2ovS1_tdrift")){
        	  FitElectronlifetime((TH2*) h);
          } else if(hname.EqualTo("hS2_Am_tdrift")){
        	  FitElectronlifetime((TH2*) h);
          }else if(hname.EqualTo("hS1_tdrift")){
        	  FitS1Zdependency((TH2*) h);
          }else if(hname.EqualTo("hS1zcorr_tdrift")){
        	  FitS1Zdependency((TH2*) h);
          }


          delete h;
          nhists++;
      } else {
          cout<<"unknown class: "<<obj ->GetName()<<endl;
      }

  }

  cout<<"Now customize histograms........"<<endl;
  PlotdT_Si_LSci_LSciPSD(file);
  PlotS2Efficiency(file);
  PlotS2vsS1(file);
  PlotF90vsS1(file);
  PlotS1S2Projection(file);
  DrawS2Fraction(file);
  CompareS1_S1zcorr(file);
  CompareS2corrections(file);

  if (outDirname.CompareTo("")!=0) gSystem->cd("-");


//fInName.Replace(fInName.Length()-5, 5, "");
//fInName.Replace(0,fInName.Last('/')+1, "");// find the last'/' and replace w/ '/Hist'
//TPDF *mPDF = new TPDF(Form("QA_%s.pdf", fInName.Data()));
//
//  mPDF->Close();

#ifndef __CINT__
//  theApp->Run();
//  delete theApp;
#endif /* __CINT __ */

  delete c1;

}

const Int_t NLSciChan = 6;
const TString LSciName[NLSciChan] = {"LSci0", "LSci1", "LSci3_7", "LSci4_6", "LSci5", "LSci8"};
const TString LSciAngle[NLSciChan] = {"4.3", "90", "40", "20", "0", "90"};
const Int_t LSci_ChID[NLSciChan] = {0, 1, 3, 4, 5, 8};
const Double_t LSci_PSD_threshold[NLSciChan] = {0.16, 0.14, 0.14, 0.12, 0.15, 0.15};

void PlotdT_Si_LSci_LSciPSD (TFile *file) {
	TLegend* leg = new TLegend(0.6, 0.55, 0.85, 0.75);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
	leg->SetFillColor(0);

	TLegend* leg2 = new TLegend(0.6, 0.5, 0.85, 0.85);
	leg2->SetBorderSize(0);
	leg2->SetTextSize(0.03);
	leg2->SetFillColor(0);


	TString hname = "hdT_Si_LSci_LSciPSD";
	vector<TH2D*> hdT_Si_LSci_LSciPSD;
	Bool_t isStack(false);
	THStack *stack = new THStack("Si_LSci_TimeCoinci", "Time Coincidence");
	for(Int_t ilsci=NLSciChan-1; ilsci>=0; ilsci--){
        hname = Form("hdT_Si_LSci_LSciPSD_%d", ilsci);
        TH2D* htmp2D = (TH2D*) file->Get(hname.Data());
        TH1D* htmp = htmp2D->ProjectionX(Form("%s_all", htmp2D->GetName()));
        htmp2D->SetAxisRange(0, LSci_PSD_threshold[ilsci],"Y");
        TH1D* htmp_1 = htmp2D->ProjectionX(Form("%s_gamma", htmp2D->GetName()));
        htmp2D->SetAxisRange(LSci_PSD_threshold[ilsci], 0.5,"Y");
        TH1D* htmp_2 = htmp2D->ProjectionX(Form("%s_neutron", htmp2D->GetName()));
        stack->Add(htmp);

    	int binmax = htmp_1->GetMaximumBin();
    	Double_t y_max = TMath::Max(htmp_1->GetMaximum(), htmp_2->GetMaximum());
    	htmp_1->SetAxisRange(0,1.2*y_max,"Y");
        htmp_1->SetLineColor(kBlack);
        htmp_1->DrawClone("HIST");
        htmp_2->SetLineColor(kRed);
        htmp_2->DrawClone("HISTsame");
    	leg->AddEntry(htmp_1,Form("Gamma (PSD_{LSci}<%3.2f)", LSci_PSD_threshold[ilsci]), "L");
    	leg->AddEntry(htmp_2,Form("Neutron (PSD_{LSci}>%3.2f)", LSci_PSD_threshold[ilsci]), "L");
    	if (isStack) htmp->SetLineWidth(0);
    	if(!isStack && ilsci==0) htmp->Scale(0.05); // scale channel 0
        isStack? leg2->AddEntry(htmp, Form("%s at %s#circ", LSciName[ilsci].Data(), LSciAngle[ilsci].Data()), "F"):
        		 leg2->AddEntry(htmp, ilsci!=0? Form("%s at %s#circ", LSciName[ilsci].Data(), LSciAngle[ilsci].Data()): Form("%s at %s#circ (x0.05)", LSciName[ilsci].Data(), LSciAngle[ilsci].Data()), "L");
        if(ilsci==0) stack->SetHistogram(htmp);

    	leg->Draw();

    	c1->Update();
//        theApp->Run();
    	c1->Print(Form("hdT_Si_LSci_Proj_%d.pdf", ilsci));
    	leg->Clear();
	}

	isStack? stack->Draw("PFC"):stack->Draw("nostack HIST PLC");
//	c1->SetLogy();
	leg2->Draw();

	c1->Update();
//         	 theApp->Run();
	isStack? c1->Print("hdT_Si_LSci_stack.pdf") : c1->Print("hdT_Si_LSci_nostack.pdf");

	delete stack;
}


void PlotS2Efficiency(TFile *file) {
	TLegend* leg = new TLegend(0.5, 0.15, 0.85, 0.5);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
	leg->SetFillColor(0);

	TString hname = "hTBAS1All";
	TH1D* hTBAS1All = (TH1D*) file->Get(hname.Data());
	if(!hTBAS1All) cout<<"hist: "<<hname.Data()<<" is not found."<<endl;
	hTBAS1All->SetMarkerColor(kGray);
	hTBAS1All->Draw();
    leg->AddEntry(hTBAS1All, "Events with at least one cluster", "P");
	TH1D* hTBAS1 = (TH1D*) file->Get("hTBAS1");
	hTBAS1->SetMarkerStyle(20);
	hTBAS1->Draw("P same PLC");
    leg->AddEntry(hTBAS1, "Events with S2", "P");
//	TH1D* hTBAS1 = (TH1D*) file->Get("hTBAS1");
//	hTBAS1->SetMarkerStyle(20);
//	hTBAS1->Draw("P same PLC");
//    leg->AddEntry(hTBAS1, "Events with S2", "P");
	leg->Draw();

	c1->Update();
//         	 theApp->Run();
	c1->Print("S1TBA_Comp.pdf");

	delete leg;
}

void PlotS2vsS1(TFile *file) {
	TLegend* leg = new TLegend(0.5, 0.15, 0.85, 0.5);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
	leg->SetFillColor(0);

	TString hname = "hS2_S1_SiTelCoinci";
	TH1F* hSiTelCoinci = (TH1F*) file->Get(hname.Data());
	if(!hSiTelCoinci) cout<<"hist: "<<hname.Data()<<" is not found."<<endl;
	hSiTelCoinci->SetMarkerColor(kGray);
	hSiTelCoinci->Draw();

	vector<TH2D*> h3Coinci_LSci;
	Int_t total(0);
	for(Int_t ilsci=0; ilsci<NLSciChan; ilsci++){
        hname = Form("hS2_S1_3Coinci_LSci_%d", ilsci);
        h3Coinci_LSci.push_back((TH2D*) file->Get(hname.Data()));
        h3Coinci_LSci[ilsci]->SetMarkerStyle(20);
        h3Coinci_LSci[ilsci]->Draw("P same PMC");
        total+=h3Coinci_LSci[ilsci]->GetEntries();
    	leg->AddEntry(h3Coinci_LSci[ilsci], Form("%s at %s#circ (%d evt)", LSciName[ilsci].Data(), LSciAngle[ilsci].Data(),(Int_t) h3Coinci_LSci[ilsci]->GetEntries()), "P");
	}
	leg->SetHeader(Form("Total: %d evt",total));
	leg->Draw();

	c1->Update();
//         	 theApp->Run();
	c1->Print("S2vsS1_all.pdf");
	delete leg;
}
void PlotS1S2Projection(TFile *file) {
	TLegend* leg = new TLegend(0.5, 0.55, 0.85, 0.85);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
	leg->SetFillColor(0);

	Int_t nrebin(2);
	THStack *stack = new THStack("S1_3Coinci", "Triple Coincidence");
	THStack *stack2 = new THStack("S2_3Coinci", "Triple Coincidence");
	TString hname;
	vector<TH2D*> h3Coinci_LSci;
	Bool_t isStack(true);
	for(Int_t ilsci=0; ilsci<NLSciChan; ilsci++){
        hname = Form("hS2_S1_3Coinci_LSci_%d", ilsci);
        h3Coinci_LSci.push_back((TH2D*) file->Get(hname.Data()));
        TH1* htmp = h3Coinci_LSci[ilsci]->ProjectionX();
        TH1* htmp2 = h3Coinci_LSci[ilsci]->ProjectionY();

        htmp->Rebin(nrebin);
        htmp->SetAxisRange(0, 600, "X");
        htmp2->Rebin(nrebin);
        isStack? leg->AddEntry(htmp, Form("%s at %s#circ (%d evt)", LSciName[ilsci].Data(), LSciAngle[ilsci].Data(),(Int_t) h3Coinci_LSci[ilsci]->GetEntries()), "F"):
        		 leg->AddEntry(htmp, Form("%s at %s#circ (%d evt)", LSciName[ilsci].Data(), LSciAngle[ilsci].Data(),(Int_t) h3Coinci_LSci[ilsci]->GetEntries()), "L");
    	stack->Add(htmp);
    	stack2->Add(htmp2);
        if(ilsci==0) stack->SetHistogram(htmp);
        if(ilsci==0) stack2->SetHistogram(htmp2);

	}
//	stack->Draw("nostack PLC");//if(stack) hs.Draw("same ][");
	isStack? stack->Draw("PFC"):stack->Draw("nostack PLC");
	leg->Draw();

	c1->Update();
//         	 theApp->Run();
	c1->Print("S1_Proj_3Coinci.pdf");

	isStack? stack2->Draw("PFC"):stack2->Draw("nostack PLC");
	leg->Draw();
	c1->Update();
//         	 theApp->Run();
	c1->Print("S2_Proj_3Coinci.pdf");

	delete leg;
	delete stack;
	delete stack2;
}

void PlotF90vsS1(TFile *file) {
	TLegend* leg = new TLegend(0.5, 0.5, 0.85, 0.85);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
	leg->SetFillColor(0);

	TString hname = "hF90_S1_SiTelCoinci";
	TH1F* hSiTelCoinci = (TH1F*) file->Get(hname.Data());
	if(!hSiTelCoinci) cout<<"hist: "<<hname.Data()<<" is not found."<<endl;
	hSiTelCoinci->SetMarkerColor(kGray);
	hSiTelCoinci->Draw();

	vector<TH2D*> h3Coinci_LSci;
	for(Int_t ilsci=0; ilsci<NLSciChan; ilsci++){
        hname = Form("hF90_S1_3Coinci_LSci_%d", ilsci);
        h3Coinci_LSci.push_back((TH2D*) file->Get(hname.Data()));
        h3Coinci_LSci[ilsci]->SetMarkerStyle(20);
        h3Coinci_LSci[ilsci]->Draw("P same PMC");
    	leg->AddEntry(h3Coinci_LSci[ilsci], Form("%s at %s#circ (%d evt)", LSciName[ilsci].Data(), LSciAngle[ilsci].Data(),(Int_t) h3Coinci_LSci[ilsci]->GetEntries()), "P");
	}
	leg->Draw();

	c1->Update();
//         	 theApp->Run();
	c1->Print("F90vsS1_all.pdf");
	delete leg;
}

void CompareS1_S1zcorr(TFile *file) {
	TLegend* leg = new TLegend(0.15, 0.35, 0.5, 0.5);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
	leg->SetFillColor(0);
	leg->SetHeader("relative resolution #sigma/#mu");

	TString hname = "hS1";
	TH1D* hS1 = (TH1D*) file->Get(hname.Data());
	if(!hS1) cout<<"hist: "<<hname.Data()<<" is not found."<<endl;

	hname = "hS1zcorr";
	TH1D* hS1zcorr = (TH1D*) file->Get(hname.Data());
	if(!hS1zcorr) cout<<"hist: "<<hname.Data()<<" is not found."<<endl;

	Double_t y_max = TMath::Max(hS1->GetMaximum(), hS1zcorr->GetMaximum());
	hS1->SetAxisRange(0,1.2*y_max,"Y");
	hS1->SetAxisRange(300, 900);
	hS1->SetTitle("S1 z-correction comparison (Am peak); S1 & S1zcorr [PE]; Number of Events");
	hS1->SetMarkerStyle(21);
	hS1->Sumw2();
	hS1->Draw();

	hS1->Fit("gaus", "EMRQ", "", 570, 720);
	TF1* func1 = (TF1*)hS1->GetFunction("gaus");
	hS1zcorr->SetLineColor(kRed);
	hS1zcorr->SetMarkerColor(kRed);
	hS1zcorr->SetMarkerStyle(20);
	hS1zcorr->Sumw2();
	hS1zcorr->Draw("sames");
	hS1zcorr->Fit("gaus", "EMRQ", "", 570, 720);
	TF1* func2 = (TF1*)hS1zcorr->GetFunction("gaus");
	func2->SetLineColor(kRed);

	leg->AddEntry(hS1, Form("%4.3f (%s)", func1->GetParameter(2)/func1->GetParameter(1), "no correction"), "PL");
	leg->AddEntry(hS1zcorr, Form("%4.3f (%s)", func2->GetParameter(2)/func2->GetParameter(1), "z-correction"), "PL");
	leg->Draw();

	gPad->Update();
    TPaveStats *st2 = (TPaveStats*)hS1zcorr->GetListOfFunctions()->FindObject("stats");
//    hS1zcorr->GetListOfFunctions()->Print();
    if(st2) {
    	Double_t Y1 = 2*st2->GetY1NDC()-st2->GetY2NDC();
    	st2->SetY2NDC(st2->GetY1NDC());
        st2->SetY1NDC(Y1);
        st2->Draw();
    }
	c1->Update();
//    theApp->Run();
	c1->Print("S1_zcorr_Comp.pdf");
	delete leg;
}

void CompareS2corrections(TFile *file) {
	gStyle->SetOptStat(11);
	TLegend* leg = new TLegend(0.15, 0.7, 0.5, 0.85);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
	leg->SetFillColor(0);
	leg->SetHeader("relative resolution #sigma/#mu");

	vector<TH1D*>hist;
	Int_t nhistos(4);
	TString hname[] = {"hS2", "hS2_zcorr", "hS2_xycorr" ,"hS2_xyzcorr"};
	TString legname[] = {"no correction", "z-correction", "xy-correction" ,"xyz-correction"};
	Double_t y_max(0);
	for (Int_t ih=0; ih<nhistos; ih++){
		hist.push_back((TH1D*) file->Get(hname[ih].Data()));
		if(!hist[ih]) cout<<"hist: "<<hname[ih].Data()<<" is not found."<<endl;
		y_max = TMath::Max(y_max, hist[ih]->GetMaximum());

		hist[ih]->SetMarkerStyle(20+ih);
		hist[ih]->Sumw2();
		(ih==0)? hist[ih]->Draw(): hist[ih]->Draw("sames PLC PMC");
		Double_t min(8e3), max(20e3);
		Int_t bin = hist[ih]->GetMaximumBin();
		min = hist[ih]->GetBinCenter(bin)-4e3;
		max = hist[ih]->GetBinCenter(bin)+4e3;
		hist[ih]->Fit("gaus", "EMRQ", "", min, max);
//		TF1* func = (TF1*)hist[ih]->FindObject("gaus");
		TF1* func = (TF1*)hist[ih]->GetFunction("gaus");
		(ih==0)? func->Draw("same"): func->Draw("same PLC");
		cout<<"relative resolution: "<<func->GetParameter(2)/func->GetParameter(1)<<endl;

		leg->AddEntry(hist[ih], Form("%4.3f (%s)", func->GetParameter(2)/func->GetParameter(1), legname[ih].Data()), "PL");
		if(ih==0) continue;
		gPad->Update();
//	    TPaveStats *st = (TPaveStats*)hist[ih]->GetListOfFunctions()->FindObject("stats");
	    TPaveStats *st = (TPaveStats*)hist[ih]->FindObject("stats");
	    if(st) {
	    	Double_t Y1 = st->GetY1NDC() -(Double_t)ih*(st->GetY2NDC()-st->GetY1NDC());
	    	Double_t Y2 = st->GetY2NDC() -(Double_t)ih*(st->GetY2NDC()-st->GetY1NDC());
	    	st->SetY2NDC(Y2);
	        st->SetY1NDC(Y1);
	        st->Draw();
	    }
	    else cout<<"No Statbox!!!"<<endl;
	}

	hist[0]->SetAxisRange(0,1.5*y_max,"Y");
	hist[0]->SetAxisRange(7e+3, 25e3);
	hist[0]->SetTitle("S2 correction comparison (Am peak); S2 [PE]; Number of Events");
	leg->Draw();

	c1->Update();
//    theApp->Run();
	c1->Print("S2_corr_Comp.pdf");
	delete leg;
}

TF1* FitS1Zdependency(TH2* h2) {
	gStyle->SetOptStat(0);

	TLegend* leg = new TLegend(0.15, 0.65, 0.4, 0.85);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
//	leg->SetFillColor(kWhite);
	leg->SetFillStyle(1001);
	leg->SetFillColorAlpha(kWhite, 0.6);

	Float_t x_fitmin(10), x_fitmax(70);
	TF1 *fpol = new TF1("fpol", "pol3(0)", x_fitmin, x_fitmax);
	fpol->SetParameters(640, 0.6, 0.05, 0.001);
	fpol->SetLineColor(kRed);
	fpol->SetLineWidth(3);

	TH2 *htmp = (TH2*)h2->DrawClone("colz");

    h2->SetAxisRange(400, 800,"Y");
    Int_t lowbin = h2->GetXaxis()->FindBin(x_fitmin);
    Int_t hibin = h2->GetXaxis()->FindBin(x_fitmax);

	h2->FitSlicesY(0, lowbin, hibin, 20);//f1 = 0, firstxbin = 0, lastxbin = -1, cut = 0
	TString hname_fit = Form("%s_1", h2->GetName());
	TH1D *h_mean = (TH1D*) gDirectory->Get(hname_fit.Data());
	if (!h_mean) cout << "hist: " << hname_fit.Data() << " is not found." << endl;
	h_mean->SetMarkerStyle(20);
	h_mean->Fit("fpol", "EMRQ");

	hname_fit = Form("%s_2", h2->GetName());
	TH1D *h_sigma = (TH1D*) gDirectory->Get(hname_fit.Data());
	if (!h_sigma) cout << "hist: " << hname_fit.Data() << " is not found." << endl;
	h_sigma->SetMarkerStyle(21);

	htmp->Draw("colz");
	h_mean->Draw("sames");
	fpol->Draw("same");
	h_sigma->Draw("same");

	gPad->Update();
	c1->Update();
//    TPaveStats *st2 = (TPaveStats*)h_mean->FindObject("stats");
    TPaveStats *st2 = (TPaveStats*)h_mean->GetListOfFunctions()->FindObject("stats");
//    h_mean->GetListOfFunctions()->Print();
    if(st2) {
        st2->SetX1NDC(0.52);st2->SetX2NDC(0.9);
        st2->SetY1NDC(0.6);st2->SetY2NDC(0.9);
        st2->Draw();
    }

	leg->AddEntry(h_mean, "Gauss. fit mean", "LP");
	leg->AddEntry(fpol, "Fit to mean", "L");
	leg->AddEntry(h_sigma, "Gauss. fit sigma", "LP");
	leg->Draw();

	c1->Update();
//         	 theApp->Run();
	c1->Print(Form("%s_Fit.pdf", h2->GetName()));
	gStyle->SetOptStat(1);
	delete leg;
	return fpol;
}


void FitElectronlifetime(TH2* h2) {
	gStyle->SetOptStat(0);
	TLegend* leg = new TLegend(0.15, 0.65, 0.4, 0.85);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
	leg->SetFillStyle(1001);
	leg->SetFillColorAlpha(kWhite, 0.6);

	Float_t x_fitmin(15), x_fitmax(60);
	TF1 *fexpo = new TF1("fexpo", "[0]*exp(-x/[1])", x_fitmin, x_fitmax);
	fexpo->SetParNames("constant", "e_lifetime");
	fexpo->SetParameters(20., 100.);
	fexpo->SetLineColor(kRed);
	fexpo->SetLineWidth(3);

	if(std::strcmp(h2->GetName(),"hS2ovS1_tdrift")!=0)
		h2->SetAxisRange(5e3, 25e3,"Y");
    Int_t lowbin = 0;//h2->GetXaxis()->FindBin(x_fitmin);
    Int_t hibin = -1;//h2->GetXaxis()->FindBin(x_fitmax);

	h2->FitSlicesY(0, lowbin, hibin, 20);//f1 = 0, firstxbin = 0, lastxbin = -1, cut = 0
	TString hname_fit = Form("%s_1", h2->GetName());
	TH1D *h_mean = (TH1D*) gDirectory->Get(hname_fit.Data());
	if (!h_mean)
		cout << "hist: " << hname_fit.Data() << " is not found." << endl;
	h_mean->SetMarkerStyle(20);
	h_mean->Fit("fexpo", "EMRQ+");

	h2->Draw("colz");
	h_mean->Draw("sames");
	fexpo->Draw("same");

	leg->AddEntry(h_mean, "Gauss. fit mean", "LP");
	leg->AddEntry(fexpo, "Fit to mean", "LP");
	leg->Draw();

	c1->Update();
//         	 theApp->Run();
	c1->Print(Form("%s_Fit.pdf", h2->GetName()));
	gStyle->SetOptStat(1);
	delete leg;
	delete fexpo;
}

void DrawS2Fraction(TFile *file) {
	TString hname = "hS2_frac_chnum";
//	TString hname = "hS2_Am_s2maxchan";
	TH2D* h2 = (TH2D*) file->Get(hname.Data());
	if(!h2) cout<<"hist: "<<hname.Data()<<" is not found."<<endl;
//	h2->Rebin2D(1, 2);
	THStack *stack = new THStack(h2, "y", "hstack", "S2 fraction; S2_ch/S2_tot"); // Make stacks from 2D histo. "y" for ProjectionY.

	stack->Draw("nostack PLC C");
	TLegend* leg = c1->BuildLegend(0.5, 0.5, 0.85, 0.85);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
	leg->SetFillColor(0);
	leg->SetNColumns(2);
	TLegendEntry* tmpl;
	Int_t ich(0);
	TIter next(leg->GetListOfPrimitives());
	while ((tmpl=(TLegendEntry*)next())){
	  tmpl->SetLabel(Form("Ch. %d", ich));
	  ich++;
	}

	leg->Draw();
	c1->SetLogy();
	c1->Update();
//         	 theApp->Run();
	c1->Print(Form("%s_Proj.pdf", h2->GetName()));
	c1->SetLogy(0);
	delete stack;
}


void DrawChannelMap(TH3 *h3){
  TCanvas *canvas = new TCanvas("canvas3D", "canvas3D", 800, 600);
  canvas->SetGrid();
  gStyle->SetTitleOffset(1.1, "Y");
  TLatex L;
  L.SetTextSize(0.03);

  THStack *stack = new THStack("ChanMap", "Channel Mapping");
  vector<Double_t> xposi, yposi;

  TAxis *xaxis = h3->GetXaxis();
  Int_t nbins = xaxis->GetNbins();
  Double_t minX = xaxis->GetXmin();
  Double_t maxX = xaxis->GetXmax();
  cout<<"nbins : "<<nbins<<" minX : "<<minX<<" maxX : "<<maxX<<endl;
  Int_t increment = 1;
  for(int i=0;i<nbins;i+=increment) {
        minX = xaxis->GetBinLowEdge(i+1);
        maxX = xaxis->GetBinUpEdge(i+increment);
        xaxis->SetRange(i+1,i+increment);
        TString hname(Form("%s_%4.2f_%4.2f", h3->GetName(), minX, maxX));
        TString htitle(Form("%s = [%4.2f, %4.2f]; bar_x; bar_y", xaxis->GetTitle(), minX, maxX));
        TH2 *tmph = (TH2*)h3->Project3D("zy");
        tmph->SetName(hname);
        tmph->SetTitle(htitle);
        tmph->SetMarkerSize(1);
        tmph->SetMarkerStyle(20);
        stack->Add(tmph);
        xposi.push_back(tmph->GetMean(1));
        yposi.push_back(tmph->GetMean(2));
  }

	stack->Draw("nostack PMC");
	for(int i=0;i<xposi.size();i++) L.DrawLatex(xposi[i], yposi[i], Form("%d",i));
    canvas->Update();
//    theApp->Run();
    canvas->Print(Form("%s.pdf", h3->GetName()));

 delete canvas;
}


void Fit_s1_tdrift_with_s2maxchan_slice(TH3 *h3){
  gSystem->Exec("mkdir Outer_chan; mkdir Inner_chan;");
  TAxis *xaxis = h3->GetXaxis();
  Int_t nbins = xaxis->GetNbins();
  Double_t minX = xaxis->GetXmin();
  Double_t maxX = xaxis->GetXmax();
  cout<<"nbins : "<<nbins<<" minX : "<<minX<<" maxX : "<<maxX<<endl;
  Int_t increment = 1;
//  vector<TF1*> fits;// somehow TF1 doesn't work with PLC option
  vector<TGraph*> gr_fits;
  TH2D *hInnerSum(NULL), *hOuterSum(NULL);
  for(int i=0;i<nbins;i+=increment) {
        minX = xaxis->GetBinLowEdge(i+1);
        maxX = xaxis->GetBinUpEdge(i+increment);
        xaxis->SetRange(i+1,i+increment);
        TString hname(Form("%s_ch%d", h3->GetName(), (int)maxX));
        TString htitle(Form("%s = %d", xaxis->GetTitle(), (int)maxX));
        TH2D *tmph = (TH2D*)h3->Project3D("zy");
        tmph->SetName(hname);
        tmph->SetTitle(htitle);
        TF1* tmpf = FitS1Zdependency(tmph);
        tmpf->SetName(Form("ch_%d", (int)maxX));
        tmpf->SetTitle(Form("Ch. %d", (int)maxX));
//        tmpf->Draw("PLC");
//        fits.push_back(tmpf);
        gr_fits.push_back(new TGraph(tmpf));

        if(IsInnerChannel((int)maxX)) {
        	if(!hInnerSum) hInnerSum = tmph; else hInnerSum->Add(tmph);
            gSystem->Exec(Form("mv hS1_tdrift_s2maxchan_ch%d_Fit.pdf Inner_chan/.", (int)maxX));
        } else {
        	if(!hOuterSum) hOuterSum = tmph; else hOuterSum->Add(tmph);
            gSystem->Exec(Form("mv hS1_tdrift_s2maxchan_ch%d_Fit.pdf Outer_chan/.", (int)maxX));
        }
  }

	TLegend* leg = new TLegend(0.7, 0.3, 0.95, 0.65);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
	leg->SetFillColor(0);


  hInnerSum->SetMarkerColor(kGray);
  hInnerSum->SetAxisRange(400, 800,"Y");
  hInnerSum->SetTitle("Inner Channel");
  hInnerSum->Draw();//Draw("colz");

  leg->SetHeader("Inner Channels");
  for (Int_t ich=0; ich<gr_fits.size(); ich++){
	  if(!IsInnerChannel(ich)) continue;
	  gr_fits[ich]->SetLineStyle(ich%10+1);
	  gr_fits[ich]->Draw("same PLC C");
	  leg->AddEntry(gr_fits[ich], Form("Ch. %d", ich), "L");
  }
  leg->Draw();
  c1->Update();
//    theApp->Run();
  c1->Print("s1_t_drift_Fit_innerChan.pdf");
  leg->Clear();

  hOuterSum->SetMarkerColor(kGray);
  hOuterSum->SetAxisRange(400, 800,"Y");
  hOuterSum->SetTitle("Outer Channel");
  hOuterSum->Draw();//Draw("colz");
  leg->SetHeader("Outer Channels");
  for (Int_t ich=0; ich<gr_fits.size(); ich++){
	  if(IsInnerChannel(ich)) continue;
	  gr_fits[ich]->SetLineStyle(ich%10+1);
	  gr_fits[ich]->Draw("same PLC C");
	  leg->AddEntry(gr_fits[ich], Form("Ch. %d", ich), "L");
  }
  leg->Draw();
  c1->Update();
//    theApp->Run();
  c1->Print("s1_t_drift_Fit_outerChan.pdf");

  delete leg;
}

//                      Ch: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21
const Int_t chan_group[] = {0, 2, 2, 0, 4, 0, 5, 2, 3, 3,  3,  4,  1,  3,  0,  1,  2,  3,  0,  1,  0,  2}; // group number for channels
const Int_t ngroups = 6;

void Fit_s1_tdrift_with_group(TH3 *h3){
  TAxis *xaxis = h3->GetXaxis();
  Int_t nbins = xaxis->GetNbins();
  Double_t maxX = xaxis->GetXmax();
  Int_t increment = 1;
  TString group_names[ngroups];
//  TH2D *hGroups[ngroups];
  TH2D **hGroups = new TH2D*[ngroups];
  for (Int_t i=0; i<ngroups; i++) hGroups[i] = NULL;
  for(int i=0;i<nbins;i+=increment) {
        Int_t ich = (Int_t) xaxis->GetBinUpEdge(i+increment);
        xaxis->SetRange(i+1,i+increment);
        TString hname(Form("%s_ch%d", h3->GetName(), ich));
        TString htitle(Form("%s = %d", xaxis->GetTitle(), ich));
        TH2D *tmph = (TH2D*)h3->Project3D("zy");
        tmph->SetName(hname);
        tmph->SetTitle(htitle);
        if(!hGroups[chan_group[ich]]) hGroups[chan_group[ich]] = tmph; else hGroups[chan_group[ich]]->Add(tmph);
        group_names[chan_group[ich]].Append(Form("_%d", ich));
  }

  vector<TGraph*> gr_group;
  Double_t mean_Am_s1_at_center(0.), t_drift_max(68.);// in PE and us
  for (Int_t i=0; i<ngroups; i++) {
	  hGroups[i]->SetName(Form("%s_ch%s", h3->GetName(), group_names[i].Data()));
	  hGroups[i]->SetTitle(Form("%s%s", xaxis->GetTitle(), group_names[i].Data()));

      TF1* tmpf = FitS1Zdependency(hGroups[i]);
      const Int_t npar(4);
      Double_t par[npar];
      tmpf->GetParameters(par);
      if(i==0) mean_Am_s1_at_center = tmpf->Eval(t_drift_max*0.5); // the middle of t_drift of Group 0 (inner channels) is reference point
      Double_t scale(1.);
      cout<<"Group "<<i<<": ";
      for(Int_t ipar=0; ipar<npar; ipar++) {
    	  par[ipar]/=mean_Am_s1_at_center; // normalize to center of t_drift of Group 0
    	  par[ipar]*=scale;                // scale to 0 to 1 for t_drift
    	  scale*=t_drift_max;           //
          cout<<par[ipar]<<", ";
      }
      cout<<endl;
      TF1 *fpol_norm = new TF1(Form("fpol_norm_%d", i), "pol3(0)", 0, 1);
      fpol_norm->SetParameters(par);
//      fpol_norm->Draw("same PLC");
      gr_group.push_back(new TGraph(fpol_norm));
  }

  TLegend* leg = new TLegend(0.3, 0.2, 0.6, 0.45);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillStyle(1001);
  leg->SetFillColorAlpha(kWhite, 1);

  TH1F *frame = (TH1F*) c1->DrawFrame(0, 0.6, 1, 1.2, "S1 correction factor with group; z position; Reciprocal of S1 correction factor");
  frame->Draw();
  for (Int_t ich=0; ich<gr_group.size(); ich++){
	  gr_group[ich]->SetLineStyle(ich%10+1);
	  gr_group[ich]->Draw("same PLC C");
	  leg->AddEntry(gr_group[ich], Form("ch%s", group_names[ich].Data()), "L");
  }
  leg->Draw();
  c1->Update();
//    theApp->Run();
  c1->Print("s1_t_drift_with_group.pdf");
  delete leg;
}

void Get_S2_MaxS2Chan_Correction(TH2* h2) {
	gStyle->SetOptStat(0);
	TLegend* leg = new TLegend(0.15, 0.65, 0.4, 0.85);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
	leg->SetFillStyle(1001);
	leg->SetFillColorAlpha(kWhite, 0.6);

	h2->SetAxisRange(5e3, 25e3,"Y");

	h2->FitSlicesY();//f1 = 0, firstxbin = 0, lastxbin = -1, cut = 0
	TString hname_fit = Form("%s_1", h2->GetName());
	TH1D *h_mean = (TH1D*) gDirectory->Get(hname_fit.Data());
	if (!h_mean)
		cout << "hist: " << hname_fit.Data() << " is not found." << endl;
	h_mean->SetMarkerStyle(20);

	hname_fit = Form("%s_2", h2->GetName());
	TH1D *h_sigma = (TH1D*) gDirectory->Get(hname_fit.Data());
	if (!h_sigma) cout << "hist: " << hname_fit.Data() << " is not found." << endl;
	h_sigma->SetMarkerStyle(21);
	h_sigma->SetMarkerColor(kRed);
	h_sigma->SetLineColor(kRed);

	h2->SetAxisRange(0, 25e3,"Y");

	h2->Draw("colz");
	h_mean->Draw("sames");
	h_sigma->Draw("same");

	leg->AddEntry(h_mean, "Gauss. fit mean", "LP");
	leg->AddEntry(h_sigma, "Gauss. fit sigma", "LP");
	leg->Draw();

	c1->SetLogz();

	c1->Update();
//         	 theApp->Run();
	c1->Print(Form("%s_Fit.pdf", h2->GetName()));
	gStyle->SetOptStat(1);
	c1->SetLogz(0);

	Double_t ref = h_mean->GetBinContent(1); // get channel 0 as a reference point
	h_mean->Scale(1./ref);
	cout<<"S2 maxS2Chan correction: "<<endl;
	for(Int_t ibin=1; ibin<=h_mean->GetNbinsX(); ibin++){
		cout<<1./h_mean->GetBinContent(ibin)<<", ";
	}
	cout<<endl;
	delete leg;
}


void DrawSlice3D(TH3 *h3){
//  cout <<key->GetName()<<endl;
//  TH3 *h3 = (TH3*)key->ReadObj();
  TCanvas *canvas = new TCanvas("canvas3D", "canvas3D", 800, 600);
//  canvas->SetLogz();
  canvas->SetGrid();
  gStyle->SetTitleOffset(1.1, "Y");

  TAxis *xaxis = h3->GetXaxis();
  Int_t nbins = xaxis->GetNbins();
  Double_t minX,maxX;
        minX = xaxis->GetXmin();
        maxX = xaxis->GetXmax();
        cout<<"nbins : "<<nbins<<" minX : "<<minX<<" maxX : "<<maxX<<endl;
//  Int_t nslice = 10;
  Int_t increment = 1;//10;//nbins/10;
//  for(int i=1;i<=nbins;i++) {
  for(int i=0;i<nbins;i+=increment) {
//              minX = xaxis->GetBinLowEdge(i);
//              maxX = xaxis->GetBinUpEdge(i);
//              xaxis->SetRange(i,i);
        minX = xaxis->GetBinLowEdge(i+1);
        maxX = xaxis->GetBinUpEdge(i+increment);
        xaxis->SetRange(i+1,i+increment);
        TString hname(Form("%s_%4.2f_%4.2f", h3->GetName(), minX, maxX));
        TString htitle(Form("%s = [%4.2f, %4.2f]", xaxis->GetTitle(), minX, maxX));
//        cout<<"hname : "<<hname.Data()<<endl;
        TH2 *tmph = (TH2*)h3->Project3D("zy");
//      TProfile2D *tmph = h3->Project3DProfile(hname);
//      hname.Prepend(histname);
        tmph->SetName(hname);
        tmph->SetTitle(htitle);
        tmph->SetAxisRange(0, 5,"Z");
        tmph->Draw("colz");
//      bin+=incre;
        canvas->Update();
#if 0
        theApp->Run();
        canvas->Print(Form("%s.pdf", hname.Data()));
#else
//        canvas->Update();
        canvas->SaveAs(Form("%s.gif+10", h3->GetName()));
#endif
  }
 delete canvas;
}


#ifndef __CINT__
int main(int argc, char **argv) {
	theApp = new TRint("App", &argc, argv, NULL, 0);
    theApp->Connect("KeyPressed(Int_t)","TSystem",gSystem,"ExitLoop()");
	if ( theApp->Argc() == 2 ) {
		std::cout << "\n==========> PlotHistograms <=============" << std::endl;
		std::cout << "==> Application start." << std::endl;
		PlotHistograms(theApp->Argv(1));

	} else if ( theApp->Argc() == 3 ) {
		std::cout << "\n==========> PlotHistograms <=============" << std::endl;
		PlotHistograms(theApp->Argv(1), theApp->Argv(2));
	} else {
		std::cout << "Usage:" <<theApp->Argc()<< std::endl;
		std::cout << "./PlotQA" << std::endl;
		return 0;
	}


	std::cout << "==> Application finished." << std::endl;
	return 0;
}
#endif /* __CINT __ */
