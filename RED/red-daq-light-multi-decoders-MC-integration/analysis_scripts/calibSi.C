#include "Utils.C"
#include "Math/SpecFunc.h"
#include "Math/Functor.h"
#include "Fit/Fitter.h"
#include "TF2.h"

#define UNDEF -100


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


double GetEdep(double E0) // Edep in the 20 um Silicon detector from Lise++ (corrected so adc/keV agrees with alpha )
{
  double A=3.00, alpha=-0.77;
  return A*pow(E0/20.,alpha);  
}


void PlotEdep()
{
  SetStyle();
  //
  const double t=20.; //um
  const double mLi=7.0160034366; // u
  double E[5]   ={ 2.250*mLi, 2.50*mLi, 3.00*mLi, 4.03*mLi, 5.4*mLi  };// MeV/u
  double Edep[5]={ 0.254*t  , 0.228*t , 0.200*t,  0.161*t,  0.1287*t };// MeV/um
  
  TGraph tg(5,E,Edep);
  tg.GetXaxis()->SetTitle("E_{Li} [MeV]");
  tg.GetYaxis()->SetTitle("E_{dep} [MeV]");

  tg.DrawClone("AP");
  gPad->SetLogx();
  gPad->SetLogy();
  
  TF1 *f=new TF1("f","[0]*pow(x/20.,[1])",1.,50.);
  tg.Fit(f,"QN","");
  f->DrawCopy("SAME");
  l.DrawLatex(0.2,0.2,Form("%5.2lf #mu m",t));
  l.DrawLatex(0.2,0.3,Form("%5.2lf %5.2lf",f->GetParameter(0),f->GetParameter(1) ));
    
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


double BiGaus(double *v,double *par)
{
  double x=v[0];
  double y=v[1];

  double meanX=par[0];
  double sigmaX=par[1];

  double meanY=par[2];
  double sigmaY=par[3];
  
  double rho=par[4];
  double norm=par[5];

  double z=pow(x-meanX,2.)/sigmaX/sigmaX + pow(y-meanY,2.)/sigmaY/sigmaY + (x-meanX)*(y-meanY)*2.*rho/sigmaX/sigmaY;

  double val=1./(2*TMath::Pi()*sigmaX*sigmaY*sqrt(1.-rho*rho))*exp(-z/2/(1-rho*rho));
  return val*norm;
  
}


double AsymGaus(double *v,double *par)
{
  double x=v[0];

  double mean=par[0];
  double sigma=par[1];
  double norm=par[2];
  double alpha=par[3];
  double cte=par[4];

  double x0=(x-mean)/sigma;
  double val=norm*TMath::Gaus(x0,0,1.)*(TMath::Erf(alpha*x0)+1.) +cte;
  return val;
  
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void FitCalib(vector<Double_t> MeanE, vector<Double_t> eMeanE, vector<Double_t> MeanDE, vector<Double_t> eMeanDE, vector<Double_t> E0 )
{

  auto chi2Function = [&](const Double_t *par)
    {
      //minimisation function computing the sum of squares of residuals
      // looping at the data points
      Double_t f=0;
      for (size_t i=0;i<MeanE.size();i++) {
	Double_t expectedEnergy = E0[i];
       
	Double_t p1 = (par[0] + par[1]*MeanDE[i]+par[2]*MeanE[i] -
		       expectedEnergy);
	Double_t errorsquare = (  (eMeanDE[i]/MeanDE[i])*(eMeanDE[i]/MeanDE[i]) +
				  (eMeanE[i]/MeanE[i])*(eMeanE[i]/MeanE[i]) )  *expectedEnergy*expectedEnergy;
	f += p1*p1/errorsquare; //minimization!
      }   
      return f;
    };
  
  ROOT::Math::Functor fcn(chi2Function,3);
  ROOT::Fit::Fitter  fitter;
  double pStart[3] = {0,1,1};
  fitter.SetFCN(fcn, pStart);
  fitter.Config().ParSettings(0).SetName("a");
  fitter.Config().ParSettings(1).SetName("b");
  fitter.Config().ParSettings(2).SetName("c");
  fitter.Config().ParSettings(0).Fix(); // Parameter Ill defined due to energy sampling

  bool ok = fitter.FitFCN();
  if (!ok) {
    Error("line3Dfit","Line3D Fit failed");
  }   
  const ROOT::Fit::FitResult & result = fitter.Result();
  result.Print(std::cout);
  
  Double_t a = result.Parameter(0);
  Double_t b = result.Parameter(1);
  Double_t c = result.Parameter(2);

  for (size_t i=0;i<MeanE.size();i++) {
    Double_t expectedEnergy = E0[i];
    Double_t p1 = a+b*MeanDE[i]+c*MeanE[i];
    printf("calibrated energy %5.2lf  (%5.2lf) \n",p1,expectedEnergy);
  }



}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Derive energy calibration using Li-Au/Li-C elastic scattering at different energies (weak dependence on thetaSi ) 
//                                Collimator 0
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

void CalibGold() 
{
  SetStyle();

  const int nF_in=6;
  const double E0[nF_in]={ 18.,26.,28.,17.92, 25.88,27.88 }; // Energy of the Scattered Li7 at theta=5 deg
  std::string NameSet[nF_in]={ "Li+Au 18 [MeV]", "Li+Au 26 [MeV]", "Li+Au 28 [MeV]", "Li+C 18 [MeV]", "Li+C 26 [MeV]", "Li+C 28 [MeV]" };
  const int nPeaks=nF_in;
  const int Collimator=0;

  int runRef[nPeaks]={635,628,621 /*Gold*/, 632,630,624};  //CH2    
  const double  Erange[nPeaks][2]={ {5600., 5900.}, {9050., 9350.}, { 9885, 10185.} /*Gold*/,    {5425.,5725.}, {8900.,9200.}, {9750.,10050.} /*CH2*/};
  const double dErange[nPeaks][2]={ {3800., 4800.}, {2600., 3600.}, {2500.,  3500.} /*Gold*/,    {3750.,4750.}, {2700.,3700.}, {2500.,3500.} /*CH2*/};
  
  //----------------------------------------------------------------------------------
  TH2F   *h2D[nF_in];
  TH1F   *h1D[nF_in];
  TGraph *t2D[nF_in];

  int color[6]={kBlue, kRed, kGreen+2, kBlue, kRed, kGreen+2};
  
  nF=nF_in;
  for (int iF=0;iF<nF;iF++)
    {
      iF_plot.push_back(iF);
      frun_i[iF]=runRef[iF]; ftime_i[iF]=10.;  IsTPC=false;

      TTree *t;
      std::string FileName=Form("reco/run_%d.root",runRef[iF]);
      TFile *f=new TFile(FileName.c_str());
      f->GetObject("reco",t);
      Summary(t,iF);

      //----------------
      const char *hname=Form("h2D%d",iF);
      int nbins=(iF/3==0?100:50);
      h2D[iF]=new TH2F(hname,hname,nbins,Erange[iF][0],Erange[iF][1],nbins,dErange[iF][0],dErange[iF][1] );
      h2D[iF]->Sumw2();
      
      std::string x=  SiESlow;
      std::string y=  SidESlow;
      
      std::string Sel=Form("%s > %5.2lf && %s < %5.2lf && %s>%5.2lf && %s <%5.2lf ", 
			   SiESlow.c_str() , Erange[iF][0] , SiESlow.c_str(), Erange[iF][1],
			   SidESlow.c_str(), dErange[iF][0], SidESlow.c_str(), dErange[iF][1] );

      int nn=t->Project(hname,Form("%s:%s", y.c_str(),x.c_str()), Sel.c_str() );

      double *v2=t->GetV2();
      double *v1=t->GetV1(); 
      t2D[iF]=new TGraph(nn,v2,v1);

      h2D[iF]->SetDirectory(0);	    
      h2D[iF]->Scale(100./h2D[iF]->Integral());

      //-------------------
      const char *hname1D=Form("h1D_%d",iF);
      h1D[iF]=new TH1F(hname1D,hname1D,100,-1.,1.);//0.98,1.02);
      h1D[iF]->Sumw2();

      x=Form("(%s*%e+%s*%e)-%e",SidESlow.c_str(), B_Conv[Collimator][0],  SiESlow.c_str(), B_Conv[Collimator][1],E0[iF] );
      nn=t->Project(hname1D,x.c_str(), Sel.c_str() );
      cout << x << endl;
      h1D[iF]->SetDirectory(0);	    
      h1D[iF]->Scale(1./h1D[iF]->Integral());
      

      f->Close();
    }
  

  //----------------------------------------------------------------------------------

  TGraphErrors tg[2][4]; //[Gold/CH2][mean E/rms E/mean dE/rms dE]
  vector<Double_t> MeanE,  eMeanE,  MeanDE,  eMeanDE,  fE0, SigmaE;
	 
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,2.*wh);
  cc[canvas]->Divide(3,2);
  cc[canvas]->Draw();  

  for (int iF=0;iF<nF;iF++)
    {
      cc[canvas]->cd(iF+1);
      
      TH2F *h0=h2D[iF];
      h0->GetXaxis()->SetTitle("E [adc]");
      h0->GetYaxis()->SetTitle("#DeltaE [adc]");
      h0->DrawCopy("COLZ");
      gPad->SetLogz();

      const Int_t npar = 6;
      Double_t par[npar] ={ (Erange[iF][0]+Erange[iF][1])/2. ,  30., (dErange[iF][0]+dErange[iF][1])/2., 100., 0.7, 100.};
      TF2 *f2 = new TF2("f2",BiGaus,Erange[iF][0],Erange[iF][1],dErange[iF][0],dErange[iF][1] , npar);
      f2->SetParameters(par);
      h0->Fit("f2","Q0");

      f2->SetLineColor(kBlack);
      double cont[2]={0.01,0.1};
      f2->SetContour(2,cont);
      f2->Draw("cont3 same");
      
      MeanE.push_back(f2->GetParameter(0));   eMeanE.push_back(10.);//f2->GetParError(0));
      MeanDE.push_back(f2->GetParameter(2));  eMeanDE.push_back(10.);//f2->GetParError(2));
      SigmaE.push_back(f2->GetParameter(1)/f2->GetParameter(0));
      fE0.push_back(E0[iF]);
      printf(" E0=%5.2lf E=%5.2lf %5.2lf DE=%5.2lf %5.2lf  ", fE0[iF], MeanE[iF],eMeanE[iF], MeanDE[iF],eMeanDE[iF] );

      l.DrawLatex(0.2,0.25,NameSet[iF].c_str());
      int iTarget=iF/3;
      for (int ipar=0;ipar<npar;ipar++) 
	{
	  if ( ipar>3 ) continue;
	  TGraphErrors *t0=&tg[iTarget][ipar];
	  t0->SetPoint(iF%3,E0[iF],(ipar%2==0?f2->GetParameter(ipar): f2->GetParameter(ipar)/f2->GetParameter(ipar-1)) );
	  t0->SetPointError(iF%3,0.,(ipar%2==0?f2->GetParError(ipar) : f2->GetParError(ipar)/f2->GetParameter(ipar-1)) );
	}
      printf("\n");
    }
  canvas++;
  
  //-------------------------------------------------------------

  std::string name[4]={ "<E> [adc]", "#sigma_{E}/<E> ", "<#Delta E> [adc]", "#sigma_{#Delta E}/<#Delta E> "};
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",1.5*ww,1.5*wh);
  cc[canvas]->Divide(2,2);
  cc[canvas]->Draw();  

    for (int ipar=0;ipar<4;ipar++) 
    {
      double x0=0.65,y0=0.25;
      TLegend *leg=new TLegend(x0,y0,x0+0.2,y0+0.15);  

      cc[canvas]->cd(ipar+1);
      for (int iTarget=0;iTarget<2;iTarget++)
	{
	  int color_i=color[iTarget];
	  TGraphErrors *t0=&tg[iTarget][ipar];
	  t0->SetLineColor(color_i);
	  t0->SetMarkerColor(color_i);
	  t0->GetYaxis()->SetTitle(name[ipar].c_str());
	  t0->GetXaxis()->SetTitle("E_{^{7}Li} [MeV]");
	  if ( ipar==1 ) { t0->SetMaximum(0.01); t0->SetMinimum(0.002); }
	  t0->DrawClone((iTarget==0?"ALP":"LP"));      
	  leg->AddEntry(t0,(iTarget==0?"Li+Au":"Li+C"),"P");
	}
      if ( ipar==0 ) leg->DrawClone();
    }
  canvas++;
  
  FitCalib(MeanE,eMeanE,MeanDE,eMeanDE,fE0);


  //-------------------------------------------------------------
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,1.5*wh);
  cc[canvas]->Draw();  

  TH2D *hFrame=new TH2D("hXS","",1000, -1., 1.,1000,0.,0.2);
  hFrame->GetXaxis()->SetTitle("E_{rec}-E_{true} [MeV]");
  hFrame->GetYaxis()->SetTitle("A.U");
  hFrame->DrawCopy();

  double x0=0.65,y0=0.55;
  TLegend *leg=new TLegend(x0,y0,x0+0.2,y0+0.35);  
    
  for (int iF=0;iF<nF_in;iF++) 
    {
      TH1F *h0=h1D[iF];
      h0->SetLineColor(color[iF]);
      h0->SetLineStyle(iF/3==0?1:2);
      h0->SetMarkerStyle(iF/3==0?20:24);
      h0->SetMarkerColor(color[iF]);
      h0->DrawCopy("E1SAME");
      leg->AddEntry(h0,NameSet[iF].c_str(),"P");
      
      h0->Fit("gaus","Q0");
      h0->GetFunction("gaus")->SetLineColor(color[iF]);
      h0->GetFunction("gaus")->SetLineStyle(iF/3==0?1:2);
      h0->GetFunction("gaus")->DrawCopy("SAME");

      double sigma=h0->GetFunction("gaus")->GetParameter(2);
      printf("E0=%5.2f | Sigma %e  (comb)  %e (E) \n",E0[iF],sigma,SigmaE[iF]);
    }
  leg->DrawClone();
  canvas++;
  
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Cross calibration of Collimator 1 using Collimator 0  (E=28 MeV)  , Spatial dependence of SiTel response
//          Using Li-Au/Li-C elastic scattering at 28 MeV
//
//  IMPORTANT: the Be/Li bands are shifted for Coll 0/1. A change in the beam energy cannot cause that. Therefore,
//            it is assumed that is due to a different detector gain. To cross calibrate and due to lack of different energies,
//             we take two runs and assume that the beam energy is the same in both cases.
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/*
 E=28 MeV   --> Derive Relative Calibration using LiC LiAu peaks  (Different E/DE gains due to spatial uniformity, different hole position )
  624 / 676   Coll 0/1  (Low thresh CH2)
  621 / 691   Coll 0/1  ( Gold )
*/
void CrossCalibGold() 
{
  SetStyle();
  const int nF_in=4;
  const double E0[nF_in]={ 28., 28., 27.88, 27.88  }; // Energy of the Scattered Li7 at theta=5 deg
  std::string NameSet[nF_in]={ "Li+Au 28 [MeV] C0", "Li+Au 28 [MeV] C1",  "Li+C 28 [MeV] C0",  "Li+C 28 [MeV] C1" };
  const int nPeaks=nF_in;


  int runRef[nPeaks]={ 621, 691 /*Gold*/, 624,676 };  //CH2    
  const double  Erange[nPeaks][2]={  { 9885, 10185.},{ 9600, 9800.} /*Gold*/,     {9750.,10050.},{9400.,9800.} /*CH2*/};
  const double dErange[nPeaks][2]={  {2500.,  3500.},{2500.,  3500.} /*Gold*/,     {2500., 3500.},{2500., 3500.} /*CH2*/};
  
  //----------------------------------------------------------------------------------
  TH2F   *h2D[nF_in];
  int color[6]={kBlue, kRed, kGreen+2, kBlue, kRed, kGreen+2};
  
  nF=nF_in;

  for (int iF=0;iF<nF;iF++)
    {
      iF_plot.push_back(iF);
      frun_i[iF]=runRef[iF]; ftime_i[iF]=10.;  IsTPC=false;

      TTree *t;
      std::string FileName=Form("reco/run_%d.root",runRef[iF]);
      TFile *f=new TFile(FileName.c_str());
      f->GetObject("reco",t);
      Summary(t,iF);

      //----------------
      const char *hname=Form("h2D%d",iF);
      int nbins=50;
      h2D[iF]=new TH2F(hname,hname,nbins,Erange[iF][0],Erange[iF][1],nbins,dErange[iF][0],dErange[iF][1] );
      h2D[iF]->Sumw2();
      
      std::string x=  SiESlow;
      std::string y=  SidESlow;
      
      std::string Sel=Form("%s > %5.2lf && %s < %5.2lf && %s>%5.2lf && %s <%5.2lf ", 
			   SiESlow.c_str() , Erange[iF][0] , SiESlow.c_str(), Erange[iF][1],
			   SidESlow.c_str(), dErange[iF][0], SidESlow.c_str(), dErange[iF][1] );

      int nn=t->Project(hname,Form("%s:%s", y.c_str(),x.c_str()), Sel.c_str() );

      h2D[iF]->SetDirectory(0);	    
      h2D[iF]->Scale(100./h2D[iF]->Integral());
      

      f->Close();
    }
  

  //----------------------------------------------------------------------------------

  vector<Double_t> MeanE,  eMeanE,  MeanDE,  eMeanDE,  fE0, SigmaE;
	 
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,2.*wh);
  cc[canvas]->Divide(2,2);
  cc[canvas]->Draw();  

  for (int iF=0;iF<nF;iF++)
    {
      cc[canvas]->cd(iF+1);
      
      TH2F *h0=h2D[iF];
      h0->GetXaxis()->SetTitle("E [adc]");
      h0->GetYaxis()->SetTitle("#DeltaE [adc]");
      h0->DrawCopy("COLZ");
      gPad->SetLogz();

      const Int_t npar = 6;
      Double_t par[npar] ={ (Erange[iF][0]+Erange[iF][1])/2. ,  30., (dErange[iF][0]+dErange[iF][1])/2., 100., 0.7, 100.};
      TF2 *f2 = new TF2("f2",BiGaus,Erange[iF][0],Erange[iF][1],dErange[iF][0],dErange[iF][1] , npar);
      f2->SetParameters(par);
      h0->Fit("f2","Q0");

      f2->SetLineColor(kBlack);
      double cont[2]={0.01,0.1};
      f2->SetContour(2,cont);
      if ( iF!=1 ) f2->Draw("cont3 same");


      MeanE.push_back(f2->GetParameter(0));   eMeanE.push_back(10.);//f2->GetParError(0));
      MeanDE.push_back(f2->GetParameter(2));  eMeanDE.push_back(10.);//f2->GetParError(2));
      SigmaE.push_back(f2->GetParameter(1)/f2->GetParameter(0));
      fE0.push_back(E0[iF]);
      printf(" E0=%5.2lf E=%5.2lf %5.2lf DE=%5.2lf %5.2lf \n ", fE0[iF], MeanE[iF],eMeanE[iF], MeanDE[iF],eMeanDE[iF] );
      l.DrawLatex(0.2,0.25,NameSet[iF].c_str());
    }
  canvas++;


}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Derive Li/Be band cuts and check Theta_Be
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

void CalibCH2(  bool FitBand=false ) 
{
  SetStyle();
  const int nF_in=3;
  int runRef[nF_in]={632,630,624};  //CH2 
  const double E0[nF_in]={18.,26.,28.};
  double thetaSi[nF_in]={ 5.5, 5.5, 5.5 };

  TH1F   *h1D[2][nF_in];
  TGraph *t2D[nF_in];
  TH2F   *h2D[nF_in];
  TProfile *hP[2][nF_in];

  int color[4]={kBlue, kRed, kGreen+2,kBlack};
  
  nF=nF_in;

  for (int iF=0;iF<nF;iF++)
    {
      iF_plot.push_back(iF);
      frun_i[iF]=runRef[iF]; ftime_i[iF]=10.;  IsTPC=false;

      TTree *t;
      std::string FileName=Form("reco/run_%d.root",runRef[iF]);
      TFile *f=new TFile(FileName.c_str());
      f->GetObject("reco",t);
      Summary(t,iF);

      //----------------
      const char *hname=Form("h2D%d",iF);
      h2D[iF]=new TH2F(hname,hname,300,5.,28.,300,1.5,7.5);
      h2D[iF]->Sumw2();
      
      std::string Sel="1";

      int nn=t->Project(hname,Form("%s:%s", SidE_dep.c_str(),SiE_dep.c_str()), Sel.c_str() );

      double *v2=t->GetV2();
      double *v1=t->GetV1(); 
      t2D[iF]=new TGraph(nn,v2,v1);

      h2D[iF]->SetDirectory(0);	    
      h2D[iF]->Scale(100./h2D[iF]->Integral());

      //----------------------
      for (int iBand=0;iBand<2;iBand++)
	{
	  //------------------
	  const char *hnameP=Form("hprof%d%d",iF,iBand);
	  hP[iBand][iF]=new TProfile(hnameP,hnameP,100,(iBand==0?5.:10.),(iBand==0?30.:20.) );
	  hP[iBand][iF]->Sumw2();
	  hP[iBand][iF]->SetDirectory(0);	    
	  for (int ii=0;ii<nn;ii++)
	    {
	      double E=v2[ii];
	      double dE=v1[ii];
	      double dE_p=5.2-4.34*log10(E/10.);
	      if ( dE>7.) continue;
	      if ( dE>dE_p && iBand==0 ) continue;
	      if ( dE<dE_p && iBand==1 ) continue;
	      hP[iBand][iF]->Fill(E,dE);
	    }
	}

      for (int iBand=0;iBand<2;iBand++)
	{
	  //-------------------
	  const char *hname1D=Form("h1D_%d%d",iF,iBand);
	  h1D[iBand][iF]=new TH1F(hname1D,hname1D,(iBand==0?300:100.),15.,(iBand==0?30.:25.));
	  h1D[iBand][iF]->Sumw2();
	  
	  nn=t->Project(hname1D,SiEnergy.c_str(), (iBand==0?IsLiBand.c_str():IsBeBand.c_str()) );
	  h1D[iBand][iF]->SetDirectory(0);	    
	  h1D[iBand][iF]->Scale(1./h1D[iBand][iF]->Integral());
	}
      f->Close();
      

    }
  

  //----------------------------------------------------------------------------------
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,1.*wh);
  cc[canvas]->Divide(3,1);
  cc[canvas]->Draw();  


  for (int iF=0;iF<nF;iF++)
    {
      cc[canvas]->cd(iF+1);
      
      TH2F *h0=h2D[iF];
      h0->GetXaxis()->SetTitle("E [MeV]");
      h0->GetYaxis()->SetTitle("#DeltaE [MeV]");
      h0->DrawCopy("COLZ");
      gPad->SetLogz();
      
      l.DrawLatex(0.5,0.85,Form("E_{beam}=%5.2lf",E0[iF]));

      for (int iBand=0;iBand<2;iBand++)
	{
	  TF1 *f=new TF1("f","([0]+[1]*(log10(x/10.)))*[2]",0.,30.);
	  f->SetLineWidth(2);
	  if ( FitBand )
	    { 
	      TProfile *h1=hP[iBand][iF];
	      h1->DrawCopy("SAME");
	      f->FixParameter(2,1.); f->SetParameter(0,4.);
	      h1->Fit(f,"QN","",0.,30.);
	      f->DrawCopy("SAME");
	      printf("Band=%d  iF=%d  %5.2lf %5.2lf  %5.2lf  N=%f \n",iBand,iF,f->GetParameter(0),f->GetParameter(1),f->GetParameter(2),h1->GetEntries());
	    }
	  else
	    {
	      for (int ii=0;ii<3;ii++)
		{
		  double par[3]={ (iBand==0?pLiBand[0]:pBeBand[0]), (iBand==0?pLiBand[1]:pBeBand[1]), (ii==0?1.:(ii==1?0.9:1.1)) };
		  f->SetParameters(par);
		  f->SetLineStyle(ii==0?1:2);
		  f->DrawCopy("SAME");
		}
	      
	    }
	}
    }
  canvas++;


  //-------------------------------------------------------------
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,1.5*wh);
  cc[canvas]->Divide(2,1);
  cc[canvas]->Draw();  


  double RangeLowBe [nF_in][2] ={ { UNDEF, UNDEF }, {17.5, 19.0}, {18.5,20.5} };
  double RangeHighBe[nF_in][2] ={ { UNDEF, UNDEF }, {20.5, 22.0}, {22.0,24.0} };

  double RangeLowLiP[nF_in][2]   ={ { UNDEF, UNDEF }, {15.0, 16.0}, {16.5,17.5} };
  double RangeHighLiP[nF_in][2]   ={ { 16.0, 17.0  }, {23.8, 24.8}, {25.2,26.5} };
  double RangeLiC[nF_in][2]   ={ { 17.4, 18.4 }, {25.0, 26.0}, {27.0,28.0} };
  
  TGraphErrors tTheta[nF_in][3];

  for (int iBand=0;iBand<2;iBand++)
    {
      cc[canvas]->cd(iBand+1);

      TH2D *hFrame=new TH2D(Form("hX%d",iBand),"",1000, 15., (iBand==0?30.:25.),1000,0.,(iBand==0?0.2:0.1));
      hFrame->GetXaxis()->SetTitle("E_{rec} [MeV]");
      hFrame->GetYaxis()->SetTitle("A.U");
      hFrame->DrawCopy();
      
      for (int iF=0;iF<nF_in;iF++) 
	{
	  if ( iBand==1 && iF==0 ) continue;
	  TH1F *h0=h1D[iBand][iF];
	  h0->SetLineColor(color[iF]);
	  h0->SetMarkerColor(color[iF]);
	  h0->DrawCopy("HISTSAME");
	  l.SetTextColor(color[iF]);
	  l.SetTextSize(0.03);
	  l.DrawLatex(0.4,0.85-0.05*double(iF),  Form("E_{beam}=%3.0lf [MeV]",E0[iF]) );
	  l.SetTextSize(0.035);
	  l.SetTextColor(kBlack);


	  for (int ipeak=0;ipeak<3;ipeak++)
	    {
	      double *range=NULL;

	      if ( iBand==1 && ipeak==0 ) continue;
	      if ( iBand==1 && ipeak==1 ) range=RangeLowBe[iF];
	      if ( iBand==1 && ipeak==2 ) range=RangeHighBe[iF];

	      if ( iBand==0 && ipeak==0 && iF==0 ) continue;
	      if ( iBand==0 && ipeak==0 ) range=RangeLowLiP[iF];
	      if ( iBand==0 && ipeak==1 ) range=RangeHighLiP[iF];
	      if ( iBand==0 && ipeak==2 ) range=RangeLiC[iF];

	      TF1 *f=new TF1("f",AsymGaus,range[0],range[1],5);
	      f->SetNpx(1000);
	      double par[5]={ (range[0]+range[1])/2., 0.5, 0.01, 0., 0.};
	      f->SetParameters(par);
	      h0->Fit(f,"QN","",range[0],range[1]);
	      f->SetLineColor(color[iF]);
	      f->SetLineWidth(1);
	      f->DrawCopy("SAME");

	      double mean=f->GetMaximumX(), err=0.2; //In the asymmetric gaussian the mean is not the peak position.
	      printf("%d %d peak=%5.2lf (p0=%5.2lf) alpha=%5.2lf\n",iBand,ipeak, mean,f->GetParameter(0),f->GetParameter(3));	      

	      int iSet;
	      if (iBand==0 && ipeak==2 ) iSet=0; // LiC
	      else if (iBand==0        ) iSet=1; // LiP
	      if (iBand==1 )             iSet=2; // Be
	      
	      TGraphErrors *t0=&tTheta[iF][iSet];
	      int ipoint=t0->GetN();
	      t0->SetPoint(ipoint,thetaSi[iF],mean);
	      t0->SetPointError(ipoint,0.,0.1);
	    }
	  l.DrawLatex(0.2,0.85,(iBand==0?"Li Band":"Be Band"));
	}
    }
  canvas++;

  //-------------------------------------------------------------

  int colorC[nChannels-1]={kGreen, kOrange, kRed, kBlue  };
  
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,1.*wh);
  cc[canvas]->Divide(3,1);
  cc[canvas]->Draw();  

  for (int iF=0;iF<nF_in;iF++) 
    {
      cc[canvas]->cd(iF+1);
      double beamEne=E0[iF];
      for (int ichan=0;ichan<nChannels-1;ichan++)
	{
	  for (int SelForward=1;SelForward>=0;SelForward--)
	    {
	      TGraph *tg=new TGraph();
	      GetEnergy(beamEne,ichan,SelForward,tg,10.);
	      if ( tg->GetN()==0 ) continue;
	      //printf("%s  | Forward=%d | theta=%5.2lf ke=%5.2lf \n", ChannelName[ichan].c_str(),SelForward,5.,tg->Eval(5.));
	      
	      tg->SetMaximum(30.);
	      tg->SetMinimum(10.);
	      tg->SetLineColor(colorC[ichan]);
	      tg->SetLineStyle(SelForward==1?20:24);
	      tg->GetXaxis()->SetTitle("#theta_{Si Detector}");
	      tg->GetYaxis()->SetTitle("KE [MeV]");
	      tg->DrawClone((ichan==0&&SelForward==1?"AL":"L"));
	      if ( SelForward==1 )
		{
		  l.SetTextColor(colorC[ichan]);
		  l.DrawLatex(0.6,(iF==0?0.85:0.35)-double(ichan)*0.05,ChannelName[ichan].c_str());
		  l.SetTextColor(kBlack);
		}
	    }
	  if ( ichan!=1 )
	    {
	      int Set=UNDEF;
	      if ( ichan==0     ) Set=0;
	      else if (ichan==2 ) Set=2;
	      else if (ichan==3 ) Set=1;	     	      
	      TGraphErrors *t0=&tTheta[iF][Set];
	      t0->SetMarkerColor(colorC[ichan]);
	      t0->DrawClone("P");
	    }
	}
      l.DrawLatex(0.2,0.85,Form("E_{beam}=%5.2lf",beamEne));
    }
  canvas++;




}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Compare collimators
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

void CompareCollimators() 
{
  SetStyle();
  const int nF_in=2;
  int runRef[nF_in]={624,676};  //CH2 
  double time_i[nF_in]={581.,110.};
  const double E0[nF_in]={28.,28.};
  std::string NameSet[2]={ "Collimator 0  28 MeV","Collimator 0 Inv  28 MeV"};
  double thetaSi[nF_in]={ 5.5, 5.0 };


  TH1F   *h1D[2][nF_in],*h1D_SiE[2][nF_in];
  TGraph *t2D[nF_in];
  TH2F   *h2D[nF_in];
  TProfile *hP[2][nF_in];

  int color[4]={kBlue, kRed, kGreen+2,kBlack};
  
  nF=nF_in;

  for (int iF=0;iF<nF;iF++)
    {
      iF_plot.push_back(iF);
      frun_i[iF]=runRef[iF]; ftime_i[iF]=time_i[iF];  IsTPC=false;

      TTree *t;
      std::string FileName=Form("reco/run_%d.root",runRef[iF]);
      TFile *f=new TFile(FileName.c_str());
      f->GetObject("reco",t);
      Summary(t,iF);

      //----------------
      const char *hname=Form("h2D%d",iF);
      h2D[iF]=new TH2F(hname,hname,300,5.,28.,300,1.5,7.5);
      h2D[iF]->Sumw2();
      
      std::string Sel="1";

      int nn=t->Project(hname,Form("%s:%s", SidE_dep.c_str(),SiE_dep.c_str()), Sel.c_str() );

      double *v2=t->GetV2();
      double *v1=t->GetV1(); 
      t2D[iF]=new TGraph(nn,v2,v1);

      h2D[iF]->SetDirectory(0);	    
      h2D[iF]->Scale(100./h2D[iF]->Integral());

      //----------------------
      for (int iBand=0;iBand<2;iBand++)
	{
	  //------------------
	  const char *hnameP=Form("hprof%d%d",iF,iBand);
	  hP[iBand][iF]=new TProfile(hnameP,hnameP,100,(iBand==0?5.:13.),(iBand==0?30.:20.) );
	  hP[iBand][iF]->Sumw2();
	  hP[iBand][iF]->SetDirectory(0);	    
	  for (int ii=0;ii<nn;ii++)
	    {
	      double E=v2[ii];
	      double dE=v1[ii];
	      double dE_p;
	      if ( iBand==0 ) dE_p=pLiBand[0]+pLiBand[1]*log10(E/10.);
	      else            dE_p=pBeBand[0]+pBeBand[1]*log10(E/10.);
	      if ( dE/dE_p<0.85 || dE/dE_p>1.15) continue;
	      hP[iBand][iF]->Fill(E,dE);
	    }
	
	  //-------------------
	  const char *hname1D=Form("h1D_%d%d",iF,iBand);
	  h1D[iBand][iF]=new TH1F(hname1D,hname1D,(iBand==0?300:100.),15.,(iBand==0?30.:25.));
	  h1D[iBand][iF]->Sumw2();
	  
	  nn=t->Project(hname1D,SiEnergy.c_str(), (iBand==0?IsLiBand.c_str():IsBeBand.c_str()) );
	  h1D[iBand][iF]->SetDirectory(0);	    
	  h1D[iBand][iF]->Scale(1./h1D[iBand][iF]->Integral());

	  //-------------------
	  const char *hname1D_SiE=Form("h1D_SiE_%d%d",iF,iBand);
	  h1D_SiE[iBand][iF]=new TH1F(hname1D_SiE,hname1D_SiE,(iBand==0?100:50),SiE_Range[0],SiE_Range[1]);
	  h1D_SiE[iBand][iF]->Sumw2();
	  
	  nn=t->Project(hname1D_SiE,SiE_dep.c_str(), (iBand==0?IsLiBand.c_str():IsBeBand.c_str()) );
	  h1D_SiE[iBand][iF]->SetDirectory(0);	    
	  h1D_SiE[iBand][iF]->Scale(1./h1D_SiE[iBand][iF]->Integral());

	}
      f->Close();
      
    }


  //----------------------------------------------------------------------------------
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,1.*wh);
  cc[canvas]->Divide(2,1);
  cc[canvas]->Draw();  

  bool Fit=false;

  for (int iF=0;iF<nF;iF++)
    {
      cc[canvas]->cd(iF+1);
      
      TH2F *h0=h2D[iF];
      h0->GetXaxis()->SetTitle("E [MeV]");
      h0->GetYaxis()->SetTitle("#DeltaE [MeV]");
      h0->DrawCopy("COLZ");
      gPad->SetLogz();
      l.DrawLatex(0.5,0.85,NameSet[iF].c_str());

      for (int iBand=0;iBand<2;iBand++)
	{
	  TF1 *f=new TF1("f","([0]+[1]*(log10(x/10.)))*[2]",0.,30.);
	  f->SetLineWidth(2);
	  for (int ii=0;ii<3;ii++)
	    {
	      double par[3]={ (iBand==0?pLiBand[0]:pBeBand[0]), (iBand==0?pLiBand[1]:pBeBand[1]), (ii==0?1.:(ii==1?0.9:1.1)) };
	      f->SetParameters(par);
	      f->SetLineStyle(ii==0?1:2);
	      f->DrawCopy("SAME");
	    }
	}
    }
  canvas++;


  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,1.5*wh);
  cc[canvas]->Divide(2,1);
  cc[canvas]->Draw();  
  
  double RangeLowBe [nF_in][2] ={  {18.5,20.5},   {18.0,19.5} };
  double RangeHighBe[nF_in][2] ={  {22.0,24.0}, {23.0,24.5} };

  double RangeLowLiP[nF_in][2]   ={  {16.5,17.5}, {16.2,17.0} };
  double RangeHighLiP[nF_in][2]   ={  {25.2,26.5},{26.,26.5} };
  double RangeLiC[nF_in][2]   ={  {27.0,28.0}, {27.0,28.0} };
  
  TGraphErrors tTheta[nF_in][3];

  for (int iBand=0;iBand<2;iBand++)
    {
      cc[canvas]->cd(iBand+1);

      TH2D *hFrame=new TH2D(Form("hX%d",iBand),"",1000, 15., (iBand==0?30.:25.),1000,0.,(iBand==0?0.2:0.1));
      hFrame->GetXaxis()->SetTitle("E_{rec} [MeV]");
      hFrame->GetYaxis()->SetTitle("A.U");
      hFrame->DrawCopy();
      
      for (int iF=0;iF<nF_in;iF++) 
	{
	  TH1F *h0=h1D[iBand][iF];
	  h0->SetLineColor(color[iF]);
	  h0->SetMarkerColor(color[iF]);
	  h0->DrawCopy("HISTSAME");

	  for (int ipeak=0;ipeak<3;ipeak++)
	    {
	      double *range=NULL;

	      if ( iBand==1 && ipeak==0 ) continue;
	      if ( iBand==1 && ipeak==1 ) range=RangeLowBe[iF];
	      if ( iBand==1 && ipeak==2 ) range=RangeHighBe[iF];

	      if ( iBand==0 && ipeak==0 ) range=RangeLowLiP[iF];
	      if ( iBand==0 && ipeak==1 ) range=RangeHighLiP[iF];
	      if ( iBand==0 && ipeak==2 ) range=RangeLiC[iF];


	      TF1 *f=new TF1("f",AsymGaus,range[0],range[1],5);
	      f->SetNpx(1000);
	      double par[5]={ (range[0]+range[1])/2., 0.5, 0.01, 0., 0.};
	      f->SetParameters(par);
	      h0->Fit(f,"QN","",range[0],range[1]);
	      f->SetLineColor(color[iF]);
	      f->SetLineWidth(1);
	      f->DrawCopy("SAME");	      	      
	      double mean=f->GetMaximumX(), err=0.2; //In the asymmetric gaussian the mean is not the peak position.
	      printf("%d %d peak=%5.2lf (p0=%5.2lf) alpha=%5.2lf\n",iBand,ipeak, mean,f->GetParameter(0),f->GetParameter(3));

	      int iSet;
	      if (iBand==0 && ipeak==2 ) iSet=0; // LiC
	      else if (iBand==0        ) iSet=1; // LiP
	      if (iBand==1 )             iSet=2; // Be
	      
	      TGraphErrors *t0=&tTheta[iF][iSet];
	      int ipoint=t0->GetN();
	      t0->SetPoint(ipoint,thetaSi[iF],mean);
	      t0->SetPointError(ipoint,0.,err);
	    }
	  l.SetTextSize(0.035);
	  l.SetTextColor(kBlack);
	  l.DrawLatex(0.2,0.85,(iBand==0?"Li Band":"Be Band"));
	  l.SetTextSize(0.025);
	  l.SetTextColor(color[iF]);
	  l.DrawLatex(0.2,0.75-0.05*double(iF),NameSet[iF].c_str());
	}
    }
  canvas++;

  //----------------------------------------------------------------

  int colorC[nChannels-1]={kGreen, kOrange, kRed, kBlue  };
  
  cc[canvas]=new TCanvas(Form("c%d",canvas),"",1.5*ww,1.5*wh);
  cc[canvas]->Draw();  

  double beamEne=28.;


  for (int ichan=0;ichan<nChannels-1;ichan++)
    {
      for (int SelForward=1;SelForward>=0;SelForward--)
	{
	  TGraph *tg=new TGraph();
	  GetEnergy(beamEne,ichan,SelForward,tg,10.);
	  if ( tg->GetN()==0 ) continue;
	  printf("%s  | Forward=%d | theta=%5.2lf ke=%5.2lf \n", ChannelName[ichan].c_str(),SelForward,5.,tg->Eval(5.));
	  
	  tg->SetMaximum(30.);
	  tg->SetMinimum(10.);
	  tg->SetLineColor(colorC[ichan]);	  
	  tg->SetLineStyle(SelForward==1?20:24);
	  tg->GetXaxis()->SetTitle("#theta_{Si Detector}");
	  tg->GetYaxis()->SetTitle("KE [MeV]");
	  tg->DrawClone((ichan==0&&SelForward==1?"AL":"L"));

	  if ( SelForward==1 )
	    {
	      l.SetTextColor(colorC[ichan]);
	      l.DrawLatex(0.6,0.35-double(ichan)*0.05,ChannelName[ichan].c_str());
	      l.SetTextColor(kBlack);
	    }
	}
      if ( ichan==1 ) continue;
      for (int iF=0;iF<nF_in;iF++) 
	{
	  int Set;
	  if ( ichan==0     ) Set=0;
	  else if (ichan==2 ) Set=2;
	  else if (ichan==3 ) Set=1;	     	      
	  TGraphErrors *t0=&tTheta[iF][Set];
	  t0->SetMarkerColor(colorC[ichan]);
	  t0->SetMarkerStyle(iF==0?20:24);
	  t0->SetLineColor(colorC[ichan]);
	  t0->DrawClone("P");
	}
    }
  l.DrawLatex(0.2,0.85,Form("E_{beam}=%5.2lf",beamEne));
  canvas++;

  //----------------------------------------------------------------------------------
  // Write histos to file
  TFile *f1 = new TFile( "bananas.root"  ,  "RECREATE" );
  for (int iF=0;iF<nF;iF++)  
    {
      h2D[iF]->Write();
      for (int iBand=0;iBand<2;iBand++)
	{
	  h1D[iBand][iF]->Write();
	  h1D_SiE[iBand][iF]->Write();
	}
    }
  f1->Close();

}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Calib Monitor
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

void CalibMonitor() 
{
  SetStyle();
  
  const int nF_in=7;
  const double E0[nF_in]={18.,26.,28.,18.,26.,28.,28.};
  const double thetaSi=23.;

  int runRef[nF_in]   ={636,627, 622, /*Gold */ 634, 629,623, 715 /*CH2*/};
  double time_i[nF_in]={701,4285,702,           891, 5149,1214. ,1481    }; 
  std::string NameSet[nF_in]={ "Li+Au 18 [MeV]", "Li+Au 26 [MeV]", "Li+Au 28 [MeV]", "Li+C 18 [MeV]", "Li+C 26 [MeV]", "Li+C 28 [MeV]", "Li+C 28 [MeV]" };

  //----------------------------------------------------------------------------------
  TH1F   *h1D[nF_in],*h1D_E[nF_in];
 
  int color[7]={kBlue, kRed, kGreen+2,kBlue, kRed, kGreen+2,kOrange};
  
  nF=nF_in;

  for (int iF=0;iF<nF;iF++)
    {
      iF_plot.push_back(iF);
      frun_i[iF]=runRef[iF]; ftime_i[iF]=time_i[iF];  IsTPC=false;

      TTree *t;
      std::string FileName=Form("reco/run_%d.root",runRef[iF]);
      TFile *f=new TFile(FileName.c_str());
      f->GetObject("reco",t);
      Summary(t,iF);

      //-------------------
      const char *hname1D=Form("h1D_%d",iF);
      h1D[iF]=new TH1F(hname1D,hname1D,(iF/3==0?300:50),395,430.);
      h1D[iF]->Sumw2();
      
      int channel=(iF/3==0?1:0);
      double ELi=GetEnergyAtSi( E0[iF],  channel, true, thetaSi );	    

      std::string Sel="1";//Form("%s>%e && %s<%e",monit.c_str(),range[iF][0],monit.c_str(),range[iF][1]);
      int nn=t->Project(hname1D,Form("%s/%e",monit.c_str(),ELi),Sel.c_str() );
      h1D[iF]->SetDirectory(0);	    
      h1D[iF]->Scale(1./h1D[iF]->Integral());


      //-------------------
      const char *hname1D_E=Form("h1D_E_%d",iF);
      h1D_E[iF]=new TH1F(hname1D_E,hname1D_E,150,-2.,1.);
      h1D_E[iF]->Sumw2();

      nn=t->Project(hname1D_E,Form("%s-%e",MonitEnergy.c_str(),ELi));

      h1D_E[iF]->SetDirectory(0);	    
      h1D_E[iF]->Scale(1./h1D_E[iF]->Integral());

      f->Close();
    }
  
  //----------------------------------------------------------------------------------

  TGraphErrors tgMean[3]; //[mean][Gold/CH2/CH2 post]

  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,1.*wh);
  cc[canvas]->Divide(2,1);
  cc[canvas]->Draw();  

  for (int ipad=0;ipad<2;ipad++)
    {
      cc[canvas]->cd(ipad+1);
      if ( ipad==0)
	{
	  TH2D *hFrame=new TH2D("hXS","",1000,  410., 430.,1000,0.01,0.25);
	  hFrame->GetXaxis()->SetTitle("Si Monitor [adc]/ E_{^{7}Li} [MeV]");
	  hFrame->GetYaxis()->SetTitle("A.U");
	  hFrame->DrawCopy();

	  double x0=0.35,y0=0.55;
	  TLegend *leg=new TLegend(x0,y0,x0+0.2,y0+0.35);  
      
	  for (int iF=0;iF<nF_in;iF++) 
	    {
	      if ( iF==nF-1 ) continue;

	      int color_i=color[iF];
	      int style_i=iF/3+1, mstyle_i=(iF/3==0?20:24);
	      TH1F *h0=h1D[iF];
	      h0->SetLineColor(color_i);   h0->SetLineStyle(style_i);	      
	      h0->SetMarkerColor(color_i); h0->SetMarkerStyle(mstyle_i);	      
	      h0->DrawCopy("E1SAME");
      
	      double peak=h0->GetBinCenter(h0->GetMaximumBin());
	      if ( iF==3 ) peak-=2;
	      double dx=(iF/3==0?0.5:5.);
	      double rate=0.;
	      for (int ibin=0;ibin<h0->GetNbinsX();ibin++)
		{
		  double x=h0->GetBinLowEdge(ibin);
		  if ( x<peak-dx || x>peak+dx ) continue;
		  double y=h0->GetBinContent(ibin);
		  rate+=y;
		}
	      
	      TF1 *f=new TF1("f",AsymGaus,peak-dx, peak+dx,5);
	      f->SetNpx(1000);
	      double par[5]={ peak , 100. , h0->GetMaximum(), 0., 0.};
	      f->SetParameters(par);
	      f->FixParameter(4,0.); f->FixParameter(3,0.);      
	      h0->Fit(f,"QN","",peak-dx, peak+dx);
	      if ( iF/3==0 ) f->ReleaseParameter(3);      
	      h0->Fit(f,"QN","",peak-dx, peak+dx);
	      f->SetLineColor(color_i);  	      f->SetLineWidth(2);
	      f->DrawCopy("SAME");	      	      
	      peak=f->GetMaximumX(); //In the asymmetric gaussian the mean is not the peak position.
	      double err=10.; 
	      printf("iF=%d peak=%5.2lf (mean AG=%5.2lf) %5.2lf rate=%5.2lf \n",iF, peak,f->GetParameter(0),err,rate);
	      leg->AddEntry(h0,NameSet[iF].c_str(),"LP");
	      
	      int channel=(iF<3?1:0);
	      double ELi=GetEnergyAtSi( E0[iF],  channel, true, thetaSi );	    
	      TGraphErrors *t0=&tgMean[(iF==(nF_in-1)?2:channel)];
	      int ipoint=t0->GetN();
	      t0->SetPoint(ipoint, peak*ELi, ELi );
	      t0->SetPointError(ipoint, err, 0.2 );	      
	    }
	  leg->DrawClone();
	  l.DrawLatex(0.65,0.85,Form("#theta_{Si}=%5.2lf",thetaSi));
	}
      else
	{
	  TH2D *hFrame=new TH2D(Form("hXS%d",ipad),"",1000,6000.,12000.,1000, 15., 30.);
	  hFrame->GetYaxis()->SetTitle("E_{^{7}Li} [MeV]");
	  hFrame->GetXaxis()->SetTitle("Si Monitor [adc]");
	  hFrame->DrawCopy();
	  
	  for (int channel=0;channel<3;channel++)
	    {
	      TGraphErrors *t0=&tgMean[channel];
	      t0->SetMarkerColor(color[channel]);
	      t0->SetLineColor(color[channel]);
	      t0->DrawClone("P");

	      if ( channel==2 ) continue;

	      TF1 *f=new TF1("f","[0]+[1]*x",6000,12000.);
	      t0->Fit(f,"QN","",6000,12000.);
	      f->SetLineColor(color[channel]);
	      f->DrawCopy("SAME");

	      printf("%s  E= %e + [adc]* %e \n",ChannelName[channel].c_str(),f->GetParameter(0),f->GetParameter(1));

	      l.SetTextColor(color[channel]); l.SetTextSize(0.03);
	      l.DrawLatex(0.4,0.85-0.05*double(channel),  ChannelName[channel].c_str() );
	      l.SetTextColor(kBlack);
	    }
	}
    }
  canvas++;

 
  //-------------------------------------------------------------

  cc[canvas]=new TCanvas(Form("c%d",canvas),"",2.5*ww,1.5*wh);
  cc[canvas]->Draw();  

  TH2D *hFrame=new TH2D("hXSS","",1000, -1.75, 0.5,1000,0.0,0.4);
  hFrame->GetXaxis()->SetTitle("E_{rec}-E_{true} [MeV]");
  hFrame->GetYaxis()->SetTitle("A.U");
  hFrame->DrawCopy();

  double x0=0.25,y0=0.55;
  TLegend *leg=new TLegend(x0,y0,x0+0.2,y0+0.35);  
    
  for (int iF=0;iF<nF_in;iF++) 
    {
      TH1F *h0=h1D_E[iF];
      h0->SetLineColor(color[iF]);
      h0->SetLineStyle(iF/3==0?1:2);
      h0->SetMarkerStyle(iF/3==0?20:24);
      h0->SetMarkerColor(color[iF]);
      h0->DrawCopy("E1SAME");
      leg->AddEntry(h0,NameSet[iF].c_str(),"LP");

    }
  
  leg->DrawClone();
  canvas++;
  
}
