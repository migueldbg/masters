#include "Utils.cc"



//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//     Utils to analyze the root file of a reconstructed normal run:
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------


void Load(int run)
{
  SetStyle();
  SetHisto(run);


  std::string FileName=Form("run_%d.root", run);
  std::ofstream *shapefile=new std::ofstream(Form("s2shape.fits.%d.%d",run,1),std::ofstream::out);

  //---------------------------
  // Main
  //---------------------------

  TFile *f=new TFile(FileName.c_str());
  TTree* t = (TTree*) f->Get("reco");
  
  EvRec0* e=new EvRec0();
  t->SetBranchAddress("recoevent",&e);

  //----
  std::vector<double> fs2,fs1,fs1p,ftd;
  std::vector<int>    fxy;
  std::vector<bool>   finner,frep1;
  //----

  printf("N=%d \n",int(t->GetEntries()));
  int nEvtFit=0,nEvtSel[2]={0,0};

  for (int i=0;i<t->GetEntries();i++)
    {
      t->GetEntry(i);
            
      int EvtNumber=e->GetEvNumber();
      int nclusters=e->GetNClusters();

      //--nclusters
      if ( nclusters<2 || nclusters>3 ) continue;
      RDCluster *c1=e->GetCluster(0), *c2=e->GetCluster(1), *c3=(nclusters==3?e->GetCluster(2):NULL);

      //-Select events with an S1+S2 + (delayed electron)   , S1 non saturated / S2 rep=1
      if ( c1->f90 < 0.1 || c2->f90>0.1 ) continue;
      //if ( c1->isSatTop  ||  c1->isSatBot || c2->isSatTop  ||  c2->isSatBot  ) continue;
      if ( c2->rep !=1 )                      continue;
      if ( ! ( nclusters==2 || ( nclusters==3 && c3->charge <s2maxEcho && (c3->start-c2->start)*FADCWidth/1.e3>tdriftEcho[0] && (c3->start-c2->start)*FADCWidth/1.e3<tdriftEcho[1] ) ) ) continue;
 
      //-Flags
      bool Rep1= ( c1->rep ==1  );                              // S1 Rep
      bool Chg  =( c1->charge>s1min && c1->charge<s1max );      // S1 charge for S2 studies
      bool ChgKr =( c1->charge>s1min_kr && c1->charge<s1max_kr );// S1 charge for Kr selection
           
      
      //--Set XY , tdrift , S2
      double tdrift=(c2->start -c1->start)*FADCWidth/1.e3;

      int XY=UNDEF; 
      double max=0.;
      for (int sipm=0;sipm<c2->charge_top.size();sipm++)
	if ( c2->charge_top[sipm]>max ) { max=c2->charge_top[sipm]; XY=sipm;}
	

      bool IsInner=false;
      for (int ixy=0;ixy<8;ixy++)
	if ( XYCentral[ixy]==XY ) IsInner=true;
      if ( XY<0 || XY>=24 ) continue;

      double S1p= c1->charge*c1->f90/0.3;
      double S1= c1->charge;
      double S2= c2->charge;
      double f100= c1->f90;
      double tba1=(c1->tot_charge_top-c1->tot_charge_bottom)/(c1->tot_charge_top+c1->tot_charge_bottom);
      double tba2=(c2->tot_charge_top-c2->tot_charge_bottom)/(c2->tot_charge_top+c2->tot_charge_bottom);
      double tba3=(nclusters==3 ? (c3->tot_charge_top-c3->tot_charge_bottom)/(c3->tot_charge_top+c3->tot_charge_bottom):UNDEF );

      nEvtSel[IsInner]++;

      //-------------------------------
      //  Save events
      //-------------------------------

      if ( Chg )
	{
	  fs2.push_back(S2); 
	  fs1.push_back(S1);
	  fs1p.push_back(S1p);
	  ftd.push_back(tdrift);
	  fxy.push_back(XY);
	  finner.push_back(IsInner);
	  frep1.push_back(Rep1);
	}


      //-------------------------------
      //  Histogramming
      //-------------------------------

      //-----
      hS1Z[IsInner]->Fill(tdrift,S1);
      if (  ChgKr )
	{
	  hS2Z[IsInner]->Fill(tdrift,S2/S1);
	  hdrift[0][IsInner]->Fill(tdrift);
	  hdrift[nclusters==2?1:2][IsInner]->Fill(tdrift);
	}


      //------
      if ( Rep1 ) 
	{
	  //--S1,f100,f1000
	  hS1[IsInner]->Fill(S1);
	  hS1F90->Fill(S1,c1->f90);	  
	}
	  
      //--tba
      if ( c1->charge>s1min_tba && c1->charge<s1max_tba )
	{  hTBA_Z[0][IsInner]->Fill(tdrift,tba1);  if ( tdrift>tdrift_min_tba &&tdrift<tdrift_max_tba ) hTBA[0][IsInner]->Fill(tba1); }
      if ( c2->charge>s2min_tba  )
	{  hTBA_Z[1][IsInner]->Fill(tdrift,tba2);  hTBA[1][IsInner]->Fill(tba2); }
	

      //--Echo
      if ( nclusters==3  )
	{
	  hTBA_Z[2][IsInner]->Fill(tdrift,tba3);  hTBA[2][IsInner]->Fill(tba3);

	  htEcho[0]->Fill(double(c3->start-c2->start)*FADCWidth/1.e3);
	  htEcho[1]->Fill(double(c3->start-c2->start-c2->mean_time)*FADCWidth/1.e3);
	  hEcho->Fill(c3->charge);
	  hEcho1->Fill(c3->charge,c3->f90);
	  hEcho2->Fill(S2,c3->charge);
	}
      

      //-------------------------------
      //  S2 shape analysis (Raw trace fits)
      //-------------------------------
      S2Shape::FitResult s2fit;
      s2fit.EvtNumber=EvtNumber;
      s2fit.td=tdrift; s2fit.s2=S2; s2fit.xy=XY; s2fit.region=IsInner; s2fit.rep=Rep1; s2fit.s1=S1;
      s2fit.tauS1=UNDEF; s2fit.etauS1=UNDEF;
      s2fit.p=UNDEF; s2fit.tau2=UNDEF; s2fit.T=UNDEF;
      s2fit.type=0;

      const int nfits=e->GetNfits();
      for (int j=0;j<nfits;j++)
	{
	  RDPulseFit *fit=e->GetFits(j);

	  int sipm=fit->sipm;
	  if (fit->type!=1 ) continue;

	  s2fit.type=fit->type;
	  s2fit.tau1=fit->par[0];      s2fit.etau1=fit->epar[0];
	  s2fit.tau2=fit->par[1];      s2fit.etau2=fit->epar[1];
	  s2fit.p=fit->par[2];         s2fit.ep=fit->epar[2];
	  s2fit.T=fit->par[3];         s2fit.eT=fit->epar[3];
	  s2fit.sigma=fit->par[4];     s2fit.esigma=fit->epar[4];
	  s2fit.norm=fit->par[5];      s2fit.enorm=fit->epar[5];
	  s2fit.t0=fit->par[6];        s2fit.et0=fit->epar[6];
	  s2fit.start=fit->par[7];     s2fit.estart=fit->epar[7];
	  s2fit.status=fit->status;    s2fit.covstatus=fit->covstatus; 
	  s2fit.chi2=fit->chi2;        s2fit.ndf=fit->ndf; 

	  *shapefile  << s2fit.type<< " " << s2fit.EvtNumber<< " "<< s2fit.td <<" "<< s2fit.s2<<" "<< s2fit.xy <<" "<<s2fit.region<<" "<< s2fit.rep
		      <<" "<< s2fit.tau1 <<" "<<s2fit.tau2 <<" "<< s2fit.p <<" "<< s2fit.T <<" "<< s2fit.sigma <<" "<< s2fit.norm <<" "<< s2fit.t0 <<" "<< s2fit.start
		      <<" "<< s2fit.chi2 <<" "<< s2fit.ndf <<" "<< s2fit.status <<" "<< s2fit.covstatus 
		      <<" "<< s2fit.etau1 <<" "<<s2fit.etau2 <<" "<< s2fit.ep <<" "<< s2fit.eT <<" "<< s2fit.esigma <<" "<< s2fit.enorm <<" "<< s2fit.et0 <<" "<< s2fit.estart
		      <<" "<< s2fit.s1    <<" "<< s2fit.tauS1 <<" "<<s2fit.etauS1 << std::endl;
	  
	
	}
    }
  printf("Finished n_sel=%d/%d  n_fit=%d \n",nEvtSel[0],nEvtSel[1],nEvtFit);
  shapefile->close();

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // (tdrift0,tdrift1)
  tdrift0=0.;
  {
    int region=1;
    TH1F *h0=hdrift[0][region];

    TF1 *fdrift = new TF1("fdrif", FuncDrift, 25., 65.,3);
    fdrift->SetNpx(10000);
    double parD[3]={50.,10.,1.};
    fdrift->SetParameters(parD);

    h0->Fit(fdrift,"QL0N","",25.,tdriftMax);
    tdrift1=fdrift->GetParameter(0);
  }

  // Lifetime 
  {
    TH2F *h0=hS2Z[1];
    TGraphErrors *t0=new TGraphErrors();
    GetProfile(h0,t0);
    TF1 *f=new TF1("f","[0]*exp(-(x-15.)/[1])",15.,tdrift1);
    f->SetParameter(0,t0->Eval(15.));
    f->SetParameter(1,100.);
    t0->Fit(f,"QN","",15.,tdrift1);
    eLifetime=f->GetParameter(1);
  }

  
  //---(S2corr,S1)
  for (int ii=0;ii<int(fs1.size());ii++)
    {
      double tdrift=ftd.at(ii);
      double S1=fs1.at(ii),S1p=fs1p.at(ii),S2=fs2.at(ii);
      double S2corr=S2/exp(-tdrift/eLifetime );
      bool Inner=finner.at(ii),  Rep=frep1.at(ii);
      if ( !Rep ) continue;
      hS1S2[Inner]->Fill(log10(S1),log10(S2corr)-log10(S1) );
      hS1S2p[Inner]->Fill(log10(S1p),log10(S2corr)-log10(S1p) );	  
      if ( S1>s1min_kr && S1<s1max_kr )
	hS2[Inner]->Fill(S2corr);
    }


  //----Fit log(S2/S1) vs log(S1) and residuals w.r.t fit
  std::vector<double> ratio[24];

  for (int type=0;type<2;type++)
    {      
      GetProfile((type==0?hS1S2[1]:hS1S2p[1]),tS1S2_param[type]);
      
      TF1 *f=new TF1("f","[0]+(x-3.)*[1]",log10(s1min_S2),log10(s1max_S2));         
      f->SetParameter(0,tS1S2_param[type]->Eval(3.));
      if ( IsAm ) f->FixParameter(1,0.35);
      else        f->SetParameter(1,1.);
      tS1S2_param[type]->Fit(f,"QN","",log10(s1min_S2),log10(s1max_S2));
      for (int ipar=0;ipar<2;ipar++) parS2S1[type][ipar]=f->GetParameter(ipar);

      for (int ii=0;ii<int(fs1.size());ii++)
	{	  
	  double tdrift=ftd.at(ii);
	  double S1=fs1.at(ii),S1p=fs1p.at(ii),S2=fs2.at(ii);
	  double S2corr=S2/exp(-tdrift/eLifetime );
	  bool Inner=finner.at(ii),  Rep1=frep1.at(ii);
	  int xy=fxy.at(ii);
	  TGraph *t0=tS1S2[Inner][type];
	  double y;
	  if ( type==0 ) y=log10(S2corr/S1)-f->Eval(log10(S1));
	  if ( type==1 ) y=log10(S2corr/S1p)-f->Eval(log10(S1p));

	  if ( type==1 ) t0->SetPoint(t0->GetN(),tdrift,y);	  
	  if ( type==0  && Rep1 && S1>s1min_S2 &&S1<s1max_S2    )
	    {
	      t0->SetPoint(t0->GetN(),tdrift,y);
	      ratio[xy].push_back(pow(10.,y));
	    }
	}
    }

  //--- S2 XY corrections
  for (int xy=0;xy<24;xy++)
    {
      if ( ratio[xy].size()==0 ) continue;
      double N=double(ratio[xy].size());
      double mean=TMath::Mean((unsigned short)N,&ratio[xy].at(0));
      double rms =TMath::RMS((unsigned short)N,&ratio[xy].at(0));
      printf(" xy=%d  mean=%5.2lf rms=%5.2lf n=%zu \n",xy,mean,rms,ratio[xy].size());
      tS2_XY[0]->SetPoint(xy,double(xy),mean);   	tS2_XY[1]->SetPoint(xy,double(xy),rms);
      tS2_XY[0]->SetPointError(xy,0.,rms/sqrt(N));      tS2_XY[1]->SetPointError(xy,0.,rms/sqrt(2.*N));      
      ratio[xy].clear();
    }

  //--- S1 XY corrections
  for (int xy=0;xy<24;xy++)
    {
      std::vector<double> s1_xy;
      for (int ii=0;ii<int(fs1.size());ii++)
	{	  
	  if ( xy !=fxy.at(ii) ) continue;
	  if ( frep1.at(ii)!=1 ) continue;
	  double S1=fs1.at(ii);
	  if ( S1<s1min_kr || S1>s1max_kr ) continue;
	  
	  s1_xy.push_back(S1);
	}
      if ( s1_xy.size()==0 ) continue;
      double N=double(s1_xy.size());
      double mean=TMath::Mean((unsigned short)N,&s1_xy.at(0));
      double rms =TMath::RMS((unsigned short)N,&s1_xy.at(0));
      printf(" xy=%d  mean=%5.2lf rms=%5.2lf n=%zu \n",xy,mean,rms,s1_xy.size());      
      tS1_XY[0]->SetPoint(xy,double(xy),mean);   	tS1_XY[1]->SetPoint(xy,double(xy),rms/mean);
      tS1_XY[0]->SetPointError(xy,0.,rms/sqrt(N));      tS1_XY[1]->SetPointError(xy,0.,rms/sqrt(2.*N)/mean);      
      s1_xy.clear();
    }
  //----------------------------------------------------------------------------
  WriteHisto(run);

  //----------------------------------------------------------------------------
}




//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//    Routine to display results after processing one run 
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void Show(int run )
{
  std::cout << "The run is taken as: " << (IsAm ? "Am241" : "Kr83m") << 
    std::endl;
  SetStyle();
  ReadHisto(run);
  Plot();
  PlotS2Shape();
  FitGasPocketAndGrid(S2Shape::tPar_XY[0][2][0], tS2_XY[0], 5.211,true);
  SaveCanvasesOnPDF(run); 
}

