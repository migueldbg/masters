#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TTimeStamp.h"
#include "red-daq/EvRec0.hh"

#define UNDEF -100
//-----------------------------------------------------------------------------------
//------------------------------Mapping             ---------------------------------
//-----------------------------------------------------------------------------------

vector<int> topID, botID;
vector<double> topMaxChanCorr;
const int nx=4,ny=6,nxy=24;

const std::string NameChannel[nxy+4]={ "F2", "F3", "F4", "F5",
				    "A1", "A2", "A3", "A4", "A5",
				    "B1", "B2", "B3", "B4", "B5",
				    "C1", "C2", "C3", "C4", "C5",
				    "D1", "D2", "D3", "D4", "D5",
				    "E2", "E3", "E4", "E5" };

void SetMapping( std::string ChanMappingFile, bool debug=false )
{
  topID.clear();   botID.clear();  topMaxChanCorr.clear();
  if ( debug >0 ) printf("Channel Mapping:  %s \n",Form("%s",ChanMappingFile.c_str()));
  std::ifstream mappingfile(Form("%s",ChanMappingFile.c_str()),std::ofstream::in);
  while (true)
    {
      int index;
      std::string name_i,name;
      double corr=1.;
      mappingfile >> index >> name_i ;
      if ( mappingfile.eof() ) break;
      
      if      ( name_i == "A4+B4" ) { name="A4"; corr=0.5; }
      else if ( name_i == "B4+A4" ) { name="A4"; corr=0.5; }
      else if ( name_i == "E3+D5" ) { name="E3"; corr=0.5; }
      else if ( name_i == "D5+E3" ) { name="E3"; corr=0.5; }
      else if ( name_i == "F2+F3" ) { name="F2"; corr=0.5; }
      else if ( name_i == "F3+F2" ) { name="F2"; corr=0.5; }
      else if ( name_i == "F4+F5" ) { name="F4"; corr=0.5; }
      else if ( name_i == "F5+F4" ) { name="F4"; corr=0.5; }
      else      name=name_i;

      int id=UNDEF;
      for ( int j=0;j<28;j++) if ( NameChannel[j]==name ) { id=j; break; }
      if ( id==UNDEF ) continue;


      bool IsBottom=(id<4);

      if ( IsBottom )   botID.push_back(id); 
      else            { topID.push_back(id); topMaxChanCorr.push_back(corr); }

      if ( debug ) std::cout << index<< "  " << name_i<< " "<< name  <<" ("<<id<<") " <<std::endl;

    }
  mappingfile.close();

}

//-----------------------------------------------------------------------------------
//------------------------------XY rec              ---------------------------------
//-----------------------------------------------------------------------------------

TGraph *tasym_xy[3][24];

const double sideTPC=5.;

const std::string  Top_XY_s[ny][nx] = {{"A1", "A2", "A3", "A4"},
				       {"B1", "B2", "A5", "B4"},
				       {"C1", "C2", "B3", "B5"},
				       {"C5", "D1", "C4", "C3"},
				       {"D5", "E2", "D4", "D3"},
				       {"E3", "E4", "D2", "E5"}};

const double Top_XY[ny][nx] = {{ 0, 1,   2,  3},
			       { 5, 6,   4,  8},
			       {10, 11,  7,  9},
			       {14, 15, 13, 12},
			       {19, 20, 18, 17},
			       {21, 22, 16, 23} };


void GetXY(int maxchan, double &xcm, double &ycm, int &x, int &y)
{
  xcm=UNDEF;
  ycm=UNDEF;

  const double x0=sideTPC/double(nx)/2.,y0=sideTPC/double(ny)/2.;

  for (int iy = 0; iy < ny; iy++) 
    for (int ix = 0; ix < nx; ix++)   
      {
	if ( Top_XY[iy][ix] == maxchan ) 
	  { 
	    xcm=(x0-sideTPC/2.) + double(ix) *sideTPC/double(nx); 
	    ycm=(y0-sideTPC/2.) + double(ny-1-iy) *sideTPC/double(ny);
	    x=ix;
	    y=iy;
	    break;
	  }
      }
}


void GetLeftRightChannel(int maxchan, int &chanL, int &chanR )
{
  chanL=UNDEF; chanR=UNDEF;

  // max chan (x,y)
  int ix0=UNDEF,iy0=UNDEF;
  for (int ix=0;ix<nx;ix++)
    for (int iy=0;iy<ny;iy++)
      if ( Top_XY[iy][ix]== maxchan ){ ix0=ix; iy0=iy;}

  //
  if ( ix0==UNDEF || iy0==UNDEF ) return;
  
  if ( ix0==0 )         { chanL=Top_XY[iy0][1];     chanR=Top_XY[iy0][0];}
  else if ( ix0==nx-1 ) { chanL=Top_XY[iy0][ix0-1]; chanR=Top_XY[iy0][ix0];}
  else                  { chanL=Top_XY[iy0][ix0-1]; chanR=Top_XY[iy0][ix0+1];}

  return;
}

void GetUpDownChannel(int maxchan, int &chanU, int &chanD )
{
  chanU=UNDEF; chanD=UNDEF;

  // max chan (x,y)
  int ix0=UNDEF,iy0=UNDEF;
  for (int ix=0;ix<nx;ix++)
    for (int iy=0;iy<ny;iy++)
      if ( Top_XY[iy][ix]== maxchan ){ ix0=ix; iy0=iy;}

  //
  if ( ix0==UNDEF || iy0==UNDEF ) return;
  
  if ( iy0==0 )         { chanU=Top_XY[1][ix0];     chanD=Top_XY[0][ix0];}
  else if ( iy0==ny-1 ) { chanU=Top_XY[iy0-1][ix0]; chanD=Top_XY[iy0][ix0];}
  else                  { chanU=Top_XY[iy0-1][ix0]; chanD=Top_XY[iy0+1][ix0];}

}

void GetXY_Rec( double lra,double uda, int maxchan, double &Xrec, double &Yrec)
{
  Xrec=UNDEF;
  Yrec=UNDEF;
  
  int x,y;
  double Xrec0,Yrec0;
  GetXY(maxchan,Xrec0,Yrec0,x,y);

  double cdfX=0.,cdfY=0.;
  double sideX=sideTPC/double(nx);
  double sideY=sideTPC/double(ny);

  for (int asym=0;asym<2;asym++)
    {
      TGraph *t0=tasym_xy[asym][maxchan];
      if ( asym==0 ) cdfX=1.-t0->Eval(lra);
      else           cdfY=t0->Eval(uda);
    }

  if ( x==0 ) cdfX=(1.-cdfX);
  if ( y==0 ) cdfY=(1.-cdfY);

  Xrec=Xrec0+cdfX*sideX-sideX/2.;
  Yrec=Yrec0+cdfY*sideY-sideY/2.;


}
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------



void Load(int run, std::string MappingFile )
{
  

  TH2F *hXY =new TH2F(Form("hXY_%d",run ),"",200, -6., 6.,200, -6., 6. );
  hXY->Sumw2();

  //--------------
  //  XY mapping 
  //--------------
  SetMapping( MappingFile, true );
  int RunXY=1169;
  TFile *fxy = new TFile(  Form("asym_xy.%d.root",RunXY)  ,  "READ" );  

  for (int asym=0;asym<3;asym++)
    for (int xy=0;xy<24;xy++)
      {
	if ( tasym_xy[asym][xy] !=NULL ) tasym_xy[asym][xy]->Delete();
	std::string name=Form("Mapping_XY_%d_%d",asym,xy);
	tasym_xy[asym][xy]= ( fxy->Get(name.c_str())==NULL ? new TGraph() : (TGraph*) fxy->Get(name.c_str()) );
	tasym_xy[asym][xy]->SetName(Form("Mapping_XY_%d_%d_%d",asym,xy,RunXY));
      }
  fxy->Close();
			   
 
  //--------------
  //  Read tree
  //--------------

  std::string FileName=Form("run_%d.root", run);

  TFile *fin=new TFile(FileName.c_str());
  if ( !fin->IsOpen() ) return;

  TTree* t = (TTree*) fin->Get("reco");
  
  EvRec0* e=new EvRec0();
  t->SetBranchAddress("recoevent",&e);

  printf("Reading tree\n");
  for (int i=0;i<t->GetEntries();i++)
    {
      t->GetEntry(i);

      // Minimal selection
      int nclusters=e->GetNClusters();
      if ( nclusters!=2 ) continue;
      RDCluster *c1=e->GetCluster(0), *c2=e->GetCluster(1);
      if ( c1->f90 < 0.2 || c1->f90>1.0 ) continue; // ER+NR
      if ( c2->f90>0.2 || c2->f90<0. ) continue;
      if ( c2->rep !=1 )                  continue;	  


      // Set XY and MaxFrac
      int XY=UNDEF;     
      double Xrec=UNDEF,Yrec=UNDEF;

      //-- Max channel
      double S2top=S2top=c2->tot_charge_top;
      double max=0.;
      for (int sipm=0;sipm<c2->charge_top.size();sipm++)
	if ( c2->charge_top[sipm]*topMaxChanCorr[sipm]/S2top > max ) { max=c2->charge_top[sipm]*topMaxChanCorr[sipm]/S2top; XY=topID[sipm]-4; }	 	  
  
      //--Left/Down and Up/Down asymmetry
      double lra2,uda2;
      int ileft,iright,iup,idown;
      GetLeftRightChannel(XY,ileft,iright);
      GetUpDownChannel(XY,iup,idown);
      double Sleft=UNDEF,Sright=UNDEF,Sup=UNDEF,Sdown=UNDEF;
      for (int sipm=0;sipm<c2->charge_top.size();sipm++)
	{
	  if ( topID[sipm]-4 == ileft  ) Sleft =c2->charge_top[sipm];
	  if ( topID[sipm]-4 == iright ) Sright=c2->charge_top[sipm];
	  if ( topID[sipm]-4 == iup    ) Sup   =c2->charge_top[sipm];
	  if ( topID[sipm]-4 == idown  ) Sdown =c2->charge_top[sipm];
	}
      if ( Sleft>0 && Sright>0 ) lra2=(Sleft-Sright)/(Sleft+Sright);
      if ( Sup>0 && Sdown>0 )    uda2=(Sup-Sdown)/(Sup+Sdown);
      
      //-- Xrec,Yrec
      GetXY_Rec(lra2,uda2,XY, Xrec,Yrec);
      hXY->Fill(Xrec,Yrec);
    }
  
  TH2F *h0=hXY;
  h0->GetYaxis()->SetTitle("y_{rec}");  
  h0->GetYaxis()->SetTitle("y_{rec}");  
  h0->GetXaxis()->SetRange(h0->GetXaxis()->FindBin(-2.8),h0->GetXaxis()->FindBin(2.8));
  h0->GetYaxis()->SetRange(h0->GetYaxis()->FindBin(-2.8),h0->GetYaxis()->FindBin(2.8));
  h0->DrawCopy("COLZ");

  for (int ix=0;ix<nx+1;ix++)
    {
      double x0=-sideTPC/2.+sideTPC/double(nx)*double(ix);
      TLine l1(x0,-sideTPC/2.,x0,sideTPC/2.);      
      l1.SetLineWidth(2.);
	  l1.DrawClone();
    }
  for (int iy=0;iy<ny+1;iy++)
    {
      double y0=-sideTPC/2.+sideTPC/double(ny)*double(iy);
      TLine l1(-sideTPC/2.,y0,sideTPC/2.,y0);      
      l1.SetLineWidth(2.);
      l1.DrawClone();
    }


  
}
