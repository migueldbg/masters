#include "RDPulseFit.hh"

//============================================================
RDPulseFit::RDPulseFit()
{
  status=-1; 
  covstatus=-1;
  chi2=-1;
  ndf=0;
  start=0; end=0;
  par.clear(); epar.clear();
  sipm=-1;
  type=-1;
}

//============================================================
RDPulseFit::~RDPulseFit()
{;}

//============================================================
RDPulseFit::RDPulseFit(const RDPulseFit& right)
{
  status=right.status;
  covstatus=right.covstatus;
  chi2=right.chi2;
  ndf=right.ndf;
  start=right.start; end=right.end;
  par=right.par;
  epar=right.epar;
  sipm=right.sipm;
  type=right.type;
}

//============================================================
const RDPulseFit& RDPulseFit::operator=(const RDPulseFit& right)
{
  status=right.status;
  covstatus=right.covstatus;
  chi2=right.chi2;
  ndf=right.ndf;
  start=right.start; end=right.end;
  par=right.par;
  epar=right.epar;
  sipm=right.sipm;
  type=right.type;

  return *this;
}

//============================================================
int RDPulseFit::operator==(const RDPulseFit& right) const
{
 
  return (status != right.status || 
	  covstatus != right.covstatus || 
	  chi2 != right.chi2 || 
	  ndf != right.ndf || 
	  start != right.start || 
	  end   != right.end || 
	  type != right.type || 
	  sipm != right.sipm)
    ? 0 : 1;
}

ClassImp ( RDPulseFit)
