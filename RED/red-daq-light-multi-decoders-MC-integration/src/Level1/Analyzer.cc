#include "Analyzer.hh"
#include "EvRaw0.hh"

#include "RDconfig.h"
#include "RDBaseline.hh"
#include "RDFind_min.hh"
#include "RDIntegral.hh"
#include "RDFiltering.hh"
#include "RDClusterization.hh"
#include "RDPulseFitter.hh"
#include "RDPulseFit.hh"

using namespace std;

//=============================================================================
Analyzer::Analyzer(bool islaser, const string& customfile) :
  fVerbosity(0),fCustomFile(customfile),fFitS2(false)
{
  fRecEvent = new EvRec0();
  // read configuration file
  //cfg_file = "cfg/" + cfg_file;
  string configfile = (islaser) ?  "cfg/laser.cfg" : "cfg/global.cfg";
  if (customfile != "")   configfile = customfile;

  RDconfig::GetInstance()->SetIsLaser(islaser);
  cfg = RDconfig::GetInstance()->get_cfg(configfile);

  fIntegral = new RDIntegral();
  fBaseline = new RDBaseline();
  fFindMin = new RDFind_min();
  fPFitter = new RDPulseFitter();  //instance of pulse fitter

  int value = cfg["fit_s2"];
  int value2 = cfg["use_baker"];
  if (value)
    fFitS2 = true;
  if (fFitS2)
    {
      fPFitter->SetUseBaker(value2);
      cout << "--> Fitting S2 waveform ";
      cout << (fPFitter->IsUsingBaker() ? "with Baker" :
	       "without Baker") << " and with DS50binning=" <<
	(fPFitter->IsDS50Binning() ? "true" : " false ") << endl;
    }
  else
    cout << "--> NOT fitting S2 waveform " << endl;

}

//=============================================================================

Analyzer::~Analyzer()
{
  if (fRecEvent)
    delete fRecEvent;
  delete fIntegral;
  delete fBaseline;
  delete fFindMin;
  delete fPFitter;
}

//=============================================================================
EvRec0* const Analyzer::GetRecEvent() const
{
  return fRecEvent;
}

//=============================================================================
EvRec0* const  Analyzer::ProcessEvent(const EvRaw0* evt)
{
  //clear
  fRecEvent->Clear();
  //Start the game: read out from raw
  EvHeader* evhead = evt->GetEvHeader();
  int evnumber = evhead->GetEventNum();
  int n_channels = evt->GetChannelNumber();
  map<int,vector<double>* > wfs = evt->GetWFs();
  //Use the first valid channel, to avoid problems
  size_t nsamples = evt->GetSampleNumber(wfs.begin()->first);

  bool IsSat=false;

  //start doing things
  vector<double> wf_all(nsamples,0);
  //Find what's the max channel number
  int maxCh = 0;
  for (map<int,vector<double>* >::iterator it=wfs.begin(); it!=wfs.end(); ++it) {
    int channumber = it->first;
    if (channumber > maxCh)
	  maxCh = channumber;
  }
  //
  //cout << "Maximum channel number: " << maxCh+1 << endl;
  vector<double> ser = RDconfig::GetInstance()->GetSER();
  //vector<double> calib = RDconfig::GetInstance()->GetLSciCalib();
  vector<double> calib(maxCh+1,1.); // = RDconfig::GetInstance()->GetLSciCalib();

  //Now loop on the channels
  vector<double> baseline_mean(maxCh+1,0);
  vector<double> baseline_rms(maxCh+1,0);
  vector<double> min_y(maxCh+1,0);
  vector<double> min_x(maxCh+1,0);
  vector<double> max_y(maxCh+1,0);
  vector<double> max_x(maxCh+1,0);
  vector<double> charge(maxCh+1,0);
  vector<double> f90(maxCh+1,0);
  vector<double> start_time(maxCh+1,0);
  double charge_tot=0.;
  double charge_tot_lsci=0;
  double f90_tot=0.;
  double fpsd_lsci=0;

  for (map<int,vector<double>* >::iterator it=wfs.begin(); it!=wfs.end(); ++it) {
    int channumber = it->first;
    vector<double> *wf=it->second;
    //cout << "Processing channel: " << it->first << " " <<
    //	RDconfig::GetInstance()->GetChannelType(channumber) << endl;

    //baseline
    fBaseline->DoIt(wf,  (int) cfg["baseline_start"], (int) cfg["baseline_stop"]);
    baseline_mean[channumber] = fBaseline -> mean();
    baseline_rms[channumber]  = fBaseline -> rms();

    if (RDconfig::GetInstance()->GetChannelType(channumber) == RDconfig::kSi) {

	    fFindMin->DoIt(wf, (int) 1, (int) wf->size() -1 , (double) fBaseline->mean(),
			              (double) cfg["cdf_threshold"], (double) cfg["fixed_thr"]);
	    min_x[channumber] = fFindMin -> x_min();
	    min_y[channumber] = fFindMin -> y_min();
	    max_x[channumber] = fFindMin -> x_max();
	    max_y[channumber] = fFindMin -> y_max();

	  } else {
	    //minimum & start time
	    fFindMin->DoIt(wf, (int) cfg["charge_start"], (int) cfg["charge_stop"], (double) fBaseline->mean() ,
			              (double) cfg["cdf_threshold"], (double) cfg["fixed_thr"]);
	    min_x[channumber] = fFindMin->x_min();
	    min_y[channumber] = fFindMin->y_min();
	    max_x[channumber] = fFindMin->x_max();
	    max_y[channumber] = fFindMin->y_max();
	  }

    start_time[channumber] = fFindMin->s_time();

    if (RDconfig::GetInstance()->IsSiPM(channumber)) {
	    // charge & f90
	    fIntegral->DoIt(wf, 0 , baseline_mean[channumber],
			                (int) cfg["charge_start"],
			                (int) cfg["charge_stop"],
			                start_time[channumber] + (int) cfg["f90_stop"]);

	    charge[channumber] = fIntegral->charge()/ser[channumber];

	  } else if (RDconfig::GetInstance()->GetChannelType(channumber) == RDconfig::kLSci) {
	    // charge & f90
	    fIntegral->DoIt(wf, 0, baseline_mean[channumber],
			                start_time[channumber] - (int) cfg["charge_prestart_lsci"],
			                start_time[channumber] + (int) cfg["charge_stop_lsci"],
			                start_time[channumber] + (int) cfg["f90_stop_lsci"]);

	    charge[channumber] = fIntegral->charge()/calib[channumber]; // only because LSci are from 0 to 8 FIXME

       /*
	     cout << "Sono qui: " << channumber << " " <<
	     cfg["charge_start"] << " " << cfg["charge_stop_lsci"] << endl;
	     cout << fIntegral->charge() << " " << ser[channumber] << endl;
	     */
	  } else /*non-SiPM, non-LSci*/ {
	    fIntegral->DoIt(wf, 0 , baseline_mean[channumber],
			                (int) cfg["charge_start"],
			                (int) cfg["charge_stop"],
			                start_time[channumber] + (int) cfg["f90_stop"]);

	    charge[channumber] = fIntegral->charge()/ser[channumber];
	  }


    if (fVerbosity) printf ("b %f q %f integral %f ser %f channel %d\n", fBaseline->m(),fBaseline->q(),
                            fIntegral->charge(), ser[channumber], channumber);

    // Exclude channels having a fake ser
    if (RDconfig::GetInstance()->IsSiPM(channumber)) {
	    f90[channumber] = fIntegral->f90();

	    if (ser[channumber] > 0) {
	      charge_tot += charge[channumber];
	      f90_tot += charge[channumber]*f90[channumber];
	    }
	  } else if (RDconfig::GetInstance()->GetChannelType(channumber) == RDconfig::kLSci) {
	    f90[channumber] = 1. - fIntegral->f90();
	    charge_tot_lsci += charge[channumber];
	    fpsd_lsci += charge[channumber]*f90[channumber];
	  }


    if (fVerbosity)	printf("charge %f f90 %f f90_tot %f\n", charge[channumber], f90[channumber], f90_tot);

    // sum top + bottom wfs
    if (RDconfig::GetInstance()->IsSiPM(channumber) && ser[channumber] > 0) {
	    for (size_t k=0;k<wf_all.size();k++) {
	      double adc = wf->at(k) - baseline_mean[channumber];
	      if (-adc > 4500 )	IsSat=true;
 	      wf_all[k] +=  adc/ser[channumber];
	    }
	  }
  }

  if (charge_tot)
    f90_tot /=  charge_tot;
  if (charge_tot_lsci)
    fpsd_lsci /= charge_tot_lsci;


  //Fill
  //take a copy of the EvHeader. Warning: memory leak prone
  fRecEvent->SetEvHeader(new EvHeader(*evhead));
  fRecEvent->SetEvNumber(evnumber);
  fRecEvent->SetBaseMean(baseline_mean);
  fRecEvent->SetBaseRMS(baseline_rms);
  fRecEvent->SetCharge(charge);
  fRecEvent->SetF90(f90);
  fRecEvent->SetChargeTot(charge_tot);
  fRecEvent->SetChargeTotLSci(charge_tot_lsci);
  fRecEvent->SetLSciPSDTot(fpsd_lsci);
  fRecEvent->SetXmin(min_x);
  fRecEvent->SetYmin(min_y);
  fRecEvent->SetXmax(max_x);
  fRecEvent->SetYmax(max_y);
  fRecEvent->SetF90Tot(f90_tot);
  fRecEvent->SetStartTime(start_time);


  //cout << "evnumber " << evnumber << " nsamples " << nsamples <<endl;

  // Sum Analysis (top + bottom)
  // moving average
  RDFiltering *ff = new RDFiltering(&wf_all, nsamples);
  vector<double> *wa_all = ff->mavg(cfg["moving_avg"]);

  // baseline
  fBaseline->DoIt(&wf_all, cfg["baseline_start"], cfg["baseline_stop"]);
  fRecEvent->SetBaselineMeanAll(fBaseline->mean());
  fRecEvent->SetBaselineRMSAll(fBaseline->rms());

  // clustering
  RDClusterization *cc = new RDClusterization(wa_all, (int) cfg["baseline_stop"], nsamples, cfg["cluster_thr"],
                                              (int) cfg["peak_sigma"], (int) cfg["cluster_start"], (int) cfg["cluster_stop"]);

  const vector<int>& vc_start_time = cc->start_time();
  const vector<int>& vc_start      = cc->start();
  const vector<int>& vc_stop       = cc->stop();
  const vector<int>& vc_rep        = cc->rep();

  const size_t n_cluster = vc_start_time.size();
  //cout << "Found:" << n_cluster <<  " clusters" << endl;

  if (n_cluster) {
    for (size_t n = 0; n < n_cluster; n++) {
	    // fill cluster
	    RDCluster* aCluster = new RDCluster();

	    // cluster info
	    aCluster->start_time = vc_start_time[n];
	    aCluster->start      = vc_start[n];
	    aCluster->stop       = vc_stop[n];
	    aCluster->rep        = vc_rep[n];

	    // find min
	    fFindMin->DoIt(&wf_all, aCluster->start, aCluster->stop, (double) fBaseline->mean(), (double) cfg["cdf_threshold"], (double) cfg["fixed_thr"]);
	    aCluster->min_x = fFindMin->x_min();
	    aCluster->min_y = fFindMin->y_min();
	    aCluster->max_x = fFindMin->x_max();
	    aCluster->max_y = fFindMin->y_max();
      aCluster->cdf_time = fFindMin->s_time();
      aCluster->fixed_time = fFindMin->f_time();

	    // charge and shape parameters

      const int ch_start = ((int) aCluster->fixed_time - (int) cfg["cluster_prestart"] < (int) aCluster->start) ?
              (int) aCluster->start : (int) aCluster->fixed_time - (int) cfg["cluster_prestart"];

	    fIntegral->DoIt(&wf_all, 0, (double) fBaseline->mean(), ch_start, aCluster->stop, ch_start + (int) cfg["f90_stop"]);

	    aCluster->charge = fIntegral->charge();

      //cout << "Charge=" << aCluster->charge << " - Start=" <<aCluster->start_time << endl;;
	    aCluster->f90 = fIntegral->f90();
	    aCluster->mean_time = fIntegral->mean_time();
	    aCluster->rms_time = fIntegral->rms_time();

	    //Integrate separately on the channels and find the maximum over the top channels
	    //(used to evaluate the position)
	    int chanMax = -1;
	    double maxValue = -DBL_MAX;

	    for (map<int,vector<double>* >::iterator it=wfs.begin(); it!=wfs.end(); ++it) {
	      int channumber = it->first;
	      vector<double> *wf = it->second;
	      fIntegral->DoIt(wf, 0 , baseline_mean[channumber],
			                  (int) aCluster->start , aCluster->stop,
			                  aCluster->start_time + (int) cfg["f90_stop"]);

	      if (RDconfig::GetInstance()->GetChannelType(channumber) == RDconfig::kSiPMTop) {
		      if (ser[channumber] > 0) {
		        aCluster->tot_charge_top += fIntegral->charge()/ser[channumber];

		        std::pair<int,int> position = RDconfig::GetInstance()->FindCoordinates(channumber);
		        int ix = position.first;
		        int iy = position.second;
		        aCluster->bar_x += ((0.5 + ix)*5/4)*fIntegral->charge()/ser[channumber]; // in cm
		        aCluster->bar_y += ((0.5 + iy)*5/6)*fIntegral->charge()/ser[channumber];

		        if (fIntegral->charge() > maxValue) {
			        maxValue = fIntegral->charge();
			        chanMax = channumber;
			      }
		      }
		      aCluster->charge_top.push_back(fIntegral->charge()/ser[channumber]);
		      /* cout << "Cluster " << n << " channel " << channumber << endl;
		      cout << "It is a TOP, area: " << fIntegral->charge() << endl;
		      cout << "There are: " << aCluster->charge_top.size() << " entries" << endl;
		      */
		    }
	      //Pad with zeroes in order to keep the ordering
	      else if (RDconfig::GetInstance()->GetChannelType(channumber) == RDconfig::kSiPMTopOFF)
		      aCluster->charge_top.push_back(0.);
	      else if (RDconfig::GetInstance()->GetChannelType(channumber) == RDconfig::kSiPMBottom)
		{
		  if (ser[channumber] > 0)
		    aCluster->tot_charge_bottom += fIntegral->charge()/ser[channumber];
		  aCluster->charge_bottom.push_back(fIntegral->charge()/ser[channumber]);
		  /*
		    cout << "Cluster " << n << " channel " << channumber << endl;
		    cout << "It is a BOTTOM, area: " << fIntegral->charge() << endl;
		    cout << "There are: " << aCluster->charge_bottom.size() << " entries" << endl;
		    cout << "Sum is: " << aCluster->tot_charge_bottom << endl;
		  */
		}
	      else if (RDconfig::GetInstance()->GetChannelType(channumber) ==
		       RDconfig::kSiPMBottomOFF)
		aCluster->charge_bottom.push_back(0.);
	    }

	  //Normalize barycenters
	  if (aCluster->tot_charge_top)
	    {
	      aCluster->bar_x /= aCluster->tot_charge_top;
	      aCluster->bar_y /= aCluster->tot_charge_top;
	    }
	  //Save value of the XY with the maximum amplitude
	  if (chanMax > 0)
	    {
	      //cout << "The top channel with maximum charge is " << chanMax << endl;
	      std::pair<int,int> pos = RDconfig::GetInstance()->FindCoordinates(chanMax);
	      aCluster->pos_x = ((0.5 + pos.first)*5/4); // in cm
	      aCluster->pos_y = ((0.5 + pos.second)*5/6);
	    }

	 //The EvRec0 will take care of releasing the memory
         fRecEvent->AddCluster(aCluster);
      }
    }


  // Fit S2
  if ( fFitS2 && fRecEvent->GetNClusters()==2 && !IsSat)
    {
      RDCluster *cs1=fRecEvent->GetCluster(0);
      RDCluster *cs2=fRecEvent->GetCluster(1);
      //Clear fitter and start fitting!
      fPFitter->Clear();

      bool f90= ( cs1->f90>0.1 && cs2->f90 < 0.1 );
      if (  f90  && cs2->charge>10000. )
      	fPFitter->FitS2(&wf_all,cs2->fixed_time,cs2->stop,fBaseline->rms() );

      //---Save fits
      if ( fPFitter->GetnFits()> 0 )
	{
	  for (int ifit=0;ifit<fPFitter->GetnFits();ifit++)
	    {
	      RDPulseFit* aFit = new RDPulseFit();
	      aFit->status      = fPFitter->fit_status(ifit);
	      aFit->covstatus   = fPFitter->fit_covstatus(ifit);
	      aFit->ndf         = fPFitter->fit_ndf(ifit);
	      aFit->start       = fPFitter->fit_start(ifit);
	      aFit->end         = fPFitter->fit_end(ifit);
	      aFit->sipm        = fPFitter->fit_sipm(ifit);
	      aFit->type        = fPFitter->fit_type(ifit);
	      aFit->chi2        = fPFitter->fit_chi2(ifit);
	      for (int ipar=0;ipar<nParPulseFitter;ipar++)
		{
		  aFit->par.push_back(fPFitter->fit_par(ifit,ipar));
		  aFit->epar.push_back(fPFitter->fit_epar(ifit,ipar));
		}

	      fRecEvent->AddFit(aFit);
	    }
	}
    }

  delete ff;
  delete cc;

  return fRecEvent;
}
