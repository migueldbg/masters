#include "MCTraces.hh"
#include "../Level1/RDconfig.h"

#include <TGraph.h>
#include <TMath.h>
#include <TRandom3.h>

const double FADCWidth = 2.0/*ns/sa*/;
const int    sa_center = 3000*int(2./FADCWidth);

MCTraces::MCTraces(std::string configfile, const std::string lArFile)
{
    InitLArModel(lArFile);
    rnd = new TRandom3(0);

    //------
    if ( configfile.empty() )
        configfile = "cfg/MCconfig.cfg";
    cfg = RDconfig::GetInstance()->get_cfg(configfile,true);

    CheckParameters();

    RecoilEnergy_min = cfg["RecoilEnergy_min"];
    RecoilEnergy_max = cfg["RecoilEnergy_max"];
    RecoilType       = (int) cfg["RecoilType"];
    SourceType       = (int) cfg["SourceType"];
    nsamples         = (int) cfg["nsamples"];
    tdrift_min_sa    = (int) cfg["tdrift_min_sa"];
    tdrift_max_sa    = (int) cfg["tdrift_max_sa"];
    nclusters        = (int) cfg["nclusters"];

    s1_start         = (int) cfg["s1_start"];
    eps1             = cfg["eps1"];
    tau1_s1          = cfg["tau1_s1"];
    tau2_s1          = cfg["tau2_s1"];
    eps2             = cfg["eps2"];
    tau1_s2          = cfg["tau1_s2"];
    tau2_s2          = cfg["tau2_s2"];
    p_s2             = cfg["p_s2"];
    transitTime_s2   = cfg["transitTime_s2"];
    sigma0           = cfg["sigma0"];
    DL               = cfg["DL"];
    Vdrift           = cfg["Vdrift"];

    //-----

    SetDriftTime();
    SetRecoilEnergy();

    ser_sipm     = RDconfig::GetInstance()->GetSER();
    serRMS_sipm  = RDconfig::GetInstance()->GetSER_RMS();
    tau_sipm     = RDconfig::GetInstance()->GetSER_Tau();
    p_sipm       = RDconfig::GetInstance()->GetSER_p();
    sigma_sipm   = RDconfig::GetInstance()->GetSER_Sigma();
    pVino_sipm   = RDconfig::GetInstance()->GetSER_pVino();
    baseline_rms = RDconfig::GetInstance()->GetBaselineRMS();

    for (int isipm=0; isipm<28; isipm++)
    {
        SetSPE(isipm);
        trace[isipm] = std::vector<double>(nsamples, 0.);
    }
}

void MCTraces::CheckParameters()
{
    if ( cfg["RecoilEnergy_min"] <= 0. )                     cfg["RecoilEnergy_min"] = 59.4;
    if ( cfg["RecoilEnergy_max"] <= 0. )                     cfg["RecoilEnergy_max"] = 59.4;
    if ( cfg["RecoilEnergy_max"] < cfg["RecoilEnergy_min"] ) cfg["RecoilEnergy_max"] = cfg["RecoilEnergy_min"];
    if ( cfg["RecoilType"] < 0 || cfg["RecoilType"] > 1 )    cfg["RecoilType"]       = 0;
    if ( cfg["SourceType"] < 0 || cfg["SourceType"] > 1 )    cfg["SourceType"]       = 0;
    if ( cfg["nsamples"] <= 0 )                              cfg["nsamples"]         = 50000;
    if ( cfg["tdrift_min_sa"] < 0 )                          cfg["tdrift_min_sa"]    = 0;
    if ( cfg["tdrift_max_sa"] < cfg["tdrift_min_sa"] )       cfg["tdrift_max_sa"]    = cfg["tdrift_min_sa"];
    if ( cfg["nclusters"] <= 0 || cfg["nclusters"] > 3 )     cfg["nclusters"]        = 3;
    if ( cfg["s1_start"] <= 0 )                              cfg["s1_start"]         = 5000;
    if ( cfg["eps1"] <= 0 )                                  cfg["eps1"]             = 0.1172;
    if ( cfg["tau1_s1"] <= 0 )                               cfg["tau1_s1"]          = 6.7;
    if ( cfg["tau2_s1"] <= 0 )                               cfg["tau2_s1"]          = 1.6e3;
    if ( cfg["eps2"] <= 0 )                                  cfg["eps2"]             = 30.;
    if ( cfg["tau1_s2"] <= 0 )                               cfg["tau1_s2"]          = 11.;
    if ( cfg["tau2_s2"] <= 0 )                               cfg["tau2_s2"]          = 3.43e3;
    if ( cfg["p_s2"] <= 0 )                                  cfg["p_s2"]             = 0.1;
    if ( cfg["transitTime_s2"] <= 0 )                        cfg["transitTime_s2"]   = 1.2e3;
    if ( cfg["sigma0"] <= 0 )                                cfg["sigma0"]           = 0.05;
    if ( cfg["DL"] <= 0 )                                    cfg["DL"]               = 4.12e-4;
    if ( cfg["Vdrift"] <= 0 )                                cfg["Vdrift"]           = 0.93;
}

MCTraces::~MCTraces()
{
    delete rnd;

    for (int iloop=0; iloop<2; iloop++)
    {
        delete tRecomb[iloop];
        delete tVisible[iloop];
        delete tQuenching[iloop];
        delete tp_s1[iloop];
    }
}

void MCTraces::Reset()
{
  if (SourceType == 1){
    // For DD gun generate 6.5% evts NR and 93.5% ER
    double p = 0.065;
    if (rnd->Rndm() < p)
        SetRecoilType(1);
    else
        SetRecoilType(0);
  }

    SetDriftTime();
    SetRecoilEnergy();

    for (int isipm=0; isipm<28; isipm++)
        trace[isipm] = std::vector<double>(nsamples, 0.);
}

void MCTraces::SetDriftTime()
{
    tdrift_sa = rnd->Integer(tdrift_max_sa-tdrift_min_sa+1) + tdrift_min_sa;
    sigma_s2  = TMath::Sqrt(sigma0*sigma0 + 2.*DL*tdrift_sa*FADCWidth/1.e3)/Vdrift *1.e3;

    cluster_start[0] = s1_start;
    cluster_start[1] = s1_start + tdrift_sa;
    cluster_start[2] = s1_start + tdrift_sa + tdrift_max_sa;
}

void MCTraces::SetRecoilEnergy()
{
  if (SourceType == 1){
    // For DD gun energy distribution
    if (RecoilType==0)
    {
        RecoilEnergy = rnd->Exp(1525);  //Exp() returns an exponential deviate: exp(-t/tau)
        while (RecoilEnergy<RecoilEnergy_min || RecoilEnergy>RecoilEnergy_max)
            RecoilEnergy = rnd->Exp(1525);
    }
    else
    {
        RecoilEnergy = rnd->Exp(1525);
        while (RecoilEnergy<RecoilEnergy_min || RecoilEnergy>RecoilEnergy_max/5.)
            RecoilEnergy = rnd->Exp(1525);
    }
  }
  else
    RecoilEnergy = rnd->Uniform(RecoilEnergy_min, RecoilEnergy_max);

  p_s1 = tp_s1[RecoilType]->Eval(RecoilEnergy);
}

void MCTraces::SetRecoilType(int type)
{
    RecoilType = type;
    p_s1 = tp_s1[RecoilType]->Eval(RecoilEnergy);
}

void MCTraces::InitLArModel(std::string filename)
{
    tVisible[0]   = new TGraph(filename.data(), "%lg %lg %*lg %*lg %*lg %*lg %*lg");
    tVisible[1]   = new TGraph(filename.data(), "%lg %*lg %lg %*lg %*lg %*lg %*lg");
    tRecomb[0]    = new TGraph(filename.data(), "%lg %*lg %*lg %lg %*lg %*lg %*lg");
    tRecomb[1]    = new TGraph(filename.data(), "%lg %*lg %*lg %*lg %lg %*lg %*lg");
    tQuenching[0] = new TGraph(filename.data(), "%lg %*lg %*lg %*lg %*lg %lg %*lg");
    tQuenching[1] = new TGraph(filename.data(), "%lg %*lg %*lg %*lg %*lg %*lg %lg");

    for (int recoil=0; recoil<2; recoil++)
        tp_s1[recoil] = new TGraph();

    // Scene ( Nuclear Recoil )
    const int n_nr = 11;
    const double E_nr[n_nr] = {16.9, 25.4, 36.1, 57.2, 75, 100, 125, 150, 175, 200, 10000};                                                             // NR energy [keV]
    const double p_s1_nr[n_nr] = {0.583182, 0.642045, 0.672097, 0.720227, 0.727032, 0.749594, 0.757822, 0.763382, 0.771537, 0.770174, 0.770174 };       // at 200V/cm

    for (int i=0;i<n_nr;i++)
        tp_s1[1]->SetPoint(i, E_nr[i], p_s1_nr[i]);

    // AAr DS-50   ( Electron Recoil )
    const int n_er=90;
    const double p_s1_er[n_er]={  0.390 ,  0.379 ,  0.368 ,  0.360 ,  0.352 ,  0.345 ,  0.340 ,  0.334 ,  0.329 ,  0.325 ,
                                  0.321 ,  0.317 ,  0.315 ,  0.312 ,  0.309 ,  0.307 ,  0.305 ,  0.303 ,  0.301 ,  0.300 ,
                                  0.298 ,  0.297 ,  0.296 ,  0.295 ,  0.294 ,  0.293 ,  0.292 ,  0.291 ,  0.291 ,  0.290 ,
                                  0.289 ,  0.289 ,  0.288 ,  0.288 ,  0.287 ,  0.287 ,  0.286 ,  0.286 ,  0.286 ,  0.285 ,
                                  0.285 ,  0.284 ,  0.284 ,  0.284 ,  0.284 ,  0.284 ,  0.284 ,  0.283 ,  0.283 ,  0.283 ,
                                  0.283 ,  0.283 ,  0.282 ,  0.282 ,  0.282 ,  0.282 ,  0.282 ,  0.282 ,  0.282 ,  0.281 ,
                                  0.282 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,
                                  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,
                                  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281 ,  0.281  };

    const double E_er[n_er]={  3.96 ,   5.77 ,   6.69 ,   7.47 ,   8.01 ,   7.84 ,   9.04 ,  10.49 ,   9.83 ,  10.32 ,
                               10.51 ,  12.77 ,  14.93 ,  12.94 ,  12.94 ,  13.71 ,  15.45 ,  16.57 ,  16.90 ,  16.04 ,
                               16.38 ,  18.49 ,  18.53 ,  19.38 ,  20.23 ,  20.89 ,  21.61 ,  22.39 ,  23.36 ,  23.70 ,
                               24.89 ,  26.11 ,  26.03 ,  27.77 ,  28.14 ,  31.37 ,  29.16 ,  28.10 ,  30.14 ,  32.87 ,
                               31.67 ,  30.46 ,  33.79 ,  36.05 ,  38.17 ,  37.29 ,  38.05 ,  38.81 ,  39.58 ,  40.00 ,
                               40.39 ,  38.76 ,  42.52 ,  43.68 ,  41.99 ,  43.64 ,  45.32 ,  46.30 ,  46.50 ,  44.89 ,
                               45.06 ,  47.01 ,  48.63 ,  50.14 ,  51.75 ,  50.61 ,  51.71 ,  52.80 ,  53.96 ,  53.52 ,
                               57.06 ,  57.74 ,  58.42 ,  59.11 ,  58.22 ,  64.38 ,  63.78 ,  63.17 ,  62.57 ,  61.97 ,
                               61.37 ,  60.76 ,  60.85 ,  60.85 ,  60.85 ,  60.85 ,  60.85 ,  60.85 ,  60.85 ,  60.85  };

    for (int i=0;i<n_er;i++)
        tp_s1[0]->SetPoint(i, E_er[i], p_s1_er[i]);
}

void MCTraces::SetSPE(int isipm)
{
    spe[isipm].clear();

    double tau   = tau_sipm.at(isipm);
    double p     = p_sipm.at(isipm);
    double sigma = sigma_sipm.at(isipm);

    for (int sa=-sa_center+1; sa<sa_center; sa++)
    {
        double low  = FADCWidth*(sa-0.5);
        double high = FADCWidth*(sa+0.5);
        double area = 0.;
        double dt   = (high-low)/100.;
        for (int ii=0; ii<100; ii++)
            area += dt*GetImpulse((low+high)/2., tau, p, sigma);
        spe[isipm].push_back(ser_sipm.at(isipm)*area);
    }
}

double MCTraces::GetImpulse(double x, double tau, double p, double sigma)
{
    double fslow = 1.-p;
    double alpha = 1./tau;
    double sigmasq = sigma*sigma;
    double sigmasq_alpha = sigmasq*alpha;

    double slow = 0.5 * alpha * TMath::Exp(-alpha*(x - 0.5*sigmasq_alpha)) * (1. + TMath::Erf((x - sigmasq_alpha)/(TMath::Sqrt2()*sigma)));
    double fast = 1./(TMath::Sqrt(2.*TMath::Pi())*sigma) * TMath::Exp(-x*x/(2.*sigmasq));

    return fast*(1.-fslow)+slow*fslow;
}

void MCTraces::GetSignal(double &Nph, int &Ne)
{
    Nph=0.; Ne=0;
    if ( RecoilType<0 || RecoilType>1 ) return;

    const double W_Ar=23.6e-3, f_Ar_ER=0.755, f_Ar_NR=0.800;
    const double I_Ar_eff=17.3e-3;   // Eff Ionization potential (to calculate fluctuations)

    double E    = RecoilEnergy;
    double f_Ar = (RecoilType==0) ? f_Ar_ER : f_Ar_NR;

    double Ev   = tVisible[RecoilType]->Eval(E)*E;

    double meanNi = Ev/W_Ar;
    int    maxNi  = int(Ev/I_Ar_eff + 0.5);
    int    Ni     = rnd->Binomial(maxNi, meanNi/double(maxNi));

    if (Ni==0 ) return;
    double recomb = tRecomb[RecoilType]->Eval(E);
    double fq     = tQuenching[RecoilType]->Eval(E);

    Ne  = rnd->Binomial(Ni, 1-recomb);
    Nph = (Ni - Ne + Ni*f_Ar)*fq;
}

void MCTraces::GetS1(int nPE, std::vector<double> &tPE)
{
    tPE.clear();

    for (int ipe=0; ipe<nPE; ipe++)
    {
        double tau = (rnd->Rndm()<p_s1) ? tau1_s1 : tau2_s1;
        double ti  = rnd->Exp(tau);
        tPE.push_back(ti);
    }
}

void MCTraces::GetS2(int Ne, std::vector<double> &tPE, double &s2)
{
    tPE.clear();

    s2=0.;
    int ie=0;
    while (true )
    {
        if ( Ne>0 && ie>=Ne ) break;

        ie++;

        double t0 = rnd->Rndm()*transitTime_s2;
        t0 += rnd->Gaus(0., sigma_s2);

        const int npe = rnd->Poisson(eps2);

        bool stop=false;
        for (int ipe=0; ipe<npe; ipe++)
        {
            double tau = (rnd->Rndm()<p_s2) ? tau1_s2 : tau2_s2;
            double ti  = t0 + rnd->Exp(tau);
            tPE.push_back(ti);
            s2+=1.;
        }
        if ( stop ) break;
    }

    // Create a tail that mimmics data (assuming the tail comes in bunches)
    /*TODO option to activate this (S2tail_option: bunches or SPE)*/
    if ( false )
    {
        int NeTail = rnd->Poisson(0.02*s2/eps2);
        printf("Number of electrons in tail=%d S2=%5.2lf (%5.2lf) eps2=%5.2lf\n",NeTail,s2,s2*0.02,eps2);
        for (int ie=0; ie<NeTail; ie++)
        {
            double t0=rnd->Exp(20.e3);
            const int npe=rnd->Poisson(eps2);
            for (int ipe=0; ipe<npe; ipe++)
            {
                double tau = (rnd->Rndm()<p_s2) ? tau1_s2 : tau2_s2;
                double ti  = t0 + rnd->Exp(tau);
                tPE.push_back(ti);
                s2+=1.;
            }
        }
    }
    // Create a tail that mimmics data (assuming the tail comes in single PE)
    if ( true )
    {
        int nPeTail=rnd->Poisson(0.02*s2);
        for (int ipe=0; ipe<nPeTail; ipe++)
        {
            double t0 = rnd->Exp(20.e3);
            tPE.push_back(t0);
            s2+=1.;
        }
    }
}

int MCTraces::GetNClusters()
{
    return nclusters;
}

std::vector<double> MCTraces::GetTrace(int isipm)
{
    return trace[isipm];
}

int MCTraces::GetSiPM_Id()
{
    // Here a parameterization of light pattern is necessary
    const double fi[2]     = { 0.5/4., 0.5/24. };
    const double fsipm[28] = {          fi[0] , fsipm[0]+fi[0] ,  fsipm[1]+fi[0] , fsipm[2]+fi[0] ,
                               fsipm[3]+fi[1] , fsipm[4]+fi[1] ,  fsipm[5]+fi[1] , fsipm[6]+fi[1] ,
                               fsipm[7]+fi[1] , fsipm[8]+fi[1] ,  fsipm[9]+fi[1] , fsipm[10]+fi[1],
                               fsipm[11]+fi[1], fsipm[12]+fi[1],  fsipm[13]+fi[1], fsipm[14]+fi[1],
                               fsipm[15]+fi[1], fsipm[16]+fi[1],  fsipm[17]+fi[1], fsipm[18]+fi[1],
                               fsipm[19]+fi[1], fsipm[20]+fi[1],  fsipm[21]+fi[1], fsipm[22]+fi[1],
                               fsipm[23]+fi[1], fsipm[24]+fi[1],  fsipm[25]+fi[1], fsipm[26]+fi[1]
                             };
    int isipm=0;
    double seed=rnd->Rndm();
    for (int ii=0;ii<28;ii++)
        if ( seed < fsipm[ii] ) { isipm=ii; break; }

    return isipm;
}

void MCTraces::GetEvent()
{
    double Nph = 0.;
    int     Ne = 0;
    GetSignal(Nph, Ne);

    // S1
    int S1 = rnd->Poisson(Nph*eps1);
    GetS1(S1, tSignal[0]);

    // S2
    if ( nclusters < 2 ) return;
    double S2=0;
    if ( Ne>0 ) GetS2(Ne, tSignal[1], S2);

    // S3
    if ( nclusters < 3 ) return;
    double S3=0;
    GetS2(1, tSignal[2], S3);
}

void MCTraces::AddToTrace(int icluster)
{
    // To be done: fluctuation in PE charge, jittering, Cross Talk and Afterpulsing

    int npe = tSignal[icluster].size();
    int start_cl = GetClusterStart(icluster);
    for (int ipe=0; ipe<npe; ipe++)
    {
        double tPE = tSignal[icluster].at(ipe);
        int isipm  = GetSiPM_Id();

        int low  = sa_center - int( 100./FADCWidth);
        int high = std::min(sa_center + int(4000./FADCWidth), int(spe[isipm].size()));
        int start_i = start_cl + int(tPE/FADCWidth);

        for (int i=low; i<high; i++)
        {
            int index = start_i - sa_center + i;
            if ( index<0 || index>=nsamples ) continue;
            trace[isipm].at(index) += (1+GetDuplicatedNpe(isipm)) * spe[isipm].at(i);
        }
    }
}

void MCTraces::AddFluctuations(int isipm)
{
    for (int i=0; i<nsamples; i++)
        trace[isipm].at(i) += (1+GetDuplicatedNpe(isipm)) * rnd->Gaus(0., baseline_rms.at(isipm));
}

int MCTraces::GetDuplicatedNpe(int isipm)
{
    int npe = 1;
    for (int counter=1; counter<=npe; counter++)
        if (rnd->Rndm() < pVino_sipm.at(isipm))
            npe++;

    //return npe-1;
    return 0;
}

int MCTraces::GetClusterStart(int icluster)
{
    return cluster_start[icluster];
}
