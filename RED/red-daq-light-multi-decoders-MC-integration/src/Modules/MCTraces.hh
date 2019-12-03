#ifndef _MCTraces_hh_
#define _MCTraces_hh_

#include <vector>
#include <string>
#include <map>

class TGraph;
class TRandom3;

class MCTraces
{
public:
     MCTraces(std::string configfile, const std::string lArFile="../Level1/MC/LArModel.txt");
    ~MCTraces();
    void Reset();
    void SetRecoilType(int type);
    int  GetNClusters();
    std::vector<double> GetTrace(int isipm);
    void GetEvent();
    void AddToTrace(int icluster);
    void AddFluctuations(int isipm);

private:
    // LAr Model
    TGraph*   tRecomb[2];
    TGraph*   tVisible[2];
    TGraph*   tQuenching[2];
    TGraph*   tp_s1[2];
    TRandom3* rnd;

    // Config
    std::map<std::string, double> cfg;

    // General parameters
    double RecoilEnergy;     /* Recoil Energy in keV => rnd->Uniform(RecoilEnergy_min, RecoilEnergy_max) */
    double RecoilEnergy_min; /* Minimum Recoil Energy in keV */
    double RecoilEnergy_max; /* Maximum Recoil Energy in keV */
    int    RecoilType;       /* Recoil Type: 0=>ER/1=>NE */
    int    SourceType;       /* 0 if gamma/beta, 1 if neutron*/
    int    nsamples;
    int    tdrift_sa;        /* drift time in samples => rnd->Integer(tdrift_max_sa-tdrift_min_sa+1) + tdrift_min_sa */
    int    tdrift_min_sa;    /* drift time in samples */
    int    tdrift_max_sa;    /* drift time in samples */

    // S1 parameters
    int    s1_start;         /* start time in samples of the S1 cluster */
    double eps1;             /* PE/ph */
    double tau1_s1;          /* ns */
    double tau2_s1;          /* ns */
    double p_s1;

    // S2 parameters
    double eps2;             /* PE/e- */
    double tau1_s2;          /* ns */
    double tau2_s2;          /* ns */
    double p_s2;
    double transitTime_s2;   /* ns */
    double sigma_s2;         /* ns */
    double sigma0;           /* mm */
    double DL;               /* mm^2/us */
    double Vdrift;           /* mm/us */

    // SiPM parameters
    std::vector<double> ser_sipm;
    std::vector<double> serRMS_sipm;
    std::vector<double> tau_sipm;
    std::vector<double> p_sipm;
    std::vector<double> sigma_sipm;
    std::vector<double> pVino_sipm;
    std::vector<double> baseline_rms;

    // Trace stotage
    int nclusters;
    int cluster_start[3];           /* start time in samples of the cluster */
    std::vector<double> spe[28];
    std::vector<double> tSignal[3]; /* times of signal cluster (0=>S1, 1=>S2, 2=>S3) */
    std::vector<double> trace[28];

    void   CheckParameters();
    void   InitLArModel(std::string filename);
    void   SetSPE(int isipm);
    void   SetRecoilEnergy();
    void   SetDriftTime();
    double GetImpulse(double x, double tau, double p, double sigma);
    void   GetSignal(double &Nph, int &Ne);
    void   GetS1(int nPE, std::vector<double> &tPE);
    void   GetS2(int Ne,  std::vector<double> &tPE, double &s2);
    int    GetSiPM_Id();
    int    GetDuplicatedNpe(int isipm);
    int    GetClusterStart(int icluster);
};

#endif
