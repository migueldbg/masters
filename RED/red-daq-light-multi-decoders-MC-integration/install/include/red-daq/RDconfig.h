#ifndef _RDCONFIG_H 
#define _RDCONFIG_H 

#include <iostream>
#include <cstring>
#include <vector>
#include <map>
#include <string>

#include <TString.h>
#include <TTree.h>

class RDconfig
{
public:
    static RDconfig* GetInstance();
    std::map<std::string, double> get_cfg(std::string fname, bool isMC=false);
    size_t SERsize(){return fSer.size();};
    size_t LScisize(){return fLSciCalib.size();};

    std::vector<double> GetSER();
    std::vector<double> GetSER_RMS();
    std::vector<double> GetBaselineRMS();
    std::vector<double> GetSER_pVino();
    std::vector<double> GetSER_Tau();
    std::vector<double> GetSER_p();
    std::vector<double> GetSER_Sigma();
    std::vector<double> GetLSciCalib();

    std::string GetSERFile(){return fSERFile;};
    void SetSERFile(std::string aname){fSERFile=aname;};
    void SetIsExtendedSER(bool val){fIsExtendedSER = val;};
    bool IsLaser(){return fIsLaser;};
    void SetIsLaser(bool val){fIsLaser = val;};
    void SetConfigFile(TString val);

    void GetSERDefault(int chanID, double &ser, double &ser_rms, double &pVino, double &tau, double &p, double &sigma, double &brms );

    enum ChannelType {kSiPMTop=1, kSiPMBottom=2, kLSci=3, kSi=4,
                      kUndefined=5, kSiPMTopOFF=6, kSiPMBottomOFF=7};
    ChannelType GetChannelType(int chan);
    TString GetChannelName(int chanID);

    TTree* GetConfigurationTree();
    TTree* GetConfigurationTree_MC();
    //int FindChannelID(int row,int column);
    std::pair<int,int> FindCoordinates(int channumber);

    bool IsTop(int chan);
    bool IsBottom(int chan);
    bool IsSiPM(int chan);
    bool IsLSci(int chan);
    bool IsSi(int chan);

private:
    RDconfig();
    static RDconfig* fInstance;
    std::vector<double> fSer, fSer_RMS, fSer_pVino, fSer_Tau, fSer_p, fSer_Sigma, fBaselineRMS;
    std::vector<double> fLSciCalib;
    std::string fSERFile;
    std::map<int,TString> fMap; //key=channel number
    std::map<TString,int> fInvMap; //key=detector name
    std::map<TString, std::pair<int,int> > fCoordinates;
    TString fConfigFile;
    bool fIsSERread;
    bool fIsExtendedSER;

    bool fIsLaser;

    void InitializeSER();
    void InitializeChannelMapping();
    bool fIsConfigRead;

    TTree* fCfg_tree;
    TTree* fCfg_MC_tree;
};

#endif
