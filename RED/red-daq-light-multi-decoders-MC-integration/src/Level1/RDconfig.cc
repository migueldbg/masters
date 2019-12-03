#include "RDconfig.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>

using namespace std;

RDconfig* RDconfig::fInstance = 0;

//=============================================================================
RDconfig::RDconfig() : 
    fIsSERread(false),fIsLaser(false),fIsConfigRead(false)
{
    //Defaults
    fSERFile = "cfg/ser.cfg";
    fConfigFile = "cfg/channelmapping.cfg";

    TString A[6][4] = {{"A1", "A2", "A3", "A4"},
                       {"B1", "B2", "A5", "B4"},
                       {"C1", "C2", "B3", "B5"},
                       {"C5", "D1", "C4", "C3"},
                       {"D5", "E2", "D4", "D3"},
                       {"E3", "E4", "D2", "E5"}};

    for (int iy = 0; iy < 6; iy++) {
        for (int ix = 0; ix < 4; ix++) {
            std::pair<int,int> pos(ix,iy);
            fCoordinates.insert(std::pair<TString,std::pair<int,int> >
                                (A[iy][ix],pos));
        }
    }
}
//=============================================================================

void RDconfig::SetConfigFile(TString val)
{
    fConfigFile = val;
}
//=============================================================================
void RDconfig::InitializeChannelMapping()
{  
    ifstream mapfile(fConfigFile.Data());

    //No file: use default values
    if (!mapfile.is_open())
    {
        cout << endl;
        fprintf (stderr,"No channel map file %s was found. Exit now!",
                 fConfigFile.Data());
        std::exit(1);
    }
    else //Read from the file
    {
        cout << endl;
        cout << "CHANNEL MAP from file " << fConfigFile << ": " << endl;

        for (int i=0; mapfile.good(); i++)
        {
            int val;
            TString key;
            mapfile >> val >> key;
            // Avoids double-reading of last line
            if (!mapfile.eof())
            {
                fMap.insert ( std::pair<int,TString>(val,key) );
                fInvMap.insert ( std::pair<TString,int>(key,val) );
                cout << "Channel: " << val << "-->" << key << endl;
            }
        }
        mapfile.close();
        if (fMap.size() > 0)
            cout << "Read " << fMap.size() << " channels from file " << endl;
        else
        {
            cerr << "No chans could be read from the file " <<  fConfigFile
                 << endl;
            cerr << " Please check" << endl;
            std::exit(1);
        }
    }
    fIsConfigRead = true;
}


//=============================================================================

void RDconfig::GetSERDefault(int chanID, double &ser, double &ser_rms, double &ser_pVino, double &tau, double &p, double &sigma, double &bRMS)
{
    const TString sipm[28]={ "F2", "F3", "F4", "F5",
                             "A1", "A2", "A3", "A4", "A5",
                             "B1", "B2", "B3", "B4", "B5",
                             "C1", "C2", "C3", "C4", "C5",
                             "D1", "D2", "D3", "D4", "D5",
                             "E2", "E3", "E4", "E5" };

    TString aName=GetChannelName(chanID);
    int sipmId=-1;
    for ( int isipm=0;isipm<28;isipm++)
        if ( sipm[isipm]==aName ) { sipmId=isipm; break; }

    //cout << aName << "  Id="<< sipmId<<std::endl;

    ser = -1.;
    tau = p = sigma = ser_rms = ser_pVino = 0.;
    if ( sipmId<0 ) return;

    double ser_i[28]={ 3370.750 , 3332.220 , 3341.620 , 3397.610 , 2493.650 , 2358.780 , 2358.190 ,
                       2484.500 , 2321.650 , 2454.240 , 2347.520 , 2374.780 , 2419.720 , 2401.280 ,
                       2409.340 , 2411.960 , 2360.520 , 2384.920 , 2630.340 , 2421.060 , 2376.230 ,
                       2435.770 , 2374.110 , 2457.540 , 2380.990 , 2420.620 , 2414.500 , 2389.650 };

    double tau_i[28]={ 602.181 , 573.133 , 592.838 , 594.734 , 545.763 , 553.927 , 563.887 ,
                       609.605 , 559.536 , 589.422 , 567.671 , 564.576 , 593.781 , 583.507 ,
                       566.833 , 562.865 , 569.298 , 595.457 , 560.520 , 567.966 , 555.742 ,
                       584.732 , 572.121 , 559.271 , 572.986 , 579.645 , 572.872 , 561.390 };

    double p_i[28]={ 0.935 , 0.934 , 0.936 , 0.942 , 0.930 , 0.929 , 0.930 ,
                     0.930 , 0.930 , 0.928 , 0.930 , 0.928 , 0.929 , 0.930 ,
                     0.929 , 0.929 , 0.930 , 0.926 , 0.928 , 0.930 , 0.926 ,
                     0.928 , 0.927 , 0.930 , 0.933 , 0.928 , 0.928 , 0.931 };

    double sigma_i[28]={ 10.318 , 10.442 , 10.361 , 9.146 , 7.186 , 6.900 , 7.246 ,
                          7.018 ,  7.161 ,  7.630 , 7.170 , 7.240 , 7.251 , 6.733 ,
                          7.771 ,  7.243 ,  7.202 , 7.858 , 7.252 , 7.163 , 7.717 ,
                          7.702 ,  7.747 ,  7.078 , 6.984 , 7.704 , 7.633 , 7.217 };

    ser=ser_i[sipmId];
    ser_rms=0.1*ser;
    ser_pVino=0.26;
    tau=tau_i[sipmId];
    p=p_i[sipmId];
    sigma=sigma_i[sipmId];
    bRMS=5;
}

void RDconfig::InitializeSER()
{
    InitializeChannelMapping();

    ifstream serfile(fSERFile.c_str());
    bool ok=serfile.is_open();

    //Laser run or no ser file
    if (fIsLaser || !ok )
    {
        for (int chanID=0;chanID<48;chanID++)
        {
            double ser,ser_rms,ser_pvino,tau,p,sigma,brms;
            GetSERDefault(chanID, ser, ser_rms, ser_pvino, tau, p, sigma, brms);
            if ( fIsLaser ) ser=1.;
            fSer.push_back(ser);
            fSer_RMS.push_back(ser_rms);
            fSer_pVino.push_back(ser_pvino);
            fSer_Tau.push_back(tau);
            fSer_p.push_back(1.-p);
            fSer_Sigma.push_back(sigma);
            fBaselineRMS.push_back(brms);
        }

        cout << endl;
        cout << "SER CONFIGURATION from default " ;
        if ( fIsLaser )  cout << " ( Laser)  All SER to 1. "<<std::endl;
        else             cout << " (!Laser)  All SER to values in run 997 "<<std::endl;

    }
    else
    {
        cout << endl;
        cout << "SER CONFIGURATION from file " << fSERFile << ": " << endl;
        if ( !fIsExtendedSER )
            cout << " No extended SER file, setting (Tau,p,Sigma) to default in run 997 "  << endl;

        for (int chanID=0; serfile.good(); chanID++)
        {
            double ser,ser_rms,ser_pvino,tau,p,sigma,brms;
            if ( !fIsExtendedSER )
            {
                GetSERDefault(chanID, ser, ser_rms, ser_pvino, tau, p, sigma, brms);
                serfile >> ser;
            }
            else
                serfile >> ser >> ser_rms >> ser_pvino >> tau >> p >> sigma >> brms;

            if (!serfile.eof())
            {
                fSer.push_back(ser);
                fSer_RMS.push_back(ser_rms);
                fSer_pVino.push_back(ser_pvino);
                fSer_Tau.push_back(tau);
                fSer_p.push_back(1.-p);
                fSer_Sigma.push_back(sigma);
                fBaselineRMS.push_back(brms);
            }
        }

        if (fSer.size() > 0)
            cout << "Read " << fSer.size() << " SER values from file " << endl;
        else
        {
            cerr << "No SER could be read from the file " <<  fSERFile  << endl;
            cerr << " Please check" << endl;
            std::exit(1);
        }
        serfile.close();
    }


    for (int i=0;i<fSer.size();i++)
        cout << "Chan "<< i << "  ser=" << fSer.at(i)  << "  ser_rms=" << fSer_RMS.at(i) << "  ser_pvino=" << fSer_pVino.at(i) << "   tau=" << fSer_Tau.at(i) << " p=" << fSer_p.at(i) << "  sigma="<< fSer_Sigma.at(i) << " Baseline RMS [adc]=" << fBaselineRMS.at(i)<< std::endl;

    fIsSERread = true;
    return;
}  


//=============================================================================

RDconfig* RDconfig::GetInstance()
{
    if (!fInstance)
        fInstance = new RDconfig();
    return fInstance;
}

//=============================================================================

std::map<string, double> RDconfig::get_cfg(string fname, bool isMC) 
{
    char cfg_name[100];
    double cfg_val;
    map<string, double> cfg;

    //First, check that file exists
    ifstream testfile(fname.c_str());
    if (!testfile.is_open())
    {
        cerr << "The file " << fname << " does not exist!" << endl;
        abort();
    }
    testfile.close();

    //Ok, now read the TTree
    TTree *fCfg = new TTree((isMC?"cfg_mc_tree":"cfg_tree"), "Configuration file tree");
    fCfg->ReadFile(fname.c_str(), "cfg_name/C:cfg_val/D");

    fCfg->SetBranchAddress("cfg_name", &cfg_name);
    fCfg->SetBranchAddress("cfg_val", &cfg_val);

    if ( !isMC ) cout << endl << "RECONSTRUCTION CONFIGURATION from file " << fname << ": " << endl;
    else         cout << endl << "MC CONFIGURATION from file " << fname << ": " << endl;
    for (int i = 0; i < fCfg->GetEntries(); i++) {
        fCfg->GetEntry(i);
        cout << left << setw(20) << cfg_name << ((cfg_val<0) ? " " : "  ") << cfg_val << endl;
        string s_cfg_name = cfg_name;
        cfg[s_cfg_name] = cfg_val;
    }

    if ( isMC ) fCfg_MC_tree=fCfg->CloneTree();
    else        fCfg_tree=fCfg->CloneTree();
    delete fCfg;

    return cfg;
}

//=============================================================================
bool RDconfig::IsTop(int chanID) 
{
    return (GetChannelType(chanID) == kSiPMTop);
}

//=============================================================================
bool RDconfig::IsBottom(int chanID) 
{
    return (GetChannelType(chanID) == kSiPMBottom);
}

//=============================================================================
bool RDconfig::IsSiPM(int chanID) 
{
    return ((GetChannelType(chanID) == kSiPMTop) ||
            (GetChannelType(chanID) == kSiPMBottom));
}

//=============================================================================
bool RDconfig::IsLSci(int chanID) 
{
    return (GetChannelType(chanID) == kLSci);
}

//=============================================================================
bool RDconfig::IsSi(int chanID) 
{
    return (GetChannelType(chanID) == kSi);
}

//=============================================================================
TTree* RDconfig::GetConfigurationTree()
{
    return fCfg_tree;
}
//=============================================================================
TTree* RDconfig::GetConfigurationTree_MC()
{
    return fCfg_MC_tree;
}

//=============================================================================
std::vector<double> RDconfig::GetSER()
{
    if (!fIsSERread)
        InitializeSER();
    return fSer;
}

//=============================================================================
std::vector<double> RDconfig::GetSER_RMS()
{
    if (!fIsSERread)
        InitializeSER();
    return fSer_RMS;
}

//=============================================================================
std::vector<double> RDconfig::GetSER_pVino()
{
    if (!fIsSERread)
        InitializeSER();
    return fSer_pVino;
}

//=============================================================================
std::vector<double> RDconfig::GetSER_Tau()
{
    if (!fIsSERread)
        InitializeSER();
    return fSer_Tau;
}

//=============================================================================
std::vector<double> RDconfig::GetSER_p()
{
    if (!fIsSERread)
        InitializeSER();
    return fSer_p;
}

//=============================================================================
std::vector<double> RDconfig::GetSER_Sigma()
{
    if (!fIsSERread)
        InitializeSER();
    return fSer_Sigma;
}

//=============================================================================
std::vector<double> RDconfig::GetBaselineRMS()
{
    if (!fIsSERread)
        InitializeSER();
    return fBaselineRMS;
}

//=============================================================================
std::vector<double> RDconfig::GetLSciCalib()
{
    // hard coded LSCi calibration FIXME
    // this calib is based on runs 566-580, see Run-DB Campaign VI (excel)
    double calib[9] = {7520/59.5, 15130/59.5, 6972/59.5, 6510/59.5, 5995/59.5, 6716/59.5, 7361/59.5, 6005/59.5, 6722/59.5};
    for (int i = 0; i < 9; i++ ) fLSciCalib.push_back(calib[i]);

    return fLSciCalib;
}

TString RDconfig::GetChannelName(int chanID)
{
    if (!fIsConfigRead)
        InitializeChannelMapping();
    if (!fMap.count(chanID))
        return "null";
    return fMap.find(chanID)->second;
}

RDconfig::ChannelType RDconfig::GetChannelType(int chanID)
{
    if (!fIsConfigRead)
        InitializeChannelMapping();
    if (!fMap.count(chanID))
        return kUndefined;
    TString val = fMap.find(chanID)->second;
    if (val.BeginsWith("A") || val.BeginsWith("B") ||
            val.BeginsWith("C") || val.BeginsWith("D") ||
            val.BeginsWith("E"))
    {
        return (val.Contains("OFF")) ?  kSiPMTopOFF : kSiPMTop;
    }
    else if (val.BeginsWith("F"))
    {
        return (val.Contains("OFF")) ?  kSiPMBottomOFF : kSiPMBottom;
    }
    else if (val.BeginsWith("L"))
        return kLSci;
    else if (val.BeginsWith("S"))
        return kSi;
    else
        return kUndefined;

}

//=============================================================================

std::pair<int,int> RDconfig::FindCoordinates(int channumber)
{
    std::pair<int,int> result(-1,-1);

    if (!IsTop(channumber))
        return result;
    if (!fMap.count(channumber))
        return result;
    TString channame = fMap.find(channumber)->second;
    if (!fCoordinates.count(channame))
        return result;
    /*cout << "det: " << channumber << " --" << channame << " " <<
    fCoordinates.find(channame)->second.first << " " <<
    fCoordinates.find(channame)->second.second << endl;
  */
    return fCoordinates.find(channame)->second;
}
