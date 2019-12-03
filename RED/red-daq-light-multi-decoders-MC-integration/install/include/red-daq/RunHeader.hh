#ifndef RunHeader_HH
#define RunHeader_HH
#include <stdlib.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>


#include<list>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<stdlib.h>
#include<stdio.h>


#include "TROOT.h"
#include "TObject.h"

using namespace std;


class RunHeader : public TObject {
    
protected:
    string   date;
    Int_t    run_number;
    string   run_type;
    Int_t    n_events;     // number of events
    Double_t fbk_bias;
    Double_t sensl_bias;
    Double_t anode;
    Double_t cathode;
    Double_t grid;
    string   comments;
/*
    Double_t    var_1;       // to be implemented
    Double_t    var_2;       // to be implemented
    Double_t    var_3;       // to be implemented
*/
    

public:

    RunHeader();
    ~RunHeader(){};

    // set method
    void SetDate(string d)            { date = d; };
    void SetRunNumber(Int_t num)      { run_number = num; };
    void SetRunType(string rt)        { run_type = rt; };
    void SetNumEvt(int en)            { n_events = en; };
    void SetFbkBias(double fbk)       { fbk_bias = fbk; };
    void SetSenslBias(double sensl)   { sensl_bias = sensl; };
    void SetAnode(double an)          { anode = an; };
    void SetCathode(double cath)      { cathode = cath; };
    void SetGrid(double gri)          { grid = gri; };
    void SetComments(string com)      { comments = com; };
/*
    void SetUserVar1(double v1)       { var_1 = v1; };  // to be implemented
    void SetUserVar2(double v2)       { var_2 = v2; };  // to be implemented
    void SetUserVar3(double v3)       { var_3 = v3; };  // to be implemented
*/
    
    
    // get method
    string   GetDate()          { return date; };
    Int_t    GetRunNumber()     { return run_number; };
    string   GetRunType()       { return run_type; };
    Int_t    GetNumEvt()        { return n_events; };
    Double_t GetFbkBias()       { return fbk_bias; };
    Double_t GetSenslBias()     { return sensl_bias; };
    Double_t GetAnode()         { return anode; };
    Double_t GetCathode()       { return cathode; };
    Double_t GetGrid()          { return grid; };
    string   GetComments()      { return comments; };
/*
    Double_t GetUserVar1()      { return var_1; };  //to be implemented
    Double_t GetUserVar2()      { return var_2; };  //to be implemented
    Double_t GetUserVar3()      { return var_3; };  //to be implemented
*/
    
    
    void Fill(string date, int run_num, string run_type,int num_ev, double fbk, double sensl, double anode, double cathode, double grid, string comments);
    
    
    
    // Dump method
    
    void Dump();

    ClassDef( RunHeader, 2 )

};

#endif
