#include "RunHeader.hh"

RunHeader::RunHeader() : 
  date("NULL"),run_number(0),run_type("NULL"),
  n_events(0),fbk_bias(0.), sensl_bias(0.),
  anode(0.), cathode(0.), grid(0.), 
  comments("NULL")
{;}

//===========================================================================

void RunHeader::Dump()
{
    cout << "Empty method" << endl;

};

//===========================================================================

void RunHeader::Fill(string s1, int i1, string s2,int i2, double d1, double d2, double d3, double d4, double d5, string s3){
    date       = s1;
     run_number = i1;
    run_type   = s2;
    n_events   = i2;     // number of events
    fbk_bias   = d1;
    sensl_bias = d2;
    anode      = d3;
    cathode    = d4;
    grid       = d5;
    comments   = s3;
/*
    var_1      = d6;    // to be implemented
    var_2      = d7;    // to be implemented
    var_3      = d8;   // to be implemented
*/
    
}

ClassImp(RunHeader)

