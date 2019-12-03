#include "EvHeader.hh"


EvHeader::EvHeader() :
 event_num(0),date("NULL"),run_type("undetermined"),
 rec_ver("NULL"),input_file("NULL"),run_number(0),
 n_channels(0),peak_buffer(0),time(0)
{;}

void EvHeader::SetNumberOfBoards(size_t nb)
{
  boardtimes.resize(nb,0);
}

void EvHeader::SetBoardTime(int id,size_t time)
{
  boardtimes.at(id) = time;
}


void EvHeader::Dump()
{
    cout << "Empty method" << endl;

};

void EvHeader::Fill(string s1, string s2, string s3, string s4, int i1, int i2, int i3, int i4, int i5){
    event_num  = i1;     // event number
    date       = s1;
    run_type   = s2;
    rec_ver    = s3;
    input_file = s4;
    run_number  = i2;
    //    n_samples   = i3;
    n_channels  = i4;
    peak_buffer = i5;
}

ClassImp(EvHeader)

