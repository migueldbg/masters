#include "Utils.C"

// TBD:  1) TPC scan in horizontal
//       2) Scrutiny of the runs in which PMT0 was placed behind the TPC.

void coincidence()
{

  SetStyle();

  double E0=28.;

  //   E=28 MeV       
  //
  //    Trigger :  ( SiTel  OR SiMon ) TPC slave     
  //                                   661 (Collimator 0)  672/680/715 (Collimator 1)   --> TPC Nominal Position
  //                                   717 (Collimator 1)                               --> 4 spaces up 
  //                                   719 (Collimator 1)                               --> 8 spaces up
  //                                   720 (Collimator 1)                               --> 6 spaces up
  //               ( SiTel AND TPC ) OR SiMon       
  //                                   683/705/716 (Collimator 1)                           --> TPC Nominal Position
  //                                   721/722                                              --> 2 spaces up      --> All channels in
  
  // 

  /*
  const int nF_in=11;
  int    run_i[nF_in]={ 661  ,672  ,680     , 715,  717, 719, 720,           683   ,705 , 716, 721   ,722   };
  double time[nF_in]={  2054., 591., 2243. , 3421, 4012, 936,1819,           9995.,2628.,4151, 10000.,7848. };
  int    plot_i[nF_in]={   1,   0  ,  0    ,   0 ,    0,   0,   1,               1,    0,   0,      1,   0  };
  */

  const int nF_in=7;
  int    run_i[nF_in]={ 672, 676, 680, 715, 717  , 719, 720  };
  double time[nF_in]={  591.,110, 2243.,3421, 4012., 936, 1819.  };
  int    plot_i[nF_in]={   1,  1  , 0,     1, 1, 1, 1  };

  nF=nF_in; IsTPC=true; iF_plot.clear();
  for (int iF=0;iF<nF;iF++) { frun_i[iF]=run_i[iF]; ftime_i[iF]=time[iF]; if ( plot_i[iF]) iF_plot.push_back(iF); }
  Analyze();

  

}
