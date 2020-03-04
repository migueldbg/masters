This folder contains a bunch of macros for the TPC-calibration analysis.

########## ly_v2_ph1.C######################

  This macro is meant to extract informations about S1 for the single-phase 
  runs taken with Am241 source.
  
  It accepts the name of the root file as input together with some other
  parameters such as, in order: the run number (used in case of an autmatic
  save of spectra), the value of drift field and fitting ranges for plots.

  S1 is corrected starting from the s1 vs TBA distribution. An order two polinomial
  is used as fitting function, and then an evaluation of the TBA in the middle of
  the chamber is done. Finally, S1 is normalized by the ratio bewtween the above mentioned
  TBA over the total value. This is done per each value of drift field.  

  Finally, corrected S1 is fitted by using a Monte Carlo function provided by Davide Franco
  that accounts also for the gaussin smearing of the detector response.

  Error bars (when available) are purely statistical.
  
  It works under ROOT, as:
  .L ly_v2_ph1.C+
  ly_v2_ph1("root_filename",run_number,200,-0.38,-0.05,400,1000)
  in case of a 200 V/cm run taken in Catania.
  
########## ly_v2_ph2.C ########################

 This macro is meant to extract informations about S1 and S2 yields together
  with the electron lifetime for the double-phase runs taken with Am241 source.
  
  It accepts the name of the root file as input together with some other
  parameters such as, in order: the run number (used in case of an autmatic
  save of spectra), the value of drift field, number of bins per each histogram,
  axis ranges of the s1 vs tdrift, s2/s1 vs tdrift, s1 and s2 corrected plots and finally
  fitting range for the above plots.

  S1 is corrected starting from the s1 vs tdrift spectrum. An order three polinomial
  is used as fitting function, and then an evaluation of the drift time in the middle of
  the chamber is done. Finally, S1 is normalized by the ratio bewtween the above mentioned
  drift time over the total value. This is done per each value of drift field. 

  S1 corrected is then used to plot S2/S1 ratio versus the drift time. So from this
  an estimation of the electron lifetime is extracted (tau).
  
  Once the above tau value is known, a correction in S2 is also performed as TMath::Exp(-tdrift/tau). 

  Finally, corrected S1 and S2 are fitted by using a Monte Carlo function provided by Davide Franco
  that accounts also for the gaussin smearing of the detector response.

  Error bars (when available) are purely statistical.
  
  It works under ROOT, as:
  .L ly_v2_ph2.C+
  ly_v2("root_filename",run_number,200,300,500,800,7,25,0,1400,0,25000,400,850,5000,16000,0,80,16,66,14,66)
  in case of a 200 V/cm run taken in Catania.
