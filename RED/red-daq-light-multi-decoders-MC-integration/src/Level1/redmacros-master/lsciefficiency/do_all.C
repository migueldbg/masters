{
  /* 
   Wrapper macro which re-calculates all "light" files 

   .x do_all.C  
  */
  gROOT->ProcessLine(".L lscicalibrator.C+");
  lscicalibrator("$HOME/red-data/run_48_int.root",true,-1); 
  lscicalibrator("$HOME/red-data/run_56_timint.root",true,1);
  lscicalibrator("$HOME/red-data/run_57_int.root",false,1);
  lscicalibrator("$HOME/red-data/run_58_int.root",false,1);
  lscicalibrator("$HOME/red-data/run_74_int.root",false);
  lscicalibrator("$HOME/red-data/run_75_int.root",true,0);
  lscicalibrator("$HOME/red-data/run_77_int.root",true,0);
  lscicalibrator("$HOME/red-data/run_78_int.root",true,0);
  lscicalibrator("$HOME/red-data/run_79_int.root",true,0);
  lscicalibrator("$HOME/red-data/run_80_int.root",true,0);
  lscicalibrator("$HOME/red-data/run_82_int.root",false,2); 

  //Checks
  //efficiency("run_48_int.root.light",9885682,100,"run48.out.root")
  //efficiency("run_58_int.root.light",9885682,100,"run58.out.root")
  //efficiency(" ",13883722,100,"runs75-79.out.root")
  //efficiency("run_74_int.root.light",13883722,100,"run74.out.root")
  //checktiming("run48.out.root","run58.out.root",true)
  //checktiming("run48.out.root","run58.out.root",false)
  //checktiming("runs75-79.out.root","run74.out.root",true)
  //checktiming("runs75-79.out.root","run74.out.root",false)
}

