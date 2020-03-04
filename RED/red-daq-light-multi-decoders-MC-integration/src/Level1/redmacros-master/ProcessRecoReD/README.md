$git clone git@baltig.infn.it:pandola/redmacros.git  
$cd redmacros  
$cd ProcessRecoReD  
$mkdir build  
$cd build  
$source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.sh  
$cmake ..  
$make  
$make install 
  
$cd ../  
$./build/ProcessRecoReD <input_file>.root

This will create <input_file>_Hist.root in the same location of the input file. So, please copy your input file to somewhere you have a write permission.  

Compile the PlotHistograms.C  
``$g++ -stdlib=libstdc++ `root-config --cflags --libs` PlotHistograms.C -o PlotHistograms``  
If your compiler complains -stdlib=libstdc++, just take it out.

Run PlotHistograms.C  
$./PlotHistograms <input_file>_Hist.root  
And press return in terminal for every plot in order to move to next plots.  

CompareRuns.C compares 1D histograms created by ProcessRecoReD.C.  
Compile the CompareRuns.C  
``$g++ -stdlib=libstdc++ `root-config --cflags --libs` CompareRuns.C -o CompareRuns``  

Run CompareRuns.C  
$./CompareRuns 1525 1526 1527 1528  
This macro assume you have the reconstruction files in "data" directory in a form of run_<run_number>.root   
