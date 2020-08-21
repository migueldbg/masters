#!/bin/sh
rm -rf red_reco_*.*

i=0
while IFS=" " read -r run_number cfg_file laser_number ch_mapping
do
   echo " 
#!/bin/sh -f
# 	   Job name (default is name of pbs script file)
#PBS -N red_reco_$i
#
# 	   Resource limits: max. wall clock time during which job can be running
#PBS -l walltime=$2:00:00
#
#          The directive below directs that the standard output and
#          error streams are to be merged, intermixed, as standard
#          output. 
#PBS -j oe
#
#PBS -d $HOME/red-daq-light/src/Level1/script/
##########################################

echo ------------------------------------------------------
echo 'Job is running on node(s): '
cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: executing queue is $PBS_QUEUE
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

# COMMAND to EXECUTE:
i=0;
cd $HOME/red-daq-light/src/Level1/script/ 
./reco_c7.sh $run_number $cfg_file $laser_number $ch_mapping

exit
" > red_reco_$i.sh
  let "i++"
done < $1

let "i--"
for j in $( eval echo {0..$i} )
  do 
    echo $j 
    qsub -N red_reco_$j -j oe -l walltime=$2:00:00 red_reco_$j.sh
done
