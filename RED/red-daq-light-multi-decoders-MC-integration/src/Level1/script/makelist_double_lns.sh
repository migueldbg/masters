ls -d -1 /storage/DATA-02/darkside/red/rawdata/lns/run_$1/run_* > run.tmp

filename="run.tmp"

rm -rf run_$1.lst
touch run_$1.lst

while read -r line
do
    name="$line"
    if [[ $name = *"b01"* ]]; then    
      echo "1 $name"
      echo "1 $name" >> run_$1.lst
    else
     echo "0 $name"
     echo "0 $name" >> run_$1.lst
    fi
done < "$filename"

rm -f run.tmp
