ls -d -1 /storage/DATA-02/darkside/red/rawdata/lns/run_$1/run_*$2* > run.tmp

filename="run.tmp"

while read -r line
do
    name="$line"
    echo "0 $name"
    echo "0 $name" >> run_$1_$2.lst
done < "$filename"

rm -f run.tmp
