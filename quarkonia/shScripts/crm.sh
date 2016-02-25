for ((i=$1 ; i<=$2;i++)) do
condor_rm $i
done