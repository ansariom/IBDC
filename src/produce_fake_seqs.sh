#!/bin/bash

# produce almost empty files with the exact number of files as upcoming sequences
#peakfile=$1
#ntotal=`grep -c ">" $peakfile`
#echo $ntotal
ntotal=$1
nseq_pfile=$2

fakefile=$3
seq_outdir=$4
seq_outbase=$5

n=$(($ntotal/100))
n=$(($n+1))
p=$(($nseq_pfile/100))

rm -f $fakefile

for i in `seq 1 $n`; do
	echo ">$i" >> $fakefile
	echo "$i" >> $fakefile
done

software/split_fasta_file.sh $fakefile $p $seq_outdir $seq_outbase

