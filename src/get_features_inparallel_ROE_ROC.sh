#!/bin/bash

ncpu=50

input_fasta=$1
outdir=$2
final_outfile=$3
roe_fwd=$4
roe_rev=$5
pwm=$6
root_overlap_dir=$7
leaf_overlap_dir=$8

basename=seq
seqs_outdir=$outdir/seqs_tilefeatures

if [ ! -d $seqs_outdir ]; then
	mkdir -p $seqs_outdir
fi

software/split_fasta_file.sh $input_fasta 200 $seqs_outdir $basename

outdir_tmp=$outdir/features
#rm -f -r $outdir_tmp

if [ ! -d $outdir_tmp ]; then
	mkdir -p $outdir_tmp
fi

count=1
for seqfile in `ls $seqs_outdir/$basename.*`; do
	echo $seqfile $outfile $maps_out
	outfile=$outdir_tmp/`basename $seqfile`.Rdat
	java -Xms10G -Xmx20G -jar software/tfbs_scan.jar GenROEROCFeatures $roe_fwd $roe_rev $seqfile $pwm $outfile --leafOC $leaf_overlap_dir --rootOC $root_overlap_dir &	

        if [ $count -lt $ncpu ]; then
		let "count=count+1"
	else
		echo "wait for batch to complete! (submitted = $count)"
                wait
                count=1
	fi
done
wait

echo "All finished!!!"

rm -f $final_outfile

cat $outdir_tmp/*.Rdat > $outdir_tmp/all.Rdat.tmp
head -n 1 $outdir_tmp/all.Rdat.tmp > $final_outfile
grep -v "_FWD_" $outdir_tmp/all.Rdat.tmp >> $final_outfile

# Cleaning up
