#!/bin/bash

inseqs=$1
nseqs_pfile=$2
base=$3
outdir=$4
pwm=$5
score_outdir=$6
psudo_count=0.01

rm -f -r $score_outdir
if [ ! -d $score_outdir ];then
        mkdir $score_outdir
fi

rm -f -r $outdir
if [ ! -d $outdir ];then
	mkdir $outdir
fi

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $inseqs | tail -n +2 > tmp.fa

split -a 3 -d -l $nseqs_pfile tmp.fa $outdir/$base.

count=1
ncpu=50
for seq_file in `ls $outdir/$base.*`; do
	if [ $count -lt $ncpu ]; then
		echo $count
		score_outbase=$score_outdir/`basename $seq_file`
		java -Xms5G -Xmx10G -jar software/tfbs_scan.jar CumScore --pseudoCounts $psudo_count -n 1 $pwm $seq_file $score_outbase &
		let "count=count+1"
	else
		wait
		count=1
	fi
done

wait


