#!/bin/bash

inseqs=$1
nseqs_pfile=$2
base=$3
outdir=$4
pwm=$5
score_outdir=$6
psudo_count=0.001

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
ncpu=5
for seq_file in `ls $outdir/$base.*`; do
	score_outbase=$score_outdir/`basename $seq_file`
	echo java -jar software/tfbs_scan.jar CumScore --pseudoCounts $psudo_count -n 1 $pwm $seq_file $score_outbase | SGE_Array -q megraw,bpp -r j8_$score_outbase.logs -m 15G 
done



