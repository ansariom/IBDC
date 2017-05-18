#!/bin/bash

indir=`readlink -f $1`
fragment_size=$2
bin_size=$3
m=$4 # normal vs. narrow modes
inbam=`readlink -f $5`
chromsize_file=`readlink -f $6`
covg_outfile=`readlink -f $7`

outbase=out_"f$fragment_size"."b$bin_size".$m
outdir=$indir/$outbase

final_peaks=$outdir/$outbase.peaks.bed

revdir=$outdir/rev
fwddir=$outdir/fwd

echo `pwd`

if [ ! -d $revdir ]; then
        mkdir -p $revdir
fi
if [ ! -d $fwddir ]; then
        mkdir -p $fwddir
fi

echo "Start finding JAMM peaks"
# Get rev and fwd reads
#######################
echo "Get rev and fwd reads"
samtools view -f 16 -b $inbam | bamToBed > $revdir/rev.bed 
samtools view -F 16 -b $inbam | bamToBed > $fwddir/fwd.bed 


# Run JAMM
#############
echo "Run JAMM .."
tmp_rev=$revdir/tmp
rev_peaks_outdir=$revdir/out
bash software/JAMM.sh -T $tmp_rev -s $revdir -g $chromsize_file -o $rev_peaks_outdir -p 40 -f $fragment_size -b $bin_size 

tmp_fwd=$fwddir/tmp
fwd_peaks_outdir=$fwddir/out
bash software/JAMM.sh -T $tmp_fwd -s $fwddir -g $chromsize_file -o $fwd_peaks_outdir -p 40 -f $fragment_size -b $bin_size 


# Combine peaks add the strand info
echo "Combine peaks add the strand info"
rev_rawpeaks=$rev_peaks_outdir/peaks/filtered.peaks.narrowPeak
rev_peaks=$rev_peaks_outdir/peaks/rev_peaks.bed
grep -v ChrC $rev_rawpeaks | grep -v ChrM | awk '{split($0,a,"\t"); split(a[7],b,"."); print a[1]"\t"a[2]"\t"a[3]"\t"a[4]"\t"b[1]"\t""-"}' | bedtools sort > $rev_peaks

fwd_rawpeaks=$fwd_peaks_outdir/peaks/filtered.peaks.narrowPeak
fwd_peaks=$fwd_peaks_outdir/peaks/fwd_peaks.bed
grep -v ChrC $fwd_rawpeaks | grep -v ChrM | awk '{split($0,a,"\t"); split(a[7],b,"."); print a[1]"\t"a[2]"\t"a[3]"\t"a[4]"\t"b[1]"\t""+"}' | bedtools sort > $fwd_peaks

tmp=tmp.txt
#final_peaks=$outdir/$outbase.peaks.bed
cat $rev_peaks $fwd_peaks > $tmp
bedtools sort -i $tmp > $final_peaks
rm -f $tmp

echo "Compute peak's coverage .."
bedtools coverage -s -d -abam $inbam -b $final_peaks > $covg_outfile

echo "Raw peaks Created!!"
