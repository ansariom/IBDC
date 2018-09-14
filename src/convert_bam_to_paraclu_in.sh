#!/bin/bash

in_bam=$1
outfile=$2
min_density=$3

out_bed=$in_bam.bed
my_awk='{OFS="\t"; if($6=="+") print $1,$6,$2,1; else print $1,$6,$3-1,1}'
bamToBed -i $in_bam | awk '{OFS="\t"; if($6=="+") print $1,$6,$2,1; else print $1,$6,$3-1,1}' > $out_bed

raw_peaks_out=$in_bam.tmp1
software/paraclu-9/paraclu 5 $out_bed > $raw_peaks_out

raw_peaks_out_sum=$in_bam.tmp2
software/paraclu-9/paraclu-cut.sh $raw_peaks_out  | awk 'OFS="\t" {print $1, $3, $4, $1"_"$3"_"$4"_"$2, $6, $2}' | sort -k1,1V -k2,2n -k6,6 > $raw_peaks_out_sum

merged_peaks=$in_bam.tmp3
bedtools merge -i $raw_peaks_out_sum -s -d 10 > $merged_peaks

tmp=$in_bam.tmp4
bedtools sort -i $merged_peaks > $tmp
cat $tmp |  awk '{print $1"\t"$2"\t"$3"\t"$1"."$2"\t.\t"$4}' > $outfile


chromsize_file=~/mitra/data/tair10/tair10_chromSize.txt
bedToBigBed $tmp $chromsize_file $outfile.bb


rm -f $in_bam.tmp*
