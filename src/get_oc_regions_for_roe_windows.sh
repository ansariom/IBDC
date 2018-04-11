#!/bin/bash


if [ $# -lt 3 ]; then
        echo "./compute_oc_overlap.sh [roe_features_wins.bed] [oc_peaks.bed] [outfile] "
        exit
fi

feature_bed=$1
oc_file=`readlink -f $2`
outbasename=$3
outfile=$4
outdir=`readlink -f $5`/oc_$outbasename
nfeatures=$6
ncpu=$7

nfeatures=5439

if [ ! -d $outdir ]; then
        mkdir -p $outdir
fi

outbase=`readlink -f $outdir/$outbasename`

fprefix="$outbase".

# This was initial methid to get each TSS map in separate file but it runs slow so we get 
# number of features then spli the file!

#awk '{print >> $1; close($1)}' $feature_bed

## This requires knowing nfeatures per TSS (compute by hand)
tail -n +2 $feature_bed > $feature_bed.tmp && mv "$feature_bed.tmp" "$feature_bed"

split -a 5 -d -l $nfeatures $feature_bed $fprefix

echo "Start computing OC overlaps ... "
count=0
total=0
for feature_file in `ls $fprefix*`; do
        if [ $count -gt $ncpu ]; then
                echo "$total files compeleted!"
                wait
                let "total=total+count"
                count=0
        else
                let "count=count+1"
        fi
	feature_bed="$feature_file".bed
	cat $feature_file | awk '{split($1, a, "_"); print $8"\t"$9"\t"$10"\t"$1"%"$2"\t.\t"a[4]}' > $feature_bed
	
	bedtools intersect -wao -a $feature_bed -b $oc_file | awk -v pad=$halfWidth '{if ($8 != -1) { max = $7; min = $8; if ($2 > max) {max = $2}; if ($3 < min){min = $3}; start = max; end = min; split($4, a, "_"); modLoc = a[3] + pad; if (a[4] == "-") {relLeft = modLoc - end; relRight= modLoc - start} else {relLeft = start - modLoc; relRight = end - modLoc}; print a[1]"_"a[2]"_"a[3]"_"a[4]"\t"$9"\t"$1"\t"start"\t"end"\t"relLeft"\t"relRight}}'> $outbase.overlap &
done


halfWidth=3000

cat $feature_bed | awk '{split($1, a, "_"); print $8"\t"$9"\t"$10"\t"$1"%"$2"\t.\t"a[4]}' > $feature_bed.tmp

bedtools intersect -wao -a $feature_bed.tmp  -b $oc_file | awk -v pad=$halfWidth '{if ($8 != -1) { max = $7; min = $8; if ($2 > max) {max = $2}; if ($3 < min){min = $3}; start = max; end = min; split($4, a, "_"); modLoc = a[3] + pad; if (a[4] == "-") {relLeft = modLoc - end; relRight= modLoc - start} else {relLeft = start - modLoc; relRight = end - modLoc}; print a[1]"_"a[2]"_"a[3]"_"a[4]"\t"$9"\t"$1"\t"start"\t"end"\t"relLeft"\t"relRight}}'> $outfile


