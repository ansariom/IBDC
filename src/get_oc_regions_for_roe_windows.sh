#!/bin/bash


if [ $# -lt 3 ]; then
        echo "./compute_oc_overlap.sh [roe_features_wins.bed] [oc_peaks.bed] [outfile] "
        exit
fi

feature_bed=$1
oc_file=`readlink -f $2`
outfile=$3

halfWidth=3000

cat $feature_bed | awk '{split($1, a, "_"); print $8"\t"$9"\t"$10"\t"$1"%"$2"\t.\t"a[4]}' > $feature_bed.tmp

bedtools intersect -wao -a $feature_bed.tmp  -b $oc_file | awk -v pad=$halfWidth '{if ($8 != -1) { max = $7; min = $8; if ($2 > max) {max = $2}; if ($3 < min){min = $3}; start = max; end = min; split($4, a, "_"); modLoc = a[3] + pad; if (a[4] == "-") {relLeft = modLoc - end; relRight= modLoc - start} else {relLeft = start - modLoc; relRight = end - modLoc}; print a[1]"_"a[2]"_"a[3]"_"a[4]"\t"$9"\t"$1"\t"start"\t"end"\t"relLeft"\t"relRight}}'> $outfile


