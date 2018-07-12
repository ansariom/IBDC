#!/bin/bash

nums=("2" "2.5" "3" "3.5" "4")
exp_levels=("all" "med_low" "med_high")

#indir=ibdc_roe-only
indir=ibdc_tile-only
outfile=$indir.performance.txt
rm -f $outfile

for e in ${exp_levels[@]}
do
	for i in ${nums[@]}
	do
		head $indir/$e/$i/[1-9]*/*held* | grep -v auroc | grep -v '==' | awk -v l=$e -v f=$i '{print $0"\t"l"\t"f}' >> $outfile
		echo $i
	done
done
