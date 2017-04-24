#!/bin/bash

# input is the peak file and output is a bed file contatining the coordonates around tss mode (upstream to downstream)

peaks=`readlink -f $1`
up=$2
down=$3
outfile=$4

awk -F ',' '{if (NR==1) {print $0} else if ($2 == "+"){start=($6-100); end=($6+100) } else {start=$6+100; end=$6-100} print $1"\t"$2"\t"$6"\t"start"\t"end}' aligned.peaks.annotated.capped_leaf.filtered

