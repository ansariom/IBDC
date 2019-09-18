#!/bin/bash

up=$1
down=$2
in_fasta=$3
outdir=$4

if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

promoter_bed=$outdir/root_promoter.bed
outfile=$outdir/root_promoter_open_regions_$up-$down".txt"

oc_peaks_root=$5
oc_peaks_leaf=$6

software/get_tss_region.sh $up $down $up $down $in_fasta > $promoter_bed 
software/get_oc_promoter_intersects.sh $promoter_bed $oc_peaks_root $outfile

promoter_bed=$outdir/leaf_promoter.bed
outfile=$outdir/leaf_promoter_open_regions_$up-$down".txt"

software/get_tss_region.sh $up $down $up $down $in_fasta > $promoter_bed
software/get_oc_promoter_intersects.sh $promoter_bed $oc_peaks_leaf $outfile
