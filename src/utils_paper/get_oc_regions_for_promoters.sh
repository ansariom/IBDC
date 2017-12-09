#!/bin/bash

up=3000
down=3000
in_fasta=ibdc_roc-model_PWMs-0.09_up1000_down500/peaks_3000_region.fa
outdir=utils_out/oc

if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

promoter_bed=$outdir/root_promoter.bed
outfile=$outdir/root_promoter_open_regions_$up-$down".txt"

oc_peaks_root=ibdc_roc-model_PWMs-0.09_up1000_down500/oc_peaks_root.bed
oc_peaks_leaf=ibdc_roc-model_PWMs-0.09_up1000_down500/oc_peaks_leaf.bed

software/get_tss_region.sh $up $down $up $down $in_fasta > $promoter_bed 
software/get_oc_promoter_intersects.sh $promoter_bed $oc_peaks_root $outfile

promoter_bed=$outdir/leaf_promoter.bed
outfile=$outdir/leaf_promoter_open_regions_$up-$down".txt"

software/get_tss_region.sh $up $down $up $down $in_fasta > $promoter_bed
software/get_oc_promoter_intersects.sh $promoter_bed $oc_peaks_leaf $outfile
