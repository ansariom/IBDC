#!/bin/bash

roe_feature_map=ibdc_roe-only/features_map.txt
foldchange=2.5
outdir=tile_vs_roe/current_oct2018/$foldchange
tile_coef_basedir=ibdc_tile-only/med_high/$foldchange
roe_coef_basedir=ibdc_roe-only/med_high/$foldchange

if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

roe_coords_out=$outdir/roe_win_coords.txt
tiles_dir=$outdir/tiles
roes_dir=$outdir/roes

if [ ! -d $tiles_dir ]; then
	mkdir -p $tiles_dir
fi

if [ ! -d $roes_dir ]; then
	mkdir -p $roes_dir
fi

# pick one tss name from features_map file to narrow down the search 
#ex_tss="AT1G01010.1_Chr1_650_+_0"
cat $roe_feature_map | head -n 30000 | grep "AT1G01010.1_Chr1_650_+_0" | cut -f2-4 > $roe_coords_out

#grep $ex_tss $roe_feature_map | cut -f2-4 > $roe_coords_out
cp -r --backup=t  $tile_coef_basedir/*/tile*coef* $tiles_dir
cp -r --backup=t $roe_coef_basedir/*/roe*coef* $roes_dir

