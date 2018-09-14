#!/bin/bash

roe_feature_map=$1
outdir=$2
tile_coef_basedir=$3
roe_coef_basedir=$4

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
ex_tss="AT1G01010.1_Chr1_650_+_0"

#grep $ex_tss $roe_feature_map | cut -f2-4 > $roe_coords_out
cp -r --backup=t  $tile_coef_basedir/*/*coef* $tiles_dir
cp -r --backup=t $roe_coef_basedir/*/*coef* $roes_dir

