#!/bin/bash

src_dir=../src/
nums=("2" "2.5" "3"  "3.5")
#nums=( "2.5" "3"  "3.5")
exp_levels=( "med_low" "med_high")
#exp_levels=( "med_high")
model_type=main
#model_type="tfbs_only"
#model_type="oc_only"

indirs=(ibdc_roe-only ibdc_tile-only)
name=(ROE_Model Tiled_Model)


outfile1=${indirs[0]}/$model_type.performance.txt
outfile2=${indirs[1]}/$model_type.performance.txt
rm -f $outfile1
rm -f $outfile2

pattern1=$model_type
pattern2=$model_type
if [ $model_type == "main" ]; then
	pattern1="roe_only"
	pattern2="tile_only"
fi

for e in ${exp_levels[@]}
do
	for i in ${nums[@]}
	do
		head ${indirs[0]}/$e/$i/[1-9]*/$pattern1*held* | grep -v auroc | grep -v '==' | awk -v l=$e -v f=$i '{print $0"\t"l"\t"f}' >> $outfile1
		head ${indirs[1]}/$e/$i/[1-9]*/$pattern2*held* | grep -v auroc | grep -v '==' | awk -v l=$e -v f=$i '{print $0"\t"l"\t"f}' >> $outfile2
		echo $i
	done
done
$src_dir/utils_paper/performance/plot_performances_models.R $outfile1 $outfile2 ROE_$model_type Tiled_$model_type $model_type
