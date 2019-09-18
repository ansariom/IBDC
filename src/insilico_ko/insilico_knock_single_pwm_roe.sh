#!/usr/bin/env bash

# main dir
indir=$1

feature_packdir=$1
model_indir=$2
outdir=$feature_packdir
all_tss_diff_info=$3


spec_outfile=$outdir/single_pwm_knock_spec_roe.txt
outfile=$outdir/single_pwm_knock_result_roe.txt

if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

tmp_roe=$feature_packdir/tmp_roe.txt
cat $feature_packdir/feature_info.tsv | grep -v OC_ | tail -n +2 | awk '{if ($7 != 0) print $1}' | uniq > $tmp_roe

cat  $tmp_roe | awk '{print $1}' | uniq | sort | awk 'BEGIN{print "experiment_id\tfeature\tset_to"}{print $1"\t^"$1"$\t0.0"}' > $spec_outfile

/local/cluster/bin/python ../src/insilico_knockdown.py --scaled_features $feature_packdir/all_tss_features_scaled.tsv \
  --diff_info $all_tss_diff_info \
  --feature_info $feature_packdir/feature_info.tsv \
  --model $model_indir/roe_only_model.sav \
  --tss_percentile 0.05 \
  --knock_spec $spec_outfile     >  $outfile

rm -f $tmp_roe
