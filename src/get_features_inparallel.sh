#!/bin/bash

ncpu=50

input_fasta=$1
outdir=$2
final_outfile=$3
final_mapfile=$4
tile_nucsUp=$5
tile_nucsDown=$6
tileWin=$7
pwm=$8

basename=seq
seqs_outdir=$outdir/seqs_tilefeatures

if [ ! -d $seqs_outdir ]; then
	mkdir -p $seqs_outdir
fi

software/split_fasta_file.sh $input_fasta 200 $seqs_outdir $basename

outdir_tmp=$outdir/features
#rm -f -r $outdir_tmp

if [ ! -d $outdir_tmp ]; then
	mkdir -p $outdir_tmp
fi

outdir_maps_tmp=$outdir/maps
if [ ! -d $outdir_maps_tmp ]; then
        mkdir -p $outdir_maps_tmp
fi

count=1
for seqfile in `ls $seqs_outdir/$basename.*`; do
	echo $seqfile $outfile $maps_out
	outfile=$outdir_tmp/`basename $seqfile`.Rdat
        maps_out=$outdir_maps_tmp/`basename $seqfile`.maps.txt
	java -Xms10G -Xmx20G -jar software/tfbs_scan.jar GenFeaturesTiledWins --nucsUp $tile_nucsUp --nucsDown $tile_nucsDown --winWidth $tileWin --mapFile $maps_out -n 5 $seqfile $pwm $outfile &	

        if [ $count -lt $ncpu ]; then
		let "count=count+1"
	else
		echo "wait for batch to complete! (submitted = $count)"
                wait
                count=1
	fi
done
wait

echo "All finished!!!"

rm -f $final_outfile

cat $outdir_tmp/*.Rdat > $outdir_tmp/all.Rdat.tmp
head -n 1 $outdir_tmp/all.Rdat.tmp > $final_outfile
grep -v "_FWD_" $outdir_tmp/all.Rdat.tmp >> $final_outfile

cat $outdir_maps_tmp/*.maps.txt > $outdir_maps_tmp/all_maps.tmp
head -n 1 $outdir_maps_tmp/all_maps.tmp > $final_mapfile
grep -v "seq_id" $outdir_maps_tmp/all_maps.tmp >> $final_mapfile

# Cleaning up
