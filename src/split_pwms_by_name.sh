#!/bin/bash

in_pwms=$1
npwms_pfile=$2
outdir=$3

rm -f -r $outdir
if [ ! -d $outdir ];then
	mkdir $outdir
fi

count=1
batch_no=1

while IFS='' read -r line || [[ -n "$line" ]]; do
    if [[ $line = \>* ]]; then
		echo $line
		SUBSTRING=$(echo $line| cut -d'>' -f 2 | cut -d' ' -f 2)
		outfile=$outdir/$SUBSTRING.txt
		echo $outfile
		count=1
    else
	l=$(echo $line | cut -d'=' -f 2)
	echo "$l" >> $outfile
    fi
    #echo "Text read from file: $line"
done < "$in_pwms"
