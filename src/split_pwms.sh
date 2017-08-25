#!/bin/bash

in_pwms=$1
npwms_pfile=$2
base=$3
outdir=$4

rm -f -r $outdir
if [ ! -d $outdir ];then
	mkdir $outdir
fi

count=1
batch_no=1
outfile=$outdir/$base.$batch_no.txt

while IFS='' read -r line || [[ -n "$line" ]]; do
    if [[ $line == \>* ]]; then
	if [ $count -gt $npwms_pfile ]; then
		batch_no=$(($batch_no+1))
		outfile=$outdir/$base.$batch_no.txt
		count=1
	else
		count=$(($count+1))
	fi
	echo "$line" >> $outfile
    else
	echo "$line" >> $outfile
    fi
    #echo "Text read from file: $line"
done < "$in_pwms"
