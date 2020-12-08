#!/bin/sh

cd .

wdir='.'
ifile='test_data.tsv'
ofile='test.pred'

# run MCCP on test data using 5-fold cross-validation
Rscript MCCP_PRS_kFold.R \
$wdir \
$ifile \
2 \
3 \
'c(4:8)' \
$ofile

###add covariates for step2
awk 'BEGIN{FS="\t";OFS="\t"} NR==FNR{a[$1]=$0;next} {if($1 in a) print $0,a[$1]}' \
$ifile $ofile | cut -f 1-5,9- > $ofile"_for_step2"
