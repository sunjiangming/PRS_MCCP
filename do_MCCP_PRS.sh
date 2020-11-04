#!/bin/sh
# module load R/4.0.2 or set R enviroment
# 
# set working directory, e.g. dir contains PRS profile
wdir=$1
# input file for training and clibration -label-known set.
tr_file=$2
# input file for predicting -label-unknown set.
te_file=$3
# column number for phenotype
y_col=$4
# column number for PRS
prs_col=$5
# column number for covariates
covar_col=$6
# impute missing values or not, impute: 1
impute=$7
# output file name
ofile=$8

## without any covariates
# Rscript MCCP_PRS.R $wdir $tr_file $te_file $y_col $prs_col 'c()' 0 $ofile
##with covariates
Rscript MCCP_PRS.R $wdir $tr_file $te_file $y_col $prs_col $covar_col $impute $ofile
