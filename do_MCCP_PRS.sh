module load R/3.4.2
#set working directory, dir contains PRS profile
wdir=$1
#input file for training and clibration -label-known set.
tr_file=$2
#input file for predicting -label-unknown set.
te_file=$3
# column number for phenotype
y_col=$4
# column number for PRS
prs_col=$5
# column number for covariates
covar_col=$6
# output file name
ofile=$7

#without any covariates
Rscript MCCP_PRS.R $wdir $tr_file $te_file $y_col $prs_col 'c()' $ofile
#with covariates
Rscript MCCP_PRS.R $wdir $tr_file $te_file $y_col $prs_col 'c("$covar_col")' $ofile
