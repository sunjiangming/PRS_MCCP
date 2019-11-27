module load R/3.4.2
#set working directory, dir contains PRS profile
wdir=$1
#input file name
ifile=$2
# column number for phenotype
y_col=$3
# column number for PRS
prs_col=$4
# column number for covariates
covar_col=$5
# output file name
ofile=$6

#without any covariates
Rscript MCCP_PRS.R $wdir $ifile  $y_col $prs_col 'c()' $ofile
#single with covariates
Rscript MCCP_PRS.R $wdir $ifile  $y_col $prs_col 'c("$covar_col")' $ofile
