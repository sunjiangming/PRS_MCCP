# PRS_MCP
A R script on translating PRS for clinical use by estimating the confidence of risk prediction at an individual level  
# How to use
sh do_mccp_prs.sh $wdir $ifile $y_col $prs_col 'c("$covar_col")' $ofile<\br>
\n\t$widr: working directory
\n\t$ifile: PRS file
\n\t$y_col: column number indciating the exact phenotype
\n\t$prs_col: column number indciating PRS
\n\t$covar_col: column numbers for covariates
\n\t$ofile: output
