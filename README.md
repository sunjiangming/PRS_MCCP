# PRS_MCP
A R script on translating PRS for clinical use by estimating the confidence of risk prediction at an individual level  
# How to use
sh do_mccp_prs.sh $wdir $ifile $y_col $prs_col 'c("$covar_col")' $ofile</br>
</br>$widr: working directory
</br>$ifile: PRS file
</br>$y_col: column number indciating the exact phenotype
</br>$prs_col: column number indciating PRS
</br>$covar_col: column numbers for covariates
</br>$ofile: output
