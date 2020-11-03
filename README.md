# PRS_MCP
A R script on translating PRS for clinical use by estimating the confidence bound of risk prediction at an individual level
# Getting Started
- Clone this repository using the following git command:

  `git clone https://github.com/sunjiangming/PRS_MCP`

- Dependencies:
 R packages: "glmnet", "doParallel", "foreach", "data.table","caret","impute"
 R packages (for plot):"rms", "viridis", "ggplot2", "gridExtra", "MLmetrics"

# Using deMeta
 'sh do_mccp_prs.sh \
    $wdir \
    $tr_file \
    $te_file \
    $y_col \
    $prs_col \
    c("$covar_col") \
    $ofile'
    
-  Arguments
```
    $widr: working directory
    $tr_file: PRS file for training and calibration with known disease status(1:disease; 0:control)
    $te_file: PRS file for predicting with unkown disease status as 0
    $y_col: column number indciating the exact phenotype
    $prs_col: column number indciating PRS
    $covar_col: column numbers for covariates
    $ofile: output
```
