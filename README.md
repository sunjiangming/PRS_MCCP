# PRS_MCCP
A script on translating PRS for clinical use by estimating the confidence bound of risk prediction at an individual level

# Getting Started
- Clone this repository using the following git command:

  `git clone https://github.com/sunjiangming/PRS_MCP`

- Prerequisites:

    R packages: "glmnet", "doParallel", "foreach", "data.table", "caret", "impute"
    
    R packages (for plot and comparisions):"rms", "viridis", "ggplot2", "gridExtra", "MLmetrics"

# Using MCCP

Here you can esimate an individual’s disease susceptibility by

    sh do_mccp_prs.sh $wdir $tr_file $te_file $y_col $prs_col c("$covar_col") $impute $ofile

 For example:

    sh do_mccp_prs.sh . train.prs test.prs 19 15 'c(2:7,12:13,22)' 0 pred.out


-  Arguments
```
    $widr: working directory
    
    $tr_file: PRS file for training and calibration with known disease status(1:disease; 0:control)
    
    $te_file: PRS file for predicting with unkown disease status as 0
    
    $y_col: column number indciating the phenotype of interest
    
    $prs_col: column number indciating PRS
    
    $covar_col: column numbers indciating covariates
    
    $impute: set as 1 if impute missing values
    
    $ofile: output file
```


-  Outputs

    In total, 6 columns were reported.
```
    predictStatus: Predicted disease status (1:disease; 0:control)
    
    predictCredibility: Credibility of the prediction
    
    predictConfidence: Confidence of the prediction
    
    prob_control: Probability to be control
    
    prob_case: Probability to be case
    
    PRS: Input Polygenic risk scores
```
# Using MCCP for comparision

  This is script used in the paper to evaluate peformance of MCCP in group when compared to classical approach (top x% and bottom x% of PRS).
  
