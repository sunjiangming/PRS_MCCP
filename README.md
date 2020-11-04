# PRS_MCCP
A script on translating PRS for clinical use by estimating the confidence bound of risk prediction at an individual level

# Getting Started
- Clone this repository using the following git command:

  `git clone https://github.com/sunjiangming/PRS_MCP`

- Prerequisites:

    R packages: "glmnet", "doParallel", "foreach", "data.table", "caret" and "impute"
    
    R packages (for comparisions): "rms"

# Using MCCP

Here you can esimate an individual’s disease susceptibility by

    sh do_mccp_prs.sh $wdir $tr_file $te_file $y_col $prs_col c("$covar_col") $impute $ofile

    For example:

    sh do_mccp_prs.sh . train.prs test.prs 19 15 'c(2:7,12:13,22)' 0 pred.out


-  Arguments
```
    $wdir: working directory, e.g. places for inputs and outputs
    
    $tr_file: PRS file for training and calibration with known disease status(1:disease; 0:control)
    
    $te_file: PRS file for predicting with unkown disease status as 0
    
    $y_col: column number indciating the phenotype of interest
    
    $prs_col: column number indciating PRS
    
    $covar_col: column numbers indciating covariates
    
    $impute: set as 1 if impute missing values
    
    $ofile: output file
```


-  Outputs

    In total, 6 columns would be reported. One can filter the predictions by setting a confidence level, as well as comparing the prob_control and prob_case.
```
    predictStatus: Predicted disease status (1:disease; 0:control)
    
    predictCredibility: Credibility of the prediction
    
    predictConfidence: Confidence of the prediction
    
    prob_control: Probability to be control
    
    prob_case: Probability to be case
    
    PRS: Input Polygenic risk scores
```
    
# Using MCCP for performance evaluation

  This is the script used in our study to evaluate peformance of MCCP when disease status are known for both training and test set. Performances in group can be compared with classical approach (top x% and bottom x% of PRS). Nevertheless, the classical approach can't report confidence of prediction.
  
-  step1: run MCCP on respective datasets.

-  step2: Compare performances with classical approach using AUC, PPV and NPV when same coverage are achieved.
  
