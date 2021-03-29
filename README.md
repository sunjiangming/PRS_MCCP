# PRS-MCCP
A script on translating PRS for clinical use by estimating the confidence bound of risk prediction at an individual level

# Getting Started
- Clone this repository using the following git command:

  `git clone https://github.com/sunjiangming/PRS_MCCP`

- Prerequisites:

    R packages: "glmnet", "doParallel", "foreach", "data.table", "caret" and "impute"
    
    R packages (for comparisions only): "rms", "pROC", "ggplot2", "viridis", "Zlig"
- OS:

    Linux
    
    Mac OSX
    
    The package has been tested on macOS Mojave (10.14.6) and Linux (Red Hat Enterprise Linux Server release 7.6). Essentially, R codes can also run under Windows enviroments.

# Using PRS-MCCP

Here you can esimate an individual’s disease susceptibility by

    sh do_mccp_prs.sh $wdir $tr_file $te_file $y_col $prs_col c("$covar_col") $impute $ofile

    For example:

    sh do_mccp_prs.sh . train.prs test.prs 19 15 'c(2:7,12:13,22)' 0 pred.out


-  Arguments
```
    $wdir: working directory, e.g. places for inputs and outputs
    
    $tr_file: input file including PRS for training and calibration with known disease status(1:disease; 0:control)
    
    $te_file: input file including PRS for predicting with unkown disease status as 0
    
    $y_col: column number indciating the phenotype of interest
    
    $prs_col: column number indciating PRS
    
    $covar_col: column numbers indciating covariates
    
    $impute: set as 1 if you want to impute missing values
    
    $ofile: output file
```


-  Outputs

    In total, 6 columns would be reported. One can filter the predictions by setting a confidence level, as well as comparing the prob_control and prob_case.
```
    predictStatus: Predicted disease status (1:disease; 0:control)
    
    predictCredibility: Credibility of the prediction
    
    predictConfidence: Confidence of the prediction
    
    prob_control: Probability value being control
    
    prob_case: Probability value being case
    
    PRS: Input Polygenic risk scores
```
    
# Example

  Under directory "papers", you can find scripts used in the manuscript to evaluate peformances of PRS-MCCP when disease status are known for both training and test set. Performances in group can be compared with classical approach (top x% and bottom x% of PRS).
 
 ``` 
     cd papers 
 ``` 
 
-  step1: Run PRS-MCCP to a random generated data set (not linked with any real world data, just an example showing how to run the script). This only gives probablitiy values to be case (p1) and control (p0) for each individual, respectively. Less than 1 minute on a Mac Pro 2020, result file "test.pred" would be generated.
```
    sh step1.sh
```
-  step2: Compare performances with top and bottom x%PRS approach using AUC, PPV and NPV when same coverage are achieved. May take long time when estimate 95% confidence intervals of AUC, PPV and NPV. This will save thee Rdata files and produce figures on AUC, PPV and NPV.
```
    sh step2.sh
```
