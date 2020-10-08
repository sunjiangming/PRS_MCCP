#!/bin/sh

cd .

#CAD
Rscript MCCP_PRS_kFold.R \
. \
UKB_PRS05_meta_4MCCP.noMissing.txt \
19 \
15 \
'c(2:7,12:13,22)' \
CAD_mccp_p05_PRS_pheno.pred

#T2D
Rscript MCCP_PRS_kFold.R \
. \
UKB_PRS05_meta_4MCCP.noMissing.txt \
21 \
17 \
'c(2:7,12:13,22)' \
T2D_mccp_p05_PRS_pheno.pred

#IBD
Rscript MCCP_PRS_kFold.R \
. \
UKB_PRS05_meta_4MCCP.noMissing.txt \
20 \
16 \
'c(2:7,12:13,22)' \
IBD_mccp_p05_PRS_pheno.pred

#BCAC
Rscript MCCP_PRS_kFold.R \
. \
UKB_PRS05_meta_4MCCP.noMissing.females.txt \
18 \
14 \
'c(2:7,12:13,22)' \
BCAC_mccp_p05_PRS_pheno.pred
