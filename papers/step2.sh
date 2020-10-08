#!/bin/sh

traits=( CAD T2D IBD BRCA )

for i in ${traits[@]}
do
 Rscript performance.R $i"mccp_p05_PRS_pheno.pred.addPCs.tsv" $i
done
