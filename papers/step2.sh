#!/bin/sh

wdir='.'
ifile='test.pred_for_step2'
output='test_output'

Rscript performance.R $wdir $ifile "score" 'c(6:10)' $output
