#!/bin/sh

wdir='.'
ifile='test.pred_for_step2'
output='test_output'

Rscript performance.R $wdir $ifile "score" 'c(7:10)' $output
