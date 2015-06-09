#!/bin/bash

echo "Running for R=0.2"
root -b -q RAA_plot_matched_unmatched_eMaxSumcand.C\+\(2\)
echo "Running for R=0.3"
root -b -q RAA_plot_matched_unmatched_eMaxSumcand.C\+\(3\)
echo "Running for R=0.4"
root -b -q RAA_plot_matched_unmatched_eMaxSumcand.C\+\(4\)

