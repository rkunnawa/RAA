#!/bin/bash

echo "Running for R=0.2"
root -b -q RAA_dataDrivenUnfoldingErrorCheck.C\+\(2\)
root -b -q RAA_analyze.C\+\(2\)

echo "Running for R=0.3"
root -b -q RAA_dataDrivenUnfoldingErrorCheck.C\+\(3\)
root -b -q RAA_analyze.C\+\(3\)

echo "Running for R=0.4"
root -b -q RAA_dataDrivenUnfoldingErrorCheck.C\+\(4\)
root -b -q RAA_analyze.C\+\(4\)

