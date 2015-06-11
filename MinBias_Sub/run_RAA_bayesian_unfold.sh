#!/bin/bash

#
# This script runs the unfolding and plotting macros for R = 0.2, 0.3, and 0.4. The pT bins are defined in the RAA/Headers/utilities.h macro along with the unfolding cut definition. so only change to be made is in that utilities macro and this should run the whole way making all the plots (right now its coded to run for the entires pT range from 65 to 500).  
# 
#

# compile the macros
root -b -l <<EOF
.L RAA_dataDrivenUnfoldingErrorCheck_new.C+
.q 
EOF

root -b -l <<EOF
.L RAA_analyze.C+
.q
EOF

echo "Running for R=0.2"
root -b -q RAA_dataDrivenUnfoldingErrorCheck_new.C\+\(2\)
root -b -q RAA_analyze.C\+\(2\)
cd ../
root -b -q RAA_plot.C\(2\)

cd MinBias_Sub/
echo "Running for R=0.3"
root -b -q RAA_dataDrivenUnfoldingErrorCheck_new.C\+\(3\)
root -b -q RAA_analyze.C\+\(3\)
cd ../
root -b -q RAA_plot.C\(3\)

cd MinBias_Sub/
echo "Running for R=0.4"
root -b -q RAA_dataDrivenUnfoldingErrorCheck_new.C\+\(4\)
root -b -q RAA_analyze.C\+\(4\)
cd ../
root -b -q RAA_plot.C\(4\)

cd MinBias_Sub/
root -b -l <<EOF
.L RAA_plot_finalpaper.C+
.q
EOF

echo "Running the final paper plotting macro "
root -b -q RAA_plot_finalpaper.C

root -b -l <<EOF
.L RAA_plot_Systematics.C+
.q
EOF

echo "Running the systematics plotting macro"
root -b -q RAA_plot_Systematics.C


#echo "Running for R=0.2"
#root -b -q RAA_dataDrivenUnfoldingErrorCheck.C\+\(2,1,2\)
#root -b -q RAA_analyze.C\+\(2,1,2\)

#echo "Running for R=0.3"
#root -b -q RAA_dataDrivenUnfoldingErrorCheck.C\+\(3,1,2\)
#root -b -q RAA_analyze.C\+\(3,1,2\)

#echo "Running for R=0.4"
#root -b -q RAA_dataDrivenUnfoldingErrorCheck.C\+\(4,1,2\)
#root -b -q RAA_analyze.C\+\(4,1,2\)

#echo "Running for R=0.2"
#root -b -q RAA_dataDrivenUnfoldingErrorCheck.C\+\(2,1,3\)
#root -b -q RAA_analyze.C\+\(2,1,3\)

#echo "Running for R=0.3"
#root -b -q RAA_dataDrivenUnfoldingErrorCheck.C\+\(3,1,3\)
#root -b -q RAA_analyze.C\+\(3,1,3\)

#echo "Running for R=0.4"
#root -b -q RAA_dataDrivenUnfoldingErrorCheck.C\+\(4,1,3\)
#root -b -q RAA_analyze.C\+\(4,1,3\)

#echo "Running for R=0.2"
#root -b -q RAA_dataDrivenUnfoldingErrorCheck.C\+\(2,1,4\)
#root -b -q RAA_analyze.C\+\(2,1,4\)

#echo "Running for R=0.3"
#root -b -q RAA_dataDrivenUnfoldingErrorCheck.C\+\(3,1,4\)
#root -b -q RAA_analyze.C\+\(3,1,4\)

#echo "Running for R=0.4"
#root -b -q RAA_dataDrivenUnfoldingErrorCheck.C\+\(4,1,4\)
#root -b -q RAA_analyze.C\+\(4,1,4\)

#echo "Running for R=0.2"
#root -b -q RAA_dataDrivenUnfoldingErrorCheck.C\+\(2,1,5\)
#root -b -q RAA_analyze.C\+\(2,1,5\)

#echo "Running for R=0.3"
#root -b -q RAA_dataDrivenUnfoldingErrorCheck.C\+\(3,1,5\)
#root -b -q RAA_analyze.C\+\(3,1,5\)

#echo "Running for R=0.4"
#root -b -q RAA_dataDrivenUnfoldingErrorCheck.C\+\(4,1,5\)
#root -b -q RAA_analyze.C\+\(4,1,5\)

