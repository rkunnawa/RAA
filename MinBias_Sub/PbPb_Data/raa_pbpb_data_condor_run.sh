#!/bin/bash

tar -xvf jetRAA_run_PbPb_Data.tar > /dev/null

export SCRAM_ARCH=slc6_amd64_gcc472
source /osg/app/cmssoft/cms/cmsset_default.sh

echo ""
echo "----------------------------------------------------"
echo "Job started on `date` at WN: `hostname` "
echo "Job is running on `uname -a`"

#cd /net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/
#eval `scram list CMSSW`
#eval `scramv1 runtime -sh`
#cd /net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/MinBias_Sub/

#echo "root directory: $ROOTSYS"

df -h

gcc --version

startfile=$1
endfile=$2
radius=$3
outfile=$4
echo "Processing..."

root -b -l <<EOF
.x RAA_read_data_pbpb.C+($startfile,$endfile,$radius,"$outfile")
.q
EOF

mv $outfile /mnt/hadoop/cms/store/user/rkunnawa/rootfiles/JetRAA/June29/

echo "Done!"

echo "Copied output files to hadoop rootfiles/JetRAA/June22" 
