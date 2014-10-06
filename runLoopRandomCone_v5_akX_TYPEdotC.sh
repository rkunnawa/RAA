#!/bin/bash



run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x LoopRandomCone_v5rev_akX_TYPE.C(0.3,5000,"DATA","pA","/net/hisrv0001/home/yenjie/scratch/production/CMSSW_5_3_3_patch3/test/PAPhysics-PF-part-v8.root","randomCones_TkpTCut0_ak3_pA_DATA.root")
.q
EOF
cp randomCones_TkpTCut0_ak3_pA_DATA.root /mnt/hadoop/cms/store/user/jrobles/PAanalysis/randomCone/v3
echo "Done with pA DATA!  "
rm randomCones_TkpTCut0_ak3_pA_DATA.root
fi

run="no"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x LoopRandomCone_v5rev_akX_TYPE.C(0.3,5000,"DATA","PbPb","/net/hisrv0001/home/icali/hadoop/HIMinBiasUPC_skimmed/MinBias-reTracking-merged/MinBias_Merged_tracking_all.root","randomCones_TkpTCut0_ak3_PbPb_DATA.root")
.q
EOF
cp randomCones_TkpTCut0_ak3_PbPb_DATA.root /mnt/hadoop/cms/store/user/jrobles/PAanalysis/randomCone/v3
echo "Done with PbPb DATA!  "
rm randomCones_TkpTCut0_ak3_PbPb_DATA.root
fi

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x LoopRandomCone_v5rev_akX_TYPE.C(0.3,5000,"HIJING","pA","/mnt/hadoop/cms/store/user/tuos/pPb/HiForest/5_3_3_patch3/pPbHijing_5_3_3_v0_HiForest2_v07.root","randomCones_TkpTCut0_ak3_pA_HIJING.root")
.q
EOF
cp randomCones_TkpTCut0_ak3_pA_HIJING.root /mnt/hadoop/cms/store/user/jrobles/PAanalysis/randomCone/v3
echo "Done with pA HIJING!  "
rm randomCones_TkpTCut0_ak3_pA_HIJING.root
fi

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x LoopRandomCone_v5rev_akX_TYPE.C(0.3,5000,"HYDJET","pA","/mnt/hadoop/cms/store/user/vzhukova/Hydjet_merged/mergedFile.root","randomCones_TkpTCut0_ak3_pA_HYDJET.root")
.q
EOF
cp randomCones_TkpTCut0_ak3_pA_HYDJET.root /mnt/hadoop/cms/store/user/jrobles/PAanalysis/randomCone/v3
echo "Done with pA HYDJET!  "
rm randomCones_TkpTCut0_ak3_pA_HYDJET.root
fi



run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x LoopRandomCone_v5rev_akX_TYPE.C(0.3,5000,"AMPT","pA","/mnt/hadoop/cms/store/user/frankmalocal/forest/HiForest_AMPT_Evening_5_3_3_v1.root","randomCones_TkpTCut0_ak3_pA_AMPT.root")
.q
EOF
cp randomCones_TkpTCut0_ak3_pA_AMPT.root /mnt/hadoop/cms/store/user/jrobles/PAanalysis/randomCone/v3
echo "Done with pA AMPT!  "
rm randomCones_TkpTCut0_ak3_pA_AMPT.root
fi
