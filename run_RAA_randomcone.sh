#!/bin/bash



run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_randomcone.C(3,"PF","Vs","MC")
.q
EOF
echo "Done with PbPb MC for akVs3PF  "
fi

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_randomcone.C(4,"PF","Vs","MC")
.q
EOF
echo "Done with PbPb MC for akVs4PF  "
fi

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_randomcone.C(5,"PF","Vs","MC")
.q
EOF
echo "Done with PbPb MC for akVs5PF  "
fi

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_randomcone.C(3,"PF","Vs","data")
.q
EOF
echo "Done with PbPb data for akVs3PF  "
fi

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_randomcone.C(4,"PF","Vs","data")
.q
EOF
echo "Done with PbPb data for akVs4PF  "
fi

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_randomcone.C(5,"PF","Vs","data")
.q
EOF
echo "Done with PbPb data for akVs5PF  "
fi
