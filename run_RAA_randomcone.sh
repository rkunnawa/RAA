#!/bin/bash

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_randomcone.C(3,"Calo","Pu","MC")
.q
EOF
echo "Done with PbPb MC for akPu3Calo  "
fi

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_randomcone.C(4,"Calo","Pu","MC")
.q
EOF
echo "Done with PbPb MC for akPu4Calo  "
fi

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_randomcone.C(5,"Calo","Pu","MC")
.q
EOF
echo "Done with PbPb MC for akPu5Calo  "
fi
