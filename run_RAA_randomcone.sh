#!/bin/bash

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_randomcone.C(2,"PF","Pu","MC")
.q
EOF
echo "Done with PbPb MC for akPu3PF  "
fi

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_randomcone.C(3,"PF","Pu","MC")
.q
EOF
echo "Done with PbPb MC for akPu4PF  "
fi

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_randomcone.C(4,"PF","Pu","MC")
.q
EOF
echo "Done with PbPb MC for akPu5PF  "
fi
