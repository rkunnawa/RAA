#!/bin/bash



run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_randomcone.C+(2,"PF","Pu","data")
.q
EOF
echo "Done with PbPb data for akPu2PF  "
fi

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_randomcone.C+(3,"PF","Pu","data")
.q
EOF
echo "Done with PbPb data for akPu3PF  "
fi

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_randomcone.C+(4,"PF","Pu","data")
.q
EOF
echo "Done with PbPb data for akPu4PF  "
fi
