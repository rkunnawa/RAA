#!/bin/bash

run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_analyze.C(3, 3,"Pu","PF",40,"0_absEta_05",1.0)
.q
EOF
echo "Done with pp R=0.3, eta Width 0 < |eta| < 0.5"
fi


run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_analyze.C(3, 3,"Pu","PF",40,"05_absEta_10",1.0)
.q
EOF
echo "Done with pp R=0.3, eta Width 0.5 < |eta| < 1.0"
fi


run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_analyze.C(3, 3,"Pu","PF",40,"10_absEta_15",1.0)
.q
EOF
echo "Done with pp R=0.3, eta Width 1.0 < |eta| < 1.5"
fi


run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_analyze.C(3, 3,"Pu","PF",40,"15_absEta_20",1.0)
.q
EOF
echo "Done with pp R=0.3, eta Width 1.5 < |eta| < 2.0"
fi


run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_analyze.C(3, 3,"Pu","PF",40,"20_absEta_25",1.0)
.q
EOF
echo "Done with pp R=0.3, eta Width 2.0 < |eta| < 2.5"
fi



run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_analyze.C(3, 5,"Pu","PF",40,"0_absEta_05",1.0)
.q
EOF
echo "Done with pp R=0.5, eta Width 0 < |eta| < 0.5"
fi


run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_analyze.C(3, 5,"Pu","PF",40,"05_absEta_10",1.0)
.q
EOF
echo "Done with pp R=0.5, eta Width 0.5 < |eta| < 1.0"
fi


run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_analyze.C(3, 5,"Pu","PF",40,"10_absEta_15",1.0)
.q
EOF
echo "Done with pp R=0.5, eta Width 1.0 < |eta| < 1.5"
fi


run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_analyze.C(3, 5,"Pu","PF",40,"15_absEta_20",1.0)
.q
EOF
echo "Done with pp R=0.5, eta Width 1.5 < |eta| < 2.0"
fi


run="yes"
echo $run
if [ $run = "yes" ]; then
root -l <<EOF
.x RAA_analyze.C(3, 5,"Pu","PF",40,"20_absEta_25",1.0)
.q
EOF
echo "Done with pp R=0.5, eta Width 2.0 < |eta| < 2.5"
fi
