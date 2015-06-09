#!/bin/bash

echo "running for PbPb R = 0.2, between eta +/- 2"
root -b -q 'RAA_plot_yetkinCutEfficiency.C+("20_eta_20",2,20,20)'

echo "running for PbPb R = 0.3, between eta +/- 2"
root -b -q 'RAA_plot_yetkinCutEfficiency.C+("20_eta_20",3,20,20)'

echo "running for PbPb R = 0.4, between eta +/- 2"
root -b -q 'RAA_plot_yetkinCutEfficiency.C+("20_eta_20",4,20,20)'
echo "running for PbPb R = 0.2, between eta +/- 1"
root -b -q 'RAA_plot_yetkinCutEfficiency.C+("10_eta_10",2,10,10)'

echo "running for PbPb R = 0.3, between eta +/- 1"
root -b -q 'RAA_plot_yetkinCutEfficiency.C+("10_eta_10",3,10,10)'

echo "running for PbPb R = 0.4, between eta +/- 1"
root -b -q 'RAA_plot_yetkinCutEfficiency.C+("10_eta_10",4,10,10)'

echo "running for PbPb R = 0.2, between eta +/- 1 to 1.8"
root -b -q 'RAA_plot_yetkinCutEfficiency.C+("10_eta_18",2,10,18)'

echo "running for PbPb R = 0.3, between eta +/- 1 to 1.8"
root -b -q 'RAA_plot_yetkinCutEfficiency.C+("10_eta_18",3,10,18)'

echo "running for PbPb R = 0.4, between eta +/- 1 to 1.8"
root -b -q 'RAA_plot_yetkinCutEfficiency.C+("10_eta_18",4,10,18)'


