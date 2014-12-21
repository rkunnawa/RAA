export SCRAM_ARCH=slc6_amd64_gcc472
#export SCRAM_ARCH=slc5_amd64_gcc462
#scram arch for slc6 is defined above
source /osg/app/cmssoft/cms/cmsset_default.sh
#source /apps/02/cmssoft/cms/cmsset_default.sh
#source /osg/app/glite/etc/profile.d/grid_env.sh

echo ""
echo "----------------------------------------------------"
echo "Job started on `date` at WN: `hostname` "
echo "Job is running on `uname -a`"

cd /net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/
#eval `scram list CMSSW`
eval `scramv1 runtime -sh`
cd /net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Macros/RAA/

echo "root directory: $ROOTSYS"

gcc --version

startfile=$1
endfile=$2
destination=$3
process=$4

echo "Processing..."

#root -b -q RAA_read_data_pbpb.C\+\($startfile,$endfile\)

root -b -q RAA_read_MinBias.C\+\($startfile,$endfile\)

#mv PbPb_data_ak3_Vs_20140820_$endfile.root ../../Output/.

#root -b -q RAA_fakecheck.C\+\($startfile,$endfile\)
#mv pbpb_ak3_Vs_fakejet_histos_$endfile.root ../../Output/.

#root -b -q RAA_duplicateEventsCheck.C\+\($startfile,$endfile\)
#mv pbpb_jet55or65_duplicate_events_loop_run_lumi_event_$endfile.txt ../../Output/.

echo "Done!"

echo "Copying output files to " $destination
