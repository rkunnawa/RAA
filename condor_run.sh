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

df -h

gcc --version

startfile=$1
endfile=$2
destination=$3
process=$4

echo "Processing..."

#root -b -q RAA_read_data_pbpb.C\+\($startfile,$endfile\)
#mv /net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_data_vz_cent_akPuPF_20150324_$endfile.root /mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/data/.

#root -b -q RAA_fragBiasCheck.C\+\($startfile,$endfile\)
#root -b -q RAA_fragBiasCheck_pp_data.C\+\($startfile,$endfile\)

#root -b -q RAA_read_MinBias.C\+\($startfile,$endfile\)
#mv /net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_MinBiasUPC_ntuple_SuperNovaRejected_akPuCalo_20150311_1.root /mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/.

#root -b -q RAA_read_mc_pp.C\+\($startfile,$endfile\)

#root -b -q RAA_read_mc.C\+\($startfile,$endfile\)

#root -b -q RAA_read_data_pp.C\+\($startfile,$endfile\)
#mv /net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/pp_data_spectra_trgObj_chMaxjtpt0p02_eMaxjtpt0p6_ak35PF_20150323_$endfile.root /mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/data/.

#root -b -q RAA_calo_pf_JetCorrelation_v2.C\+\($startfile,$endfile\)
#mv /net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/PbPb_MC_calo_pf_jet_correlation_mcclosure_histograms_deltaR_0p2_akPu2_20150407_$endfile.root /mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PbPb/2011/mc/.

root -b -q RAA_calo_pf_JetCorrelation_pp.C\+\($startfile,$endfile\)
mv /net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/pp_Data_calo_pf_jet_correlation_deltaR_0p2_ak3_20150409_$endfile.root /mnt/hadoop/cms/store/user/rkunnawa/rootfiles/PP/2013/.

#mv PbPb_data_ak3_Vs_20140820_$endfile.root ../../Output/.

#root -b -q RAA_fakecheck.C\+\($startfile,$endfile\)
#mv pbpb_ak3_Vs_fakejet_histos_$endfile.root ../../Output/.

#root -b -q RAA_duplicateEventsCheck.C\+\($startfile,$endfile\)

echo "Done!"

echo "Copying output files to " $destination
