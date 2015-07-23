#!/bin/sh

counter=0
incrementer=1

destination=/net/hisrv0001/home/rkunnawa/WORK/RAA/CMSSW_5_3_20/src/Output/
filelist=jetRAA_pp_mc_forest.txt

#filelist=jetRAA_pp_data_forest.txt
#filelist=jetRAA_MinBiasUPC_forest.txt
#filelist=14010_MinBiasUPC_forest.txt
#filelist=jetRAA_PbPb_data_forest.txt
#filelist=jetRAA_PbPb_mc_forest.txt
#filelist=PbPb_HydjetMinBias_forest.txt
#filelist=PbPb_MinBiasUPC_forest.txt
#filelist2=jet80_filelist.txt
#filelist=wxie_MinBiasUPC_all.txtfor 
#filelist=pPbMCBForestList.txt

#mkdir -p $destination

nFiles=`wc -l < $filelist`

#export X509_USER_PROXY=/net/hisrv0001/home/rkunnawa/myproxy/

echo "nFiles in list: $nFiles"
while [ $counter -lt $1 ]
   do
	  echo $counter >> Submitted

	  Error=`echo "err/$counter" | sed "s/root/err/g"`
	  Output=`echo "out/$counter" | sed "s/root/out/g"`
	  Log=`echo "log/$counter" | sed "s/root/log/g"`        

          startfile=$(( $counter * $2 ))
          endfile=$(( ($counter + 1) * $2 ))
        if [ $endfile -gt $nFiles ]; then
            let endfile=$nFiles
            let counter=$1
        fi
	# Condor submit file
	cat > subfile <<EOF
Universe       = vanilla
# files will be copied back to this dir
Initialdir     = .
#tell condor where my grid certificate it
#x509userproxy=/tmp/x509up_u2142
# run my script
Executable     = condor_run.sh
+AccountingGroup = "group_cmshi.rkunnawa"
#+IsMadgraph = 1
Arguments      = $startfile $endfile $destination \$(Process)
# input files. in this case, there are none.
Input          = /dev/null
# log files
Error          = $Error
Output         = $Output
Log            = $Log
# get the environment (path, etc.)
Getenv         = True
# prefer to run on fast computers
Rank           = kflops
# only run on 64 bit computers
Requirements   = Arch == "X86_64"
# should write all output & logs to a local directory
# and then transfer it back to Initialdir on completion
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
# specify any extra input files (for example, an orcarc file)
transfer_input_files = $filelist

Queue
EOF

	# submit the job
	echo "submitting condor_run.sh $startfile $endfile to $destination ..." 
	condor_submit subfile
	counter=$(($counter + 1))
done
