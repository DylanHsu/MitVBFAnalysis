#!/bin/bash
function splitPath {
  local IFS=/
  local pieces
  # Disable glob expansion
  set -f
  pieces=( $@ ) 
  set +f
  #printf '%d\n' "${#pieces[@]}"
  #printf '%s\n' "${pieces[@]}"
  echo ${pieces[${#pieces[@]}-1]}
}

voms-proxy-init -voms cms
#cp -v /tmp/x509up_u$UID $WD/x509up

# name of the job
jobName=$1
# text file with list of file locations
listOfDataFiles=$2
# (electrons | muons)
flavor=$3
typeOne=$4 #(true|false)

startDir=`pwd`
receptacle=/data/t3home000/${USER}/vbfTriggerPlots
#receptacle=/mnt/hadoop/scratch/${USER}/vbfTriggerPlots

# Tar up necessary CMSSW areas
mkdir -p ${receptacle}/${jobName}
echo "Making tarball of CMSSW: ${receptacle}/CMSSW_${jobName}.tgz"
cd $CMSSW_BASE
tar  --exclude-vcs -chzf "${receptacle}/CMSSW_${jobName}.tgz" src/PandaAnalysis src/PandaCore src/PandaTree src/MitVBFAnalysis/macros src/.rootlogon.C python biglib bin lib objs test external

# Set up libs that the job will use to run
# Need the shell script and the CMSSW code
to_run="$CMSSW_BASE/src/MitVBFAnalysis/T3/condor-run ./batch_vbfTriggers.sh"
libs="${receptacle}/CMSSW_${jobName}.tgz $CMSSW_BASE/src/MitVBFAnalysis/T3/batch_vbfTriggers.sh"
jobNumber=0
cd $CMSSW_BASE/src

# Take the list of input files and chop it up into many smaller ones,
# then write the list of the smaller pieces to jobFileSets.txt
split -a 6 -l 4 $listOfDataFiles -d ${receptacle}/${jobName}/jobFiles_
ls ${receptacle}/${jobName}/jobFiles_* > ${receptacle}/${jobName}/jobFileSets.txt

# Loop over jobFileSets and submit a job for each set of N files
while read line
do
    cd ${receptacle}/${jobName}
    options=( $line )
    batchFile=${options[0]}
    batchFileBase=`splitPath ${batchFile}`
    outputFile="vbf_batchTree_triggerEff_${jobNumber}.root"
    echo "Setting up job with input files: ${libs}; output files ${outputFile}"
    args="--auxiliary-input ${libs} ${batchFile} --output ${outputFile}"
    args="${args} --task-name ${jobName}"
    args="${args} --job-args \"${flavor} ${batchFileBase} ${jobNumber} ${typeOne} ${jobName} ${CMSSW_VERSION}\" "
    #echo "$to_run $args"
    eval "$to_run $args"
    ((jobNumber++))
done < "${receptacle}/${jobName}/jobFileSets.txt"
cd ${startDir}
