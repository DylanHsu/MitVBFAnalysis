source /cvmfs/cms.cern.ch/cmsset_default.sh
scram project CMSSW ${6}
pwd
basedir=`pwd`
#cmsrel ${6}
mkdir -p ${6}
mv CMSSW_${5}.tgz ${6}
cd ${6}
tar -zxf "CMSSW_${5}.tgz"
cd src
eval `scram runtime -sh`
mv ${basedir}/jobFiles_* .
rm  MitVBFAnalysis/macros/*.so
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:.:MitVBFAnalysis/macros
root -b -l -q MitVBFAnalysis/macros/triggerEff.C+\(\"$1\",\"$2\",\"$3\",$4\)
mv vbf_batchTree_triggerEff_${3}.root ${basedir}
