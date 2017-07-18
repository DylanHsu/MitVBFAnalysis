# MitVBFAnalysis

## Installation
Requires `CMSSW_8_0_26_patch1` or higher. Install in `CMSSW_*/src` along with Panda dependencies.

ACLiC compiling scripts require you to add the following lines to `CMSSW_*/src/.rootlogon.C`:

`gSystem->Load("libPandaTreeObjects.so");`

`gSystem->Load("libPandaCoreTools.so");`

## Trigger measurement
We measure the efficiency of MET triggers in a leptonic W selection. The selection is performed on Panda data files by `macros/triggerEff.C` which compiles via ACLiC and can be submitted to MIT T3 and flocked to T2. Panda files are listed in `data`. From your `CMSSW_*/src` directory, call:

`MitVBFAnalysis/T3/submit_vbfTriggers.sh (name of your job) (file containing Panda addresses) (electrons|muons) (want type1 MET triggers? true|false)`

When your jobs are finished, use `hadd -ff` to concatenate the output, then check out the results using the macros in `plotting`.

