# RecoBTag-BTagAnalyzerLite

## Software setup

```
cmsrel CMSSW_7_4_7
cd CMSSW_7_4_7/src
cmsenv

setenv CMSSW_GIT_REFERENCE /cvmfs/cms.cern.ch/cmssw.git.daily
git cms-init

git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git
git fetch --tags btv-cmssw

git cms-merge-topic -u cms-btv-pog:FixSoftElectronTagger-v1_from-CMSSW_7_4_1
git cms-merge-topic -u cms-btv-pog:CSVv2InCombinedMVA-v1_from-CMSSW_7_4_5

git clone -b V00-00-01 git://github.com/cms-btv-pog/cms-EventCounter.git MyAnalysis/EventCounter
git clone -b 7_4_X git@github.com:cms-btv-pog/RecoBTag-BTagAnalyzerLite.git RecoBTag/BTagAnalyzerLite

scram b -j8

cd RecoBTag/BTagAnalyzerLite/test/

cmsRun runBTagAnalyzerLite_cfg.py miniAOD=True maxEvents=100 reportEvery=1 wantSummary=True
```