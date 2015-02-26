RecoBTag-BTagAnalyzerLite
=========================
cmsrel CMSSW_7_4_0_pre7
cd CMSSW_7_4_0_pre7/src
cmsenv

git cms-merge-topic -u cms-btv-pog:PATBTaggingUpdates_from-CMSSW_7_4_0_pre7

git clone -b V00-00-01 git://github.com/cms-btv-pog/cms-EventCounter.git MyAnalysis/EventCounter
git clone -b 7_4_X_dev git@github.com:cms-btv-pog/RecoBTag-BTagAnalyzerLite.git RecoBTag/BTagAnalyzerLite

scram b -j8

cd RecoBTag/BTagAnalyzerLite/test/

cmsRun runBTagAnalyzerLite_cfg.py maxEvents=100 reportEvery=1 wantSummary=True
