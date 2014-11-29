RecoBTag-BTagAnalyzerLite
=========================
cmsrel CMSSW_7_3_0_pre3
cd CMSSW_7_3_0_pre3/src
cmsenv

git clone -b 7_3_X git@github.com:cms-btv-pog/RecoBTag-BTagAnalyzerLite.git RecoBTag/BTagAnalyzerLite

scram b -j5

cd RecoBTag/BTagAnalyzerLite/test/

cmsRun runBTagAnalyzerLite_cfg.py maxEvents=100 reportEvery=1 wantSummary=True
