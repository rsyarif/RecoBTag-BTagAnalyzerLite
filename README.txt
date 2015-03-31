RecoBTag-BTagAnalyzerLite - Implementing Shower Deconstruction (http://arxiv.org/abs/1211.3140, http://arxiv.org/abs/1102.3480)
=========================
cmsrel CMSSW_7_4_1
cd CMSSW_7_4_1/src
cmsenv
git cms-init
git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git
git fetch btv-cmssw
git cms-merge-topic -u cms-btv-pog:BoostedDoubleSVTagger-WithWeightFiles-v2_from-CMSSW_7_4_1

git clone -b V00-00-01 git://github.com/cms-btv-pog/cms-EventCounter.git MyAnalysis/EventCounter
git clone -b 7_4_X_v2.06 git@github.com:cms-btv-pog/RecoBTag-BTagAnalyzerLite.git RecoBTag/BTagAnalyzerLite

copy SD files ( ShowerDeconstruction/ ) in  working area CMSSW_X/src

make sure folder data/ exist under ShowerDeconstruction/
and has the file input_data.dat in data/ for SD configurations

scram b -j8

cd RecoBTag/BTagAnalyzerLite/test/

cmsRun runBTagAnalyzerLite_cfg.py miniAOD=True maxEvents=100 reportEvery=1 wantSummary=True
