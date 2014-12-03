RecoBTag-BTagAnalyzerLite
=========================
cmsrel CMSSW_5_3_20
cd CMSSW_5_3_20/src
cmsenv

git cms-merge-topic cms-btv-pog:Nsubjettiness_Qjets-backport-V02_from-CMSSW_5_3_13
git cms-merge-topic cms-btv-pog:IVFVertexNTracksFix-V01_from-CMSSW_5_3_13
git cms-merge-topic cms-btv-pog:SVClustering-V01_from-CMSSW_5_3_16
git cms-merge-topic cms-btv-pog:KtPruning_BDRSFiltering-V01_from-CMSSW_5_3_16
git cms-merge-topic cms-btv-pog:CSVV2-V02_from-CMSSW_5_3_20

git clone -b V00-00-01 git://github.com/cms-btv-pog/cms-EventCounter.git MyAnalysis/EventCounter
git clone -b 5_3_X_v1.01 git@github.com:cms-btv-pog/RecoBTag-BTagAnalyzerLite.git RecoBTag/BTagAnalyzerLite

scram b -j5

cd RecoBTag/BTagAnalyzerLite/test/

cmsRun runBTagAnalyzerLite_cfg.py maxEvents=100 reportEvery=1 wantSummary=True
