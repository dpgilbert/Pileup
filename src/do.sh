#!/bin/bash

INDIR=/hadoop/cms/store/group/snt/run2_moriond17/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/V08-00-16/
MAXEVT=-1

#breaks sometimes
#for i in {10..278}
#breaks sometimes

#for i in {10..60}
#for i in {61..110}
#for i in {111..160}
#for i in {161..210}
for i in {210..278}
do
    eval 'nohup nice -n 10 ./PileupCorrection.exe output/merged_ntuple_${i} ${INDIR}/merged_ntuple_${i}.root ${MAXEVT} >& logs/${i}_log.txt &'
done
