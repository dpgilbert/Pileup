#!/bin/bash

INDIR=/hadoop/cms/store/group/snt/run2_moriond17/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/V08-00-16/
MAXEVT=-1

for i in {10..278}
do
    eval 'nice -n 10 ./PileupCorrection.exe output/merged_ntuple_${i} ${INDIR}/merged_ntuple_${i}.root ${MAXEVT}'
done
