#!/bin/bash

INDIR=/hadoop/cms/store/group/snt/run2_moriond17/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/V08-00-16/
declare -a SAMPLES=(merged_ntuple_10)
MAXEVT=-1
OUTPUT=merged10x

for SAMPLE in ${SAMPLES[@]}; do
    eval 'nohup nice -n 10 ./PileupCorrection.exe ${OUTPUT} ${INDIR}/${SAMPLE}.root ${MAXEVT} >& ${SAMPLE}_log.txt &'
done
