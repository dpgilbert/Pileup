#!/bin/bash

cd ../src/output_RecoNpv

hadd -f hists.root merged*

cd -

mkdir -p RecoNpv
mkdir -p RecoNpv/pngs
mkdir -p RecoNpv/pdfs

python plotter_RecoNpv.py