#!/bin/bash

cd ../src/output

hadd -f hists.root merged*

cd -

mkdir -p pngs
mkdir -p pdfs

python plotter.py