#!/bin/sh

###########################################################################
###  This script creates base input for pseudo-data using ttbar related
###  selectionRoot files for use in plotterclass.cc (13 TeV)
###########################################################################

FILE=_run2015C
TTSig=_ttbarsignalplustau
TTBkgd=_ttbarbgviatau

for c in ee emu mumu; do
    DIR=selectionRoot/Nominal/$c/$c
    rm $DIR$FILE.root
    hadd $DIR$FILE.root $DIR$TTSig.root $DIR$TTBkgd.root
    echo " "
done

wait
echo "...Done!"


