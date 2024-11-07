#!/bin/sh

###########################################################################
###  This script creates base input for pseudo-data using ttbar related
###  selectionRoot files for use in plotterclass.cc (13 TeV)
###########################################################################

FILE=_run2015P
TTSig=_ttbarsignalplustau
TTBkgd=_ttbarbgviatau


for syst in Nominal JES_DOWN JES_UP JER_DOWN JER_UP MATCH_UP MATCH_DOWN MASS_UP MASS_DOWN PU_UP PU_DOWN TRIG_UP TRIG_DOWN SCALE_UP SCALE_DOWN MCATNLO POWHEG POWHEGHERWIG SPINCORR \
            LEPT_UP LEPT_DOWN \
            KIN_UP KIN_DOWN \
            BTAG_UP BTAG_DOWN BTAG_LJET_UP BTAG_LJET_DOWN \
            BTAG_PT_UP BTAG_PT_DOWN BTAG_ETA_UP BTAG_ETA_DOWN \
            BTAG_LJET_PT_UP BTAG_LJET_PT_DOWN BTAG_LJET_ETA_UP BTAG_LJET_ETA_DOWN \
            BTAG_BEFF_UP BTAG_BEFF_DOWN BTAG_CEFF_UP BTAG_CEFF_DOWN BTAG_LEFF_UP BTAG_LEFF_DOWN; do
    for c in ee emu mumu; do
        DIR=selectionRoot/$syst/$c/$c
        rm $DIR$FILE.root
        hadd $DIR$FILE.root $DIR$TTSig*.root $DIR$TTBkgd*.root
        echo " "
    done
    wait
done
wait
echo "...Done!"


