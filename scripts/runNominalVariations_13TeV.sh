#!/bin/sh

source $(dirname `readlink -f $0`)/parallelTools.sh

for sys in JES_UP JES_DOWN JER_UP JER_DOWN \
           JESAbsoluteStat_UP JESAbsoluteStat_DOWN \
           JESAbsoluteScale_UP JESAbsoluteScale_DOWN \
           JESAbsoluteFlavMap_UP JESAbsoluteFlavMap_DOWN \
           JESAbsoluteMPFBias_UP JESAbsoluteMPFBias_DOWN \
           JESFragmentation_UP JESFragmentation_DOWN \
           JESSinglePionECAL_UP JESSinglePionECAL_DOWN \
           JESSinglePionHCAL_UP JESSinglePionHCAL_DOWN \
           JESFlavorQCD_UP JESFlavorQCD_DOWN \
           JESTimePtEta_UP JESTimePtEta_DOWN \
           JESRelativeBal_UP JESRelativeBal_DOWN \
           JESRelativeJEREC1_UP JESRelativeJEREC1_DOWN \
           JESRelativeJEREC2_UP JESRelativeJEREC2_DOWN \
           JESRelativeJERHF_UP JESRelativeJERHF_DOWN \
           JESRelativePtBB_UP JESRelativePtBB_DOWN \
           JESRelativePtEC1_UP JESRelativePtEC1_DOWN \
           JESRelativePtEC2_UP JESRelativePtEC2_DOWN \
           JESRelativePtHF_UP JESRelativePtHF_DOWN \
           JESRelativeFSR_UP JESRelativeFSR_DOWN \
           JESRelativeStatFSR_UP JESRelativeStatFSR_DOWN \
           JESRelativeStatEC_UP JESRelativeStatEC_DOWN \
           JESRelativeStatHF_UP JESRelativeStatHF_DOWN \
           JESPileUpDataMC_UP JESPileUpDataMC_DOWN \
           JESPileUpPtRef_UP JESPileUpPtRef_DOWN \
           JESPileUpPtEC1_UP JESPileUpPtEC1_DOWN \
           JESPileUpPtEC2_UP JESPileUpPtEC2_DOWN \
           JESPileUpPtHF_UP JESPileUpPtHF_DOWN \
           JESPileUpPtBB_UP JESPileUpPtBB_DOWN \
           PU_UP PU_DOWN \
           TRIG_UP TRIG_DOWN \
           TRIG_ETA_UP TRIG_ETA_DOWN \
           LEPT_UP LEPT_DOWN \
           UNCLUSTERED_UP UNCLUSTERED_DOWN \
           KIN_UP KIN_DOWN \
           BTAG_UP BTAG_DOWN BTAG_LJET_UP BTAG_LJET_DOWN \
           BTAG_PT_UP BTAG_PT_DOWN BTAG_ETA_UP BTAG_ETA_DOWN \
           BTAG_LJET_PT_UP BTAG_LJET_PT_DOWN BTAG_LJET_ETA_UP BTAG_LJET_ETA_DOWN; do # \
#            BTAG_BEFF_UP BTAG_BEFF_DOWN BTAG_CEFF_UP BTAG_CEFF_DOWN BTAG_LEFF_UP BTAG_LEFF_DOWN; do

    for c in ee emu mumu; do
        w
        $LA -f ttbarsignalplustau.root -c $c -s $sys&
        $LA -f ttbarsignalplustau.root -c $c -s $sys --bgviatau&
        $LA -f dy -d 11 -c $c -s $sys&
        $LA -f dy -d 13 -c $c -s $sys&
        $LA -f dy -d 15 -c $c -s $sys&
        $LA -f ttbarZ -c $c -s $sys&
    done

    for i in qcd tw.root ttbarbg.root wtol ww wz zz ttbarW ttgjets; do
        w
        $LA -f $i -s $sys&
    done

done
wait

