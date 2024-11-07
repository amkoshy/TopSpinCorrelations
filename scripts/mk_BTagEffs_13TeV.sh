#!/bin/sh

#######################################################
###  This scripts creates the BTagging efficiency 
###  histograms for each systematic variation
#######################################################

source $(dirname `readlink -f $0`)/parallelTools.sh

for syst in Nominal \
            JER_UP JER_DOWN \
            JES_UP JES_DOWN \
            PU_UP PU_DOWN \
            TRIG_UP TRIG_DOWN \
            TRIG_ETA_UP TRIG_ETA_DOWN \
            LEPT_UP LEPT_DOWN \
            UNCLUSTERED_UP UNCLUSTERED_DOWN \
            KIN_UP KIN_DOWN \
            BTAG_UP BTAG_DOWN BTAG_LJET_UP BTAG_LJET_DOWN \
            BTAG_PT_UP BTAG_PT_DOWN BTAG_ETA_UP BTAG_ETA_DOWN \
            BTAG_LJET_PT_UP BTAG_LJET_PT_DOWN BTAG_LJET_ETA_UP BTAG_LJET_ETA_DOWN; do
    for c in ee emu mumu; do
        w
        $LA -f ttbarsignalplustau.root -c $c -s $syst&
    done
done

for syst in JESAbsoluteStat_UP JESAbsoluteStat_DOWN \
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
            JESPileUpPtBB_UP JESPileUpPtBB_DOWN; do
    for c in ee emu mumu; do
        w
        $LA -f ttbarsignalplustau.root -c $c -s $syst&
    done
done

#######################################################
###  Signal variations
###  PDF systematic will use the Nominal btag effciencies
#######################################################
for c in ee emu mumu; do
    w
    $LA -f ttbarsignalplustau.root -s PDF_ALPHAS_UP -c $c &
    $LA -f ttbarsignalplustau.root -s PDF_ALPHAS_DOWN -c $c &
    $LA -f ttbarsignalplustau.root -s MESCALE_UP -c $c &
    $LA -f ttbarsignalplustau.root -s MESCALE_DOWN -c $c &
    $LA -f ttbarsignalplustau.root -s MEFACSCALE_UP -c $c &
    $LA -f ttbarsignalplustau.root -s MEFACSCALE_DOWN -c $c &
    $LA -f ttbarsignalplustau.root -s MERENSCALE_UP -c $c &
    $LA -f ttbarsignalplustau.root -s MERENSCALE_DOWN -c $c &
    $LA -f ttbarsignalplustau.root -s BFRAG_UP -c $c &
    $LA -f ttbarsignalplustau.root -s BFRAG_DOWN -c $c &
    $LA -f ttbarsignalplustau.root -s BFRAG_CENTRAL -c $c &
    $LA -f ttbarsignalplustau.root -s BFRAG_PETERSON -c $c &
    $LA -f ttbarsignalplustau.root -s BSEMILEP_UP -c $c &
    $LA -f ttbarsignalplustau.root -s BSEMILEP_DOWN -c $c &
    $LA -f tau_psfsrscaleup.root -c $c &
    $LA -f tau_psfsrscaledown.root -c $c &
    $LA -f tau_psisrscaleup.root -c $c &
    $LA -f tau_psisrscaledown.root -c $c &
    $LA -f tau_175_massup.root -c $c &
    $LA -f tau_169_massdown.root -c $c &
    $LA -f tau_uetuneup.root -c $c &
    $LA -f tau_uetunedown.root -c $c &
    $LA -f tau_matchup.root -c $c &
    $LA -f tau_matchdown.root -c $c &
    $LA -f tau_powhegv2Herwig.root -c $c &
    $LA -f tau_fromDilepton_amcatnlofxfx.root -c $c &
    #$LA -f tau_madgraphmlm.root -c $c &
    $LA -f tau_erdon.root -c $c &
    $LA -f tau_erdonretune.root -c $c &
    $LA -f tau_gluonmovetune.root -c $c &
done



wait

if [ "$isNAF" = 1 ]; then
    echo "Please check your jobs with qstat -u $USER | grep load_Analysis"
else
    echo "Processing all nominal samples finished!"
fi


