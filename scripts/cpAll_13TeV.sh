#!/bin/sh

source $(dirname `readlink -f $0`)/parallelTools.sh

# in the excludeList, put all distributions that you dont want. Separate them with a |, i.e. HypLLBarDPhi|HypNeutrinopT
# excludeList='HypNeutrinopT|bcp|_step|akr|bkr|SF|#'
# plotList=`awk '{print $1}' HistoList_control_13TeV | grep Hyp| grep -Erv $excludeList`

#echo
#echo "******************************************************************"
#echo "FIXME: this macro needs an update for 2016 recommendations, use at own risk... "
#echo "Creating control plots for signal variations ONLY"
#echo "If you want to create control plots for all/any other variation please run"
#echo "install/bin/Histo13TeV -t cp -s <systematic>"
#echo "******************************************************************"
#echo
#echo

#excludeList="^$|step1|step2|step3|step4|akr|bkr|SF|#|bcp"

excludeList="#"
plotList=`awk '{print $1}' HistoList_control_13TeV | grep -Ev $excludeList`

echo "Please press any key to start unfolding the following distributions in parallel or press Ctrl-C to cancel:"
echo "$plotList" | perl -l40 -pe ''
read -n 1 -s
echo ""

########################
## draw all systematic variations for control plots
## no unceratinty band will be plotted yet
########################
for i in $plotList; do
    $HISTO13TeV -t cp -s all -p +$i&
    w
done

########################
## draw Nominal control plot including uncertainty band
##   plese notice the ' -b'
########################
#
#echo "----------------------------------------------------------------"
#echo "Now sumbitting jobs to draw the control plot with the error band"
#for i in $plotList; do 
#    $HISTO13TeV -t cp -s Nominal -p +$i  -b&
#    w
#done

#Nominal JER_UP JER_DOWN
#
# JESAbsoluteStat_UP JESAbsoluteStat_DOWN JESAbsoluteScale_UP JESAbsoluteScale_DOWN JESAbsoluteMPFBias_UP JESAbsoluteMPFBias_DOWN JESFragmentation_UP JESFragmentation_DOWN 
# JESSinglePionECAL_UP JESSinglePionECAL_DOWN  
# JESSinglePionHCAL_UP JESSinglePionHCAL_DOWN 
# JESFlavorQCD_UP JESFlavorQCD_DOWN 
# JESTimePtEta_UP JESTimePtEta_DOWN 
# JESRelativeBal_UP JESRelativeBal_DOWN 
# JESRelativeJEREC1_UP JESRelativeJEREC1_DOWN
# JESRelativePtBB_UP JESRelativePtBB_DOWN 
# JESRelativePtEC1_UP JESRelativePtEC1_DOWN 
# JESRelativeFSR_UP JESRelativeFSR_DOWN 
# JESRelativeStatFSR_UP JESRelativeStatFSR_DOWN 
# JESRelativeStatEC_UP JESRelativeStatEC_DOWN 
# JESPileUpDataMC_UP JESPileUpDataMC_DOWN 
# JESPileUpPtRef_UP JESPileUpPtRef_DOWN 
# JESPileUpPtEC1_UP JESPileUpPtEC1_DOWN 
# JESPileUpPtBB_UP JESPileUpPtBB_DOWN 
#    
# UNCLUSTERED_UP UNCLUSTERED_DOWN PU_UP PU_DOWN 
# TRIG_UP TRIG_DOWN TRIG_ETA_UP TRIG_ETA_DOWN LEPT_UP LEPT_DOWN 
# DY_UP DY_DOWN BG_UP BG_DOWN KIN_UP KIN_DOWN 
# BTAG_UP BTAG_DOWN BTAG_LJET_UP BTAG_LJET_DOWN 
# BTAG_PT_UP BTAG_PT_DOWN BTAG_ETA_UP BTAG_ETA_DOWN 
# BTAG_LJET_PT_UP BTAG_LJET_PT_DOWN BTAG_LJET_ETA_UP BTAG_LJET_ETA_DOWN 
# MESCALE_UP MESCALE_DOWN MEFACSCALE_UP MEFACSCALE_DOWN MERENSCALE_UP MERENSCALE_DOWN 
# PSISRSCALE_UP PSISRSCALE_DOWN PSFSRSCALE_UP PSFSRSCALE_DOWN 
# BFRAG_UP BFRAG_DOWN BFRAG_PETERSON 
# BSEMILEP_UP BSEMILEP_DOWN 
# ERDON ERDONRETUNE GLUONMOVETUNE 
# UETUNE_UP UETUNE_DOWN 
# MATCH_UP MATCH_DOWN 
# MASS_UP MASS_DOWN 
# PDF_ALPHAS_UP PDF_ALPHAS_DOWN






