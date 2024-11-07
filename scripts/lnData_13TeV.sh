#!/bin/zsh

BASEDIR=`pwd`

cd selectionRoot
echo "Creating links for alternative ttbar MC samples ... " 
ln -s AMCATNLOFXFX POWHEG
ln -s POWHEGV2HERWIG MCATNLO
ln -s MADGRAPHMLM POWHEGHERWIG
cd $BASEDIR

foreach channel (ee emu mumu)
   
   foreach syst (JES_DOWN JES_UP JER_DOWN JER_UP MATCH_UP MATCH_DOWN MASS_UP MASS_DOWN PU_UP PU_DOWN TRIG_UP TRIG_DOWN TRIG_ETA_UP TRIG_ETA_DOWN UNCLUSTERED_UP UNCLUSTERED_DOWN \
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
                 MESCALE_UP MESCALE_DOWN MEFACSCALE_UP MEFACSCALE_DOWN MERENSCALE_UP MERENSCALE_DOWN PDF_ALPHAS_UP PDF_ALPHAS_DOWN \
                 BFRAG_UP BFRAG_DOWN BFRAG_CENTRAL BFRAG_PETERSON BSEMILEP_UP BSEMILEP_DOWN \
                 ERDON ERDONRETUNE GLUONMOVETUNE \
                 PSISRSCALE_UP PSISRSCALE_DOWN PSFSRSCALE_UP PSFSRSCALE_DOWN UETUNE_UP UETUNE_DOWN MCATNLO POWHEG POWHEGHERWIG SPINCORR \
                 LEPT_UP LEPT_DOWN \
                 KIN_UP KIN_DOWN \
                 BTAG_UP BTAG_DOWN BTAG_LJET_UP BTAG_LJET_DOWN \
                 BTAG_PT_UP BTAG_PT_DOWN BTAG_ETA_UP BTAG_ETA_DOWN \
                 BTAG_LJET_PT_UP BTAG_LJET_PT_DOWN BTAG_LJET_ETA_UP BTAG_LJET_ETA_DOWN \
                 BTAG_BEFF_UP BTAG_BEFF_DOWN BTAG_CEFF_UP BTAG_CEFF_DOWN BTAG_LEFF_UP BTAG_LEFF_DOWN
        )
        if [ -d selectionRoot/$syst/$channel ] ; then
            cd selectionRoot/$syst/$channel/
            echo
            echo "Creating Links in ... " 
            pwd
            echo 
            ln -s ../../Nominal/$channel/*.root .
            cd $BASEDIR
        fi
	if [ -d ddaInput/$syst/$channel ] ; then
            cd ddaInput/$syst/$channel/
            echo
            echo "Creating Links in ... "
            pwd
            echo 
            ln -s ../../Nominal/$channel/*.root .
            cd $BASEDIR
        fi
   end

end

#Don't add here the systematics 'name' corresponding to the absence of a tag in the "ttbarsignalplustau.root" root file name
#(eg, MESCALE (like JES_UP) is created via the running over the nominal sample "ttbarsignalplustau.root", while MASS_DOWN is using the "ttbarsignalplustau_massdown.root") 
foreach SignalSyst (MATCH_UP MATCH_DOWN MASS_UP MASS_DOWN PSISRSCALE_UP PSISRSCALE_DOWN PSFSRSCALE_UP PSFSRSCALE_DOWN UETUNE_UP UETUNE_DOWN ERDON ERDONRETUNE GLUONMOVETUNE MCATNLO POWHEG POWHEGHERWIG SPINCORR)

    if [ -d selectionRoot/$SignalSyst ] ; then
        rm selectionRoot/$SignalSyst/*/*_ttbarbg.root
        rm selectionRoot/$SignalSyst/*/*_ttbarbgviatau.root
        rm selectionRoot/$SignalSyst/*/*_ttbarsignalplustau.root
    fi
end

#foreach SignalSyst (MASS_UP MASS_DOWN)
#
#    if [ -d selectionRoot/$SignalSyst ] ; then
#        rm selectionRoot/$SignalSyst/*/*_tw.root
#    fi
#end

# mkdir selectionRoot/HADUP selectionRoot/HADDOWN
# 
# foreach channel (ee emu mumu)
# 
#    mkdir selectionRoot/HADUP/$channel
#    mkdir selectionRoot/HADDOWN/$channel
#    
#    cp selectionRoot/MCATNLO/$channel/* selectionRoot/HADUP/$channel
#    cp selectionRoot/POWHEG/$channel/* selectionRoot/HADDOWN/$channel
# 
#    cd selectionRoot/HADUP/$channel/
#    echo
#    echo "Creating Links in ... " 
#    pwd
#    echo 
#    ln -s ../../Nominal/$channel/* .
#    cd $BASEDIR
# 
#    cd selectionRoot/HADDOWN/$channel/
#    echo
#    echo "Creating Links in ... " 
#    pwd
#    echo 
#    ln -s ../../Nominal/$channel/* .
#    cd $BASEDIR
# 
#    rm selectionRoot/HADUP/$channel/ttbarbg.root
#    rm selectionRoot/HADUP/$channel/ttbarsignalplustau.root
#    rm selectionRoot/HADDOWN/$channel/ttbarbg.root
#    rm selectionRoot/HADDOWN/$channel/ttbarsignalplustau.root
# 
# end

