#!/bin/zsh
foreach channel (ee emu mumu)

  foreach Syst  (JES_UP JES_DOWN JER_UP JER_DOWN \
                 #JESAbsoluteStat_UP JESAbsoluteStat_DOWN \
                 #JESAbsoluteScale_UP JESAbsoluteScale_DOWN \
                 #JESAbsoluteFlavMap_UP JESAbsoluteFlavMap_DOWN \
                 #JESAbsoluteMPFBias_UP JESAbsoluteMPFBias_DOWN \
                 #JESFragmentation_UP JESFragmentation_DOWN \
                 #JESSinglePionECAL_UP JESSinglePionECAL_DOWN \
                 #JESSinglePionHCAL_UP JESSinglePionHCAL_DOWN \
                 #JESFlavorQCD_UP JESFlavorQCD_DOWN \
                 #JESTimePtEta_UP JESTimePtEta_DOWN \
                 #JESRelativeJEREC1_UP JESRelativeJEREC1_DOWN \
                 #JESRelativeJEREC2_UP JESRelativeJEREC2_DOWN \
                 #JESRelativeJERHF_UP JESRelativeJERHF_DOWN \
                 #JESRelativePtBB_UP JESRelativePtBB_DOWN \
                 #JESRelativePtEC1_UP JESRelativePtEC1_DOWN \
                 #JESRelativePtEC2_UP JESRelativePtEC2_DOWN \
                 #JESRelativePtHF_UP JESRelativePtHF_DOWN \
                 #JESRelativeFSR_UP JESRelativeFSR_DOWN \
                 #JESRelativeStatFSR_UP JESRelativeStatFSR_DOWN \
                 #JESRelativeStatEC_UP JESRelativeStatEC_DOWN \
                 #JESRelativeStatHF_UP JESRelativeStatHF_DOWN \
                 #JESPileUpDataMC_UP JESPileUpDataMC_DOWN \
                 #JESPileUpPtRef_UP JESPileUpPtRef_DOWN \
                 #JESPileUpPtEC1_UP JESPileUpPtEC1_DOWN \
                 #JESPileUpPtEC2_UP JESPileUpPtEC2_DOWN \
                 #JESPileUpPtHF_UP JESPileUpPtHF_DOWN \
                 #JESPileUpPtBB_UP JESPileUpPtBB_DOWN \
                 PU_UP PU_DOWN TRIG_UP TRIG_DOWN TRIG_ETA_UP TRIG_ETA_DOWN \
                 LEPT_UP LEPT_DOWN \
                 UNCLUSTERED_UP UNCLUSTERED_DOWN \
                 KIN_UP KIN_DOWN \
                 PDF_ALPHAS_UP PDF_ALPHAS_DOWN \
                 MESCALE_UP MESCALE_DOWN MEFACSCALE_UP MEFACSCALE_DOWN MERENSCALE_UP MERENSCALE_DOWN PSISRSCALE_UP PSISRSCALE_DOWN PSFSRSCALE_UP PSFSRSCALE_DOWN UETUNE_UP UETUNE_DOWN MASS_UP MASS_DOWN \
                 BFRAG_UP BFRAG_DOWN BFRAG_CENTRAL BFRAG_PETERSON BSEMILEP_UP BSEMILEP_DOWN ERDON ERDONRETUNE GLUONMOVETUNE \ 
                 MATCH_UP MATCH_DOWN \
                 POWHEG POWHEGHERWIG POWHEGV2HERWIG MCATNLO \
                 BTAG_UP BTAG_DOWN BTAG_LJET_UP BTAG_LJET_DOWN \
                 BTAG_PT_UP BTAG_PT_DOWN BTAG_ETA_UP BTAG_ETA_DOWN \
                 BTAG_LJET_PT_UP BTAG_LJET_PT_DOWN BTAG_LJET_ETA_UP BTAG_LJET_ETA_DOWN
)
    ./install/bin/loadTUnfoldHisto ${Syst} ${channel}
  end

end
