#!/bin/zsh

mkdir -p TUnfoldHistoFileLists
rm TUnfoldHistoFileLists/TUnfoldHistoFileList*

foreach sample (run smu_run se_run dyee dymumu dytautau ww wz zz wtolnu ttbarW ttbarZ single ttbarbg ttbarsignal dy1050 dy50inf)

   foreach channel (ee emu mumu combined)

      foreach Syst  (
		     JES_UP JES_DOWN JER_UP JER_DOWN \
                JESAbsoluteStat_UP JESAbsoluteStat_DOWN \
                JESAbsoluteScale_UP JESAbsoluteScale_DOWN \
                JESAbsoluteFlavMap_UP JESAbsoluteFlavMap_DOWN \
                JESAbsoluteMPFBias_UP JESAbsoluteMPFBias_DOWN \
                JESFragmentation_UP JESFragmentation_DOWN \
                JESSinglePionECAL_UP JESSinglePionECAL_DOWN \
                JESSinglePionHCAL_UP JESSinglePionHCAL_DOWN \
                JESFlavorQCD_UP JESFlavorQCD_DOWN \
                JESTimePtEta_UP JESTimePtEta_DOWN \
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
                PDF_0_CENTRAL \
                PDF_10_DOWN PDF_10_UP \
                PDF_11_DOWN PDF_11_UP \
                PDF_12_DOWN PDF_12_UP \
                PDF_13_DOWN PDF_13_UP \
                PDF_14_DOWN PDF_14_UP \
                PDF_15_DOWN PDF_15_UP \
                PDF_16_DOWN PDF_16_UP \
                PDF_17_DOWN PDF_17_UP \
                PDF_18_DOWN PDF_18_UP \
                PDF_19_DOWN PDF_19_UP \
                PDF_1_DOWN PDF_1_UP \
                PDF_20_DOWN PDF_20_UP \
                PDF_21_DOWN PDF_21_UP \
                PDF_22_DOWN PDF_22_UP \
                PDF_23_DOWN PDF_23_UP \
                PDF_24_DOWN PDF_24_UP \
                PDF_25_DOWN PDF_25_UP \
                PDF_26_DOWN PDF_26_UP \
                PDF_27_DOWN PDF_27_UP \
                PDF_28_DOWN PDF_28_UP \
                PDF_29_DOWN PDF_29_UP \
                PDF_2_DOWN PDF_2_UP \
                PDF_30_DOWN PDF_30_UP \
                PDF_31_DOWN PDF_31_UP \
                PDF_32_DOWN PDF_32_UP \
                PDF_33_DOWN PDF_33_UP \
                PDF_34_DOWN PDF_34_UP \
                PDF_35_DOWN PDF_35_UP \
                PDF_36_DOWN PDF_36_UP \
                PDF_37_DOWN PDF_37_UP \
                PDF_38_DOWN PDF_38_UP \
                PDF_39_DOWN PDF_39_UP \
                PDF_3_DOWN PDF_3_UP \
                PDF_40_DOWN PDF_40_UP \
                PDF_41_DOWN PDF_41_UP \
                PDF_42_DOWN PDF_42_UP \
                PDF_43_DOWN PDF_43_UP \
                PDF_44_DOWN PDF_44_UP \
                PDF_45_DOWN PDF_45_UP \
                PDF_46_DOWN PDF_46_UP \
                PDF_47_DOWN PDF_47_UP \
                PDF_48_DOWN PDF_48_UP \
                PDF_49_DOWN PDF_49_UP \
                PDF_4_DOWN PDF_4_UP \
                PDF_50_DOWN PDF_50_UP \
                PDF_5_DOWN PDF_5_UP \
                PDF_6_DOWN PDF_6_UP \
                PDF_7_DOWN PDF_7_UP \
                PDF_8_DOWN PDF_8_UP \
                PDF_9_DOWN PDF_9_UP \
                PU_UP PU_DOWN TRIG_UP TRIG_DOWN TRIG_ETA_UP TRIG_ETA_DOWN \
                LEPT_UP LEPT_DOWN \
                UNCLUSTERED_UP UNCLUSTERED_DOWN \
                KIN_UP KIN_DOWN \
                PDF_ALPHAS_UP PDF_ALPHAS_DOWN \
                MESCALE_UP MESCALE_DOWN MEFACSCALE_UP MEFACSCALE_DOWN MERENSCALE_UP MERENSCALE_DOWN PSISRSCALE_UP PSISRSCALE_DOWN PSFSRSCALE_UP PSFSRSCALE_DOWN UETUNE_UP UETUNE_DOWN MASS_UP MASS_DOWN \
                BFRAG_UP BFRAG_DOWN BFRAG_CENTRAL BFRAG_PETERSON BSEMILEP_UP BSEMILEP_DOWN ERDON ERDONRETUNE GLUONMOVETUNE \ 
                MATCH_UP MATCH_DOWN \
                POWHEGV2HERWIG \
                BTAG_UP BTAG_DOWN BTAG_LJET_UP BTAG_LJET_DOWN \
                BTAG_PT_UP BTAG_PT_DOWN BTAG_ETA_UP BTAG_ETA_DOWN \
                BTAG_LJET_PT_UP BTAG_LJET_PT_DOWN BTAG_LJET_ETA_UP BTAG_LJET_ETA_DOWN \
                AMCATNLOFXFX
)

        if [ -d UnfoldingHistos/$Syst/$channel ] ; then
            ls -1 UnfoldingHistos/$Syst/$channel/histosTUnfold_${channel}_$sample*.root >> TUnfoldHistoFileLists/TUnfoldHistoFileList_$Syst\_$channel.txt
        fi

      end

      #ls -1 UnfoldingHistos/Nominal/$channel/histosTUnfold_${channel}_${sample}*.root | grep -v _amcatnlofxfx >> TUnfoldHistoFileLists/TUnfoldHistoFileList_Nominal\_$channel.txt
      #ls -1 UnfoldingHistos/Nominal/$channel/histosTUnfold_${channel}_${sample}*.root | grep -v _amcatnlofxfx >> TUnfoldHistoFileLists/TUnfoldHistoFileList_Nominal\_combined.txt
      ls -1 UnfoldingHistos/Nominal/$channel/histosTUnfold_${channel}_${sample}*.root | grep -v 1050.root | grep -v 50inf.root  >> TUnfoldHistoFileLists/TUnfoldHistoFileList_Nominal\_$channel.txt

   end

end
