#!/usr/bin/env python
import os
import ROOT


eras = [
"2016ULpreVFP",
"2016ULpostVFP",
"2017UL",
"2018UL",
]

eras2016 = [
"2016ULpreVFP",
"2016ULpostVFP",
]

channels = [
"ee",
"emu",
"mumu",
]

nominal = [
"Nominal",
]

expersysts = [
"JES_UP", "JES_DOWN",
"JER_UP","JER_DOWN",
"PU_UP", "PU_DOWN",
"TRIG_UP", "TRIG_DOWN",
"ELE_ID_UP", "ELE_ID_DOWN",
"ELE_RECO_UP", "ELE_RECO_DOWN",
"ELE_SCALESMEARING_UP", "ELE_SCALESMEARING_DOWN",
"ELE_SCALE_SYST_UP", "ELE_SCALE_SYST_DOWN",
"ELE_SCALE_GAIN_UP", "ELE_SCALE_GAIN_DOWN",
"ELE_SCALE_STAT_UP", "ELE_SCALE_STAT_DOWN",
"MUON_ID_UP", "MUON_ID_DOWN",
"MUON_ISO_UP", "MUON_ISO_DOWN",
"MUON_SCALE_UP", "MUON_SCALE_DOWN",
"L1PREFIRING_UP", "L1PREFIRING_DOWN",
"UNCLUSTERED_UP", "UNCLUSTERED_DOWN",
"BTAG_CORR_UP", "BTAG_CORR_DOWN",
"BTAG_LJET_CORR_UP", "BTAG_LJET_CORR_DOWN",
"BTAG_UNCORR_UP", "BTAG_UNCORR_DOWN",
"BTAG_LJET_UNCORR_UP", "BTAG_LJET_UNCORR_DOWN",
"TOP_PT",
]

JESsysts = [
"JESAbsoluteMPFBias_UP", "JESAbsoluteMPFBias_DOWN",
"JESAbsoluteScale_UP", "JESAbsoluteScale_DOWN",
"JESAbsoluteStat_UP", "JESAbsoluteStat_DOWN",
"JESFlavorQCD_UP", "JESFlavorQCD_DOWN",
"JESFlavorZJet_UP", "JESFlavorZJet_DOWN",
"JESFlavorRealistic_UP", "JESFlavorRealistic_DOWN",
"JESFlavorPhotonJet_UP", "JESFlavorPhotonJet_DOWN",
"JESFlavorPureGluon_UP", "JESFlavorPureGluon_DOWN",
"JESFlavorPureQuark_UP", "JESFlavorPureQuark_DOWN",
"JESFlavorPureCharm_UP", "JESFlavorPureCharm_DOWN",
"JESFlavorPureBottom_UP", "JESFlavorPureBottom_DOWN",
"JESFragmentation_UP", "JESFragmentation_DOWN",
"JESPileUpDataMC_UP", "JESPileUpDataMC_DOWN",
"JESPileUpPtBB_UP", "JESPileUpPtBB_DOWN",
"JESPileUpPtEC1_UP", "JESPileUpPtEC1_DOWN",
"JESPileUpPtEC2_UP", "JESPileUpPtEC2_DOWN",
"JESPileUpPtHF_UP", "JESPileUpPtHF_DOWN",
"JESPileUpPtRef_UP", "JESPileUpPtRef_DOWN",
"JESRelativeFSR_UP", "JESRelativeFSR_DOWN",
"JESRelativeJEREC1_UP", "JESRelativeJEREC1_DOWN",
"JESRelativeJEREC2_UP", "JESRelativeJEREC2_DOWN",
"JESRelativeJERHF_UP", "JESRelativeJERHF_DOWN",
"JESRelativePtBB_UP", "JESRelativePtBB_DOWN",
"JESRelativePtEC1_UP", "JESRelativePtEC1_DOWN",
"JESRelativePtEC2_UP", "JESRelativePtEC2_DOWN",
"JESRelativePtHF_UP", "JESRelativePtHF_DOWN",
"JESRelativeBal_UP", "JESRelativeBal_DOWN",
"JESRelativeSample_UP", "JESRelativeSample_DOWN",
"JESRelativeStatEC_UP", "JESRelativeStatEC_DOWN",
"JESRelativeStatFSR_UP", "JESRelativeStatFSR_DOWN",
"JESRelativeStatHF_UP", "JESRelativeStatHF_DOWN",
"JESSinglePionECAL_UP", "JESSinglePionECAL_DOWN",
"JESSinglePionHCAL_UP", "JESSinglePionHCAL_DOWN",
"JESTimePtEta_UP", "JESTimePtEta_DOWN",
]

theorysysts = [
"PDF_ALPHAS_UP", "PDF_ALPHAS_DOWN",
"MESCALE_UP", "MESCALE_DOWN",
"MEFACSCALE_UP", "MEFACSCALE_DOWN",
"MERENSCALE_UP", "MERENSCALE_DOWN",
"BFRAG_UP", "BFRAG_DOWN",
"BFRAG_PETERSON", "BFRAG_CENTRAL",
"BSEMILEP_UP", "BSEMILEP_DOWN",

"MASS_UP","MASS_DOWN",
"UETUNE_UP","UETUNE_DOWN",
"MATCH_UP","MATCH_DOWN",
"ERDON","ERDONRETUNE","GLUONMOVETUNE",
"AMCATNLOFXFX", "MADGRAPHMLM", "POWHEGV2HERWIG",

#"PSSCALE_WEIGHT_2_UP","PSSCALE_WEIGHT_2_DOWN",
#"PSSCALE_WEIGHT_3_UP","PSSCALE_WEIGHT_3_DOWN",
#"PSSCALE_WEIGHT_4_UP","PSSCALE_WEIGHT_4_DOWN",
#"PSSCALE_WEIGHT_5_UP","PSSCALE_WEIGHT_5_DOWN",
"PSSCALE_WEIGHT_6_UP","PSSCALE_WEIGHT_6_DOWN",
"PSSCALE_WEIGHT_7_UP","PSSCALE_WEIGHT_7_DOWN",

#"PSSCALE_WEIGHT_8_UP","PSSCALE_WEIGHT_8_DOWN",
#"PSSCALE_WEIGHT_9_UP","PSSCALE_WEIGHT_9_DOWN",
#"PSSCALE_WEIGHT_10_UP","PSSCALE_WEIGHT_10_DOWN",
#"PSSCALE_WEIGHT_11_UP","PSSCALE_WEIGHT_11_DOWN",
#"PSSCALE_WEIGHT_12_UP","PSSCALE_WEIGHT_12_DOWN",
#"PSSCALE_WEIGHT_13_UP","PSSCALE_WEIGHT_13_DOWN",
#"PSSCALE_WEIGHT_14_UP","PSSCALE_WEIGHT_14_DOWN",
#"PSSCALE_WEIGHT_15_UP","PSSCALE_WEIGHT_15_DOWN",
#"PSSCALE_WEIGHT_16_UP","PSSCALE_WEIGHT_16_DOWN",
#"PSSCALE_WEIGHT_17_UP","PSSCALE_WEIGHT_17_DOWN",
#"PSSCALE_WEIGHT_18_UP","PSSCALE_WEIGHT_18_DOWN",
#"PSSCALE_WEIGHT_19_UP","PSSCALE_WEIGHT_19_DOWN",
#"PSSCALE_WEIGHT_20_UP","PSSCALE_WEIGHT_20_DOWN",
#"PSSCALE_WEIGHT_21_UP","PSSCALE_WEIGHT_21_DOWN",
#"PSSCALE_WEIGHT_22_UP","PSSCALE_WEIGHT_22_DOWN",
#"PSSCALE_WEIGHT_23_UP","PSSCALE_WEIGHT_23_DOWN",

"PDF_0_CENTRAL",
"PDF_1_UP", "PDF_1_DOWN",
"PDF_2_UP", "PDF_2_DOWN",
"PDF_3_UP", "PDF_3_DOWN",
"PDF_4_UP", "PDF_4_DOWN",
"PDF_5_UP", "PDF_5_DOWN",
"PDF_6_UP", "PDF_6_DOWN",
"PDF_7_UP", "PDF_7_DOWN",
"PDF_8_UP", "PDF_8_DOWN",
"PDF_9_UP", "PDF_9_DOWN",
"PDF_10_UP", "PDF_10_DOWN",
"PDF_11_UP", "PDF_11_DOWN",
"PDF_12_UP", "PDF_12_DOWN",
"PDF_13_UP", "PDF_13_DOWN",
"PDF_14_UP", "PDF_14_DOWN",
"PDF_15_UP", "PDF_15_DOWN",
"PDF_16_UP", "PDF_16_DOWN",
"PDF_17_UP", "PDF_17_DOWN",
"PDF_18_UP", "PDF_18_DOWN",
"PDF_19_UP", "PDF_19_DOWN",
"PDF_20_UP", "PDF_20_DOWN",
"PDF_21_UP", "PDF_21_DOWN",
"PDF_22_UP", "PDF_22_DOWN",
"PDF_23_UP", "PDF_23_DOWN",
"PDF_24_UP", "PDF_24_DOWN",
"PDF_25_UP", "PDF_25_DOWN",
"PDF_26_UP", "PDF_26_DOWN",
"PDF_27_UP", "PDF_27_DOWN",
"PDF_28_UP", "PDF_28_DOWN",
"PDF_29_UP", "PDF_29_DOWN",
"PDF_30_UP", "PDF_30_DOWN",
"PDF_31_UP", "PDF_31_DOWN",
"PDF_32_UP", "PDF_32_DOWN",
"PDF_33_UP", "PDF_33_DOWN",
"PDF_34_UP", "PDF_34_DOWN",
"PDF_35_UP", "PDF_35_DOWN",
"PDF_36_UP", "PDF_36_DOWN",
"PDF_37_UP", "PDF_37_DOWN",
"PDF_38_UP", "PDF_38_DOWN",
"PDF_39_UP", "PDF_39_DOWN",
"PDF_40_UP", "PDF_40_DOWN",
"PDF_41_UP", "PDF_41_DOWN",
"PDF_42_UP", "PDF_42_DOWN",
"PDF_43_UP", "PDF_43_DOWN",
"PDF_44_UP", "PDF_44_DOWN",
"PDF_45_UP", "PDF_45_DOWN",
"PDF_46_UP", "PDF_46_DOWN",
"PDF_47_UP", "PDF_47_DOWN",
"PDF_48_UP", "PDF_48_DOWN",
"PDF_49_UP", "PDF_49_DOWN",
"PDF_50_UP", "PDF_50_DOWN",
]

virtualsysts = [
"TOT_PDF_UP","TOT_PDF_DOWN",
"TOT_SCALE_UP","TOT_SCALE_DOWN",
"TOT_BFRAG_UP", "TOT_BFRAG_DOWN",
"TOT_COLORREC_UP", "TOT_COLORREC_DOWN",
"DY_UP","DY_DOWN",
"BG_UP","BG_DOWN",
#"TW_UP", "TW_DOWN",
]

completedatasamples = [
[
"run2016B",
"run2016C",
"run2016D",
"run2016E",
"run2016F1",
"se_run2016B",
"se_run2016C",
"se_run2016D",
"se_run2016E",
"se_run2016F1",
"smu_run2016B",
"smu_run2016C",
"smu_run2016D",
"smu_run2016E",
"smu_run2016F1",
],
[
"run2016F2",
"run2016G",
"run2016H",
"se_run2016F2",
"se_run2016G",
"se_run2016H",
"smu_run2016F2",
"smu_run2016G",
"smu_run2016H",
],
[
"run2017B",
"run2017C",
"run2017D",
"run2017E",
"run2017F",
"se_run2017B",
"se_run2017C",
"se_run2017D",
"se_run2017E",
"se_run2017F",
"smu_run2017B",
"smu_run2017C",
"smu_run2017D",
"smu_run2017E",
"smu_run2017F",
],
[
"run2018A",
"run2018B",
"run2018C",
"run2018D",
"se_run2018A",
"se_run2018B",
"se_run2018C",
"se_run2018D",
"smu_run2018A",
"smu_run2018B",
"smu_run2018C",
"smu_run2018D",
]
]

signalMCsamples = [
"ttbarsignalplustau_fromDilepton",
"ttbarsignalviatau_fromDilepton",
]

backgroundMCsamples  = [
"dyee1050",
"dyee50inf_amcatnlofxfx",
"dyee50inf_madgraphmlm",
"dymumu1050",
"dymumu50inf_amcatnlofxfx",
"dymumu50inf_madgraphmlm",
"dytautau1050",
"dytautau50inf_amcatnlofxfx",
"dytautau50inf_madgraphmlm",
"wtolnu",
"wwtoall",
"wztoall",
"zztoall",
"ttbarWjetstolnu",
"ttbarWjetstoqq",
"ttbarZtollnunu",
"ttbarZtoqq",
"singleantitop_tw",
"singletop_tw",
"ttbarbg_fromDilepton",
"ttbarbg_fromHadronic",
"ttbarbg_fromLjets",
]


# Make File Lists for each era/channel

for era in eras:

    checkfile =  open("check_FileLists_"+era+".txt","w")

    if not os.path.exists("FileLists_"+era.replace("UL","")): os.makedirs("FileLists_"+era.replace("UL",""))
    
    for channel in channels:

        for syst in nominal + expersysts + JESsysts + theorysysts + virtualsysts:

            file = open("FileLists_"+era.replace("UL","")+"/HistoFileList_"+syst+"_"+channel+".txt","w")

            if (era.replace("UL","") == "2016preVFP"): datasamples = completedatasamples[0]
            if (era.replace("UL","") == "2016postVFP"): datasamples = completedatasamples[1]
            if (era.replace("UL","") == "2017"): datasamples = completedatasamples[2]
            if (era.replace("UL","") == "2018"): datasamples = completedatasamples[3]

            for sample in datasamples:

                file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+".root\n")
                if not os.path.exists("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+".root"): checkfile.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+".root\n")

            for sample in backgroundMCsamples + signalMCsamples:

                if (syst in nominal + expersysts + JESsysts) and (sample in signalMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")
                    if not os.path.exists("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root"): checkfile.write("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in theorysysts) and (sample in signalMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")
                    if not os.path.exists("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root"): checkfile.write("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in nominal + expersysts + JESsysts + theorysysts) and (sample in backgroundMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")
                    if not os.path.exists("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root"): checkfile.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")
                
                elif (syst in virtualsysts):

                    file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")
                    if not os.path.exists("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root"): checkfile.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

            file.close()

    checkfile.close()

    # combined channel
    for syst in nominal + expersysts + JESsysts + theorysysts + virtualsysts:

        file = open("FileLists_"+era.replace("UL","")+"/HistoFileList_"+syst+"_combined.txt","w")

        if (era.replace("UL","") == "2016preVFP"): datasamples = completedatasamples[0]
        if (era.replace("UL","") == "2016postVFP"): datasamples = completedatasamples[1]
        if (era.replace("UL","") == "2017"): datasamples = completedatasamples[2]
        if (era.replace("UL","") == "2018"): datasamples = completedatasamples[3]
        
        for sample in datasamples:

            for channel in channels:
                            
                file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+".root\n")
            
        for sample in backgroundMCsamples + signalMCsamples:

            for channel in channels:

                if (syst in nominal + expersysts + JESsysts) and (sample in signalMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in theorysysts) and (sample in signalMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in nominal + expersysts + JESsysts + theorysysts) and (sample in backgroundMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in virtualsysts):

                    file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")
                        
        file.close()


# Make File Lists for 2016 for each channel
if not os.path.exists("FileLists_2016"): os.makedirs("FileLists_2016")
    
for channel in channels:

    for syst in nominal + expersysts + JESsysts + theorysysts + virtualsysts:

        file = open("FileLists_2016/HistoFileList_"+syst+"_"+channel+".txt","w")

        for era in eras2016:

            if (era.replace("UL","") == "2016preVFP"): datasamples = completedatasamples[0]
            if (era.replace("UL","") == "2016postVFP"): datasamples = completedatasamples[1]

            for sample in datasamples:
            
                file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+".root\n")

        for sample in backgroundMCsamples + signalMCsamples:

            for era in eras2016:

                if (syst in nominal + expersysts + JESsysts) and (sample in signalMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in theorysysts) and (sample in signalMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in nominal + expersysts + JESsysts + theorysysts) and (sample in backgroundMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in virtualsysts):

                    file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

        file.close()


# Make File Lists for 2016 for combined channel
for syst in nominal + expersysts + JESsysts + theorysysts + virtualsysts:

    file = open("FileLists_2016/HistoFileList_"+syst+"_combined.txt","w")

    for era in eras2016:

        if (era.replace("UL","") == "2016preVFP"): datasamples = completedatasamples[0]
        if (era.replace("UL","") == "2016postVFP"): datasamples = completedatasamples[1]
        
        for channel in channels:

            for sample in datasamples:
            
                file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+".root\n")
            
    for sample in backgroundMCsamples + signalMCsamples:
            
        for era in eras2016:

            for channel in channels:

                if (syst in nominal + expersysts + JESsysts) and (sample in signalMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in theorysysts) and (sample in signalMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in nominal + expersysts + JESsysts + theorysysts) and (sample in backgroundMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in virtualsysts):

                    file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

    file.close()


# Make File Lists for full run II for each channel
if not os.path.exists("FileLists_fullRun2UL"): os.makedirs("FileLists_fullRun2UL")
    
for channel in channels:

    for syst in nominal + expersysts + JESsysts + theorysysts + virtualsysts:

        file = open("FileLists_fullRun2UL/HistoFileList_"+syst+"_"+channel+".txt","w")

        for era in eras:

            if (era.replace("UL","") == "2016preVFP"): datasamples = completedatasamples[0]
            if (era.replace("UL","") == "2016postVFP"): datasamples = completedatasamples[1]
            if (era.replace("UL","") == "2017"): datasamples = completedatasamples[2]
            if (era.replace("UL","") == "2018"): datasamples = completedatasamples[3]

            for sample in datasamples:
            
                file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+".root\n")

        for sample in backgroundMCsamples + signalMCsamples:

            for era in eras:

                if (syst in nominal + expersysts + JESsysts) and (sample in signalMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in theorysysts) and (sample in signalMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in nominal + expersysts + JESsysts + theorysysts) and (sample in backgroundMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in virtualsysts):

                    file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

        file.close()


# Make File Lists for full run II for combined channel
for syst in nominal + expersysts + JESsysts + theorysysts + virtualsysts:

    file = open("FileLists_fullRun2UL/HistoFileList_"+syst+"_combined.txt","w")

    for era in eras:

        if (era.replace("UL","") == "2016preVFP"): datasamples = completedatasamples[0]
        if (era.replace("UL","") == "2016postVFP"): datasamples = completedatasamples[1]
        if (era.replace("UL","") == "2017"): datasamples = completedatasamples[2]
        if (era.replace("UL","") == "2018"): datasamples = completedatasamples[3]
        
        for channel in channels:

            for sample in datasamples:
            
                file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+".root\n")
            
    for sample in backgroundMCsamples + signalMCsamples:
            
        for era in eras:

            for channel in channels:

                if (syst in nominal + expersysts) and (sample in signalMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in theorysysts) and (sample in signalMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in nominal + expersysts + JESsysts + theorysysts) and (sample in backgroundMCsamples):

                    file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                elif (syst in virtualsysts):

                    file.write("selectionRoot_"+era.replace("UL","")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

    file.close()
