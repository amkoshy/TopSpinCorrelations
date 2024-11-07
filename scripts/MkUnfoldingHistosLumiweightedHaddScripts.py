#!/usr/bin/env python
import os

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
"BFRAG_PETERSON","BFRAG_CENTRAL",
"BSEMILEP_UP", "BSEMILEP_DOWN",

"MASS_UP","MASS_DOWN",
"UETUNE_UP","UETUNE_DOWN",
"MATCH_UP","MATCH_DOWN",
"ERDON","ERDONRETUNE","GLUONMOVETUNE",

"PSSCALE_WEIGHT_6_UP","PSSCALE_WEIGHT_6_DOWN",
"PSSCALE_WEIGHT_7_UP","PSSCALE_WEIGHT_7_DOWN",

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
"AMCATNLOFXFX",
]

datasamples = {
"2016ULpreVFP":
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
"2016ULpostVFP":
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
"2017UL":
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
"2018UL":
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
}


signalsamples = [
"ttbarsignalplustau_fromDilepton",
"ttbarsignalviatau_fromDilepton",
]

backgroundsamples = [
"dyee1050",
"dymumu1050",
"dytautau1050",
"ttbarWjetstolnu",
"ttbarWjetstoqq",
"singleantitop_tw",
"singletop_tw",
"wtolnu",
"wwtoall",
"wztoall",
"zztoall",
"ttbarbg_fromDilepton",
"ttbarbg_fromLjets",
"ttbarbg_fromHadronic",
"dyee50inf_amcatnlofxfx",
"dyee50inf_madgraphmlm",
"dymumu50inf_amcatnlofxfx",
"dymumu50inf_madgraphmlm",
"dytautau50inf_amcatnlofxfx",
"dytautau50inf_madgraphmlm",
"ttbarZtollnunu",
"ttbarZtoqq",
]


# fullRun2UL
haddfile = open("scripts/hadd_UnfoldingHistos_Lumiweighted_fullRun2UL.sh","w")
checkfile =  open("check_UnfoldingHistos_Lumiweighted_fullRun2UL.txt","w")

if not os.path.exists("UnfoldingHistos_Lumiweighted_fullRun2"): os.makedirs("UnfoldingHistos_Lumiweighted_fullRun2")

for syst in nominal + expersysts + theorysysts:

    if syst in nominal:

        if not os.path.exists("UnfoldingHistos_Lumiweighted_fullRun2"+"/"+syst): os.makedirs("UnfoldingHistos_Lumiweighted_fullRun2"+"/"+syst)

        for channel in channels:

            if not os.path.exists("UnfoldingHistos_Lumiweighted_fullRun2"+"/"+syst+"/"+channel): os.makedirs("UnfoldingHistos_Lumiweighted_fullRun2"+"/"+syst+"/"+channel)

            for sample in signalsamples + backgroundsamples:

                haddfile.write("nohup hadd -f -j 12 ")

                haddfile.write("UnfoldingHistos_Lumiweighted_fullRun2/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_fullRun2UL.root ")

                for era in eras:

                    haddfile.write("UnfoldingHistos_Lumiweighted_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root ")

                haddfile.write("\n")

            haddfile.write("nohup hadd -f -j 12 ")

            haddfile.write("UnfoldingHistos_Lumiweighted_fullRun2/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_run_fullRun2.root ")

            for era in eras:

                for sample in datasamples[era]:

                    haddfile.write("UnfoldingHistos_Lumiweighted_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+".root ")

            haddfile.write("\n")

    else:

        if not os.path.exists("UnfoldingHistos_Lumiweighted_fullRun2"+"/"+syst): os.makedirs("UnfoldingHistos_Lumiweighted_fullRun2"+"/"+syst)

        for channel in channels:

            if not os.path.exists("UnfoldingHistos_Lumiweighted_fullRun2"+"/"+syst+"/"+channel): os.makedirs("UnfoldingHistos_Lumiweighted_fullRun2"+"/"+syst+"/"+channel)

            for sample in signalsamples:

                haddfile.write("nohup hadd -f -j 12 ")

                haddfile.write("UnfoldingHistos_Lumiweighted_fullRun2/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_fullRun2UL.root ")

                for era in eras:

                    haddfile.write("UnfoldingHistos_Lumiweighted_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root ")

                haddfile.write("\n")


                    
# 2016UL
haddfile = open("scripts/hadd_UnfoldingHistos_Lumiweighted_2016UL.sh","w")
checkfile =  open("check_UnfoldingHistos_Lumiweighted_2016UL.txt","w")

if not os.path.exists("UnfoldingHistos_Lumiweighted_2016"): os.makedirs("UnfoldingHistos_Lumiweighted_2016")

for syst in nominal + expersysts + theorysysts:

    if syst in nominal:

        if not os.path.exists("UnfoldingHistos_Lumiweighted_2016"+"/"+syst): os.makedirs("UnfoldingHistos_Lumiweighted_2016"+"/"+syst)

        for channel in channels:

            if not os.path.exists("UnfoldingHistos_Lumiweighted_2016"+"/"+syst+"/"+channel): os.makedirs("UnfoldingHistos_Lumiweighted_2016"+"/"+syst+"/"+channel)

            for sample in signalsamples + backgroundsamples:

                haddfile.write("nohup hadd -f -j 12 ")

                haddfile.write("UnfoldingHistos_Lumiweighted_2016/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_2016UL.root ")

                for era in eras2016:

                    haddfile.write("UnfoldingHistos_Lumiweighted_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root ")

                haddfile.write("\n")

            haddfile.write("nohup hadd -f -j 12 ")

            haddfile.write("UnfoldingHistos_Lumiweighted_2016/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_run_2016.root ")

            for era in eras2016:

                for sample in datasamples[era]:

                    haddfile.write("UnfoldingHistos_Lumiweighted_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+".root ")

            haddfile.write("\n")

    else:

        if not os.path.exists("UnfoldingHistos_Lumiweighted_2016"+"/"+syst): os.makedirs("UnfoldingHistos_Lumiweighted_2016"+"/"+syst)

        for channel in channels:

            if not os.path.exists("UnfoldingHistos_Lumiweighted_2016"+"/"+syst+"/"+channel): os.makedirs("UnfoldingHistos_Lumiweighted_2016"+"/"+syst+"/"+channel)

            for sample in signalsamples:

                haddfile.write("nohup hadd -f -j 12 ")

                haddfile.write("UnfoldingHistos_Lumiweighted_2016/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_2016UL.root ")

                for era in eras2016:

                    haddfile.write("UnfoldingHistos_Lumiweighted_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root ")

                haddfile.write("\n")


                    

                


