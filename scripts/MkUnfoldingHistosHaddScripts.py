#!/usr/bin/env python
import os
import ROOT

eras = [
"2016ULpreVFP",
"2016ULpostVFP",
"2017UL",
"2018UL",
]

channels = [
"ee",
"emu",
"mumu",
]

nominal = [
"Nominal"
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
"PDF_ALPHAS_UP","PDF_ALPHAS_DOWN",
"MESCALE_UP","MESCALE_DOWN",
"MEFACSCALE_UP","MEFACSCALE_DOWN",
"MERENSCALE_UP","MERENSCALE_DOWN",
"BFRAG_UP","BFRAG_DOWN",
"BFRAG_PETERSON", "BFRAG_CENTRAL",
"BSEMILEP_UP","BSEMILEP_DOWN",

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
]

theorysysts_dedicated = {
"MASS_UP":"ttbarsignalplustau_fromDilepton_173_massup", 
"MASS_DOWN":"ttbarsignalplustau_fromDilepton_171_massdown",
"UETUNE_UP":"ttbarsignalplustau_fromDilepton_uetuneup", 
"UETUNE_DOWN":"ttbarsignalplustau_fromDilepton_uetunedown",
"MATCH_UP":"ttbarsignalplustau_fromDilepton_matchup", 
"MATCH_DOWN":"ttbarsignalplustau_fromDilepton_matchdown",
"ERDON":"ttbarsignalplustau_fromDilepton_erdon", 
"ERDONRETUNE":"ttbarsignalplustau_fromDilepton_erdonretune", 
"GLUONMOVETUNE":"ttbarsignalplustau_fromDilepton_gluonmovetune",
"AMCATNLOFXFX":"ttbarsignalplustau_amcatnlofxfx",
"MADGRAPHMLM":"ttbarsignalplustau_madgraphmlm",
"POWHEGV2HERWIG":"ttbarsignalplustau_powhegv2Herwig",
}

virtualsysts = [
"DY_UP","DY_DOWN",
"BG_UP","BG_DOWN",
#"TW_UP", "TW_DOWN",
]

signalMCsamples = [
"ttbarsignalplustau_fromDilepton",
"ttbarsignalviatau_fromDilepton",
]

backgroundMCsamples  = [
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
]

dy50infMCsamples = [
"dyee50inf_amcatnlofxfx",
"dyee50inf_madgraphmlm",
"dymumu50inf_amcatnlofxfx",
"dymumu50inf_madgraphmlm",
"dytautau50inf_amcatnlofxfx",
"dytautau50inf_madgraphmlm",
]

ttbarZMCsamples = [
"ttbarZtollnunu",
"ttbarZtoqq",
]

nSignalFiles = -1
nSystFiles = -1

for era in eras:

    if (era == "2016ULpreVFP"):
        nSignalFiles = 20
        nSystFiles = 20
    elif (era == "2016ULpostVFP"):
        nSignalFiles = 20
        nSystFiles = 20
    elif (era == "2017UL"):
        nSignalFiles = 30
        nSystFiles = 20
    elif (era == "2018UL"):
        nSignalFiles = 40
        nSystFiles = 20

    haddfile = open("scripts/hadd_UnfoldingHistos_"+era+".sh","w")
    checkfile =  open("check_UnfoldingHistos_"+era+".txt","w")

    for channel in channels:

        for syst in nominal + expersysts + theorysysts + list(theorysysts_dedicated.keys()):

            for sample in signalMCsamples:

                if (syst in nominal + expersysts + theorysysts):

                    haddfile.write("nohup hadd -f -j 12 ")

                    haddfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root ")

                    if (syst in list(theorysysts_dedicated.keys()) and sample == "ttbarsignalplustau_fromDilepton"): sample = theorysysts_dedicated[syst]
                    elif (syst in list(theorysysts_dedicated.keys()) and sample == "ttbarsignalviatau_fromDilepton"): sample = theorysysts_dedicated[syst].replace("plustau","viatau")

                    for i in range(0,nSignalFiles):

                        haddfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root ")
                        if not os.path.exists("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"+" \n")

                    haddfile.write(" &> nohuplogs/hadd_"+syst+"_"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".out ")

                    haddfile.write("\n")

                elif (syst in list(theorysysts_dedicated.keys())):

                    haddfile.write("nohup hadd -f -j 12 ")

                    haddfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root ")

                    if (syst in list(theorysysts_dedicated.keys()) and sample == "ttbarsignalplustau_fromDilepton"): sample = theorysysts_dedicated[syst]
                    elif (syst in list(theorysysts_dedicated.keys()) and sample == "ttbarsignalviatau_fromDilepton"): sample = theorysysts_dedicated[syst].replace("plustau","viatau")

                    for i in range(0,nSystFiles):

                        haddfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root ")
                        if not os.path.exists("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"+" \n")

                    haddfile.write(" &> nohuplogs/hadd_"+syst+"_"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".out ")

                    haddfile.write("\n")


            if (syst in expersysts + theorysysts + list(theorysysts_dedicated.keys())): continue


            for sample in backgroundMCsamples:

                if not os.path.exists("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root"): checkfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root"+" \n")


            for sample in dy50infMCsamples:

                haddfile.write("nohup hadd -f -j 12 ")

                haddfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root ")

                for i in range(0,10):

                    haddfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root ")
                    if not os.path.exists("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"+" \n")

                haddfile.write(" &> nohuplogs/hadd_"+syst+"_"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".out ")

                haddfile.write("\n")


            for sample in ttbarZMCsamples:

                haddfile.write("nohup hadd -f -j 12 ")

                haddfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root ")

                for i in range(0,5):

                    haddfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root ")
                    if not os.path.exists("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"+" \n")

                haddfile.write(" &> nohuplogs/hadd_"+syst+"_"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".out ")

                haddfile.write("\n")


        for syst in virtualsysts:

            for sample in backgroundMCsamples:

                if ((syst == "DY_UP" or syst == "DY_DOWN") and (sample == "dyee1050" or sample == "dymumu1050" or sample == "dytautau1050")):
                
                    if not os.path.exists("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root"): checkfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root"+" \n")

            if (syst == "DY_UP" or syst == "DY_DOWN"):

                for sample in dy50infMCsamples:

                    haddfile.write("nohup hadd -f -j 12 ")

                    haddfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root ")

                    for i in range(0,10):

                        haddfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root ")
                        if not os.path.exists("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"+" \n")

                    haddfile.write(" &> nohuplogs/hadd_"+syst+"_"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".out ")

                    haddfile.write("\n")


            for sample in backgroundMCsamples:

                if ((syst == "BG_UP" or syst == "BG_DOWN") and (sample != "dyee1050" and sample != "dymumu1050" and sample != "dytautau1050")):

                    if not os.path.exists("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root"): checkfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root"+" \n")

            if (syst == "BG_UP" or syst == "BG_DOWN"):

                for sample in ttbarZMCsamples:

                    haddfile.write("nohup hadd -f -j 12 ")

                    haddfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".root ")

                    for i in range(0,5):

                        haddfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root ")
                        if not os.path.exists("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write("UnfoldingHistos_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+"histosTUnfold_"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"+" \n")

                    haddfile.write(" &> nohuplogs/hadd_"+syst+"_"+"histosTUnfold_"+channel+"_"+sample+"_"+era+".out ")

                    haddfile.write("\n")


            haddfile.write("\n")

    haddfile.close()
    checkfile.close()
