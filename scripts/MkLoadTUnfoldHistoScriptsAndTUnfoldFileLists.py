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

theorysysts = [
"PDF_ALPHAS_UP", "PDF_ALPHAS_DOWN",
"MESCALE_UP", "MESCALE_DOWN",
"MEFACSCALE_UP", "MEFACSCALE_DOWN",
"MERENSCALE_UP", "MERENSCALE_DOWN",
"BFRAG_UP", "BFRAG_DOWN",
"BFRAG_PETERSON",
"BSEMILEP_UP", "BSEMILEP_DOWN",

#"PSSCALE_WEIGHT_2_UP","PSSCALE_WEIGHT_2_DOWN",
#"PSSCALE_WEIGHT_3_UP","PSSCALE_WEIGHT_3_DOWN",
#"PSSCALE_WEIGHT_4_UP","PSSCALE_WEIGHT_4_DOWN",
#"PSSCALE_WEIGHT_5_UP","PSSCALE_WEIGHT_5_DOWN",
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
"dymumu1050",
"dytautau1050",
"wtolnu",
"wwtoall",
"wztoall",
"zztoall",
"ttbarWjetstolnu",
"ttbarWjetstoqq",
"singleantitop_tw",
"singletop_tw",
"ttbarbg_fromDilepton",
"ttbarbg_fromHadronic",
"ttbarbg_fromLjets",
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

# Make File Lists for each era/channel

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

    LoadTUnfoldscriptfile = open("scripts/LoadTUnfoldHisto_"+era+".sh","w")
    selectionRoot2Dscriptfile = open("scripts/selectionRoot2D_"+era+".sh","w")

    parallelLoadTUnfoldscriptfile = open("scripts/LoadTUnfoldHisto_"+era+"_parallel.sh","w")

    checkfile =  open("check_TUnfoldFileLists_"+era+".txt","w")

    if era == "2016ULpreVFP": location = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016preVFP"
    elif era == "2016ULpostVFP": location = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2016postVFP"
    elif era == "2017UL": location = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2017"
    elif era == "2018UL": location = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/spinCorrInput_2018"

    if not os.path.exists("TUnfoldFileLists_"+era.replace("UL","")): os.makedirs("TUnfoldFileLists_"+era.replace("UL",""))
    
    for channel in channels:

        for syst in nominal + expersysts + theorysysts + virtualsysts + list(theorysysts_dedicated.keys()):

            LoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+era+".out \n")
            selectionRoot2Dscriptfile.write("nohup python scripts/makeselectionRoot2D.py -y "+era+" -s "+syst+" -c "+channel+" &> nohuplogs/selectionRoot2D_"+syst+"_"+channel+"_"+era+".out \n")

            file = open("TUnfoldFileLists_"+era.replace("UL","")+"/TUnfoldFileList_"+syst+"_"+channel+".txt","w")

            if (era.replace("UL","") == "2016preVFP"): datasamples = completedatasamples[0]
            if (era.replace("UL","") == "2016postVFP"): datasamples = completedatasamples[1]
            if (era.replace("UL","") == "2017"): datasamples = completedatasamples[2]
            if (era.replace("UL","") == "2018"): datasamples = completedatasamples[3]

            for sample in datasamples:

                if (syst in nominal):

                    file.write(location+"/Nominal/"+channel+"/"+channel+"_"+sample+".root\n")
                    parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+".out \n")
                    if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+".root\n")

            for sample in backgroundMCsamples + signalMCsamples + dy50infMCsamples + ttbarZMCsamples:

                if (syst in nominal):

                    if sample in signalMCsamples:

                        for i in range(0,nSignalFiles):

                            file.write(location+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")
                            parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+"_"+era+"_"+str(i)+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+"_"+str(i)+".out \n")
                            if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")

                    elif sample in dy50infMCsamples:

                        for i in range(0,10):

                            file.write(location+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")
                            parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+"_"+era+"_"+str(i)+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+"_"+str(i)+".out \n")
                            if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")

                    elif sample in ttbarZMCsamples:

                        for i in range(0,5):

                            file.write(location+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")
                            parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+"_"+era+"_"+str(i)+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+"_"+str(i)+".out \n")
                            if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")

                    elif sample in backgroundMCsamples:

                        file.write(location+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")
                        parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+"_"+era+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+".out \n")
                        if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")


                elif (syst in expersysts + theorysysts + list(theorysysts_dedicated.keys())) and (sample in signalMCsamples):
                    
                    if syst in list(theorysysts_dedicated.keys()):
                        if sample == "ttbarsignalplustau_fromDilepton": sample = theorysysts_dedicated[syst]
                        elif sample == "ttbarsignalviatau_fromDilepton": sample = theorysysts_dedicated[syst].replace("plustau","viatau")

                    if (syst in expersysts + theorysysts):
                        
                        for i in range(0,nSignalFiles):

                            file.write(location+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")
                            parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+"_"+era+"_"+str(i)+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+"_"+str(i)+".out \n")
                            if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")


                    else:

                        for i in range(0,nSystFiles):

                            file.write(location+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")
                            parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+"_"+era+"_"+str(i)+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+"_"+str(i)+".out \n")
                            if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/"+syst+"/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")

                #elif (syst in theorysysts + list(theorysysts_dedicated.keys())) and (sample in backgroundMCsamples + dy50infMCsamples + ttbarZMCsamples):

                    #if sample in dy50infMCsamples:

                        #for i in range(0,10):

                            #file.write(location+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")
                            #parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+"_"+era+"_"+str(i)+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+"_"+str(i)+".out \n")
                            #if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")

                    #elif sample in ttbarZMCsamples:

                        #for i in range(0,5):

                            #file.write(location+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")
                            #parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+"_"+era+"_"+str(i)+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+"_"+str(i)+".out \n")
                            #if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")

                    #else:

                        #file.write(location+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")
                        #parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+"_"+era+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+".out \n")
                        #if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")


                elif (syst in virtualsysts):

                    if ((sample in dy50infMCsamples) and (syst == "DY_UP" or syst == "DY_DOWN")):

                        for i in range(0,10):

                            file.write(location+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")
                            parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+"_"+era+"_"+str(i)+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+"_"+str(i)+".out \n")
                            if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")

                    elif ((sample in ttbarZMCsamples) and (syst == "BG_UP" or syst == "BG_DOWN")):

                        for i in range(0,5):

                            file.write(location+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")
                            parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+"_"+era+"_"+str(i)+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+"_"+str(i)+".out \n")
                            if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")

                    elif ((sample in backgroundMCsamples)  and (syst == "BG_UP" or syst == "BG_DOWN") and (sample != "dyee1050" and sample != "dymumu1050" and sample != "dytautau1050")):

                        file.write(location+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")
                        parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+"_"+era+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+".out \n")
                        if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                    elif ((sample in backgroundMCsamples)  and (syst == "DY_UP" or syst == "DY_DOWN") and (sample == "dyee1050" or sample == "dymumu1050" or sample == "dytautau1050")):

                        file.write(location+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")
                        parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+"_"+era+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+".out \n")
                        if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+".root\n")

                    #elif sample in signalMCsamples:

                        #if era == "2018UL":
                        
                            #for i in range(0,40):

                                #file.write(location+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")
                                #parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+"_"+era+"_"+str(i)+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+"_"+str(i)+".out \n")
                                #if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")

                        #else:
                        
                            #for i in range(0,20):

                                #file.write(location+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")
                                #parallelLoadTUnfoldscriptfile.write("nohup ./install/bin/loadTUnfoldHisto -y "+era.replace("UL","")+" -s "+syst+" -c "+channel+" -f "+channel+"_"+sample+"_"+era+"_"+str(i)+".root &> nohuplogs/LoadTUnfoldHisto_"+syst+"_"+channel+"_"+sample+"_"+era+"_"+str(i)+".out \n")
                                #if not os.path.exists(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root"): checkfile.write(location.replace("root://xrootd.rcac.purdue.edu//","/mnt/hadoop/")+"/Nominal/"+channel+"/"+channel+"_"+sample+"_"+era+"_"+str(i)+".root\n")
                    
                
            file.close()

    LoadTUnfoldscriptfile.close()
    selectionRoot2Dscriptfile.close()

    parallelLoadTUnfoldscriptfile.close()

    checkfile.close()
