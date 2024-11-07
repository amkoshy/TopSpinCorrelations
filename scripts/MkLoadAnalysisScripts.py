#!/usr/bin/env python

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

samples  = [
"dy1050_201 -d 11",
"dy1050_201 -d 13",
"dy1050_201 -d 15",
"ttbarWjetstolnu_201",
"ttbarWjetstoqq_201",
"top_tw_201",
"wtolnu_201",
"wwtoall_201",
"wztoall_201",
"zztoall_201",
"ttbarbg_fromDilepton_201",
"ttbarbg_fromLjets_201",
"ttbarbg_fromHadronic_201",
]

nominal = [
"Nominal"
]

expersysts = [
"JES_UP", "JES_DOWN",
"JER_UP", "JER_DOWN",
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

#nominalmodes = "spinCorr cp kinReco kinRecoQualityStudies"
nominalmodes = "spinCorr kinReco kinRecoQualityStudies"
modes = "spinCorr"

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

    BTag_scriptfile = open("scripts/LoadAnalysis_"+"BTag_"+era+".sh","w")
    BTag_scriptfile.write("#!/bin/bash \n")
    
    for channel in channels:

        for syst in nominal + expersysts + JESsysts + theorysysts + theorysysts_dedicated.keys():

            if (syst in nominal + expersysts + JESsysts + theorysysts):
                
                sample = "ttbarsignalplustau_fromDilepton"
                    
                for i in range(0,nSignalFiles):

                    BTag_scriptfile.write("nohup ./install/bin/load_Analysis -f "+sample+"_"+era+"_"+str(i)+".root -c "+channel+" -s "+syst+" &> nohuplogs/BTag_"+sample+"_"+syst+"_"+channel+"_"+era+"_"+str(i)+".out \n")

            elif syst in theorysysts_dedicated.keys():

                sample = theorysysts_dedicated[syst]

                for i in range(0,nSystFiles):
        
                    BTag_scriptfile.write("nohup ./install/bin/load_Analysis -f "+sample+"_"+era+"_"+str(i)+".root -c "+channel+" &> nohuplogs/BTag_"+sample+"_"+syst+"_"+channel+"_"+era+"_"+str(i)+".out \n")

        # PS

        #for j in range(10,14):
                                     
            #if era == "2018UL":
                    
                #for i in range(0,40):

                    #sample = "ttbarsignalplustau_fromDilepton"
                
                    #BTag_scriptfile.write("nohup ./install/bin/load_Analysis -f "+sample+"_"+era+"_"+str(i)+".root -c "+channel+" -s PSSCALE_WEIGHT --ps "+str(j)+" &> nohuplogs/Btag_"+sample+"_PSSCALE_WEIGHT_"+str(j)+"_"+channel+"_"+era+"_"+str(i)+".out \n")

            #else: 

                #for i in range(0,20):

                    #sample = "ttbarsignalplustau_fromDilepton"
                
                    #BTag_scriptfile.write("nohup ./install/bin/load_Analysis -f "+sample+"_"+era+"_"+str(i)+".root -c "+channel+" -s PSSCALE_WEIGHT --ps "+str(j)+" &> nohuplogs/Btag_"+sample+"_PSSCALE_WEIGHT_"+str(j)+"_"+channel+"_"+era+"_"+str(i)+".out \n")

        BTag_scriptfile.write("\n")


    BTag_scriptfile.close()



    Nominal_scriptfile = open("scripts/LoadAnalysis_SpinCorr_"+"Nominal_MC_"+era+".sh","w")
    Nominal_scriptfile.write("#!/bin/bash \n")
    
    for channel in channels:

        for sample in samples:

            Nominal_scriptfile.write("nohup ./install/bin/load_Analysis -f "+sample+" -m "+nominalmodes+" -c "+channel+" -s Nominal &> nohuplogs/SpinCorr_Nominal_MC_"+sample.replace(" -d ","_").replace(" --signalviatau","_signalviatau").replace(".root","")+"_"+channel+"_"+era+".out \n")

        for i in range(0,nSignalFiles):
            
            Nominal_scriptfile.write("nohup ./install/bin/load_Analysis -f ttbarsignalplustau_fromDilepton_"+era+"_"+str(i)+".root -m "+nominalmodes+" -c "+channel+" -s Nominal &> nohuplogs/SpinCorr_Nominal_MC_ttbarsignalplustau_fromDilepton_"+str(i)+"_"+channel+"_"+era+".out \n")

        for i in range(0,nSignalFiles):

            Nominal_scriptfile.write("nohup ./install/bin/load_Analysis -f ttbarsignalplustau_fromDilepton_"+era+"_"+str(i)+".root --signalviatau -m "+nominalmodes+" -c "+channel+" -s Nominal &> nohuplogs/SpinCorr_Nominal_MC_ttbarsignalplustau_fromDilepton_"+str(i)+"_signalviatau_"+channel+"_"+era+".out \n")


        for i in range(0,10):

            Nominal_scriptfile.write("nohup ./install/bin/load_Analysis -f dy50inf_amcatnlofxfx_"+era+"_"+str(i)+".root -d 11 -m "+nominalmodes+" -c "+channel+" -s Nominal &> nohuplogs/SpinCorr_Nominal_MC_"+"dy50inf_amcatnlofxfx_"+str(i)+"_11_"+channel+"_"+era+".out \n")
            Nominal_scriptfile.write("nohup ./install/bin/load_Analysis -f dy50inf_amcatnlofxfx_"+era+"_"+str(i)+".root -d 13 -m "+nominalmodes+" -c "+channel+" -s Nominal &> nohuplogs/SpinCorr_Nominal_MC_"+"dy50inf_amcatnlofxfx_"+str(i)+"_13_"+channel+"_"+era+".out \n")
            Nominal_scriptfile.write("nohup ./install/bin/load_Analysis -f dy50inf_amcatnlofxfx_"+era+"_"+str(i)+".root -d 15 -m "+nominalmodes+" -c "+channel+" -s Nominal &> nohuplogs/SpinCorr_Nominal_MC_"+"dy50inf_amcatnlofxfx_"+str(i)+"_15_"+channel+"_"+era+".out \n")

        for i in range(0,10):

            Nominal_scriptfile.write("nohup ./install/bin/load_Analysis -f dy50inf_madgraphmlm_"+era+"_"+str(i)+".root -d 11 -m "+nominalmodes+" -c "+channel+" -s Nominal &> nohuplogs/SpinCorr_Nominal_MC_"+"dy50inf_madgraphmlm_"+str(i)+"_11_"+channel+"_"+era+".out \n")
            Nominal_scriptfile.write("nohup ./install/bin/load_Analysis -f dy50inf_madgraphmlm_"+era+"_"+str(i)+".root -d 13 -m "+nominalmodes+" -c "+channel+" -s Nominal &> nohuplogs/SpinCorr_Nominal_MC_"+"dy50inf_madgraphmlm_"+str(i)+"_13_"+channel+"_"+era+".out \n")
            Nominal_scriptfile.write("nohup ./install/bin/load_Analysis -f dy50inf_madgraphmlm_"+era+"_"+str(i)+".root -d 15 -m "+nominalmodes+" -c "+channel+" -s Nominal &> nohuplogs/SpinCorr_Nominal_MC_"+"dy50inf_madgraphmlm_"+str(i)+"_15_"+channel+"_"+era+".out \n")


        for i in range(0,5):

            Nominal_scriptfile.write("nohup ./install/bin/load_Analysis -f ttbarZtollnunu_"+era+"_"+str(i)+".root -m "+nominalmodes+" -c "+channel+" -s Nominal &> nohuplogs/SpinCorr_Nominal_MC_ttbarZtollnunu_"+str(i)+"_"+channel+"_"+era+".out \n")
            Nominal_scriptfile.write("nohup ./install/bin/load_Analysis -f ttbarZtoqq_"+era+"_"+str(i)+".root -m "+nominalmodes+" -c "+channel+" -s Nominal &> nohuplogs/SpinCorr_Nominal_MC_ttbarZtoqq_"+str(i)+"_"+channel+"_"+era+".out \n")


        Nominal_scriptfile.write("\n")


    Nominal_scriptfile.close()



    NominalVariations_scriptfile = open("scripts/LoadAnalysis_SpinCorr_"+"NominalVariations_"+era+".sh","w")
    NominalVariations_scriptfile.write("#!/bin/bash \n")

    for syst in expersysts + JESsysts:

        for channel in channels:
        
            #for sample in samples:

                #NominalVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f "+sample+" -m "+modes+" -c "+channel+" -s "+syst+" &> nohuplogs/SpinCorr_"+syst+"_"+sample.replace(" -d ","_").replace(" --signalviatau","_signalviatau").replace(".root","")+"_"+channel+"_"+era+".out \n")


            for i in range(0,nSignalFiles):

                NominalVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f ttbarsignalplustau_fromDilepton_"+era+"_"+str(i)+".root -m "+modes+" -c "+channel+" -s "+syst+" &> nohuplogs/SpinCorr_"+syst+"_ttbarsignalplustau_fromDilepton_"+str(i)+"_"+channel+"_"+era+".out \n")

            for i in range(0,nSignalFiles):

                NominalVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f ttbarsignalplustau_fromDilepton_"+era+"_"+str(i)+".root --signalviatau -m "+modes+" -c "+channel+" -s "+syst+" &> nohuplogs/SpinCorr_"+syst+"_ttbarsignalplustau_fromDilepton_"+str(i)+"_signalviatau_"+channel+"_"+era+".out \n")


            #for i in range(0,10):

                #NominalVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f dy50inf_amcatnlofxfx_"+era+"_"+str(i)+".root -d 11 -m "+modes+" -c "+channel+" -s "+syst+" &> nohuplogs/SpinCorr_"+syst+"_"+"dy50inf_amcatnlofxfx_"+str(i)+"_11_"+channel+"_"+era+".out \n")
                #NominalVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f dy50inf_amcatnlofxfx_"+era+"_"+str(i)+".root -d 13 -m "+modes+" -c "+channel+" -s "+syst+" &> nohuplogs/SpinCorr_"+syst+"_"+"dy50inf_amcatnlofxfx_"+str(i)+"_13_"+channel+"_"+era+".out \n")
                #NominalVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f dy50inf_amcatnlofxfx_"+era+"_"+str(i)+".root -d 15 -m "+modes+" -c "+channel+" -s "+syst+" &> nohuplogs/SpinCorr_"+syst+"_"+"dy50inf_amcatnlofxfx_"+str(i)+"_15_"+channel+"_"+era+".out \n")

            #for i in range(0,10):

                #NominalVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f dy50inf_madgraphmlm_"+era+"_"+str(i)+".root -d 11 -m "+modes+" -c "+channel+" -s "+syst+" &> nohuplogs/SpinCorr_"+syst+"_"+"dy50inf_madgraphmlm_"+str(i)+"_11_"+channel+"_"+era+".out \n")
                #NominalVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f dy50inf_madgraphmlm_"+era+"_"+str(i)+".root -d 13 -m "+modes+" -c "+channel+" -s "+syst+" &> nohuplogs/SpinCorr_"+syst+"_"+"dy50inf_madgraphmlm_"+str(i)+"_13_"+channel+"_"+era+".out \n")
                #NominalVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f dy50inf_madgraphmlm_"+era+"_"+str(i)+".root -d 15 -m "+modes+" -c "+channel+" -s "+syst+" &> nohuplogs/SpinCorr_"+syst+"_"+"dy50inf_madgraphmlm_"+str(i)+"_15_"+channel+"_"+era+".out \n")


            #for i in range(0,5):

                #NominalVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f ttbarZtollnunu_"+era+"_"+str(i)+".root -m "+modes+" -c "+channel+" -s "+syst+" &> nohuplogs/SpinCorr_"+syst+"_MC_ttbarZtollnunu_"+str(i)+"_"+channel+"_"+era+".out \n")
                #NominalVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f ttbarZtoqq_"+era+"_"+str(i)+".root -m "+modes+" -c "+channel+" -s "+syst+" &> nohuplogs/SpinCorr_"+syst+"_MC_ttbarZtoqq_"+str(i)+"_"+channel+"_"+era+".out \n")


    NominalVariations_scriptfile.close()


    
    TheoryVariations_scriptfile = open("scripts/LoadAnalysis_SpinCorr_"+"TheoryVariations_"+era+".sh","w")
    TheoryVariations_scriptfile.write("#!/bin/bash \n")

    for syst in theorysysts + theorysysts_dedicated.keys():

        for channel in channels:

            if syst in theorysysts_dedicated.keys():

                for i in range(0,nSystFiles):

                    sample = theorysysts_dedicated[syst]+"_"+era+"_"+str(i)+".root"

                    TheoryVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f "+sample+" -m "+modes+" -c "+channel+" &> nohuplogs/SpinCorr_"+syst+"_"+sample.replace(" -d ","_").replace(" --signalviatau","_signalviatau").replace(".root","")+"_"+channel+"_"+era+".out \n")
                    TheoryVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f "+sample+" --signalviatau -m "+modes+" -c "+channel+" &> nohuplogs/SpinCorr_"+syst+"_"+sample.replace(" -d ","_").replace(" --signalviatau","_signalviatau").replace(".root","")+"_signalviatau_"+channel+"_"+era+".out \n")
                        
            else:

                for i in range(0,nSignalFiles):

                    sample = "ttbarsignalplustau_fromDilepton_"+era+"_"+str(i)+".root"

                    TheoryVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f "+sample+" -m "+modes+" -c "+channel+" -s "+syst+" &> nohuplogs/SpinCorr_"+syst+"_"+sample.replace(" -d ","_").replace(" --signalviatau","_signalviatau").replace(".root","")+"_"+channel+"_"+era+".out \n")
                    TheoryVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f "+sample+" --signalviatau -m "+modes+" -c "+channel+" -s "+syst+" &> nohuplogs/SpinCorr_"+syst+"_"+sample.replace(" -d ","_").replace(" --signalviatau","_signalviatau").replace(".root","")+"_signalviatau_"+channel+"_"+era+".out \n")



    for channel in channels:

        # PDF

        for j in range(101):

            for i in range(0,nSignalFiles):

                sample = "ttbarsignalplustau_fromDilepton_"+era+"_"+str(i)+".root"

                TheoryVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f ttbarsignalplustau_fromDilepton_"+era+"_"+str(i)+".root -m "+modes+" -c "+channel+" -s PDF --pdf "+str(j)+" &> nohuplogs/SpinCorr_PDF_"+str(j)+"_"+sample.replace(" -d ","_").replace(" --signalviatau","_signalviatau").replace(".root","")+"_"+channel+"_"+era+".out \n")
                TheoryVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f ttbarsignalplustau_fromDilepton_"+era+"_"+str(i)+".root --signalviatau -m "+modes+" -c "+channel+" -s PDF --pdf "+str(j)+" &> nohuplogs/SpinCorr_PDF_"+str(j)+"_"+sample.replace(" -d ","_").replace(" --signalviatau","_signalviatau").replace(".root","")+"_signalviatau_"+channel+"_"+era+".out \n")
             

        # PS

        for j in range(10,14):

            for i in range(0,nSignalFiles):

                sample = "ttbarsignalplustau_fromDilepton_"+era+"_"+str(i)+".root"
                
                TheoryVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f ttbarsignalplustau_fromDilepton_"+era+"_"+str(i)+".root -m "+modes+" -c "+channel+" -s PSSCALE_WEIGHT --ps "+str(j)+" &> nohuplogs/SpinCorr_PS_"+str(j)+"_"+sample.replace(" -d ","_").replace(" --signalviatau","_signalviatau").replace(".root","")+"_"+channel+"_"+era+".out \n")
                TheoryVariations_scriptfile.write("nohup ./install/bin/load_Analysis -f ttbarsignalplustau_fromDilepton_"+era+"_"+str(i)+".root --signalviatau -m "+modes+" -c "+channel+" -s PSSCALE_WEIGHT --ps "+str(j)+" &> nohuplogs/SpinCorr_PS_"+str(j)+"_"+sample.replace(" -d ","_").replace(" --signalviatau","_signalviatau").replace(".root","")+"_signalviatau_"+channel+"_"+era+".out \n")


    TheoryVariations_scriptfile.close()
