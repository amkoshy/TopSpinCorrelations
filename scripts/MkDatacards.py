#!/usr/bin/env python 
import os
import ROOT

channellist = [
#"combined",
#"ee",
#"emu",
"mumu"
]


distlist = [

"dyScaling_Allh1_step4",
"dyScaling_Allh1_step5",
"dyScaling_Allh1_step6",
#"dyScaling_Allh1_step7",
#"dyScaling_Allh1_step8",

#"dyScaling_TTh1_step4",
#"dyScaling_TTh1_step5",
#"dyScaling_TTh1_step6",
#"dyScaling_TTh1_step7",
#"dyScaling_TTh1_step8",

"dyScaling_Zh1_step4zWindow",
"dyScaling_Zh1_step5zWindow",
"dyScaling_Zh1_step6zWindow",
#"dyScaling_Zh1_step7zWindow",
#"dyScaling_Zh1_step8zWindow",

#"DIMFull",
#"DIMFull_fullSel",
#"HypAntiLeptonBj",
#"HypAntiLeptonBk",
#"HypAntiLeptonBn",
#"HypAntiLeptonBq",
#"HypAntiLeptonBr",
#"HypLeptonBj",
#"HypLeptonBk",
#"HypLeptonBn",
#"HypLeptonBq",
#"HypLeptonBr",

#"HypLLBarBMjj",
#"HypLLBarBMkk",
#"HypLLBarBMnn",
#"HypLLBarBMqq",
#"HypLLBarBMrr",
#"HypLLBarBPjj",
#"HypLLBarBPkk",
#"HypLLBarBPnn",
#"HypLLBarBPqq",
#"HypLLBarBPrr",

#"HypLLBarCMnk",
#"HypLLBarCMnr",
#"HypLLBarCMrk",
#"HypLLBarCPnk",
#"HypLLBarCPnr",
#"HypLLBarCPrk",
#"HypLLBarCkk",
#"HypLLBarCkn",
#"HypLLBarCkr",
#"HypLLBarCnk",
#"HypLLBarCnn",
#"HypLLBarCnr",
#"HypLLBarCrk",
#"HypLLBarCrn",
#"HypLLBarCrr",

#"HypLLBarcHel",
#"HypLLBarcLab",
#"HypLLBarDPhi",

#"HypTTBarMass",
#"HypTTBarRapidity",
#"HypTTBarpT"


]


systlist = [
"Nominal",
#"ELE_RECO",
#"ELE_ID",
#"ELE_SCALESMEARING",
"MUON_ID",
"MUON_ISO",
"MUON_SCALE",
#"L1PREFIRING",
"TRIG",
#"JEREta0",
#"JEREta1",
#"JES",
#"JER",
#"PU",
#"KIN",
#"UNCLUSTERED",
#"BTAG",
#"BTAG_LJET",

]

samplist = [
"ttbarsignal",
"ttbarother",
"zjets",
"diboson",
"wjets",
"singlet",
"ttV"
]


for channel in channellist:

#    if (channel != "combined"): continue

    for dist in distlist:

#        if (dist != "HypAntiLeptonBj"): continue

        #  Make the input root files

        for syst in systlist:


            if syst == "Nominal":
                filename = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_October2021/CMSSW_10_6_27/src/TopAnalysis/Configuration/analysis/diLeptonic/Plots_2017/"+syst+"/"+channel+"/"+dist+"_source.root"

                os.system("rootcp --recreate "+filename+":"+dist+"_data combineTool_DileptonMass/"+dist+"_"+channel+"_input.root:data_obs")

                for samp in samplist:
                    os.system("rootcp "+filename+":"+dist+"_"+samp+" combineTool_DileptonMass/"+dist+"_"+channel+"_input.root:"+samp)


#                os.system("rootcp "+filename+":"+dist+"_allttbar "+dist+"_"+channel+"_input.root:ttbarother")
#                os.system("rootcp "+filename+":"+dist+"_ratio "+dist+"_"+channel+"_input.root:ratio")


            elif syst != "Nominal":

                for updown_i in ["UP", "DOWN"]:

                    if updown_i == "UP":
                        updown_o = "Up"
                    elif updown_i == "DOWN":
                        updown_o = "Down"

                    filename = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_October2021/CMSSW_10_6_27/src/TopAnalysis/Configuration/analysis/diLeptonic/Plots_2017/"+syst+"_"+updown_i+"/"+channel+"/"+dist+"_source.root"

                    for samp in samplist:                
                        os.system("rootcp "+filename+":"+dist+"_"+samp+" combineTool_DileptonMass/"+dist+"_"+channel+"_input.root:"+samp+"_"+syst+updown_o)

#                    os.system("rootcp "+filename+":"+dist+"_allttbar "+dist+"_"+channel+"_input.root:ttbarother_"+syst+updown_o)
#                    os.system("rootcp "+filename+":"+dist+"_ratio "+dist+"_"+channel+"_input.root:ratio_"+syst+updown_o)



        outputfile = ROOT.TFile.Open("combineTool_DileptonMass/"+dist+"_"+channel+"_input.root" ,"READ")


        # Make the .txt data card

        datacard = open("combineTool_DileptonMass/"+dist+"_"+channel+".txt","w")

        datacard.write("imax 1 number of bins \n")
        datacard.write("jmax 6 number of processes minus 1 \n")
        datacard.write("kmax "+str(len(systlist) - 1 + 2)+" number of nuisance parameters (sources of systematic uncertainty) \n")
        datacard.write("--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n")
        datacard.write("bin"+" \t \t "+dist+"_"+channel+"  \n")
        datacard.write("observation "+str(outputfile.Get("data_obs").Integral())+" \n")
        datacard.write("---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n")
        datacard.write("shapes"+"\t"+"*"+"\t"+"*"+"\t"+dist+"_"+channel+"_input.root \t $PROCESS \t $PROCESS_$SYSTEMATIC \n")
        datacard.write("---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n")
        datacard.write("bin"+" \t \t "+dist+"_"+channel+" \t \t "+dist+"_"+channel+" \t \t "+dist+"_"+channel+" \t \t "+dist+"_"+channel+" \t \t "+dist+"_"+channel+" \t \t "+dist+"_"+channel+" \t \t "+dist+"_"+channel+" \n")
        datacard.write("process         zjets         ttbarother         ttbarsignal         diboson         wjets         singlet         ttV \n")
        datacard.write("process        0         1         2         3         4         5         6 \n")
        datacard.write("rate"+"\t"+str(outputfile.Get("zjets").Integral())+"\t"+str(outputfile.Get("ttbarother").Integral())+"\t"+str(outputfile.Get("ttbarsignal").Integral())+"\t"+str(outputfile.Get("diboson").Integral())+"\t"+str(outputfile.Get("wjets").Integral())+"\t"+str(outputfile.Get("singlet").Integral())+"\t"+str(outputfile.Get("ttV").Integral())+" \n")
        datacard.write("---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n")


        for syst in systlist:

            if syst != "Nominal": datacard.write(syst+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")


        datacard.write("DY"+" \t \t "+"lnN"+" \t \t "+"1.2"+" \t \t "+"-"+" \t \t "+"-"+" \t \t "+"-"+" \t \t "+"-"+" \t \t "+"-"+" \t \t "+"-" +" \n")
        datacard.write("BG"+" \t \t "+"lnN"+" \t \t "+"-"+" \t \t "+"1.3"+" \t \t "+"1.3"+" \t \t "+"1.3"+" \t \t "+"1.3"+" \t \t "+"1.3"+" \t \t "+"1.3" +" \n")



#        datacard.write("MASS"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("BSEMILEP"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("UETUNE"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("MATCH"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("PDF_ALPHAS "+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("JES"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("JER"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("PU"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("LEPT"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("TOT_TRIG"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("TOT_SCALE"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("TOT_BFRAG"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("TOT_COLORREC"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("TOT_PDF"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("KIN"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("UNCLUSTERED"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("TOT_BTAG"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("TOT_BTAG_LJET"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("DY"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
#        datacard.write("BG"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")

        datacard.close()


        # Make the run script


        runscript = open("combineTool_DileptonMass/"+dist+"_"+channel+".sh","w")



        runscript.write("#!/bin/sh \n")
        runscript.write("text2workspace.py "+dist+"_"+channel+".txt -m 172.5 -o "+dist+"_"+channel+"_combine.root  \n")
        runscript.write("combine -M AsymptoticLimits "+dist+"_"+channel+"_combine.root  -m 172.5 \n")
        runscript.write("combineTool.py -M Impacts -d "+dist+"_"+channel+"_combine.root  -m 172.5 --rMin -1 --rMax 2 --robustFit 1 --doInitialFit \n")
        runscript.write("combineTool.py -M Impacts -d "+dist+"_"+channel+"_combine.root -m 172.5 --rMin -1 --rMax 2 --robustFit 1 --doFits \n")
        runscript.write("combineTool.py -M Impacts -d "+dist+"_"+channel+"_combine.root -m 172.5 --rMin -1 --rMax 2 --robustFit 1 --output "+dist+"_"+channel+"_impacts.json \n")
        runscript.write("plotImpacts.py -i "+dist+"_"+channel+"_impacts.json -o "+dist+"_"+channel+"  \n")

        runscript.close()


        os.system("mkdir -p nohuplogs")


        outputfile.Close()

        os.system("nohup bash "+dist+"_"+channel+".sh  &> nohuplogs/"+dist+"_"+channel+".out &")
