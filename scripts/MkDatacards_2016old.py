#!/usr/bin/env python 
import os
import ROOT

channellist = [
"combined",
"ee",
"emu",
"mumu"
]


distlist = [

"DIMFull",
"DIMFull_fullSel",

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
#"MASS",
"BSEMILEP",
#"UETUNE",
#"MATCH",
"PDF_ALPHAS",
"JES",
"JER",
"PU",
"LEPT",
"TOT_TRIG",
#"TOT_SCALE",
"TOT_BFRAG",
#"TOT_COLORREC",
#"TOT_PDF",
"KIN",
"UNCLUSTERED",
"TOT_BTAG",
"TOT_BTAG_LJET",
"DY",
"BG",
#"TOP_PT"
]

samplist = [
"allmc",
"signalmc",
"allttbar",
"zjets",
"diboson",
"wjets",
"singlet",
"ttbarzw"
]


for channel in channellist:

    for dist in distlist:
        
        #  Make the input root files

        for syst in systlist:

            if syst != "Nominal":

                for updown_i in ["UP", "DOWN"]:

                    if updown_i == "UP":
                        updown_o = "Up"
                    elif updown_i == "DOWN":
                        updown_o = "Down"

                    filename = "/home/jthiema/TopSpinCorr/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/Plots/"+syst+"_"+updown_i+"/"+channel+"/"+dist+"_source.root"

                    for samp in samplist:                
                        os.system("rootcp "+filename+":"+dist+"_"+samp+" combineTool_DileptonMass/"+dist+"_"+channel+"_input.root:"+samp+"_"+syst+updown_o)

                    os.system("rootcp "+filename+":"+dist+"_allttbar combineTool_DileptonMass/"+dist+"_"+channel+"_input.root:ttbarother_"+syst+updown_o)


            elif syst == "Nominal":
                filename = "/home/jthiema/TopSpinCorr/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/Plots/"+syst+"/"+channel+"/"+dist+"_source.root"
                for samp in samplist:
                    os.system("rootcp "+filename+":"+dist+"_"+samp+" combineTool_DileptonMass/"+dist+"_"+channel+"_input.root:"+samp)

                os.system("rootcp "+filename+":"+dist+"_data combineTool_DileptonMass/"+dist+"_"+channel+"_input.root:data_obs")

                os.system("rootcp "+filename+":"+dist+"_allttbar combineTool_DileptonMass/"+dist+"_"+channel+"_input.root:ttbarother")

        

        # correct the input root files for missing ttbar other


        file = ROOT.TFile.Open("combineTool_DileptonMass/"+dist+"_"+channel+"_input.root" ,"READ")

        outputfile = ROOT.TFile.Open("combineTool_DileptonMass/"+dist+"_"+channel+"_inputcorrected.root" ,"RECREATE")

        for syst in systlist:

            if syst != "Nominal":

                for updown_o in ["Up", "Down"]:

                    allttbar = file.Get("allttbar"+"_"+syst+updown_o)
                    signalttbar = file.Get("signalmc"+"_"+syst+updown_o)
                    ttbarother = file.Get("ttbarother"+"_"+syst+updown_o)

                    allmc = file.Get("allmc"+"_"+syst+updown_o)
                    zjets = file.Get("zjets"+"_"+syst+updown_o)
                    diboson = file.Get("diboson"+"_"+syst+updown_o)
                    wjets = file.Get("wjets"+"_"+syst+updown_o)
                    singlet = file.Get("singlet"+"_"+syst+updown_o)
                    ttbarzw = file.Get("ttbarzw"+"_"+syst+updown_o)

                    ttbarother.Add(ttbarother, signalttbar, 1, -1)
                    ttbarother.Write()
                    allttbar.Write()
                    signalttbar.Write()
                    allmc.Write()
                    zjets.Write()
                    diboson.Write()
                    wjets.Write()
                    singlet.Write()
                    ttbarzw.Write()

            elif syst == "Nominal":

                allttbar = file.Get("allttbar")
                signalttbar = file.Get("signalmc")
                ttbarother = file.Get("ttbarother")


                allmc = file.Get("allmc")
                zjets = file.Get("zjets")
                diboson = file.Get("diboson")
                wjets = file.Get("wjets")
                singlet = file.Get("singlet")
                ttbarzw = file.Get("ttbarzw")

                data_obs = file.Get("data_obs")

                ttbarother.Add(ttbarother, signalttbar, 1, -1)
                ttbarother.Write()
                allttbar.Write()
                signalttbar.Write()
                allmc.Write()
                zjets.Write()
                diboson.Write()
                wjets.Write()
                singlet.Write()
                ttbarzw.Write()

                data_obs.Write()

        file.Close()
        outputfile.Close()



for channel in channellist:

    for dist in distlist:

        shapesfile = ROOT.TFile.Open("combineTool_DileptonMass/"+dist+"_"+channel+"_inputcorrected.root" ,"READ")

        # Make the .txt data card

        datacard = open("combineTool_DileptonMass/"+dist+"_"+channel+".txt","w")

        datacard.write("imax 1 number of bins \n")
        datacard.write("jmax 6 number of processes minus 1 \n")
        datacard.write("kmax 14 number of nuisance parameters (sources of systematic uncertainty) \n")
        datacard.write("--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n")
        datacard.write("bin"+" \t \t "+dist+"_"+channel+"  \n")
        datacard.write("observation "+str(shapesfile.Get("data_obs").Integral())+" \n")
        datacard.write("---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n")
        datacard.write("shapes"+"\t"+"*"+"\t"+"*"+"\t"+dist+"_"+channel+"_inputcorrected.root \t $PROCESS \t $PROCESS_$SYSTEMATIC \n")
        datacard.write("---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n")
        datacard.write("bin"+" \t \t "+dist+"_"+channel+" \t \t "+dist+"_"+channel+" \t \t "+dist+"_"+channel+" \t \t "+dist+"_"+channel+" \t \t "+dist+"_"+channel+" \t \t "+dist+"_"+channel+" \t \t "+dist+"_"+channel+" \n")
        datacard.write("process         zjets         ttbarother         signalmc         diboson         wjets         singlet         ttbarzw \n")
        datacard.write("process        0         1         2         3         4         5         6 \n")
        datacard.write("rate"+"\t"+str(shapesfile.Get("zjets").Integral())+"\t"+str(shapesfile.Get("ttbarother").Integral())+"\t"+str(shapesfile.Get("signalmc").Integral())+"\t"+str(shapesfile.Get("diboson").Integral())+"\t"+str(shapesfile.Get("wjets").Integral())+"\t"+str(shapesfile.Get("singlet").Integral())+"\t"+str(shapesfile.Get("ttbarzw").Integral())+" \n")
        datacard.write("---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n")


        #datacard.write("MASS"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        datacard.write("BSEMILEP"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        #datacard.write("UETUNE"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        #datacard.write("MATCH"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        datacard.write("PDF_ALPHAS "+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        datacard.write("JES"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        datacard.write("JER"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        datacard.write("PU"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        datacard.write("LEPT"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        datacard.write("TOT_TRIG"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        #datacard.write("TOT_SCALE"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        datacard.write("TOT_BFRAG"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        #datacard.write("TOT_COLORREC"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        #datacard.write("TOT_PDF"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        datacard.write("KIN"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        datacard.write("UNCLUSTERED"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        datacard.write("TOT_BTAG"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        datacard.write("TOT_BTAG_LJET"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        datacard.write("DY"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")
        datacard.write("BG"+" \t \t "+"shape"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1"+" \t \t "+"1" +" \n")

        datacard.close()


        # Make the run script


        runscript = open("combineTool_DileptonMass/"+dist+"_"+channel+".sh","w")



        runscript.write("#!/bin/sh \n")
        runscript.write("text2workspace.py "+dist+"_"+channel+".txt -m 172.5 -o "+dist+"_"+channel+"_combine.root  \n")
        runscript.write("combine -M AsymptoticLimits "+dist+"_"+channel+"_combine.root  -m 172.5 \n")
        runscript.write("combineTool.py -M Impacts -d "+dist+"_"+channel+"_combine.root  -m 172.5 --rMin -1 --rMax 2 --robustFit 1 --doInitialFit --robustHesse 1 \n")
        runscript.write("combineTool.py -M Impacts -d "+dist+"_"+channel+"_combine.root -m 172.5 --rMin -1 --rMax 2 --robustFit 1 --doFits --robustHesse 1 \n")
        runscript.write("combineTool.py -M Impacts -d "+dist+"_"+channel+"_combine.root -m 172.5 --rMin -1 --rMax 2 --robustFit 1 --robustHesse 1 --output "+dist+"_"+channel+"_impacts.json \n")
        runscript.write("plotImpacts.py -i "+dist+"_"+channel+"_impacts.json -o "+dist+"_"+channel+"  \n")

        runscript.close()


        os.system("mkdir -p nohuplogs")


        shapesfile.Close()

        os.system("nohup bash "+dist+"_"+channel+".sh  &> nohuplogs/"+dist+"_"+channel+".out &")
