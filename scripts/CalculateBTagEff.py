#!/usr/bin/env python
import os
import subprocess

import ROOT
from ROOT import TString
from ROOT import TStyle
from ROOT import TSystem
from ROOT import TFile
from ROOT import TCanvas
from ROOT import TH2
from ROOT import TPaveText
from ROOT import TText


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

systs = [
"Nominal",

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

"PDF_ALPHAS_UP", "PDF_ALPHAS_DOWN",
"MESCALE_UP", "MESCALE_DOWN",
"MEFACSCALE_UP", "MEFACSCALE_DOWN",
"MERENSCALE_UP", "MERENSCALE_DOWN",
"BFRAG_UP", "BFRAG_DOWN",
"BFRAG_PETERSON", "BFRAG_CENTRAL",
"BSEMILEP_UP", "BSEMILEP_DOWN",

"PSSCALE_WEIGHT_6_UP","PSSCALE_WEIGHT_6_DOWN",
"PSSCALE_WEIGHT_7_UP","PSSCALE_WEIGHT_7_DOWN",

"MASS_UP","MASS_DOWN",
"UETUNE_UP","UETUNE_DOWN",
"MATCH_UP","MATCH_DOWN",
"ERDON","ERDONRETUNE","GLUONMOVETUNE",
"AMCATNLOFXFX", "MADGRAPHMLM", "POWHEGV2HERWIG",

]



for era in eras:

    print(era)

    for syst in systs:

        print(syst)

        for channel in channels:

            print(channel)

            sample = "ttbarsignalplustau_fromDilepton"

            filename = "selectionRoot_"+era.replace("UL","")+"/"+"BTagEff_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+sample+"_"+era+".root"
            outfilename = "selectionRoot_"+era.replace("UL","")+"/"+"BTagEff_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"+sample+".root"
            outfoldername = "selectionRoot_"+era.replace("UL","")+"/"+"BTagEff_"+era.replace("UL","")+"/"+syst+"/"+channel+"/"

            if not os.path.isfile(filename):
                continue

            file = TFile.Open(filename,"READ")

            beff = file.Get("beff2D").Clone()
            beff.SetDirectory(0)
            ceff = file.Get("ceff2D").Clone()
            ceff.SetDirectory(0)
            leff = file.Get("leff2D").Clone()
            leff.SetDirectory(0)
            bjets = file.Get("bjets2D").Clone()
            bjets.SetDirectory(0)
            cjets = file.Get("cjets2D").Clone()
            cjets.SetDirectory(0)
            ljets = file.Get("ljets2D").Clone()
            ljets.SetDirectory(0)
            btagged = file.Get("bjetsTagged2D").Clone()
            btagged.SetDirectory(0)
            ctagged = file.Get("cjetsTagged2D").Clone()
            ctagged.SetDirectory(0)
            ltagged = file.Get("ljetsTagged2D").Clone()
            ltagged.SetDirectory(0)
            medians = file.Get("medians").Clone()
            medians.SetDirectory(0)

            file.Close()

            outHisto1 = beff.Clone()

            for binX in range(1, beff.GetNbinsX() + 1):
                for binY in range(1, beff.GetNbinsY() + 1):
                    num1 = btagged.GetBinContent(binX, binY)
                    denum1 = bjets.GetBinContent(binX, binY)
                    content1 = num1 / denum1
                    err1 = ROOT.TMath.Sqrt(content1 * (1 - content1)/bjets.GetBinContent(binX, binY))
                    outHisto1.SetBinContent(binX, binY, content1)
                    outHisto1.SetBinError(binX, binY, err1)

            outHisto2 = ceff.Clone()

            for binX in range(1, ceff.GetNbinsX() + 1):
                for binY in range(1, ceff.GetNbinsY() + 1):
                    num2 = ctagged.GetBinContent(binX, binY)
                    denum2 = cjets.GetBinContent(binX, binY)
                    content2 = 0. if denum2==0. else num2 / denum2
                    err2 = 0. if denum2==0. else ROOT.TMath.Sqrt(content2 * (1 - content2)/cjets.GetBinContent(binX, binY))
                    outHisto2.SetBinContent(binX, binY, content2)
                    outHisto2.SetBinError(binX, binY, err2)

            outHisto3 = leff.Clone()

            for binX in range(1, leff.GetNbinsX() + 1):
                for binY in range(1, leff.GetNbinsY() + 1):
                    num3 = ltagged.GetBinContent(binX, binY)
                    denum3 = ljets.GetBinContent(binX, binY)
                    content3 = num3 / denum3
                    err3 = ROOT.TMath.Sqrt(content3 * (1 - content3)/ljets.GetBinContent(binX, binY))
                    outHisto3.SetBinContent(binX, binY, content3)
                    outHisto3.SetBinError(binX, binY, err3)

            if not os.path.exists(outfoldername):
                os.makedirs(outfoldername)

            outfile = TFile.Open(outfilename,"RECREATE")

            outHisto1.Write()
            outHisto2.Write()
            outHisto3.Write()
            bjets.Write()
            btagged.Write()
            cjets.Write()
            ctagged.Write()
            ljets.Write()
            ltagged.Write()
            medians.Write()
            outfile.Close()

