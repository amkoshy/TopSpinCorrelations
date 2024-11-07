from ROOT import TFile, TH1D, TH2D
import os
import argparse



def get_lumi_from_samplename(samplename) :

    # 2016
    if   "2016ULpreVFP"  in samplename : lumi_ = 19500.0
    elif "2016ULpostVFP" in samplename : lumi_ = 16810.0
    elif "2016UL" in samplename : lumi_ = 36310.0
    
    # 2017
    elif "2017UL" in samplename : lumi_ = 41480.0
    
    # 2018
    elif "2018UL" in samplename : lumi_ = 59830.0
    
    # full Run2
    elif "fullRun2UL" in samplename : lumi_ = 137620.0

    return lumi_


def get_dySF_from_samplename(samplename,channel) :

    dySF = 1.0 

    if  "2018" in samplename:
        if    channel == "ee"       : dySF = 1.042
        elif  channel == "emu"      : dySF = 1.044
        elif  channel == "mumu"     : dySF = 1.047
        elif  channel == "combined" : dySF = 1.045
    
    elif "2017" in samplename :
        if    channel == "ee"       : dySF = 1.078
        elif  channel == "emu"      : dySF = 1.091
        elif  channel == "mumu"     : dySF = 1.103
        elif  channel == "combined" : dySF = 1.095

    elif "2016postVFP" in samplename :
        if    channel == "ee"       : dySF = 1.112
        elif  channel == "emu"      : dySF = 1.128
        elif  channel == "mumu"     : dySF = 1.144
        elif  channel == "combined" : dySF = 1.133

    elif "2016preVFP" in samplename :
        if    channel == "ee"       : dySF = 1.081
        elif  channel == "emu"      : dySF = 1.076
        elif  channel == "mumu"     : dySF = 1.071
        elif  channel == "combined" : dySF = 1.074

    elif "fullRun2UL" in samplename :
        if    channel == "ee"       : dySF = 1.067
        elif  channel == "emu"      : dySF = 1.073
        elif  channel == "mumu"     : dySF = 1.079
        elif  channel == "combined" : dySF = 1.075

    return dySF 


def get_cross_section_from_samplename(samplename, channel) :

    topxsec = 830.91

    # alternative generators
    if ("ttbarsignal" in samplename) and ("AMCATNLOFXFX" in samplename) : xsection = topxsec
    elif ("ttbarsignal" in samplename) and ("POWHEGV2HERWIG" in samplename) : xsection = topxsec
    
    # w/o tau signal files
    elif   "ee_ttbarsignalplustau_fromDilepton"   in samplename : xsection = topxsec * 0.10706 * 0.964261576
    elif "emu_ttbarsignalplustau_fromDilepton"  in samplename : xsection = topxsec * 0.10706 * 0.957058875
    elif "mumu_ttbarsignalplustau_fromDilepton" in samplename : xsection = topxsec * 0.10706 * 0.949909976
    
    # w tau signal files
    elif "ee_ttbarsignalviatau_fromDilepton"   in samplename : xsection = topxsec * 0.10706 * 1.029827957
    elif "emu_ttbarsignalviatau_fromDilepton"  in samplename : xsection = topxsec * 0.10706 * 1.026209047
    elif "mumu_ttbarsignalviatau_fromDilepton" in samplename : xsection = topxsec * 0.10706 * 1.022670477

    # backgrounds
    
    # ttbar backgrounds
    elif "bg_fromDilepton" in samplename : xsection = topxsec * 0.10706
    elif "fromLjets"       in samplename : xsection = topxsec * 0.44113
    elif "fromHadronic"    in samplename : xsection = topxsec * 0.45441
    
    # elif ("ttbar"  in samplename) and not("ttbarW" in samplename) and not("ttbarW" in samplename) : xsection = topxsec 
    
    # Single top
    elif ("single" in samplename)    and ("tw" in samplename) : xsection = 35.85 * 0.54559
    elif ("single" in samplename)    and ("_s" in samplename) : xsection = 10.32
    elif ("singletop" in samplename) and ("_t" in samplename) : xsection = 136.02

    # Single antitop
    elif ("singleantitop" in samplename) and ("_t" in samplename) : xsection = 80.95
    
    # VV
    elif "ww" in samplename : xsection = 118.7
    elif "wz" in samplename : xsection = 47.13
    elif "zz" in samplename : xsection = 16.523
    
    # Drell-Yan
    elif "1050" in samplename               : xsection = 22635.1 * get_dySF_from_samplename(samplename, channel)
    elif "50inf_amcatnlofxfx" in samplename : xsection = 0.5 * 3.*2075.14 * get_dySF_from_samplename(samplename, channel)
    elif "50inf_madgraphmlm"  in samplename : xsection = 0.5 * 3.*2075.14 * get_dySF_from_samplename(samplename, channel)
    
    # Smaller backgrounds
    elif "wtolnu"          in samplename : xsection = 61526.7
    elif "ttbarWjetstolnu" in samplename : xsection = 0.2043
    elif "ttbarWjetstoqq"  in samplename : xsection = 0.4062
    elif "ttbarZtollnunu"  in samplename : xsection = 0.2529
    elif "ttbarZtoqq"      in samplename : xsection = 0.5297
    else : xsection = topxsec

    return xsection


def get_lumi_weight(root_fileptr, lumi, xsection) :
    nevents = root_fileptr.Get('hNrOfEvts').GetBinContent(1)
    lumi_wt = (lumi * xsection) / nevents
    return lumi_wt


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


dileptonic_dir = os.getenv("CMSSW_BASE") + "/src/TopAnalysis/Configuration/analysis/diLeptonic"


def main():

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-y', '--era', dest='era',
                            action='store', default='',
                            help='Era (example: 2018UL, 2017UL, 2016ULpreVFP, 2016ULpostVFP')
    parser.add_argument('-c', '--channel', dest='channel',
                            action='store', default='',
                            help='Channel (example: ee, emu, mumu')
    parser.add_argument('-s', '--systematic', dest='systematic',
                            action='store', default='',
                            help='Systematic (example: Nominal, JER_UP, etc.')

    
    opts, opts_unknown = parser.parse_known_args()


    for era in eras:

        if str(opts.era) not in era: continue

        if not os.path.exists("UnfoldingHistos_Lumiweighted_"+era.replace("UL","")): os.makedirs("UnfoldingHistos_Lumiweighted_"+era.replace("UL",""))

        for systematic in nominal + expersysts + theorysysts:

            if str(opts.systematic) not in systematic: continue

            if not os.path.exists("UnfoldingHistos_Lumiweighted_"+era.replace("UL","")+"/"+systematic): os.makedirs("UnfoldingHistos_Lumiweighted_"+era.replace("UL","")+"/"+systematic)

            for channel in channels:

                if str(opts.channel) not in channel: continue

                if not os.path.exists("UnfoldingHistos_Lumiweighted_"+era.replace("UL","")+"/"+systematic+"/"+channel): os.makedirs("UnfoldingHistos_Lumiweighted_"+era.replace("UL","")+"/"+systematic+"/"+channel)

                TUnfoldHistoFileLists_dir = "TUnfoldHistoFileLists_"+era.replace("UL","")
                TUnfoldHistoFileList_filename = "TUnfoldHistoFileList_"+systematic+"_"+channel+".txt"

                TUnfoldFileList = open(TUnfoldHistoFileLists_dir + "/" + TUnfoldHistoFileList_filename, "r")

                filelist = TUnfoldFileList.readlines()

                for filename in filelist:

                    filename = "UnfoldingHistos_"+era.replace("UL","")+"/"+systematic+"/"+channel+"/"+os.path.splitext(os.path.basename(filename))[0]+".root"

                    print("Processing: " + filename)

                    inHistFile = TFile.Open ( "UnfoldingHistos_"+era.replace("UL","")+"/"+systematic+"/"+channel+"/"+os.path.splitext(os.path.basename(filename))[0]+".root" , "READ" )

                    selectionRootFile = TFile.Open ( "selectionRoot_"+era.replace("UL","")+"/"+systematic+"/"+channel+"/"+os.path.splitext(os.path.basename(filename))[0].replace("histosTUnfold_","")+".root" , "READ" )
                    trueLevelNoRenormalisationWeightSum_hist = selectionRootFile.Get("trueLevelNoRenormalisationWeightSum_hist")
                    trueLevelWeightSum_hist = selectionRootFile.Get("trueLevelWeightSum_hist")
                    globalNormalisationFactor = trueLevelNoRenormalisationWeightSum_hist.GetBinContent(1)/trueLevelWeightSum_hist.GetBinContent(1)

                    outHistFile = TFile.Open ( "UnfoldingHistos_Lumiweighted_"+era.replace("UL","")+"/"+systematic+"/"+channel+"/"+os.path.splitext(os.path.basename(filename))[0]+".root" , "RECREATE" )

                    if "run201" in filename and systematic == "Nominal":

                        for key in inHistFile.GetListOfKeys():
                            inHistFile.cd()
                            name = key.GetName()
                            input_hist = inHistFile.Get(name)

                            outHistFile.cd()
                            input_hist.Write()

                    else:

                        for key in inHistFile.GetListOfKeys():
                            inHistFile.cd()
                            name = key.GetName()
                            input_hist = inHistFile.Get(name)
                            Lumi = get_lumi_from_samplename(filename)
                            XSection = get_cross_section_from_samplename(filename, channel)
                            lumiWeight = get_lumi_weight(inHistFile, Lumi, XSection)

                            outHistFile.cd()
                            input_hist.Scale(globalNormalisationFactor)
                            input_hist.Scale(lumiWeight)
                            input_hist.Write()


                    inHistFile.Close()
                    outHistFile.Close()



##############################

if __name__ == "__main__":

    main()

##############################
