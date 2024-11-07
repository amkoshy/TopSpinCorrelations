import os
import sys
import uproot
from ROOT import TFile, TObject, TH1D, TH2D, TUnfoldBinning, TUnfoldBinningXML, cout
import hist
import numpy as np
import argparse
import awkward as ak
from uncertainties import ufloat, unumpy

def remove_underoverflow(histogram, doGen = True , doReco = True):

    if len(histogram.shape) == 1:
        
        x_range = len(histogram[:]) - 1

        histogram[1] = histogram[0] + histogram[1]
        histogram[0] = histogram[0] - histogram[0]
    
        histogram[x_range - 1] = histogram[x_range] + histogram[x_range - 1]
        histogram[x_range] = histogram[x_range] - histogram[x_range]
    
    elif len(histogram.shape) == 2:
        
        x_range = len(histogram[:,0]) - 1 
        y_range = len(histogram[0,:]) - 1
        
        if doGen:
            histogram[1,:] = histogram[0,:] + histogram[1,:]
            histogram[0,:] = histogram[0,:] - histogram[0,:]
        
        if doReco:
            histogram[:,1] = histogram[:,0] + histogram[:,1]
            histogram[:,0] = histogram[:,0] - histogram[:,0]
    
        if doGen:
            histogram[x_range - 1,:] = histogram[x_range,:] + histogram[x_range - 1,:]
            histogram[x_range,:] = histogram[x_range,:] - histogram[x_range,:]

        if doReco:
            histogram[:,y_range - 1] = histogram[:,y_range] + histogram[:,y_range - 1]
            histogram[:,y_range] = histogram[:,y_range] - histogram[:,y_range]
        
        
    elif len(histogram.shape) == 3:
        
        x_range = len(histogram[:,0,0]) - 1
        y_range = len(histogram[0,:,0]) - 1
        z_range = len(histogram[0,0,:]) - 1

        histogram[1,:,:] = histogram[0,:,:] + histogram[1,:,:]
        histogram[0,:,:] = histogram[0,:,:] - histogram[0,:,:]
        histogram[:,1,:] = histogram[:,0,:] + histogram[:,1,:]
        histogram[:,0,:] = histogram[:,0,:] - histogram[:,0,:]
        histogram[:,:,1] = histogram[:,:,0] + histogram[:,:,1]
        histogram[:,:,0] = histogram[:,:,0] - histogram[:,:,0]
    
        histogram[x_range - 1,:,:] = histogram[x_range,:,:] + histogram[x_range - 1,:,:]
        histogram[x_range,:,:] = histogram[x_range,:,:] - histogram[x_range,:,:]
        histogram[:,y_range - 1,:] = histogram[:,y_range,:] + histogram[:,y_range - 1,:]
        histogram[:,y_range,:] = histogram[:,y_range,:] - histogram[:,y_range,:]
        histogram[:,:,z_range - 1] = histogram[:,:,z_range] + histogram[:,:,z_range - 1]
        histogram[:,:,z_range] = histogram[:,:,z_range] - histogram[:,:,z_range]
        
    elif len(histogram.shape) == 4:

        w_range = len(histogram[:,0,0,0]) - 1
        x_range = len(histogram[0,:,0,0]) - 1
        y_range = len(histogram[0,0,:,0]) - 1
        z_range = len(histogram[0,0,0,:]) - 1

        if doGen:
            histogram[1,:,:,:] = histogram[0,:,:,:] + histogram[1,:,:,:]
            histogram[0,:,:,:] = histogram[0,:,:,:] - histogram[0,:,:,:]
            histogram[:,1,:,:] = histogram[:,0,:,:] + histogram[:,1,:,:]
            histogram[:,0,:,:] = histogram[:,0,:,:] - histogram[:,0,:,:]
        
        if doReco:
            histogram[:,:,1,:] = histogram[:,:,0,:] + histogram[:,:,1,:]
            histogram[:,:,0,:] = histogram[:,:,0,:] - histogram[:,:,0,:]
            histogram[:,:,:,1] = histogram[:,:,:,0] + histogram[:,:,:,1]
            histogram[:,:,:,0] = histogram[:,:,:,0] - histogram[:,:,:,0]

        if doGen:
            histogram[w_range - 1,:,:,:] = histogram[w_range,:,:,:] + histogram[w_range - 1,:,:,:]
            histogram[w_range,:,:,:] = histogram[w_range,:,:,:] - histogram[w_range,:,:,:]
            histogram[:,x_range - 1,:,:] = histogram[:,x_range,:,:] + histogram[:,x_range - 1,:,:]
            histogram[:,x_range,:,:] = histogram[:,x_range,:,:] - histogram[:,x_range,:,:]

        if doReco:
            histogram[:,:,y_range - 1,:] = histogram[:,:,y_range,:] + histogram[:,:,y_range - 1,:]
            histogram[:,:,y_range,:] = histogram[:,:,y_range,:] - histogram[:,:,y_range,:]
            histogram[:,:,:,z_range - 1] = histogram[:,:,:,z_range] + histogram[:,:,:,z_range - 1]
            histogram[:,:,:,z_range] = histogram[:,:,:,z_range] - histogram[:,:,:,z_range]
    
    return histogram


def compute_reco(var_name, reco_axes, step8_data, step8_weights, symmetrize = False):

    reco = hist.Hist(*reco_axes, storage=hist.storage.Weight())

    if symmetrize == False:
        reco.fill(*step8_data, weight=step8_weights)
    elif symmetrize == True:
        reco.fill(*[1.0*np.array(step8_data[0])] + [data for data in step8_data[1:] if step8_data[1:]], weight=(0.5*np.array(step8_weights)))
        reco.fill(*[-1.0*np.array(step8_data[0])] + [data for data in step8_data[1:] if step8_data[1:]], weight=(0.5*np.array(step8_weights)))
        reco.variances(flow=True)[...] = 2.0 * reco.variances(flow=True)[...]

    remove_underoverflow(reco.values(flow=True))
    remove_underoverflow(reco.variances(flow=True))

    return reco


def Convert_to_TH1(prefix, var, hist1D):


    TH1 = TH1D(prefix+"_"+var, 
               prefix+"_"+var, 
               hist1D.values(flow=False).shape[0], 
               np.array(GetBinning(hist1D.axes)[0]),
               #hist1D.axes[0][0][0], 
               #hist1D.axes[0][len(hist1D.axes[0])-1][1]
    )
    TH1.Sumw2()

    for i in range(0, hist1D.values(flow=False).shape[0] + 1):
        TH1.SetBinContent(i,hist1D.values(flow=True)[i])
        TH1.SetBinError(i,np.sqrt(hist1D.variances(flow=True)[i]))                     

    return TH1


def Unwrap_to_TH1(prefix, var1D, var2D, hist2D):
    
    TH1 = TH1D(prefix+"_"+var1D+"_vs_"+var2D, 
                prefix+"_"+var1D+"_vs_"+var2D, 
                (hist2D.values(flow=False).shape[0])*(hist2D.values(flow=False).shape[1]), 
                hist2D.axes[0][0][0], 
                hist2D.axes[0][len(hist2D.axes[0])-1][1]
                )
    
    TH1.Sumw2()

    i_range = hist2D.values(flow=False).shape[0]+1
    j_range = hist2D.values(flow=False).shape[1]+1
                        
    for i in range(1,i_range):
        for j in range(1,j_range):
            TH1.SetBinContent((j-1)*(i_range-1)+i, hist2D.values(flow=True)[i,j])
            TH1.SetBinError((j-1)*(i_range-1)+i, np.sqrt(hist2D.variances(flow=True)[i,j]))
    
    return TH1


def BinFinely(Original_Binning, nFineBins):
    Fine_Binning = sum([[Original_Binning[j]+(Original_Binning[j+1]-Original_Binning[j])*(i/nFineBins) for i in range (0,nFineBins)] for j in range(0,len(Original_Binning)-1)] + [[Original_Binning[len(Original_Binning)-1]]],[])
    return Fine_Binning


def GetBinning(axeslist):
    return [[axeslist[i][j][0] for j in range(0,len(axeslist[i]))] + [axeslist[i][len(axeslist[i])-1][1]] for i in range(0,len(axeslist))]


# To run, example:
# python scripts/makeTUnfoldHisto.py -y 2016ULpostVFP -s Nominal -c emu

def main():

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-y', '--era', dest='era',
                            action='store', default='',
                            help='Era (example: 2018UL, 2017UL, 2016ULpreVFP, 2016ULpostVFP')
    parser.add_argument('-c', '--channel', dest='channel',
                            action='store', default='',
                            help='Channel (example: ee, emu, mumu')
    parser.add_argument('-s', '--systematic', dest='systematic',
                            action='store', default='Nominal',
                            help='Systematic (example: Nominal, JER_UP, etc.')
    parser.add_argument('-f', '--file', dest='file',
                            action='store', default='',
                            help='File (example:ttbarsignalplustau_fromDilepton')

    opts, opts_unknown = parser.parse_known_args()

    era = opts.era
    systematic = opts.systematic
    channel = opts.channel
    filepattern = opts.file

    if not os.path.exists("selectionRoot_"+era.replace("UL","")+"/"+systematic+"/"+channel): os.makedirs("selectionRoot_"+era.replace("UL","")+"/"+systematic+"/"+channel)

    dileptonic_dir = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic"
    TUnfoldFileLists_dir = "TUnfoldFileLists_"+era.replace("UL","")
    TUnfoldFileList_filename = "TUnfoldFileList_"+systematic+"_"+channel+".txt"

    TUnfoldFileList = open(dileptonic_dir + "/" + TUnfoldFileLists_dir + "/" + TUnfoldFileList_filename, "r")

    filelist = TUnfoldFileList.readlines()

    for file in filelist:

        print("Processing: " + file)

        outHistFile = TFile.Open ( "selectionRoot_"+era.replace("UL","")+"/"+systematic+"/"+channel+"/"+os.path.splitext(os.path.basename(file))[0]+".root" , "UPDATE" )
        outHistFile.cd()

        trueLevelNoRenormalisationWeightSum_hist = outHistFile.Get("trueLevelNoRenormalisationWeightSum_hist")
        trueLevelWeightSum_hist = outHistFile.Get("trueLevelWeightSum_hist")
        globalNormalisationFactor = trueLevelNoRenormalisationWeightSum_hist.GetBinContent(1)/trueLevelWeightSum_hist.GetBinContent(1)

        step8tree = uproot.open(file+':ttBar_treeVariables_step8')

        tree_vars = [

            "ll_cHel","ll_cLab","llbar_delta_phi","llbar_delta_eta",
            "b1k","b2k","b1r","b2r","b1n","b2n","b1j","b2j","b1q","b2q",
            "c_kk","c_rr","c_nn","c_kj","c_rq",
            "c_rk","c_kr","c_nr","c_rn","c_nk","c_kn","c_rj","c_jr",

            "eventWeight",

            "ttbar_mass", "top_scatteringangle_ttbarframe", "top_pt", "n_extraJets_iso08"

        ]


        step8 = ak.Array(
            ak.zip( dict( (tree_var, step8tree[tree_var].array()) for tree_var in tree_vars ) )
        )

        vars1D_dict = {

            "LLBarcHel" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["ll_cHel"],
                "can_symmetrize" : True,
            },
            "LLBarcLab" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["ll_cLab"],
                "can_symmetrize" : True,
            },
            "LLBarDPhi" : {
                "bin_edge_low" : 0.0, 
                "bin_edge_high" : np.pi, 
                "n_reco_bins" : 12, 
                "val_8" : step8["llbar_delta_phi"],
                "can_symmetrize" : False,
            },
            "LLBarDEta" : {
                "bin_edge_low" : 0.0, 
                "bin_edge_high" : 5.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["llbar_delta_eta"],
                "can_symmetrize" : False,
            },

            "AntiLeptonBk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b1k"],
                "can_symmetrize" : True,
            },
            "LeptonBk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b2k"],
                "can_symmetrize" : True,
            },
            "AntiLeptonBr" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b1r"],
                "can_symmetrize" : True,
            },
            "LeptonBr" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b2r"],
                "can_symmetrize" : True,
            },
            "AntiLeptonBn" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b1n"],
                "can_symmetrize" : True,
            },
            "LeptonBn" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b2n"],
                "can_symmetrize" : True,
            },
            "AntiLeptonBj" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b1j"],
                "can_symmetrize" : True,
            },
            "LeptonBj" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b2j"],
                "can_symmetrize" : True,
            },
            "AntiLeptonBq" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b1q"],
                "can_symmetrize" : True,
            },
            "LeptonBq" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b2q"],
                "can_symmetrize" : True,
            },

            "LLBarCkk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_kk"],
                "can_symmetrize" : True,
            },
            "LLBarCrr" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_rr"],
                "can_symmetrize" : True,
            },
            "LLBarCnn" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_nn"],
                "can_symmetrize" : True,
            },
            "LLBarCkj" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_kj"],
                "can_symmetrize" : True,
            },
            "LLBarCrq" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_rq"],
                "can_symmetrize" : True,
            },

            "LLBarCPrk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_rk"] + step8["c_kr"],
                "can_symmetrize" : True,
            },
            "LLBarCMrk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_rk"] - step8["c_kr"],
                "can_symmetrize" : True,
            },
            "LLBarCPnr" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_nr"] + step8["c_kn"],
                "can_symmetrize" : True,
            },
            "LLBarCMnr" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_nr"] - step8["c_kn"],
                "can_symmetrize" : True,
            },
            "LLBarCPnk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_nk"] + step8["c_kn"],
                "can_symmetrize" : True,
            },
            "LLBarCMnk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_nk"] - step8["c_kn"],
                "can_symmetrize" : True,
            },
            "LLBarCPrj" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_rj"] + step8["c_jr"],
                "can_symmetrize" : True,
            },
            "LLBarCMrj" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_rj"] - step8["c_jr"],
                "can_symmetrize" : True,
            },

            "LLBarChan" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : +step8["c_kk"] - step8["c_rr"] - step8["c_nn"],
                "can_symmetrize" : True,
            },
            "LLBarCsca" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : -step8["c_kk"] + step8["c_rr"] - step8["c_nn"],
                "can_symmetrize" : True,
            },
            "LLBarCtra" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : -step8["c_kk"] - step8["c_rr"] + step8["c_nn"],
                "can_symmetrize" : True,
            },
            "LLBarCkjL" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : -step8["c_kj"] - step8["c_rr"] - step8["c_nn"],
                "can_symmetrize" : True,
            },
            "LLBarCrqL" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : -step8["c_kk"] - step8["c_rq"] - step8["c_nn"],
                "can_symmetrize" : True,
            },

            "LLBarCrkP" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : -step8["c_rk"] - step8["c_kr"] - step8["c_nn"],
                "can_symmetrize" : True,
            },
            "LLBarCrkM" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : -step8["c_rk"] + step8["c_kr"] - step8["c_nn"],
                "can_symmetrize" : True,
            },
            "LLBarCnrP" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : -step8["c_nr"] - step8["c_rn"] - step8["c_kk"],
                "can_symmetrize" : True,
            },
            "LLBarCnrM" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : -step8["c_nr"] + step8["c_rn"] - step8["c_kk"],
                "can_symmetrize" : True,
            },
            "LLBarCnkP" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : -step8["c_nk"] - step8["c_kn"] - step8["c_rr"],
                "can_symmetrize" : True,
            },
            "LLBarCnkM" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : -step8["c_nk"] + step8["c_kn"] - step8["c_rr"],
                "can_symmetrize" : True,
            },

        }


        binning_ttbar_mass = [300., 450., 600., 800., 2000.]
        binning_top_scatteringangle_ttbarframe = [-1.0, -0.5, 0.0, +0.5, +1.0]
        binning_top_pt = [0., 80., 150., 250., 550.]
        binning_n_extraJets_iso08 = [-0.5, 0.5, 1.5, 2.5, 3.5]


        vars2D_dict = {

            "TTBarMass" : {
                "bin_edge_low" : 300., 
                "bin_edge_high" : 1500., 
                "reco_binning" : BinFinely(binning_ttbar_mass, 1),
                "val_8" : step8["ttbar_mass"],
                "can_symmetrize" : False,
            },

            "ScatteringAngle_TTBarFrame" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "reco_binning" : BinFinely(binning_top_scatteringangle_ttbarframe, 1),
                "val_8" : step8["top_scatteringangle_ttbarframe"],
                "can_symmetrize" : False,
            },

            "ToppT" : {
                "bin_edge_low" : 0., 
                "bin_edge_high" : 550., 
                "reco_binning" : BinFinely(binning_top_pt, 1),
                "val_8" : step8["top_pt"],
                "can_symmetrize" : False,
            },

            "ExtraJets" : {
                "bin_edge_low" : -0.5, 
                "bin_edge_high" : 3.5, 
                "reco_binning" : binning_n_extraJets_iso08,
                "val_8" : step8["n_extraJets_iso08"],
                "can_symmetrize" : False,
            },

        }

        nFineBins_1D = 2
        nFineBins_2D = 1

        do_symmetrize = True

        # 1D

        prefix = ""

        for var1D in vars1D_dict.keys():

            if(do_symmetrize*vars1D_dict[var1D]["can_symmetrize"] == True): prefix = "HypSym"
            else: prefix = "Hyp"

            reco_axes = [hist.axis.Regular(nFineBins_1D*vars1D_dict[var1D]["n_reco_bins"], vars1D_dict[var1D]["bin_edge_low"], vars1D_dict[var1D]["bin_edge_high"], name="reco_"+var1D, label="reco_"+var1D, underflow=True, overflow=True)]

            reco = compute_reco(var1D, reco_axes, [vars1D_dict[var1D]["val_8"]], [step8['eventWeight']], do_symmetrize*vars1D_dict[var1D]["can_symmetrize"])
            hreco = Convert_to_TH1(prefix, var1D, reco)
            hreco.Scale(globalNormalisationFactor)
            hreco.Write(hreco.GetName(),TObject.kOverwrite)

        # 2D

        for var2D in vars2D_dict.keys():

            reco_axes_2D = [hist.axis.Variable(BinFinely(vars2D_dict[var2D]["reco_binning"],nFineBins_2D), name="reco2D_"+var2D, label="reco2D_"+var2D, underflow=True, overflow=True)]

            for var1D in vars1D_dict.keys():

                if(do_symmetrize*vars1D_dict[var1D]["can_symmetrize"] == True): prefix = "HypSym"
                else: prefix = "Hyp"

                reco_axes = [hist.axis.Regular(nFineBins_2D*vars1D_dict[var1D]["n_reco_bins"], vars1D_dict[var1D]["bin_edge_low"], vars1D_dict[var1D]["bin_edge_high"], name="reco_"+var1D, label="reco_"+var1D, underflow=True, overflow=True)]

                reco = compute_reco(var1D+"_"+var2D, reco_axes + reco_axes_2D, [vars1D_dict[var1D]["val_8"], vars2D_dict[var2D]["val_8"]], [step8['eventWeight']], do_symmetrize*vars1D_dict[var1D]["can_symmetrize"])
                hreco = Unwrap_to_TH1(prefix, var1D, var2D, reco)
                hreco.Scale(globalNormalisationFactor)
                hreco.Write(hreco.GetName(),TObject.kOverwrite)


            if var2D == "ExtraJets":
                reco_axes_2D = [hist.axis.Variable(BinFinely(vars2D_dict[var2D]["reco_binning"],1), name="reco2D_"+var2D, label="reco2D_"+var2D, underflow=True, overflow=True)]
            else:
                reco_axes_2D = [hist.axis.Variable(BinFinely(vars2D_dict[var2D]["reco_binning"],4), name="reco2D_"+var2D, label="reco2D_"+var2D, underflow=True, overflow=True)]

            if(do_symmetrize*vars2D_dict[var2D]["can_symmetrize"] == True): prefix = "HypSym"
            else: prefix = "Hyp"

            reco = compute_reco(var2D, reco_axes_2D, [vars2D_dict[var2D]["val_8"]], [step8['eventWeight']], do_symmetrize*vars2D_dict[var2D]["can_symmetrize"])
            hreco = Convert_to_TH1(prefix, var2D, reco)
            hreco.Scale(globalNormalisationFactor)
            hreco.Write(hreco.GetName(),TObject.kOverwrite)


        outHistFile.Close()


##############################

if __name__ == "__main__":

    main()

##############################  
