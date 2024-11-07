import os
import uproot
from ROOT import TFile, TH1D, TH2D
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


def compute_reco(var_name, reco_axes, step8_data, step8_weights):
    
    reco = hist.Hist(*reco_axes, storage=hist.storage.Weight())
    
    reco.fill(*step8_data, weight=step8_weights)
    
    remove_underoverflow(reco.values(flow=True))    
    remove_underoverflow(reco.variances(flow=True))
    
    return reco


def Convert_to_TH1(prefix, PE, var, hist1D):


    TH1 = TH1D(prefix+"_"+var+str(PE), 
               prefix+"_"+var+str(PE), 
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


def Unwrap_to_TH1(prefix, PE, var1D, var2D, hist2D):
    
    TH1 = TH1D(prefix+"_"+var1D+"_"+var2D+str(PE), 
                prefix+"_"+var1D+"_"+var2D+str(PE), 
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
    parser.add_argument('-n', '--nPE', dest='nPE',
                            action='store', default='1',
                            help='Number of Pseudoexperiments (example: 50')
    parser.add_argument('-r', '--RNGseed', dest='RNGseed',
                            action='store', default='12345',
                            help='Seed for RNG (example: 12345')
    parser.add_argument('-f', '--file', dest='file',
                            action='store', default='',
                            help='File (example:ttbarsignalplustau_fromDilepton')

    opts, opts_unknown = parser.parse_known_args()

    era = opts.era
    systematic = opts.systematic
    channel = opts.channel
    nPE = int(opts.nPE)
    RNGseed = opts.RNGseed
    filepattern = opts.file

    if not os.path.exists("UnfoldingHistos_"+era.replace("UL","")+"/"+systematic+"/"+channel): os.makedirs("UnfoldingHistos_"+era.replace("UL","")+"/"+systematic+"/"+channel)

    dileptonic_dir = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_January2023/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic"
    TUnfoldFileLists_dir = "TUnfoldFileLists_"+era.replace("UL","")
    TUnfoldFileList_filename = "TUnfoldFileList_"+systematic+"_"+channel+".txt"

    TUnfoldFileList = open(dileptonic_dir + "/" + TUnfoldFileLists_dir + "/" + TUnfoldFileList_filename, "r")

    filelist = TUnfoldFileList.readlines()

    for file in filelist:
        print("Processing: " + file)

        if "run201" not in file: continue

        inTreeFile  = TFile.Open ( "spinCorrInput_"+era.replace("UL","")+"/"+systematic+"/"+channel+"/"+os.path.splitext(os.path.basename(file))[0]+".root" , "READ" )

        outHistFile = TFile.Open ( "UnfoldingHistos_"+era.replace("UL","")+"/"+systematic+"/"+channel+"/histosTUnfold_"+os.path.splitext(os.path.basename(file))[0]+".root" , "UPDATE" )
        outHistFile.cd()

        wtdEvts = inTreeFile.Get("weightedEvents")
        wtdEvts.Write("hNrOfEvts")

        inTreeFile.Close()

        step8tree = uproot.open(file+':ttBar_treeVariables_step8')

        isTTBarSignal = "ttbarsignal" in file


        tree_vars = [

            "ll_cHel","ll_cLab","llbar_delta_phi","llbar_delta_eta",
            "b1k","b2k","b1r","b2r","b1n","b2n","b1j","b2j","b1q","b2q",
            "c_kk","c_rr","c_nn","c_kj","c_rq","c_rq",   
            "c_rk","c_kr","c_nr","c_rn","c_nk","c_kn","c_rj","c_jr",

            "eventWeight", "trueLevelWeight",

            "ttbar_mass", "top_scatteringangle_ttbarframe",

        ]


        step8 = ak.Array(
            ak.zip( dict( (tree_var, step8tree[tree_var].array()) for tree_var in tree_vars ) )
        )


        vars_dict = {

            "ll_cHel" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["ll_cHel"],
                     },
            "ll_cLab" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["ll_cLab"],
                  },
            "llbar_delta_phi" : {
                "bin_edge_low" : 0.0, 
                "bin_edge_high" : np.pi, 
                "n_reco_bins" : 12, 
                "val_8" : step8["llbar_delta_phi"],
                  },
            "llbar_delta_eta" : {
                "bin_edge_low" : 0.0, 
                "bin_edge_high" : 5.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["llbar_delta_eta"],
                  },

            "b1k" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b1k"],
                  },
            "b2k" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b2k"],
                  },
            "b1r" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b1r"],
                  },
            "b2r" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b2r"],
                  },
            "b1n" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b1n"],
                  },
            "b2n" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b2n"],
                  },
            "b1j" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b1j"],
                  },
            "b2j" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b2j"],
                  },
            "b1q" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b1q"],
                  },
            "b2q" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["b2q"],
                  },

            "c_kk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_kk"],
                  },
            "c_rr" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_rr"],
                  },
            "c_nn" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_nn"],
                  },
            "c_kj" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_kj"],
                  },
            "c_rq" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_rq"],
                  },

            "c_Prk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_rk"] + step8["c_kr"],
                  },
            "c_Mrk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_rk"] - step8["c_kr"],
                  },
            "c_Pnr" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_nr"] + step8["c_kn"],
                  },
            "c_Mnr" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_nr"] - step8["c_kn"],
                  },
            "c_Pnk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_nk"] + step8["c_kn"],
                  },
            "c_Mnk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_nk"] - step8["c_kn"],
                  },
            "c_Prj" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_rj"] + step8["c_jr"],
                  },
            "c_Mrj" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_reco_bins" : 12, 
                "val_8" : step8["c_rj"] - step8["c_jr"],
                  },

            "c_han" : {
                "bin_edge_low" : -1.0,
                "bin_edge_high" : 1.0,
                "n_reco_bins" : 12,
                "val_8" : +step8["c_kk"] - step8["c_rr"] - step8["c_nn"],
            },
            "c_sca" : {
                "bin_edge_low" : -1.0,
                "bin_edge_high" : 1.0,
                "n_reco_bins" : 12,
                "val_8" : -step8["c_kk"] + step8["c_rr"] - step8["c_nn"],
            },
            "c_tra" : {
                "bin_edge_low" : -1.0,
                "bin_edge_high" : 1.0,
                "n_reco_bins" : 12,
                "val_8" : -step8["c_kk"] - step8["c_rr"] + step8["c_nn"],
            },
            "c_kjL" : {
                "bin_edge_low" : -1.0,
                "bin_edge_high" : 1.0,
                "n_reco_bins" : 12,
                "val_8" : -step8["c_kj"] - step8["c_rr"] - step8["c_nn"],
            },
            "c_rqL" : {
                "bin_edge_low" : -1.0,
                "bin_edge_high" : 1.0,
                "n_reco_bins" : 12,
                "val_8" : -step8["c_kk"] - step8["c_rq"] - step8["c_nn"],
            },

            "c_rkP" : {
                "bin_edge_low" : -1.0,
                "bin_edge_high" : 1.0,
                "n_reco_bins" : 12,
                "val_8" : -step8["c_rk"] - step8["c_kr"] - step8["c_nn"],
            },
            "c_rkM" : {
                "bin_edge_low" : -1.0,
                "bin_edge_high" : 1.0,
                "n_reco_bins" : 12,
                "val_8" : -step8["c_rk"] + step8["c_kr"] - step8["c_nn"],
            },
            "c_nrP" : {
                "bin_edge_low" : -1.0,
                "bin_edge_high" : 1.0,
                "n_reco_bins" : 12,
                "val_8" : -step8["c_nr"] - step8["c_rn"] - step8["c_kk"],
            },
            "c_nrM" : {
                "bin_edge_low" : -1.0,
                "bin_edge_high" : 1.0,
                "n_reco_bins" : 12,
                "val_8" : -step8["c_nr"] + step8["c_rn"] - step8["c_kk"],
            },
            "c_nkP" : {
                "bin_edge_low" : -1.0,
                "bin_edge_high" : 1.0,
                "n_reco_bins" : 12,
                "val_8" : -step8["c_nk"] - step8["c_kn"] - step8["c_rr"],
            },
            "c_nkM" : {
                "bin_edge_low" : -1.0,
                "bin_edge_high" : 1.0,
                "n_reco_bins" : 12,
                "val_8" : -step8["c_nk"] + step8["c_kn"] - step8["c_rr"],
            },

        }


        binning_ttbar_mass = [300, 450, 600, 800, 2000]
        binning_top_scatteringangle_ttbarframe = [-1.0, -0.5, 0.0, +0.5, +1.0]

        vars2D_dict = {

            "ttbar_mass" : {
                "reco_binning" : BinFinely(binning_ttbar_mass, 2),
                "val_8" : step8["ttbar_mass"],
                        },

            "top_scatteringangle_ttbarframe" : {
                "reco_binning" : BinFinely(binning_top_scatteringangle_ttbarframe, 2),
                "val_8" : step8["top_scatteringangle_ttbarframe"],
                                            },

        }


        rand = np.random.default_rng(RNGseed)

        for PE in range(0, nPE):

            PE_samples = rand.poisson(1.0, len(step8['eventWeight']))

            PE_weights = [np.repeat(step8['eventWeight'], PE_samples, axis=0)]

            for var1D in vars_dict.keys():

                reco_axes =[hist.axis.Regular(
                    vars_dict[var1D]["n_reco_bins"], 
                    vars_dict[var1D]["bin_edge_low"], 
                    vars_dict[var1D]["bin_edge_high"], 
                    name="hrecoBootstrap_"+var1D+str(PE), 
                    label="hrecoBootstrap_"+var1D+str(PE), 
                    underflow=True, overflow=True
                )]

                # 1D Histograms
                hist_PE = compute_reco(var1D, 
                    reco_axes, 
                    [np.repeat(vars_dict[var1D]["val_8"], PE_samples, axis=0)], 
                    PE_weights)

                h_PE = Convert_to_TH1("hrecoBootstrap", PE, var1D, hist_PE)
                h_PE.Write()

                #2D Histograms
                for var2D in vars2D_dict.keys():

                    reco_axes_2D = [hist.axis.Variable(
                                        vars2D_dict[var2D]["reco_binning"], 
                                        name="hrecoBootstrap2D_"+var2D+str(PE), 
                                        label="hrecoBootstrap2D_ttbarmass"+var2D+str(PE), 
                                        underflow=True, overflow=True)]

                    hist2D_PE = compute_reco(var1D + "_" + var2D, 
                                             reco_axes+reco_axes_2D, 
                                             [np.repeat(vars_dict[var1D]["val_8"], PE_samples, axis=0)]+[np.repeat(vars2D_dict["ttbar_mass"]["val_8"], PE_samples, axis=0)], 
                                             PE_weights)

                    h2D_PE = Unwrap_to_TH1("hrecoBootstrap", PE, var1D, var2D,  hist2D_PE)
                    h2D_PE.Write()



        outHistFile.Close()


##############################

if __name__ == "__main__":

    main()

##############################  
