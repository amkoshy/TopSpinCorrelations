import os
import sys
import uproot
from ROOT import TFile, TH1D, TH2D, TUnfoldBinning, TUnfoldBinningXML, cout
import hist
import numpy as np
import argparse
import awkward as ak
from uncertainties import ufloat, unumpy
from contextlib import contextmanager

try:
    import ctypes
    from ctypes.util import find_library
except ImportError:
    libc = None
else:
    try:
        libc = ctypes.cdll.msvcrt # Windows
    except OSError:
        libc = ctypes.cdll.LoadLibrary(find_library('c'))

def flush(stream):
    try:
        libc.fflush(None)
        stream.flush()
    except (AttributeError, ValueError, IOError):
        pass # unsupported


def fileno(file_or_fd):
    fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd


@contextmanager
def stdout_redirected(to=os.devnull, stdout=None):
    if stdout is None:
       stdout = sys.stdout

    stdout_fd = fileno(stdout)
    with os.fdopen(os.dup(stdout_fd), 'w') as copied: 
        stdout.flush()
        try:
            os.dup2(fileno(to), stdout_fd)
        except ValueError:
            with open(to, 'w') as to_file:
                os.dup2(to_file.fileno(), stdout_fd)
        try:
            yield stdout
        finally:
            stdout.flush()
            os.dup2(copied.fileno(), stdout_fd)


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


def compute_response_matrix(var_name, gen_axes, reco_axes, gen_step8_data, step8_data, gen_step0_data, gen_step8_weights, step8_weights, gen_step0_weights, symmetrize = False):

    axes = gen_axes + reco_axes

    response_matrix = hist.Hist(*axes, storage=hist.storage.Weight())
    response_matrix_genUnderflow = hist.Hist(*axes, storage=hist.storage.Weight())
    response_matrix_genUnderflowCorrection = hist.Hist(*axes, storage=hist.storage.Weight())

    # fill the reconstructed events
    if symmetrize == False:
        response_matrix.fill(*gen_step8_data + step8_data, weight=step8_weights)
    elif symmetrize == True:
        response_matrix.fill(*[1.0*np.array(gen_step8_data[0])] + [data for data in gen_step8_data[1:] if gen_step8_data[1:]] + [1.0*np.array(step8_data[0])] + [data for data in step8_data[1:] if step8_data[1:]], weight=(0.5*np.array(step8_weights)))
        response_matrix.fill(*[-1.0*np.array(gen_step8_data[0])] + [data for data in gen_step8_data[1:] if gen_step8_data[1:]] + [-1.0*np.array(step8_data[0])] + [data for data in step8_data[1:] if step8_data[1:]], weight=(0.5*np.array(step8_weights)))
        response_matrix.variances(flow=True)[...] = 2.0 * response_matrix.variances(flow=True)[...]

    remove_underoverflow(response_matrix.values(flow=True))    
    remove_underoverflow(response_matrix.variances(flow=True))

    # fill the generator level events in the reconstruction underflow bin
    # this is back at step0 so that we can unfold back to full phase space
    if symmetrize == False:
        response_matrix_genUnderflow.fill(*gen_step0_data + [(reco_axes[j][0][0] - 1)*np.ones(len(gen_step0_data[0])) for j in range(0,len(reco_axes))], weight=gen_step0_weights)
    elif symmetrize == True:
        response_matrix_genUnderflow.fill(*[1.0*np.array(gen_step0_data[0])] + [data for data in gen_step0_data[1:] if gen_step0_data[1:]] + [(reco_axes[j][0][0] - 1)*np.ones(len(gen_step0_data[0])) for j in range(0,len(reco_axes))], weight=(0.5*np.array(gen_step0_weights)))
        response_matrix_genUnderflow.fill(*[-1.0*np.array(gen_step0_data[0])] + [data for data in gen_step0_data[1:] if gen_step0_data[1:]] + [(reco_axes[j][0][0] - 1)*np.ones(len(gen_step0_data[0])) for j in range(0,len(reco_axes))], weight=(0.5*np.array(gen_step0_weights)))

        response_matrix_genUnderflow.variances(flow=True)[...] = 2.0 * response_matrix_genUnderflow.variances(flow=True)[...]

    remove_underoverflow(response_matrix_genUnderflow.values(flow=True), True, False)
    remove_underoverflow(response_matrix_genUnderflow.variances(flow=True), True, False)
    
    # next create a separate response matrix to fill the reco weights
    # this is because we need to subtract the entries AND variances and can't
    # do that with a "simple" negative weight (this would only subtract entries,
    # variances are sum of weights squared)

    # fill this new matrix
    if symmetrize == False:
        response_matrix_genUnderflowCorrection.fill(*gen_step8_data + [(reco_axes[j][0][0] - 1)*np.ones(len(gen_step8_data[0])) for j in range(0,len(reco_axes))], weight=step8_weights)
    elif symmetrize == True:
        response_matrix_genUnderflowCorrection.fill(*[1.0*np.array(gen_step8_data[0])] + [data for data in gen_step8_data[1:] if gen_step8_data[1:]] + [(reco_axes[j][0][0] - 1)*np.ones(len(gen_step8_data[0])) for j in range(0,len(reco_axes))], weight=(0.5*np.array(step8_weights)))
        response_matrix_genUnderflowCorrection.fill(*[-1.0*np.array(gen_step8_data[0])] + [data for data in gen_step8_data[1:] if gen_step8_data[1:]] + [(reco_axes[j][0][0] - 1)*np.ones(len(gen_step8_data[0])) for j in range(0,len(reco_axes))], weight=(0.5*np.array(step8_weights)))

    remove_underoverflow(response_matrix_genUnderflowCorrection.values(flow=True), True, False)
    remove_underoverflow(response_matrix_genUnderflowCorrection.variances(flow=True), True, False)
    
    # set values and variances equal to the difference between response_matrix
    # and _response_matrix to account for these events that aren't reconstructed
    response_matrix[...] = np.concatenate(
        ((response_matrix.values(flow=True) + response_matrix_genUnderflow.values(flow=True) - response_matrix_genUnderflowCorrection.values(flow=True))[..., None],
         (response_matrix.variances(flow=True) + response_matrix_genUnderflow.variances(flow=True) - response_matrix_genUnderflowCorrection.variances(flow=True))[..., None]),
        axis=len(axes)
    )

    return response_matrix


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


def compute_gen(var_name, gen_axes, gen_step0_data, gen_step0_weights, symmetrize = False):
    
    gen = hist.Hist(*gen_axes, storage=hist.storage.Weight())

    if symmetrize == False:    
        gen.fill(*gen_step0_data, weight=gen_step0_weights)
    elif symmetrize == True:
        gen.fill(*[1.0*np.array(gen_step0_data[0])] + [data for data in gen_step0_data[1:] if gen_step0_data[1:]], weight=(0.5*np.array(gen_step0_weights)))
        gen.fill(*[-1.0*np.array(gen_step0_data[0])] + [data for data in gen_step0_data[1:] if gen_step0_data[1:]], weight=(0.5*np.array(gen_step0_weights)))
        gen.variances(flow=True)[...] = 2.0 * gen.variances(flow=True)[...]
    
    remove_underoverflow(gen.values(flow=True))    
    remove_underoverflow(gen.variances(flow=True))

    return gen


def compute_resolution(var_name, residual_axes, gen_axes, step8_data, gen_step8_data, step8_weights):
    
    axes = residual_axes + gen_axes
    
    resolution = hist.Hist(*axes, storage=hist.storage.Weight())

    resolution.fill(*np.subtract(step8_data, gen_step8_data), *gen_step8_data, weight=step8_weights)
    
    remove_underoverflow(resolution.values(flow=True))    
    remove_underoverflow(resolution.variances(flow=True))
    
    return resolution


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


def Convert_to_TH2(prefix, var, hist2D):

    TH2 = TH2D(prefix+"_"+var, 
               prefix+"_"+var, 
               hist2D.values(flow=False).shape[0], 
               np.array(GetBinning(hist2D.axes)[0]),
#               hist2D.axes[0][0][0], 
#               hist2D.axes[0][len(hist2D.axes[0])-1][1],
               hist2D.values(flow=False).shape[1], 
               np.array(GetBinning(hist2D.axes)[1]),
#               hist2D.axes[1][0][0], 
#               hist2D.axes[1][len(hist2D.axes[1])-1][1]
    )
    TH2.Sumw2()
    
    i_range = hist2D.values(flow=False).shape[0]+1
    j_range = hist2D.values(flow=False).shape[1]+1

    for i in range(0,i_range):
        for j in range(0,j_range):
            TH2.SetBinContent(i,j,hist2D.values(flow=True)[i,j])
            TH2.SetBinError(i,j,np.sqrt(hist2D.variances(flow=True)[i,j]))

    return TH2


def Unwrap_to_TH1(prefix, var1D, var2D, hist2D):
    
    TH1 = TH1D(prefix+"_"+var1D+"_"+var2D, 
                prefix+"_"+var1D+"_"+var2D, 
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


def Unwrap_to_TH2(prefix, var1D, var2D, resmat2D):

    i_range = resmat2D.view(flow=True).shape[1]-2
    j_range = resmat2D.view(flow=True).shape[3]-2
    m_range = resmat2D[:,0,:,0].view(flow=True).shape[0]-2
    n_range = resmat2D[:,0,:,0].view(flow=True).shape[1]-2

    resmat2D_unwrapped_values = np.zeros(((i_range)*(m_range)+2,(j_range)*(n_range)+2))
    resmat2D_unwrapped_variances = np.zeros(((i_range)*(m_range)+2,(j_range)*(n_range)+2))

    for i in range(1,i_range+1):
        for j in range(1,j_range+1):
            for m in range(1,m_range+1):
                for n in range(1,n_range+1):
                    resmat2D_unwrapped_values[(i-1)*(m_range) + m,(j-1)*(n_range) + n] = resmat2D.values(flow=True)[m,i,n,j]
                    resmat2D_unwrapped_variances[(i-1)*(m_range) + m,(j-1)*(n_range) + n] = resmat2D.variances(flow=True)[m,i,n,j]

    # Underflow for gen axes only.  No true underflow and no overflow whatsoever; just generated events that were not reconstructed in underflow bins.
    for i in range(1,i_range+1):
        for m in range(1,m_range+1):
            resmat2D_unwrapped_values[(i-1)*(m_range) + m,0] = resmat2D.values(flow=True)[m,i,0,0]
            resmat2D_unwrapped_variances[(i-1)*(m_range) + m,0] = resmat2D.variances(flow=True)[m,i,0,0]
    
    resmat2D_unwrapped = hist.Hist(
        hist.axis.Regular((i_range)*(m_range), 0, (i_range)*(m_range), name="gen", label="gen", underflow=True, overflow=True),
        hist.axis.Regular((j_range)*(n_range), 0, (j_range)*(n_range), name="reco", label="reco", underflow=True, overflow=True),
        storage=hist.storage.Weight()
    )

    resmat2D_unwrapped[...] = np.concatenate(
        ((resmat2D_unwrapped_values)[..., None],
         (resmat2D_unwrapped_variances)[..., None]),
        axis=2
    )

    TH2 = TH2D(prefix+"_"+var1D+"_"+var2D, 
               prefix+"_"+var1D+"_"+var2D, 
               resmat2D_unwrapped.values(flow=False).shape[0], 
               resmat2D.axes[0][0][0], 
               resmat2D.axes[0][len(resmat2D.axes[0])-1][1],
               resmat2D_unwrapped.values(flow=False).shape[1], 
               resmat2D.axes[0][0][0], 
               resmat2D.axes[0][len(resmat2D.axes[0])-1][1]
              )
    TH2.Sumw2()
    
    i_range = resmat2D_unwrapped.values(flow=False).shape[0]+1
    j_range = resmat2D_unwrapped.values(flow=False).shape[1]+1

    for i in range(0,i_range):
        for j in range(0,j_range):
            TH2.SetBinContent(i,j,resmat2D_unwrapped.values(flow=True)[i,j])
            TH2.SetBinError(i,j,np.sqrt(resmat2D_unwrapped.variances(flow=True)[i,j]))

    return TH2


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

    if not os.path.exists("UnfoldingHistos_"+era.replace("UL","")+"/"+systematic+"/"+channel): os.makedirs("UnfoldingHistos_"+era.replace("UL","")+"/"+systematic+"/"+channel)

    dileptonic_dir = "/depot/cms/top/jthiema/AnalysisFramework/TopSpinCorr_FullRunIIUL_September2022/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic"
    TUnfoldFileLists_dir = "TUnfoldFileLists_"+era.replace("UL","")
    TUnfoldFileList_filename = "TUnfoldFileList_"+systematic+"_"+channel+".txt"

    TUnfoldFileList = open(dileptonic_dir + "/" + TUnfoldFileLists_dir + "/" + TUnfoldFileList_filename, "r")

    filelist = TUnfoldFileList.readlines()

    for file in filelist:
        
        if filepattern not in file: continue

        print("Processing: " + file)

        inTreeFile  = TFile.Open ( "spinCorrInput_"+era.replace("UL","")+"/"+systematic+"/"+channel+"/"+os.path.splitext(os.path.basename(file))[0]+".root" , "READ" )

        outHistFile = TFile.Open ( "UnfoldingHistos_"+era.replace("UL","")+"/"+systematic+"/"+channel+"/histosTUnfold_"+os.path.splitext(os.path.basename(file))[0]+".root" , "RECREATE" )
        outHistFile.cd()

        wtdEvts = inTreeFile.Get("weightedEvents")
        wtdEvts.Write("hNrOfEvts")

        inTreeFile.Close()

        isTTBarSignal = "ttbarsignal" in file


        step0tree = uproot.open(file+':ttBar_treeVariables_step0')
        step8tree = uproot.open(file+':ttBar_treeVariables_step8')


        tree_vars = [

            "ll_cHel","ll_cLab","llbar_delta_phi","llbar_delta_eta",
            "b1k","b2k","b1r","b2r","b1n","b2n","b1j","b2j","b1q","b2q",
            "c_kk","c_rr","c_nn","c_kj","c_rq","c_rq",   
            "c_rk","c_kr","c_nr","c_rn","c_nk","c_kn","c_rj","c_jr",

            "gen_ll_cHel","gen_ll_cLab","gen_llbar_delta_phi","gen_llbar_delta_eta",
            "gen_b1k","gen_b2k","gen_b1r","gen_b2r","gen_b1n","gen_b2n","gen_b1j","gen_b2j","gen_b1q","gen_b2q",
            "gen_c_kk","gen_c_rr","gen_c_nn","gen_c_kj","gen_c_rq","gen_c_rq",   
            "gen_c_rk","gen_c_kr","gen_c_nr","gen_c_rn","gen_c_nk","gen_c_kn","gen_c_rj","gen_c_jr",

            "eventWeight", "trueLevelWeight",

            "ttbar_mass", "top_scatteringangle_ttbarframe",

            "gen_ttbar_mass", "gen_top_scatteringangle_ttbarframe",

            "gen_l_pt", "gen_lbar_pt", "gen_l_eta", "gen_lbar_eta", 
            "gen_b_pt", "gen_bbar_pt", "gen_b_eta", "gen_bbar_eta",
            #"gen_l_pdgid", "gen_lbar_pdgid",

        ]


        step0 = ak.Array(
            ak.zip( dict( (var, step0tree[var].array()) for var in tree_vars ) )
        )

        step8 = ak.Array(
            ak.zip( dict( (tree_var, step8tree[tree_var].array()) for tree_var in tree_vars ) )
        )

        visgen_mask = (
            (step0.gen_l_pt > 20.0) & (step0.gen_lbar_pt > 20.0)
            & (abs(step0.gen_l_eta) < 2.4) & (abs(step0.gen_lbar_eta) < 2.4)
            & (step0.gen_b_pt > 30) & (step0.gen_bbar_pt > 30)
            & (abs(step0.gen_b_eta) < 2.4) & (abs(step0.gen_bbar_eta) < 2.4)
        )

        visgen_step0 = step0[visgen_mask]


        vars_dict = {

            "ll_cHel" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_ll_cHel"],
                "visgen_val_0" : visgen_step0["gen_ll_cHel"],
                "val_8" : step8["ll_cHel"],
                "gen_val_8" : step8["gen_ll_cHel"], 
                "can_symmetrize" : True,
            },
            "ll_cLab" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_ll_cLab"],
                "visgen_val_0" : visgen_step0["gen_ll_cLab"],
                "val_8" : step8["ll_cLab"],
                "gen_val_8" : step8["gen_ll_cLab"], 
                "can_symmetrize" : True,
            },
            "llbar_delta_phi" : {
                "bin_edge_low" : 0.0, 
                "bin_edge_high" : np.pi, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -1.0, 
                "residual_bin_edge_high" : 1.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_llbar_delta_phi"],
                "visgen_val_0" : visgen_step0["gen_llbar_delta_phi"],
                "val_8" : step8["llbar_delta_phi"],
                "gen_val_8" : step8["gen_llbar_delta_phi"], 
                "can_symmetrize" : False,
            },
            "llbar_delta_eta" : {
                "bin_edge_low" : 0.0, 
                "bin_edge_high" : 5.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -1.0, 
                "residual_bin_edge_high" : 1.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_llbar_delta_eta"],
                "visgen_val_0" : visgen_step0["gen_llbar_delta_eta"],
                "val_8" : step8["llbar_delta_eta"],
                "gen_val_8" : step8["gen_llbar_delta_eta"], 
                "can_symmetrize" : False,
            },

            "b1k" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_b1k"],
                "visgen_val_0" : visgen_step0["gen_b1k"],
                "val_8" : step8["b1k"],
                "gen_val_8" : step8["gen_b1k"], 
                "can_symmetrize" : True,
            },
            "b2k" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_b2k"],
                "visgen_val_0" : visgen_step0["gen_b2k"],
                "val_8" : step8["b2k"],
                "gen_val_8" : step8["gen_b2k"], 
                "can_symmetrize" : True,
            },
            "b1r" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_b1r"],
                "visgen_val_0" : visgen_step0["gen_b1r"],
                "val_8" : step8["b1r"],
                "gen_val_8" : step8["gen_b1r"], 
                "can_symmetrize" : True,
            },
            "b2r" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_b2r"],
                "visgen_val_0" : visgen_step0["gen_b2r"],
                "val_8" : step8["b2r"],
                "gen_val_8" : step8["gen_b2r"], 
                "can_symmetrize" : True,
            },
            "b1n" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_b1n"],
                "visgen_val_0" : visgen_step0["gen_b1n"],
                "val_8" : step8["b1n"],
                "gen_val_8" : step8["gen_b1n"], 
                "can_symmetrize" : True,
            },
            "b2n" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_b2n"],
                "visgen_val_0" : visgen_step0["gen_b2n"],
                "val_8" : step8["b2n"],
                "gen_val_8" : step8["gen_b2n"], 
                "can_symmetrize" : True,
            },
            "b1j" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_b1j"],
                "visgen_val_0" : visgen_step0["gen_b1j"],
                "val_8" : step8["b1j"],
                "gen_val_8" : step8["gen_b1j"], 
                "can_symmetrize" : True,
            },
            "b2j" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_b2j"],
                "visgen_val_0" : visgen_step0["gen_b2j"],
                "val_8" : step8["b2j"],
                "gen_val_8" : step8["gen_b2j"], 
                "can_symmetrize" : True,
            },
            "b1q" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_b1q"],
                "visgen_val_0" : visgen_step0["gen_b1q"],
                "val_8" : step8["b1q"],
                "gen_val_8" : step8["gen_b1q"], 
                "can_symmetrize" : True,
            },
            "b2q" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_b2q"],
                "visgen_val_0" : visgen_step0["gen_b2q"],
                "val_8" : step8["b2q"],
                "gen_val_8" : step8["gen_b2q"], 
                "can_symmetrize" : True,
            },

            "c_kk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_c_kk"],
                "visgen_val_0" : visgen_step0["gen_c_kk"],
                "val_8" : step8["c_kk"],
                "gen_val_8" : step8["gen_c_kk"], 
                "can_symmetrize" : True,
            },
            "c_rr" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_c_rr"],
                "visgen_val_0" : visgen_step0["gen_c_rr"],
                "val_8" : step8["c_rr"],
                "gen_val_8" : step8["gen_c_rr"], 
                "can_symmetrize" : True,
            },
            "c_nn" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_c_nn"],
                "visgen_val_0" : visgen_step0["gen_c_nn"],
                "val_8" : step8["c_nn"],
                "gen_val_8" : step8["gen_c_nn"], 
                "can_symmetrize" : True,
            },
            "c_kj" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_c_kj"],
                "visgen_val_0" : visgen_step0["gen_c_kj"],
                "val_8" : step8["c_kj"],
                "gen_val_8" : step8["gen_c_kj"], 
                "can_symmetrize" : True,
            },

            "c_Prk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_c_rk"] + step0["gen_c_kr"],
                "visgen_val_0" : visgen_step0["gen_c_rk"] + visgen_step0["gen_c_kr"],
                "val_8" : step8["c_rk"] + step8["c_kr"],
                "gen_val_8" : step8["gen_c_rk"] + step8["gen_c_kr"], 
                "can_symmetrize" : True,
            },
            "c_Mrk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_c_rk"] - step0["gen_c_kr"],
                "visgen_val_0" : visgen_step0["gen_c_rk"] - visgen_step0["gen_c_kr"],
                "val_8" : step8["c_rk"] - step8["c_kr"],
                "gen_val_8" : step8["gen_c_rk"] - step8["gen_c_kr"], 
                "can_symmetrize" : True,
            },
            "c_Pnr" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_c_nr"] + step0["gen_c_kn"],
                "visgen_val_0" : visgen_step0["gen_c_nr"] + visgen_step0["gen_c_kn"],
                "val_8" : step8["c_nr"] + step8["c_kn"],
                "gen_val_8" : step8["gen_c_nr"] + step8["gen_c_kn"], 
                "can_symmetrize" : True,
            },
            "c_Mnr" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_c_nr"] - step0["gen_c_kn"],
                "visgen_val_0" : visgen_step0["gen_c_nr"] - visgen_step0["gen_c_kn"],
                "val_8" : step8["c_nr"] - step8["c_kn"],
                "gen_val_8" : step8["gen_c_nr"] - step8["gen_c_kn"], 
                "can_symmetrize" : True,
            },
            "c_Pnk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_c_nk"] + step0["gen_c_kn"],
                "visgen_val_0" : visgen_step0["gen_c_nk"] + visgen_step0["gen_c_kn"],
                "val_8" : step8["c_nk"] + step8["c_kn"],
                "gen_val_8" : step8["gen_c_nk"] + step8["gen_c_kn"], 
                "can_symmetrize" : True,
            },
            "c_Mnk" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_c_nk"] - step0["gen_c_kn"],
                "visgen_val_0" : visgen_step0["gen_c_nk"] - visgen_step0["gen_c_kn"],
                "val_8" : step8["c_nk"] - step8["c_kn"],
                "gen_val_8" : step8["gen_c_nk"] - step8["gen_c_kn"], 
                "can_symmetrize" : True,
            },
            "c_Prj" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_c_rj"] + step0["gen_c_jr"],
                "visgen_val_0" : visgen_step0["gen_c_rj"] + visgen_step0["gen_c_jr"],
                "val_8" : step8["c_rj"] + step8["c_jr"],
                "gen_val_8" : step8["gen_c_rj"] + step8["gen_c_jr"], 
                "can_symmetrize" : True,
            },
            "c_Mrj" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_c_rj"] - step0["gen_c_jr"],
                "visgen_val_0" : visgen_step0["gen_c_rj"] - visgen_step0["gen_c_jr"],
                "val_8" : step8["c_rj"] - step8["c_jr"],
                "gen_val_8" : step8["gen_c_rj"] - step8["gen_c_jr"], 
                "can_symmetrize" : True,
            },

            "c_han" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : +step0["gen_c_kk"] - step0["gen_c_rr"] - step0["gen_c_nn"],
                "visgen_val_0" : +visgen_step0["gen_c_kk"] - visgen_step0["gen_c_rr"] - visgen_step0["gen_c_nn"],
                "val_8" : +step8["c_kk"] - step8["c_rr"] - step8["c_nn"],
                "gen_val_8" : +step8["gen_c_kk"] - step8["gen_c_rr"] - step8["gen_c_nn"], 
                "can_symmetrize" : True,
            },
            "c_sca" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : -step0["gen_c_kk"] + step0["gen_c_rr"] - step0["gen_c_nn"],
                "visgen_val_0" : -visgen_step0["gen_c_kk"] + visgen_step0["gen_c_rr"] - visgen_step0["gen_c_nn"],
                "val_8" : -step8["c_kk"] + step8["c_rr"] - step8["c_nn"],
                "gen_val_8" : -step8["gen_c_kk"] + step8["gen_c_rr"] - step8["gen_c_nn"], 
                "can_symmetrize" : True,
            },
            "c_tra" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : -step0["gen_c_kk"] - step0["gen_c_rr"] + step0["gen_c_nn"],
                "visgen_val_0" : -visgen_step0["gen_c_kk"] - visgen_step0["gen_c_rr"] + visgen_step0["gen_c_nn"],
                "val_8" : -step8["c_kk"] - step8["c_rr"] + step8["c_nn"],
                "gen_val_8" : -step8["gen_c_kk"] - step8["gen_c_rr"] + step8["gen_c_nn"], 
                "can_symmetrize" : True,
            },
            "c_kjL" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : -step0["gen_c_kj"] - step0["gen_c_rr"] - step0["gen_c_nn"],
                "visgen_val_0" : -visgen_step0["gen_c_kj"] - visgen_step0["gen_c_rr"] - visgen_step0["gen_c_nn"],
                "val_8" : -step8["c_kj"] - step8["c_rr"] - step8["c_nn"],
                "gen_val_8" : -step8["gen_c_kj"] - step8["gen_c_rr"] - step8["gen_c_nn"], 
                "can_symmetrize" : True,
            },
            "c_rqL" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : -step0["gen_c_kk"] - step0["gen_c_rq"] - step0["gen_c_nn"],
                "visgen_val_0" : -visgen_step0["gen_c_kk"] - visgen_step0["gen_c_rq"] - visgen_step0["gen_c_nn"],
                "val_8" : -step8["c_kk"] - step8["c_rq"] - step8["c_nn"],
                "gen_val_8" : -step8["gen_c_kk"] - step8["gen_c_rq"] - step8["gen_c_nn"], 
                "can_symmetrize" : True,
            },

            "c_rkP" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : -step0["gen_c_rk"] - step0["gen_c_kr"] - step0["gen_c_nn"],
                "visgen_val_0" : -visgen_step0["gen_c_rk"] - visgen_step0["gen_c_kr"] - visgen_step0["gen_c_nn"],
                "val_8" : -step8["c_rk"] - step8["c_kr"] - step8["c_nn"],
                "gen_val_8" : -step8["gen_c_rk"] - step8["gen_c_kr"] - step8["gen_c_nn"], 
                "can_symmetrize" : True,
            },
            "c_rkM" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : -step0["gen_c_rk"] + step0["gen_c_kr"] - step0["gen_c_nn"],
                "visgen_val_0" : -visgen_step0["gen_c_rk"] + visgen_step0["gen_c_kr"] - visgen_step0["gen_c_nn"],
                "val_8" : -step8["c_rk"] + step8["c_kr"] - step8["c_nn"],
                "gen_val_8" : -step8["gen_c_rk"] + step8["gen_c_kr"] - step8["gen_c_nn"], 
                "can_symmetrize" : True,
            },
            "c_nrP" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : -step0["gen_c_nr"] - step0["gen_c_rn"] - step0["gen_c_kk"],
                "visgen_val_0" : -visgen_step0["gen_c_nr"] - visgen_step0["gen_c_rn"] - visgen_step0["gen_c_kk"],
                "val_8" : -step8["c_nr"] - step8["c_rn"] - step8["c_kk"],
                "gen_val_8" : -step8["gen_c_nr"] - step8["gen_c_rn"] - step8["gen_c_kk"], 
                "can_symmetrize" : True,
            },
            "c_nrM" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : -step0["gen_c_nr"] + step0["gen_c_rn"] - step0["gen_c_kk"],
                "visgen_val_0" : -visgen_step0["gen_c_nr"] + visgen_step0["gen_c_rn"] - visgen_step0["gen_c_kk"],
                "val_8" : -step8["c_nr"] + step8["c_rn"] - step8["c_kk"],
                "gen_val_8" : -step8["gen_c_nr"] + step8["gen_c_rn"] - step8["gen_c_kk"], 
                "can_symmetrize" : True,
            },
            "c_nkP" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : -step0["gen_c_nk"] - step0["gen_c_kn"] - step0["gen_c_rr"],
                "visgen_val_0" : -visgen_step0["gen_c_nk"] - visgen_step0["gen_c_kn"] - visgen_step0["gen_c_rr"],
                "val_8" : -step8["c_nk"] - step8["c_kn"] - step8["c_rr"],
                "gen_val_8" : -step8["gen_c_nk"] - step8["gen_c_kn"] - step8["gen_c_rr"], 
                "can_symmetrize" : True,
            },
            "c_nkM" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "n_gen_bins" : 6, 
                "n_reco_bins" : 12, 
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : -step0["gen_c_nk"] + step0["gen_c_kn"] - step0["gen_c_rr"],
                "visgen_val_0" : -visgen_step0["gen_c_nk"] + visgen_step0["gen_c_kn"] - visgen_step0["gen_c_rr"],
                "val_8" : -step8["c_nk"] + step8["c_kn"] - step8["c_rr"],
                "gen_val_8" : -step8["gen_c_nk"] + step8["gen_c_kn"] - step8["gen_c_rr"], 
                "can_symmetrize" : True,
            },

        }


        binning_ttbar_mass = [300., 450., 600., 800., 2000.]
        binning_top_scatteringangle_ttbarframe = [-1.0, -0.5, 0.0, +0.5, +1.0]

        vars2D_dict = {

            "ttbar_mass" : {
                "bin_edge_low" : 300., 
                "bin_edge_high" : 2000., 
                "gen_binning" : binning_ttbar_mass,
                "reco_binning" : BinFinely(binning_ttbar_mass, 2),
                "residual_bin_edge_low" : -1000.0, 
                "residual_bin_edge_high" : 1000.0, 
                "n_residual_bins" : 100,
                "gen_val_0" : step0["gen_ttbar_mass"],
                "visgen_val_0" : visgen_step0["gen_ttbar_mass"],
                "val_8" : step8["ttbar_mass"],
                "gen_val_8" : step8["gen_ttbar_mass"],
                "can_symmetrize" : False,
            },

            "top_scatteringangle_ttbarframe" : {
                "bin_edge_low" : -1.0, 
                "bin_edge_high" : 1.0, 
                "gen_binning" : binning_top_scatteringangle_ttbarframe,
                "reco_binning" : BinFinely(binning_top_scatteringangle_ttbarframe, 2),
                "residual_bin_edge_low" : -2.0, 
                "residual_bin_edge_high" : 2.0, 
                "n_residual_bins" : 80,
                "gen_val_0" : step0["gen_top_scatteringangle_ttbarframe"],
                "visgen_val_0" : visgen_step0["gen_top_scatteringangle_ttbarframe"],
                "val_8" : step8["top_scatteringangle_ttbarframe"],
                "gen_val_8" : step8["gen_top_scatteringangle_ttbarframe"],
                "can_symmetrize" : False,
            },

        }

        nFineBins_1D = 4
        nFineBins_2D = 2

        do_symmetrize = True

        # Remember the original sys.stdout object to revert so that we can close files
        stdout_obj = sys.stdout

        if not os.path.exists("binning"): os.makedirs("binning")

        TUnfoldBinningXML.WriteDTD("binning/tunfoldbinning.dtd")


        # 1D

        for var1D in vars_dict.keys():

            gen_axes = [hist.axis.Regular(nFineBins_1D*vars_dict[var1D]["n_gen_bins"], vars_dict[var1D]["bin_edge_low"], vars_dict[var1D]["bin_edge_high"], name="gen_"+var1D, label="gen_"+var1D, underflow=True, overflow=True),]
            reco_axes = [hist.axis.Regular(nFineBins_1D*vars_dict[var1D]["n_reco_bins"], vars_dict[var1D]["bin_edge_low"], vars_dict[var1D]["bin_edge_high"], name="reco_"+var1D, label="reco_"+var1D, underflow=True, overflow=True)]
            residual_axes = [hist.axis.Regular(vars_dict[var1D]["n_residual_bins"], vars_dict[var1D]["residual_bin_edge_low"], vars_dict[var1D]["residual_bin_edge_high"], name="residual_"+var1D, label="residual_"+var1D, underflow=True, overflow=True)]
            visgen_axes = [hist.axis.Regular(nFineBins_1D*vars_dict[var1D]["n_gen_bins"], vars_dict[var1D]["bin_edge_low"], vars_dict[var1D]["bin_edge_high"], name="visgen_"+var1D, label="visgen_"+var1D, underflow=True, overflow=True),]

            reco = compute_reco(var1D, reco_axes, [vars_dict[var1D]["val_8"]], [step8['eventWeight']], do_symmetrize*vars_dict[var1D]["can_symmetrize"])
            hreco = Convert_to_TH1("hreco", var1D, reco)
            hreco.Write()

            if (isTTBarSignal == False): continue

            resmat = compute_response_matrix(var1D, gen_axes, reco_axes, [vars_dict[var1D]["gen_val_8"]], [vars_dict[var1D]["val_8"]], [vars_dict[var1D]["gen_val_0"]], [step8['trueLevelWeight']], [step8['eventWeight']], [step0['trueLevelWeight']], do_symmetrize*vars_dict[var1D]["can_symmetrize"])
            hresmat = Convert_to_TH2("hrecoVsgen", var1D, resmat)
            hresmat.Write()

            gen = compute_gen(var1D, gen_axes, [vars_dict[var1D]["gen_val_0"]], [step0['trueLevelWeight']], do_symmetrize*vars_dict[var1D]["can_symmetrize"])
            hgen = Convert_to_TH1("hgen", var1D, gen)
            hgen.Write()

            visgen = compute_gen(var1D, visgen_axes, [vars_dict[var1D]["visgen_val_0"]], [visgen_step0['trueLevelWeight']], do_symmetrize*vars_dict[var1D]["can_symmetrize"])
            hvisgen = Convert_to_TH1("hvisgen", var1D, visgen)
            hvisgen.Write()

            resolution = compute_resolution(var1D, residual_axes, reco_axes, [vars_dict[var1D]["val_8"]], [vars_dict[var1D]["gen_val_8"]], [step8['eventWeight']])
            hresolution = Convert_to_TH2("hresolutionbins", var1D, resolution)
            hresolution.Write()


            if systematic != "Nominal": continue

            generatorBinning = TUnfoldBinning("generator")
            detectorBinning = TUnfoldBinning("detector")
            genDistribution = generatorBinning.AddBinning("ttbargen")
            detectorDistribution = detectorBinning.AddBinning("ttbarreco")

            detectorDistribution.AddAxis(var1D, len(GetBinning(reco_axes)[0])-1, np.array(GetBinning(reco_axes)[0]), False, False)
            genDistribution.AddAxis(var1D, len(GetBinning(gen_axes)[0])-1, np.array(GetBinning(gen_axes)[0]), False, False)

            with open("binning/"+var1D+"_binning.xml","w") as binning_xml, stdout_redirected(binning_xml):
                TUnfoldBinningXML.ExportXML(detectorBinning,cout,True,False)
                TUnfoldBinningXML.ExportXML(generatorBinning,cout,False,True)
                flush(binning_xml)

            rebin_txt = open("binning/"+var1D+"_rebin.txt","w")
            sys.stdout = rebin_txt
            print(nFineBins_1D)
            sys.stdout = stdout_obj
            rebin_txt.close()


        # 2D

        for var2D in vars2D_dict.keys():

            gen_axes_2D = [hist.axis.Variable(BinFinely(vars2D_dict[var2D]["gen_binning"],nFineBins_2D), name="gen2D_"+var2D, label="gen2D_"+var2D, underflow=True, overflow=True)]
            reco_axes_2D = [hist.axis.Variable(BinFinely(vars2D_dict[var2D]["reco_binning"],nFineBins_2D), name="reco2D_"+var2D, label="reco2D_"+var2D, underflow=True, overflow=True)]
            residual_axes_2D = [hist.axis.Regular(vars2D_dict[var2D]["n_residual_bins"], vars2D_dict[var2D]["residual_bin_edge_low"], vars2D_dict[var2D]["residual_bin_edge_high"], name="residual_"+var2D, label="residual_"+var2D, underflow=True, overflow=True)]
            visgen_axes_2D = [hist.axis.Variable(BinFinely(vars2D_dict[var2D]["gen_binning"],nFineBins_2D), name="visgen2D_"+var2D, label="visgen2D_"+var2D, underflow=True, overflow=True)]


            for var1D in vars_dict.keys():

                gen_axes = [hist.axis.Regular(nFineBins_2D*vars_dict[var1D]["n_gen_bins"], vars_dict[var1D]["bin_edge_low"], vars_dict[var1D]["bin_edge_high"], name="gen_"+var1D, label="gen_"+var1D, underflow=True, overflow=True),]
                reco_axes = [hist.axis.Regular(nFineBins_2D*vars_dict[var1D]["n_reco_bins"], vars_dict[var1D]["bin_edge_low"], vars_dict[var1D]["bin_edge_high"], name="reco_"+var1D, label="reco_"+var1D, underflow=True, overflow=True)]
                visgen_axes = [hist.axis.Regular(nFineBins_2D*vars_dict[var1D]["n_gen_bins"], vars_dict[var1D]["bin_edge_low"], vars_dict[var1D]["bin_edge_high"], name="visgen_"+var1D, label="visgen_"+var1D, underflow=True, overflow=True),]

                if (isTTBarSignal == False): continue

                reco = compute_reco(var1D+"_vs_"+var2D, reco_axes + reco_axes_2D, [vars_dict[var1D]["val_8"], vars2D_dict[var2D]["val_8"]], [step8['eventWeight']], do_symmetrize*vars_dict[var1D]["can_symmetrize"])
                hreco = Unwrap_to_TH1("hreco", var1D, var2D, reco)
                hreco.Write()

                resmat = compute_response_matrix(var1D+"_vs_"+var2D, gen_axes + gen_axes_2D, reco_axes + reco_axes_2D, [vars_dict[var1D]["gen_val_8"], vars2D_dict[var2D]["gen_val_8"]], [vars_dict[var1D]["val_8"], vars2D_dict[var2D]["val_8"]], [vars_dict[var1D]["gen_val_0"], vars2D_dict[var2D]["gen_val_0"]], [step8['trueLevelWeight']], [step8['eventWeight']], [step0['trueLevelWeight']], do_symmetrize*vars_dict[var1D]["can_symmetrize"])
                hresmat = Unwrap_to_TH2("hrecoVsgen", var1D, var2D, resmat)
                hresmat.Write()

                gen = compute_gen(var1D+"_vs_"+var2D, gen_axes + gen_axes_2D, [vars_dict[var1D]["gen_val_0"], vars2D_dict[var2D]["gen_val_0"]], [step0['trueLevelWeight']], do_symmetrize*vars_dict[var1D]["can_symmetrize"])
                hgen = Unwrap_to_TH1("hgen", var1D, var2D, gen)
                hgen.Write()

                visgen = compute_gen(var1D+"_vs_"+var2D, gen_axes + gen_axes_2D, [vars_dict[var1D]["visgen_val_0"], vars2D_dict[var2D]["visgen_val_0"]], [visgen_step0['trueLevelWeight']], do_symmetrize*vars_dict[var1D]["can_symmetrize"])
                hvisgen = Unwrap_to_TH1("hvisgen", var1D, var2D, visgen)
                hvisgen.Write()


                if systematic != "Nominal": continue

                generatorBinning = TUnfoldBinning("generator")
                detectorBinning = TUnfoldBinning("detector")
                genDistribution = generatorBinning.AddBinning("ttbargen")
                detectorDistribution = detectorBinning.AddBinning("ttbarreco")

                detectorDistribution.AddAxis(var1D, len(GetBinning(reco_axes)[0])-1, np.array(GetBinning(reco_axes)[0]), False, False)
                detectorDistribution.AddAxis(var2D, len(GetBinning(reco_axes_2D)[0])-1, np.array(GetBinning(reco_axes_2D)[0]), False, False)
                genDistribution.AddAxis(var1D, len(GetBinning(gen_axes)[0])-1, np.array(GetBinning(gen_axes)[0]), False, False)
                genDistribution.AddAxis(var2D, len(GetBinning(gen_axes_2D)[0])-1, np.array(GetBinning(gen_axes_2D)[0]), False, False)


                with open("binning/"+var1D+"_"+var2D+"_binning.xml","w") as binning_xml, stdout_redirected(binning_xml):
                    TUnfoldBinningXML.ExportXML(detectorBinning,cout,True,False)
                    TUnfoldBinningXML.ExportXML(generatorBinning,cout,False,True)
                    flush(binning_xml)

                rebin_txt = open("binning/"+var1D+"_"+var2D+"_rebin.txt","w")
                sys.stdout = rebin_txt
                print(nFineBins_2D)
                sys.stdout = stdout_obj
                rebin_txt.close()


            reco = compute_reco(var2D, reco_axes_2D, [vars2D_dict[var2D]["val_8"]], [step8['eventWeight']], do_symmetrize*vars2D_dict[var2D]["can_symmetrize"])
            hreco = Convert_to_TH1("hreco", var2D, reco)
            hreco.Write()

            if (isTTBarSignal == False): continue

            resmat = compute_response_matrix(var2D, gen_axes_2D, reco_axes_2D, [vars2D_dict[var2D]["gen_val_8"]], [vars2D_dict[var2D]["val_8"]], [vars2D_dict[var2D]["gen_val_0"]], [step8['trueLevelWeight']], [step8['eventWeight']], [step0['trueLevelWeight']], do_symmetrize*vars2D_dict[var2D]["can_symmetrize"])
            hresmat = Convert_to_TH2("hrecoVsgen", var2D, resmat)
            hresmat.Write()

            gen = compute_gen(var2D, gen_axes_2D, [vars2D_dict[var2D]["gen_val_0"]], [step0['trueLevelWeight']], do_symmetrize*vars2D_dict[var2D]["can_symmetrize"])
            hgen = Convert_to_TH1("hgen", var2D, gen)
            hgen.Write()

            visgen = compute_gen(var2D, visgen_axes_2D, [vars2D_dict[var2D]["visgen_val_0"]], [visgen_step0['trueLevelWeight']], do_symmetrize*vars2D_dict[var2D]["can_symmetrize"])
            hvisgen = Convert_to_TH1("hvisgen", var2D, visgen)
            hvisgen.Write()

            resolution = compute_resolution(var2D, residual_axes_2D, reco_axes_2D, [vars2D_dict[var2D]["val_8"]], [vars2D_dict[var2D]["gen_val_8"]], [step8['eventWeight']])
            hresolution = Convert_to_TH2("hresolutionbins", var2D, resolution)
            hresolution.Write()


            if systematic != "Nominal": continue

            generatorBinning = TUnfoldBinning("generator")
            detectorBinning = TUnfoldBinning("detector")
            genDistribution = generatorBinning.AddBinning("ttbargen")
            detectorDistribution = detectorBinning.AddBinning("ttbarreco")

            detectorDistribution.AddAxis(var2D, len(GetBinning(reco_axes_2D)[0])-1, np.array(GetBinning(reco_axes_2D)[0]), False, False)
            genDistribution.AddAxis(var2D, len(GetBinning(gen_axes_2D)[0])-1, np.array(GetBinning(gen_axes_2D)[0]), False, False)

            with open("binning/"+var2D+"_binning.xml","w") as binning_xml, stdout_redirected(binning_xml):
                TUnfoldBinningXML.ExportXML(detectorBinning,cout,True,False)
                TUnfoldBinningXML.ExportXML(generatorBinning,cout,False,True)
                flush(binning_xml)

            rebin_txt = open("binning/"+var2D+"_rebin.txt","w")
            sys.stdout = rebin_txt
            print(nFineBins_2D)
            sys.stdout = stdout_obj
            rebin_txt.close()


        outHistFile.Close()


##############################

if __name__ == "__main__":

    main()

##############################  
