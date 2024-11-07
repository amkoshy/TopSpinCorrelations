import tqdm
import math
import uproot
import numpy  as np

def main() :
    # 2016 variable list
    var_list = ['b1k', 'b2k','b1r','b2r','b1n','b2n','b1j','b2j','b1q','b2q',
                'c_kk','c_rr','c_nn','c_Prk','c_Mrk','c_Pnr','c_Mnr','c_Pnk','c_Mnk',
                'll_cHel','ll_cLab','llbar_delta_phi']
    
    BASEDIR = "/scratch/brown/bakshi3/TopSpinCorr_Run2_generalized_ND/CMSSW_10_6_30/src/TopAnalysis/Configuration/analysis/diLeptonic/TUnfoldResults_2016/Nominal/combined/Bootstrap_op_root_files_1D_2016_vars_1000PE/"

    fileptr_dict = {}
    for var in var_list :
        fileptr_dict[var] = uproot.open(BASEDIR + str(var) + ".root")

    nPE      = 200
    nbinsvar = 6
    nbinstot = nbinsvar * len(var_list)
    r        = np.zeros((nbinstot,nbinstot))
    cov      = np.zeros((nbinstot,nbinstot))

    print('Begin processing matrices with nPE :: ' + str(nPE))
    # Create an array with varname repeated 6 times and flatten it
    # Some fun list comprehension stuff

    n_times_over   = [[var] * nbinsvar for var in var_list]
    flattened_vars = [item for sublist in n_times_over for item in sublist]

    from tqdm import tqdm

    for xbin in tqdm(range(nbinstot)) :
        for ybin in range(nbinstot) :
            
            mu_x = 0
            mu_y = 0
            
            # Returns a string corresponding to varname
            xvar = flattened_vars[xbin]
            yvar = flattened_vars[ybin]
            
            # ##################
            # Computing the mean
            # ##################
            
            for i in range(nPE) :
                X = fileptr_dict[xvar][xvar + '_pseudo' + str(i) + 'TUnfResult_rebinnedA'].values()
                Y = fileptr_dict[yvar][yvar + '_pseudo' + str(i) + 'TUnfResult_rebinnedA'].values()   
                
                # modulo 6 for the var
                mu_x += X[xbin%6]
                mu_y += Y[ybin%6]
                
            mu_x /= nPE
            mu_y /= nPE

            sum_xy  = 0
            diff_x2 = 0
            diff_y2 = 0
            
            for i in range(nPE) :
                X = fileptr_dict[xvar][xvar +'_pseudo' + str(i) + 'TUnfResult_rebinnedA'].values()
                Y = fileptr_dict[yvar][yvar +'_pseudo' + str(i) + 'TUnfResult_rebinnedA'].values()   
                
                diff_x2 += (X[xbin%6] - mu_x)**2
                diff_y2 += (Y[ybin%6] - mu_y)**2
                
                sum_xy  += (X[xbin%6] - mu_x) * (Y[ybin%6] - mu_y)
                
            cov[xbin][ybin] = sum_xy / nPE
            
            # Correlation is covariance normalized by the 2 variances
            r[xbin][ybin]   = sum_xy / math.sqrt(diff_x2 * diff_y2)

    np.savetxt('Correlation_matrix_200PE.txt',   r, fmt='%1.6f')
    np.savetxt('Covariance_matrix_200PE.txt' , cov, fmt='%1.6f')

if __name__ == '__main__' :
    main()
