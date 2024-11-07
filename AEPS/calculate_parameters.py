import numpy as np
import os
import time
from ROOT import TFile, TH1D, TH2D

def calculate_aeps(T0, T7, T8, variable, S_gen0, S_gen7, S_gen8, S_reco8, B_E):
    """
    Calculate acceptance, efficiency, purity, and stability.
    
    Parameters:
    T0, T7, T8 (ROOT.TTree): Input trees
    variable (str): Variable name
    S_gen0, S_gen7, S_gen8, S_reco8 (str): Branch names
    B_E (list): Bin edges
    
    Calculates:
    A, E, AE, P, S (numpy arrays): Acceptance, efficiency, acceptance-efficiency, purity, stability
    """
    
    # Initialize histograms
    fHist0 = TH1D('Step0_' + variable, 'Gen0', len(B_E) - 1, B_E)
    fHist1 = TH1D('Step7gen_' + variable, 'Gen', len(B_E) - 1, B_E)
    fHist2 = TH1D('rec8gen_' + variable, 'Rec', len(B_E) - 1, B_E)
    fHist3 = TH1D('Step8gen_' + variable, 'Gen', len(B_E) - 1, B_E)
    fHist4 = TH2D('genRec' + variable, 'Gen VS reco', len(B_E) - 1, B_E, len(B_E) - 1, B_E)
    
    # Fill histograms
    for tree, hist, branch in [(T0, fHist0, S_gen0), (T7, fHist1, S_gen7), (T8, fHist2, S_reco8), (T8, fHist3, S_gen8)]:
        for i in range(tree.GetEntries()):
            tree.GetEntry(i)
            gen = getattr(tree, branch)
            hist.Fill(gen)
            
    for i in range(T8.GetEntries()):
        T8.GetEntry(i)
        gen = getattr(T8, S_gen8)
        reco = getattr(T8, S_reco8)
        fHist4.Fill(gen, reco)
    
    # Get histogram contents
    gen0_content = [fHist0.GetBinContent(i) for i in range(1, len(B_E))]
    gen7_content = [fHist1.GetBinContent(i) for i in range(1, len(B_E))]
    reco8_content = [fHist2.GetBinContent(i) for i in range(1, len(B_E))]
    gen8_content = [fHist3.GetBinContent(i) for i in range(1, len(B_E))]
    gen8VSrec8_content = [[fHist4.GetBinContent(i, j) for j in range(1, len(B_E))] for i in range(1, len(B_E))]
    gen8VSrec8_diagonal = [fHist4.GetBinContent(i, i) for i in range(1, len(B_E))]
    
    # Calculate AEPS
    np.seterr(divide='ignore', invalid='ignore')
    acceptance = np.nan_to_num(np.array(gen7_content) / np.array(gen0_content), posinf=0)
    efficiency = np.nan_to_num(np.array(reco8_content) / np.array(gen7_content), posinf=0)
    acceptance_efficiency = np.nan_to_num(np.array(reco8_content) / np.array(gen0_content), posinf=0)
    purity = np.nan_to_num(np.array(gen8VSrec8_diagonal) / np.array(reco8_content), posinf=0)
    stability = np.nan_to_num(np.array(gen8VSrec8_diagonal) / np.array(gen7_content), posinf=0)
    
    # Save results
    variableresult = '1D/' + variable + '/AEPS/'
    if not os.path.exists(variableresult):
        os.makedirs(variableresult)
    
    np.savetxt(os.path.join(my_path, variableresult, 'gen8VSrec8_content.csv'), gen8VSrec8_content, delimiter=',')
    np.savetxt(os.path.join(my_path, variableresult, 'A.csv'), acceptance, delimiter=',')
    np.savetxt(os.path.join(my_path, variableresult, 'E.csv'), efficiency, delimiter=',')
    np.savetxt(os.path.join(my_path, variableresult, 'AE.csv'), acceptance_efficiency, delimiter=',')
    np.savetxt(os.path.join(my_path, variableresult, 'P.csv'), purity, delimiter=',')
    np.savetxt(os.path.join(my_path, variableresult, 'S.csv'), stability, delimiter=',')

def calculate_2d_ae(T0, T7, T8, X_variable, Y_variable, X_gen0, Y_gen0, X_gen7, Y_gen7, X_rec8, Y_rec8, X_B_E, Y_B_E):
    """
    Calculate 2D acceptance and efficiency.
    
    Parameters:
    T0, T7, T8 (ROOT.TTree): Input trees
    X_variable, Y_variable (str): Variable names
    X_gen0, Y_gen0, X_gen7, Y_gen7, X_rec8, Y_rec8 (str): Branch names
    X_B_E, Y_B_E (list): Bin edges
    
     Calculates:
    acceptance, efficiency, acceptance_efficiency (numpy arrays)
    """
    
    def fill_histogram(tree, gen_X, gen_Y, weight_branch, bins_x, bins_y):
        hist = np.zeros((len(bins_x) - 1, len(bins_y) - 1))
        N = tree.GetEntries()
        for i in range(N):
            if i % 100000 == 0:
                print(i)
            tree.GetEntry(i)
            x = getattr(tree, gen_X)
            y = getattr(tree, gen_Y)
            weight = getattr(tree, weight_branch)
            for x_bin in range(len(bins_x) - 1):
                for y_bin in range(len(bins_y) - 1):
                    if bins_x[x_bin] < x < bins_x[x_bin + 1] and bins_y[y_bin] < y < bins_y[y_bin + 1]:
                        hist[x_bin][y_bin] += weight
        return hist
    
    gen_step0 = fill_histogram(T0, X_gen0, Y_gen0, "trueLevelWeight", X_B_E, Y_B_E)
    gen_step7 = fill_histogram(T7, X_gen7, Y_gen7, "trueLevelWeight", X_B_E, Y_B_E)
    rec_step8 = fill_histogram(T8, X_rec8, Y_rec8, "eventWeight", X_B_E, Y_B_E)
    
    acceptance = np.nan_to_num(gen_step7 / gen_step0, posinf=0)
    efficiency = np.nan_to_num(rec_step8 / gen_step7, posinf=0)
    acceptance_efficiency = np.nan_to_num(rec_step8 / gen_step0, posinf=0)
    
    # Save results
    variableresult = f'2D/{X_variable}___{Y_variable}/AE/'
    if not os.path.exists(variableresult):
        os.makedirs(variableresult)
    
    np.savetxt(os.path.join(my_path, variableresult, 'gen0.csv'), gen_step0, delimiter=',')
    np.savetxt(os.path.join(my_path, variableresult, 'rec8.csv'), rec_step8, delimiter=',')
    np.savetxt(os.path.join(my_path, variableresult, 'gen7.csv'), gen_step7, delimiter=',')
    np.savetxt(os.path.join(my_path, variableresult, 'A.csv'), acceptance, delimiter=',')
    np.savetxt(os.path.join(my_path, variableresult, 'E.csv'), efficiency, delimiter=',')
    np.savetxt(os.path.join(my_path, variableresult, 'AE.csv'), acceptance_efficiency, delimiter=',')
    np.savetxt(os.path.join(my_path, variableresult, 'X_B_E.csv'), X_B_E, delimiter=',')
    np.savetxt(os.path.join(my_path, variableresult, 'Y_B_E.csv'), Y_B_E, delimiter=',')

def calculate_2d_ps(T8, X_variable, Y_variable, X_gen8, Y_gen8, X_rec8, Y_rec8, X_B_E, Y_B_E):
    """
    Calculate 2D purity and stability.
    
    Parameters:
    T8 (ROOT.TTree): Input tree
    X_variable, Y_variable (str): Variable names
    X_gen8, Y_gen8, X_rec8, Y_rec8 (str): Branch names
    X_B_E, Y_B_E (list): Bin edges
    
     Calculates:
    purity, stability (numpy arrays)
    """
    
    def fill_gen_vs_rec_4d(tree, X_gen, Y_gen, X_rec, Y_rec, weight_branch, X_B_E, Y_B_E):
        gen_vs_rec_4d = np.zeros(((len(Y_B_E) - 1) * (len(X_B_E) - 1), (len(X_B_E) - 1) * (len(Y_B_E) - 1)))
        N = tree.GetEntries()
        for i in range(N):
            if i % 5000 == 0:
                print(i)
            tree.GetEntry(i)
            gen_x, gen_y, rec_x, rec_y, weight = getattr(tree, X_gen), getattr(tree, Y_gen), getattr(tree, X_rec), getattr(tree, Y_rec), getattr(tree, weight_branch)
            
            for x_bin in range(len(X_B_E) - 1):
                for x_bin_rec in range(len(X_B_E) - 1):
                    for y_bin in range(len(Y_B_E) - 1):
                        for y_bin_rec in range(len(Y_B_E) - 1):
                            if (X_B_E[x_bin] < gen_x < X_B_E[x_bin + 1] and 
                                X_B_E[x_bin_rec] < rec_x < X_B_E[x_bin_rec + 1] and 
                                Y_B_E[y_bin] < gen_y < Y_B_E[y_bin + 1] and 
                                Y_B_E[y_bin_rec] < rec_y < Y_B_E[y_bin_rec + 1]):
                                gen_vs_rec_4d[x_bin_rec + (y_bin_rec * (len(X_B_E) - 1))][x_bin + (y_bin * (len(X_B_E) - 1))] += weight
        return gen_vs_rec_4d
    
    gen_vs_rec_4d = fill_gen_vs_rec_4d(T8, X_gen8, Y_gen8, X_rec8, Y_rec8, "eventWeight", X_B_E, Y_B_E)
    
    # Save results
    variableresult = f'2D/{X_variable}___{Y_variable}/PS/'
    if not os.path.exists(variableresult):
        os.makedirs(variableresult)
    np.savetxt(os.path.join(my_path, variableresult, 'GenVSRec4D.csv'), gen_vs_rec_4d, delimiter=',')
    
    # Load gen7 data
    csv_loc_ae = f'2D/{X_variable}___{Y_variable}/AE/'
    gen_7 = np.genfromtxt(os.path.join(my_path, csv_loc_ae, 'gen7.csv'), delimiter=',')
    gen7 = gen_7.transpose().reshape((len(X_B_E) - 1) * (len(Y_B_E) - 1))
    
    # Calculate purity and stability
    zero_axis_sum = np.sum(gen_vs_rec_4d, axis=0)
    one_axis_sum = np.sum(gen_vs_rec_4d, axis=1)
    np.savetxt(os.path.join(my_path, variableresult, 'ZeroAxisSum.csv'), zero_axis_sum, delimiter=',')
    np.savetxt(os.path.join(my_path, variableresult, 'OneAxisSum.csv'), one_axis_sum, delimiter=',')
    
    purity = np.nan_to_num(np.diag(gen_vs_rec_4d) / one_axis_sum, posinf=0)
    stability = np.nan_to_num(np.diag(gen_vs_rec_4d) / gen7, posinf=0)
    
    # Reshape and transpose purity and stability
    purity = purity.reshape((len(Y_B_E) - 1, len(X_B_E) - 1)).transpose()
    stability = stability.reshape((len(Y_B_E) - 1, len(X_B_E) - 1)).transpose()
    
    np.savetxt(os.path.join(my_path, variableresult, 'purity.csv'), purity, delimiter=',')
    np.savetxt(os.path.join(my_path, variableresult, 'stability.csv'), stability, delimiter=',')


def main():
    start_time = time.time()
    my_path = os.getcwd()
    print(my_path)

    # Load ROOT file and trees
    file_path = '../emu_ttbarsignalplustau.root'
    root_file = TFile(file_path)
    trees = {
        'Tstep0': root_file.Get('ttBar_treeVariables_step0;46'),
        'Tstep7': root_file.Get('ttBar_treeVariables_step7;16'),
        'Tstep8': root_file.Get('ttBar_treeVariables_step8;18')
    }

    # Define bin edges
    bin_edges = {
        'delta_phi': [i / 20 * TMath.Pi() for i in range(21)],
        'delta_eta': [i / 20 * 6.0 for i in range(21)],
        'm_ttbar': [0, 450, 550, 800, 1500],
        'pt_top': [0, 100, 200, 400, 1000]
    }

    # Calculate AEPS and PS for 1D and 2D cases
    one_d_calculations = [
        ('delta_phi_llbar', 'gen_llbar_delta_phi', 'gen_llbar_delta_phi', 'gen_llbar_delta_phi', 'llbar_delta_phi', bin_edges['delta_phi'])
    ]

    two_d_ae_calculations = [
        ('delta_phi_llbar', 'm_ttbar', 'gen_llbar_delta_phi', 'gen_ttbar_mass', 'gen_llbar_delta_phi', 'gen_ttbar_mass', 'llbar_delta_phi', 'ttbar_mass', bin_edges['delta_phi'], bin_edges['m_ttbar']),
        ('delta_phi_llbar', 'pt_top', 'gen_llbar_delta_phi', 'gen_top_pt', 'gen_llbar_delta_phi', 'gen_top_pt', 'llbar_delta_phi', 'top_pt', bin_edges['delta_phi'], bin_edges['pt_top']),
        ('delta_eta_llbar', 'm_ttbar', 'gen_llbar_delta_eta', 'gen_ttbar_mass', 'gen_llbar_delta_eta', 'gen_ttbar_mass', 'llbar_delta_eta', 'ttbar_mass', bin_edges['delta_eta'], bin_edges['m_ttbar']),
        ('delta_eta_llbar', 'pt_top', 'gen_llbar_delta_eta', 'gen_top_pt', 'gen_llbar_delta_eta', 'gen_top_pt', 'llbar_delta_eta', 'top_pt', bin_edges['delta_eta'], bin_edges['pt_top']),
        ('m_ttbar', 'delta_phi_llbar', 'gen_ttbar_mass', 'gen_llbar_delta_phi', 'gen_ttbar_mass', 'gen_llbar_delta_phi', 'ttbar_mass', 'llbar_delta_phi', bin_edges['m_ttbar'], bin_edges['delta_phi'])
    ]

    two_d_ps_calculations = [
        ('delta_phi_llbar', 'm_ttbar', 'gen_llbar_delta_phi', 'gen_ttbar_mass', 'llbar_delta_phi', 'ttbar_mass', bin_edges['delta_phi'], bin_edges['m_ttbar']),
        ('delta_phi_llbar', 'pt_top', 'gen_llbar_delta_phi', 'gen_top_pt', 'llbar_delta_phi', 'top_pt', bin_edges['delta_phi'], bin_edges['pt_top']),
        ('delta_eta_llbar', 'm_ttbar', 'gen_llbar_delta_eta', 'gen_ttbar_mass', 'llbar_delta_eta', 'ttbar_mass', bin_edges['delta_eta'], bin_edges['m_ttbar']),
        ('delta_eta_llbar', 'pt_top', 'gen_llbar_delta_eta', 'gen_top_pt', 'llbar_delta_eta', 'top_pt', bin_edges['delta_eta'], bin_edges['pt_top']),
        ('m_ttbar', 'delta_phi_llbar', 'gen_ttbar_mass', 'gen_llbar_delta_phi', 'ttbar_mass', 'llbar_delta_phi', bin_edges['m_ttbar'], bin_edges['delta_phi'])
    ]

    for calculation in one_d_calculations:
        calculate_aeps(trees['Tstep0'], trees['Tstep7'], trees['Tstep8'], *calculation)
    

    for calculation in two_d_ae_calculations:
        calculate_2d_ae(trees['Tstep0'], trees['Tstep7'], trees['Tstep8'], *calculation)

    for calculation in two_d_ps_calculations:
        calculate_2d_ps(trees['Tstep8'], *calculation)


if __name__ == "__main__":
    main()