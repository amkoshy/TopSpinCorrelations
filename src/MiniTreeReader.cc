#include <MiniTreeReader.h>


MiniTreeReader::MiniTreeReader(TTree &tree, const Systematic::Systematic& sys)
{

    tree.SetBranchStatus("*", 0);

    tree.SetBranchStatus ("weight", 1);
    tree.SetBranchAddress("weight", &weight);
    tree.SetBranchStatus ("leptonSF", 1);
    tree.SetBranchAddress("leptonSF", &leptonSF);
    tree.SetBranchStatus ("l1PrefiringWeight", 1);
    tree.SetBranchAddress("l1PrefiringWeight", &l1PrefiringWeight);
    tree.SetBranchStatus ("triggerSF", 1);
    tree.SetBranchAddress("triggerSF", &triggerSF);
    tree.SetBranchStatus ("btagSF", 1);
    tree.SetBranchAddress("btagSF", &btagSF);
    tree.SetBranchStatus ("pileupSF", 1);
    tree.SetBranchAddress("pileupSF", &pileupSF);

    // DNN_rho =  tree.SetBranchAddress("DNN_rho", &DNN_rho);
    // tree.SetBranchStatus("DNN_rho", 0);

    tree.SetBranchStatus ("passStep3", 1);
    tree.SetBranchAddress("passStep3", &passStep3);
    tree.SetBranchStatus ("eventNumber", 0);
    tree.SetBranchStatus ("passStep4", 0);
    tree.SetBranchStatus ("passStep5", 0);
    tree.SetBranchStatus ("passStep6", 0);
    tree.SetBranchStatus ("passStep7", 0);
    tree.SetBranchStatus ("passStep8", 0);
    tree.SetBranchStatus ("passStep8Loose", 0);
    tree.SetBranchStatus ("lepton1_pt", 1);
    tree.SetBranchAddress("lepton1_pt", &lepton1_pt);
    tree.SetBranchStatus ("lepton1_eta", 1);
    tree.SetBranchAddress("lepton1_eta", &lepton1_eta);
    tree.SetBranchStatus ("lepton1_phi", 1);
    tree.SetBranchAddress("lepton1_phi", &lepton1_phi);
    tree.SetBranchStatus ("lepton1_m", 1);
    tree.SetBranchAddress("lepton1_m", &lepton1_m);
    tree.SetBranchStatus ("lepton2_pt", 1);
    tree.SetBranchAddress("lepton2_pt", &lepton2_pt);
    tree.SetBranchStatus ("lepton2_eta", 1);
    tree.SetBranchAddress("lepton2_eta", &lepton2_eta);
    tree.SetBranchStatus ("lepton2_phi", 1);
    tree.SetBranchAddress("lepton2_phi", &lepton2_phi);
    tree.SetBranchStatus ("lepton2_m", 1);
    tree.SetBranchAddress("lepton2_m", &lepton2_m);
    tree.SetBranchStatus ("met_pt", 1);
    tree.SetBranchAddress("met_pt", &met_pt);
    tree.SetBranchStatus ("met_phi", 1);
    tree.SetBranchAddress("met_phi", &met_phi);
    tree.SetBranchStatus ("met_significance", 1);
    tree.SetBranchAddress("met_significance", &met_significance);
    tree.SetBranchStatus ("n_jets", 1);
    tree.SetBranchAddress("n_jets", &n_jets);
    tree.SetBranchStatus ("jets_pt", 1);
    tree.SetBranchAddress("jets_pt", &jets_pt);
    tree.SetBranchStatus ("jets_eta", 1);
    tree.SetBranchAddress("jets_eta", &jets_eta);
    tree.SetBranchStatus ("jets_phi", 1);
    tree.SetBranchAddress("jets_phi", &jets_phi);
    tree.SetBranchStatus ("jets_m", 1);
    tree.SetBranchAddress("jets_m", &jets_m);
    tree.SetBranchStatus ("jets_btag", 1);
    tree.SetBranchAddress("jets_btag", &jets_btag);
    // tree.SetBranchStatus ("jets_charge", 0);
    // tree.SetBranchAddress("jets_charge", &jets_charge);




    // tree.SetBranchStatus("kinReco_*", 0);
    tree.SetBranchStatus ("hasKinRecoSolution", 1);
    tree.SetBranchAddress("hasKinRecoSolution", &hasKinRecoSolution);
    tree.SetBranchStatus ("kinReco_top_pt", 0);
    tree.SetBranchAddress("kinReco_top_pt", &kinReco_top_pt);
    tree.SetBranchStatus ("kinReco_top_eta", 0);
    tree.SetBranchAddress("kinReco_top_eta", &kinReco_top_eta);
    tree.SetBranchStatus ("kinReco_top_phi", 0);
    tree.SetBranchAddress("kinReco_top_phi", &kinReco_top_phi);
    tree.SetBranchStatus ("kinReco_top_m", 0);
    tree.SetBranchAddress("kinReco_top_m", &kinReco_top_m);
    tree.SetBranchStatus ("kinReco_antitop_pt", 0);
    tree.SetBranchAddress("kinReco_antitop_pt", &kinReco_antitop_pt);
    tree.SetBranchStatus ("kinReco_antitop_eta", 0);
    tree.SetBranchAddress("kinReco_antitop_eta", &kinReco_antitop_eta);
    tree.SetBranchStatus ("kinReco_antitop_phi", 0);
    tree.SetBranchAddress("kinReco_antitop_phi", &kinReco_antitop_phi);
    tree.SetBranchStatus ("kinReco_antitop_m", 0);
    tree.SetBranchAddress("kinReco_antitop_m", &kinReco_antitop_m);
    tree.SetBranchStatus ("kinReco_lepton_pt", 0);
    tree.SetBranchAddress("kinReco_lepton_pt", &kinReco_lepton_pt);
    tree.SetBranchStatus ("kinReco_lepton_eta", 0);
    tree.SetBranchAddress("kinReco_lepton_eta", &kinReco_lepton_eta);
    tree.SetBranchStatus ("kinReco_lepton_phi", 0);
    tree.SetBranchAddress("kinReco_lepton_phi", &kinReco_lepton_phi);
    tree.SetBranchStatus ("kinReco_lepton_m", 0);
    tree.SetBranchAddress("kinReco_lepton_m", &kinReco_lepton_m);
    tree.SetBranchStatus ("kinReco_antilepton_pt", 0);
    tree.SetBranchAddress("kinReco_antilepton_pt", &kinReco_antilepton_pt);
    tree.SetBranchStatus ("kinReco_antilepton_eta", 0);
    tree.SetBranchAddress("kinReco_antilepton_eta", &kinReco_antilepton_eta);
    tree.SetBranchStatus ("kinReco_antilepton_phi", 0);
    tree.SetBranchAddress("kinReco_antilepton_phi", &kinReco_antilepton_phi);
    tree.SetBranchStatus ("kinReco_antilepton_m", 0);
    tree.SetBranchAddress("kinReco_antilepton_m", &kinReco_antilepton_m);
    tree.SetBranchStatus ("kinReco_b_pt", 0);
    tree.SetBranchAddress("kinReco_b_pt", &kinReco_b_pt);
    tree.SetBranchStatus ("kinReco_b_eta", 0);
    tree.SetBranchAddress("kinReco_b_eta", &kinReco_b_eta);
    tree.SetBranchStatus ("kinReco_b_phi", 0);
    tree.SetBranchAddress("kinReco_b_phi", &kinReco_b_phi);
    tree.SetBranchStatus ("kinReco_b_m", 0);
    tree.SetBranchAddress("kinReco_b_m", &kinReco_b_m);
    tree.SetBranchStatus ("kinReco_antib_pt", 0);
    tree.SetBranchAddress("kinReco_antib_pt", &kinReco_antib_pt);
    tree.SetBranchStatus ("kinReco_antib_eta", 0);
    tree.SetBranchAddress("kinReco_antib_eta", &kinReco_antib_eta);
    tree.SetBranchStatus ("kinReco_antib_phi", 0);
    tree.SetBranchAddress("kinReco_antib_phi", &kinReco_antib_phi);
    tree.SetBranchStatus ("kinReco_antib_m", 0);
    tree.SetBranchAddress("kinReco_antib_m", &kinReco_antib_m);
    tree.SetBranchStatus ("kinReco_nonbjet_pt", 1);
    tree.SetBranchAddress("kinReco_nonbjet_pt", &kinReco_nonbjet_pt);
    tree.SetBranchStatus ("kinReco_nonbjet_eta", 1);
    tree.SetBranchAddress("kinReco_nonbjet_eta", &kinReco_nonbjet_eta);
    tree.SetBranchStatus ("kinReco_nonbjet_phi", 1);
    tree.SetBranchAddress("kinReco_nonbjet_phi", &kinReco_nonbjet_phi);
    tree.SetBranchStatus ("kinReco_nonbjet_m", 1);
    tree.SetBranchAddress("kinReco_nonbjet_m", &kinReco_nonbjet_m);
    tree.SetBranchStatus ("kinReco_rho", 1);
    tree.SetBranchAddress("kinReco_rho", &kinReco_rho);

    tree.SetBranchStatus ("looseKinReco_*", 0);
    tree.SetBranchStatus ("hasLooseKinRecoSolution", 0);
    tree.SetBranchAddress("hasLooseKinRecoSolution", &hasLooseKinRecoSolution);
    tree.SetBranchAddress("looseKinReco_ttbar_pt", &looseKinReco_ttbar_pt);
    tree.SetBranchAddress("looseKinReco_ttbar_eta", &looseKinReco_ttbar_eta);
    tree.SetBranchAddress("looseKinReco_ttbar_phi", &looseKinReco_ttbar_phi);
    tree.SetBranchAddress("looseKinReco_ttbar_m", &looseKinReco_ttbar_m);
    tree.SetBranchAddress("looseKinReco_b1_pt", &looseKinReco_b1_pt);
    tree.SetBranchAddress("looseKinReco_b1_eta", &looseKinReco_b1_eta);
    tree.SetBranchAddress("looseKinReco_b1_phi", &looseKinReco_b1_phi);
    tree.SetBranchAddress("looseKinReco_b1_m", &looseKinReco_b1_m);
    tree.SetBranchAddress("looseKinReco_b2_pt", &looseKinReco_b2_pt);
    tree.SetBranchAddress("looseKinReco_b2_eta", &looseKinReco_b2_eta);
    tree.SetBranchAddress("looseKinReco_b2_phi", &looseKinReco_b2_phi);
    tree.SetBranchAddress("looseKinReco_b2_m", &looseKinReco_b2_m);
    tree.SetBranchAddress("looseKinReco_l1_pt", &looseKinReco_l1_pt);
    tree.SetBranchAddress("looseKinReco_l1_eta", &looseKinReco_l1_eta);
    tree.SetBranchAddress("looseKinReco_l1_phi", &looseKinReco_l1_phi);
    tree.SetBranchAddress("looseKinReco_l1_m", &looseKinReco_l1_m);
    tree.SetBranchAddress("looseKinReco_l2_pt", &looseKinReco_l2_pt);
    tree.SetBranchAddress("looseKinReco_l2_eta", &looseKinReco_l2_eta);
    tree.SetBranchAddress("looseKinReco_l2_phi", &looseKinReco_l2_phi);
    tree.SetBranchAddress("looseKinReco_l2_m", &looseKinReco_l2_m);
    tree.SetBranchAddress("looseKinReco_nonbjet_pt", &looseKinReco_nonbjet_pt);
    tree.SetBranchAddress("looseKinReco_nonbjet_eta", &looseKinReco_nonbjet_eta);
    tree.SetBranchAddress("looseKinReco_nonbjet_phi", &looseKinReco_nonbjet_phi);
    tree.SetBranchAddress("looseKinReco_nonbjet_m", &looseKinReco_nonbjet_m);
    tree.SetBranchAddress("looseKinReco_rho", &looseKinReco_rho);


    tree.SetBranchStatus ("gen_additional_jet_pt", 0);
    tree.SetBranchAddress("gen_additional_jet_pt", &gen_additional_jet_pt);
    tree.SetBranchStatus ("gen_additional_jet_eta", 0);
    tree.SetBranchAddress("gen_additional_jet_eta", &gen_additional_jet_eta);
    tree.SetBranchStatus ("gen_additional_jet_phi", 0);
    tree.SetBranchAddress("gen_additional_jet_phi", &gen_additional_jet_phi);
    tree.SetBranchStatus ("gen_additional_jet_m", 0);
    tree.SetBranchAddress("gen_additional_jet_m", &gen_additional_jet_m);
    tree.SetBranchStatus ("gen_additional_jet_withNu_pt", 0);
    tree.SetBranchAddress("gen_additional_jet_withNu_pt", &gen_additional_jet_withNu_pt);
    tree.SetBranchStatus ("gen_additional_jet_withNu_eta", 0);
    tree.SetBranchAddress("gen_additional_jet_withNu_eta", &gen_additional_jet_withNu_eta);
    tree.SetBranchStatus ("gen_additional_jet_withNu_phi", 0);
    tree.SetBranchAddress("gen_additional_jet_withNu_phi", &gen_additional_jet_withNu_phi);
    tree.SetBranchStatus ("gen_additional_jet_withNu_m", 0);
    tree.SetBranchAddress("gen_additional_jet_withNu_m", &gen_additional_jet_withNu_m);
    tree.SetBranchStatus ("gen_top_pt", &gen_top_pt);
    tree.SetBranchAddress("gen_top_pt", &gen_top_pt);
    tree.SetBranchStatus ("gen_top_eta", &gen_top_eta);
    tree.SetBranchAddress("gen_top_eta", &gen_top_eta);
    tree.SetBranchStatus ("gen_top_phi", &gen_top_phi);
    tree.SetBranchAddress("gen_top_phi", &gen_top_phi);
    tree.SetBranchStatus ("gen_top_m", &gen_top_m);
    tree.SetBranchAddress("gen_top_m", &gen_top_m);
    tree.SetBranchStatus ("gen_antitop_pt", &gen_antitop_pt);
    tree.SetBranchAddress("gen_antitop_pt", &gen_antitop_pt);
    tree.SetBranchStatus ("gen_antitop_eta", &gen_antitop_eta);
    tree.SetBranchAddress("gen_antitop_eta", &gen_antitop_eta);
    tree.SetBranchStatus ("gen_antitop_phi", &gen_antitop_phi);
    tree.SetBranchAddress("gen_antitop_phi", &gen_antitop_phi);
    tree.SetBranchStatus ("gen_antitop_m", &gen_antitop_m);
    tree.SetBranchAddress("gen_antitop_m", &gen_antitop_m);
    tree.SetBranchStatus ("gen_lepton_pt", 0);
    tree.SetBranchAddress("gen_lepton_pt", &gen_lepton_pt);
    tree.SetBranchStatus ("gen_lepton_eta", 0);
    tree.SetBranchAddress("gen_lepton_eta", &gen_lepton_eta);
    tree.SetBranchStatus ("gen_lepton_phi", 0);
    tree.SetBranchAddress("gen_lepton_phi", &gen_lepton_phi);
    tree.SetBranchStatus ("gen_lepton_m", 0);
    tree.SetBranchAddress("gen_lepton_m", &gen_lepton_m);
    tree.SetBranchStatus ("gen_antilepton_pt", 0);
    tree.SetBranchAddress("gen_antilepton_pt", &gen_antilepton_pt);
    tree.SetBranchStatus ("gen_antilepton_eta", 0);
    tree.SetBranchAddress("gen_antilepton_eta", &gen_antilepton_eta);
    tree.SetBranchStatus ("gen_antilepton_phi", 0);
    tree.SetBranchAddress("gen_antilepton_phi", &gen_antilepton_phi);
    tree.SetBranchStatus ("gen_antilepton_m", 0);
    tree.SetBranchAddress("gen_antilepton_m", &gen_antilepton_m);
    tree.SetBranchStatus ("gen_b_pt", 0);
    tree.SetBranchAddress("gen_b_pt", &gen_b_pt);
    tree.SetBranchStatus ("gen_b_eta", 0);
    tree.SetBranchAddress("gen_b_eta", &gen_b_eta);
    tree.SetBranchStatus ("gen_b_phi", 0);
    tree.SetBranchAddress("gen_b_phi", &gen_b_phi);
    tree.SetBranchStatus ("gen_b_m", 0);
    tree.SetBranchAddress("gen_b_m", &gen_b_m);
    tree.SetBranchStatus ("gen_antib_pt", 0);
    tree.SetBranchAddress("gen_antib_pt", &gen_antib_pt);
    tree.SetBranchStatus ("gen_antib_eta", 0);
    tree.SetBranchAddress("gen_antib_eta", &gen_antib_eta);
    tree.SetBranchStatus ("gen_antib_phi", 0);
    tree.SetBranchAddress("gen_antib_phi", &gen_antib_phi);
    tree.SetBranchStatus ("gen_antib_m", 0);
    tree.SetBranchAddress("gen_antib_m", &gen_antib_m);
    tree.SetBranchStatus ("gen_neutrino_pt", 0);
    tree.SetBranchAddress("gen_neutrino_pt", &gen_neutrino_pt);
    tree.SetBranchStatus ("gen_neutrino_eta", 0);
    tree.SetBranchAddress("gen_neutrino_eta", &gen_neutrino_eta);
    tree.SetBranchStatus ("gen_neutrino_phi", 0);
    tree.SetBranchAddress("gen_neutrino_phi", &gen_neutrino_phi);
    tree.SetBranchStatus ("gen_neutrino_m", 0);
    tree.SetBranchAddress("gen_neutrino_m", &gen_neutrino_m);
    tree.SetBranchStatus ("gen_antineutrino_pt", 0);
    tree.SetBranchAddress("gen_antineutrino_pt", &gen_antineutrino_pt);
    tree.SetBranchStatus ("gen_antineutrino_eta", 0);
    tree.SetBranchAddress("gen_antineutrino_eta", &gen_antineutrino_eta);
    tree.SetBranchStatus ("gen_antineutrino_phi", 0);
    tree.SetBranchAddress("gen_antineutrino_phi", &gen_antineutrino_phi);
    tree.SetBranchStatus ("gen_antineutrino_m", 0);
    tree.SetBranchAddress("gen_antineutrino_m", &gen_antineutrino_m);
    tree.SetBranchStatus ("gen_met_pt", 0);
    tree.SetBranchAddress("gen_met_pt", &gen_met_pt);
    tree.SetBranchStatus ("gen_met_eta", 0);
    tree.SetBranchAddress("gen_met_eta", &gen_met_eta);
    tree.SetBranchStatus ("gen_met_phi", 0);
    tree.SetBranchAddress("gen_met_phi", &gen_met_phi);
    tree.SetBranchStatus ("gen_met_m", 0);
    tree.SetBranchAddress("gen_met_m", &gen_met_m);
    tree.SetBranchStatus ("gen_rho", 0);
    tree.SetBranchAddress("gen_rho", &gen_rho);
    tree.SetBranchStatus ("gen_rhoWithNu", 0);
    tree.SetBranchAddress("gen_rhoWithNu", &gen_rhoWithNu);
    tree.SetBranchStatus ("pseudo_top_pt", 0);
    tree.SetBranchAddress("pseudo_top_pt", &pseudo_top_pt);
    tree.SetBranchStatus ("pseudo_top_eta", 0);
    tree.SetBranchAddress("pseudo_top_eta", &pseudo_top_eta);
    tree.SetBranchStatus ("pseudo_top_phi", 0);
    tree.SetBranchAddress("pseudo_top_phi", &pseudo_top_phi);
    tree.SetBranchStatus ("pseudo_top_m", 0);
    tree.SetBranchAddress("pseudo_top_m", &pseudo_top_m);
    tree.SetBranchStatus ("pseudo_antitop_pt", 0);
    tree.SetBranchAddress("pseudo_antitop_pt", &pseudo_antitop_pt);
    tree.SetBranchStatus ("pseudo_antitop_eta", 0);
    tree.SetBranchAddress("pseudo_antitop_eta", &pseudo_antitop_eta);
    tree.SetBranchStatus ("pseudo_antitop_phi", 0);
    tree.SetBranchAddress("pseudo_antitop_phi", &pseudo_antitop_phi);
    tree.SetBranchStatus ("pseudo_antitop_m", 0);
    tree.SetBranchAddress("pseudo_antitop_m", &pseudo_antitop_m);
    tree.SetBranchStatus ("pseudo_lepton_pt", 0);
    tree.SetBranchAddress("pseudo_lepton_pt", &pseudo_lepton_pt);
    tree.SetBranchStatus ("pseudo_lepton_eta", 0);
    tree.SetBranchAddress("pseudo_lepton_eta", &pseudo_lepton_eta);
    tree.SetBranchStatus ("pseudo_lepton_phi", 0);
    tree.SetBranchAddress("pseudo_lepton_phi", &pseudo_lepton_phi);
    tree.SetBranchStatus ("pseudo_lepton_m", 0);
    tree.SetBranchAddress("pseudo_lepton_m", &pseudo_lepton_m);
    tree.SetBranchStatus ("pseudo_antilepton_pt", 0);
    tree.SetBranchAddress("pseudo_antilepton_pt", &pseudo_antilepton_pt);
    tree.SetBranchStatus ("pseudo_antilepton_eta", 0);
    tree.SetBranchAddress("pseudo_antilepton_eta", &pseudo_antilepton_eta);
    tree.SetBranchStatus ("pseudo_antilepton_phi", 0);
    tree.SetBranchAddress("pseudo_antilepton_phi", &pseudo_antilepton_phi);
    tree.SetBranchStatus ("pseudo_antilepton_m", 0);
    tree.SetBranchAddress("pseudo_antilepton_m", &pseudo_antilepton_m);
    tree.SetBranchStatus ("pseudo_b_pt", 0);
    tree.SetBranchAddress("pseudo_b_pt", &pseudo_b_pt);
    tree.SetBranchStatus ("pseudo_b_eta", 0);
    tree.SetBranchAddress("pseudo_b_eta", &pseudo_b_eta);
    tree.SetBranchStatus ("pseudo_b_phi", 0);
    tree.SetBranchAddress("pseudo_b_phi", &pseudo_b_phi);
    tree.SetBranchStatus ("pseudo_b_m", 0);
    tree.SetBranchAddress("pseudo_b_m", &pseudo_b_m);
    tree.SetBranchStatus ("pseudo_antib_pt", 0);
    tree.SetBranchAddress("pseudo_antib_pt", &pseudo_antib_pt);
    tree.SetBranchStatus ("pseudo_antib_eta", 0);
    tree.SetBranchAddress("pseudo_antib_eta", &pseudo_antib_eta);
    tree.SetBranchStatus ("pseudo_antib_phi", 0);
    tree.SetBranchAddress("pseudo_antib_phi", &pseudo_antib_phi);
    tree.SetBranchStatus ("pseudo_antib_m", 0);
    tree.SetBranchAddress("pseudo_antib_m", &pseudo_antib_m);
    tree.SetBranchStatus ("pseudo_additional_jet_pt", 0);
    tree.SetBranchAddress("pseudo_additional_jet_pt", &pseudo_additional_jet_pt);
    tree.SetBranchStatus ("pseudo_additional_jet_eta", 0);
    tree.SetBranchAddress("pseudo_additional_jet_eta", &pseudo_additional_jet_eta);
    tree.SetBranchStatus ("pseudo_additional_jet_phi", 0);
    tree.SetBranchAddress("pseudo_additional_jet_phi", &pseudo_additional_jet_phi);
    tree.SetBranchStatus ("pseudo_additional_jet_m", 0);
    tree.SetBranchAddress("pseudo_additional_jet_m", &pseudo_additional_jet_m);
    tree.SetBranchStatus ("pseudo_rho", 0);
    tree.SetBranchAddress("pseudo_rho", &pseudo_rho);
    tree.SetBranchStatus ("gen_partonLevel_additional_jet_pt", 1);
    tree.SetBranchAddress("gen_partonLevel_additional_jet_pt", &gen_partonLevel_additional_jet_pt);
    tree.SetBranchStatus ("gen_partonLevel_additional_jet_eta", 1);
    tree.SetBranchAddress("gen_partonLevel_additional_jet_eta", &gen_partonLevel_additional_jet_eta);
    tree.SetBranchStatus ("gen_partonLevel_additional_jet_phi", 1);
    tree.SetBranchAddress("gen_partonLevel_additional_jet_phi", &gen_partonLevel_additional_jet_phi);
    tree.SetBranchStatus ("gen_partonLevel_additional_jet_m", 1);
    tree.SetBranchAddress("gen_partonLevel_additional_jet_m", &gen_partonLevel_additional_jet_m);
    tree.SetBranchStatus ("gen_partonLevel_rho", 1);
    tree.SetBranchAddress("gen_partonLevel_rho", &gen_partonLevel_rho);

    tree.SetBranchStatus("var_*", 0);
    switch(sys.type()){
        case Systematic::meRenScale:
            if(sys.variation()==Systematic::Variation::up){
                tree.SetBranchStatus ("var_MERENSCALE_UP", 1);
                tree.SetBranchAddress("var_MERENSCALE_UP", &systVarWeight_MERENSCALE_UP);
                tree.SetBranchStatus ("var_BTAG_MERENSCALE_UP", 1);
                tree.SetBranchAddress("var_BTAG_MERENSCALE_UP", &systVarWeight_BTAG_MERENSCALE_UP);

                }
            if(sys.variation()==Systematic::Variation::down){
                tree.SetBranchStatus ("var_MERENSCALE_DOWN", 1);
                tree.SetBranchAddress("var_MERENSCALE_DOWN", &systVarWeight_MERENSCALE_DOWN);
                tree.SetBranchStatus ("var_BTAG_MERENSCALE_DOWN", 1);
                tree.SetBranchAddress("var_BTAG_MERENSCALE_DOWN", &systVarWeight_BTAG_MERENSCALE_DOWN);
            }
            break;
        case Systematic::meFacScale:
            if(sys.variation()==Systematic::Variation::up){
                tree.SetBranchStatus ("var_MEFACSCALE_UP", 1);
                tree.SetBranchAddress("var_MEFACSCALE_UP", &systVarWeight_MEFACSCALE_UP);
                tree.SetBranchStatus ("var_BTAG_MEFACSCALE_UP", 1);
                tree.SetBranchAddress("var_BTAG_MEFACSCALE_UP", &systVarWeight_BTAG_MEFACSCALE_UP);

            }
            if(sys.variation()==Systematic::Variation::down){
                tree.SetBranchStatus ("var_MEFACSCALE_DOWN", 1);
                tree.SetBranchAddress("var_MEFACSCALE_DOWN", &systVarWeight_MEFACSCALE_DOWN);
                tree.SetBranchStatus ("var_BTAG_MEFACSCALE_DOWN", 1);
                tree.SetBranchAddress("var_BTAG_MEFACSCALE_DOWN", &systVarWeight_BTAG_MEFACSCALE_DOWN);

            }
            break;
        case Systematic::meScale:
            if(sys.variation()==Systematic::Variation::up){
                tree.SetBranchStatus ("var_MESCALE_UP", 1);
                tree.SetBranchAddress("var_MESCALE_UP", &systVarWeight_MESCALE_UP);
                tree.SetBranchStatus ("var_BTAG_MESCALE_UP", 1);
                tree.SetBranchAddress("var_BTAG_MESCALE_UP", &systVarWeight_BTAG_MESCALE_UP);
            }
            if(sys.variation()==Systematic::Variation::down){
                tree.SetBranchStatus ("var_MESCALE_DOWN", 1);
                tree.SetBranchAddress("var_MESCALE_DOWN", &systVarWeight_MESCALE_DOWN);
                tree.SetBranchStatus ("var_BTAG_MESCALE_DOWN", 1);
                tree.SetBranchAddress("var_BTAG_MESCALE_DOWN", &systVarWeight_BTAG_MESCALE_DOWN);
            }
            break;
        case Systematic::bFrag:
            if(sys.variation()==Systematic::Variation::up){
                tree.SetBranchStatus ("var_BFRAG_UP", 1);
                tree.SetBranchAddress("var_BFRAG_UP", &systVarWeight_BFRAG_UP);
                tree.SetBranchStatus ("var_BTAG_BFRAG_UP", 1);
                tree.SetBranchAddress("var_BTAG_BFRAG_UP", &systVarWeight_BTAG_BFRAG_UP);
            }
            if(sys.variation()==Systematic::Variation::down){
                tree.SetBranchStatus ("var_BFRAG_DOWN", 1);
                tree.SetBranchAddress("var_BFRAG_DOWN", &systVarWeight_BFRAG_DOWN);
                tree.SetBranchStatus ("var_BTAG_BFRAG_DOWN", 1);
                tree.SetBranchAddress("var_BTAG_BFRAG_DOWN", &systVarWeight_BTAG_BFRAG_DOWN);
            }
            break;
        case Systematic::bFrag_Peterson:
             tree.SetBranchStatus ("var_BFRAG_PETERSON", 1);
             tree.SetBranchAddress("var_BFRAG_PETERSON", &systVarWeight_BFRAG_PETERSON);
             tree.SetBranchStatus ("var_BTAG_BFRAG_PETERSON", 1);
             tree.SetBranchAddress("var_BTAG_BFRAG_PETERSON", &systVarWeight_BTAG_BFRAG_PETERSON);
            break;
        case Systematic::bSemilep:
            if(sys.variation()==Systematic::Variation::up){
                tree.SetBranchAddress("var_BSEMILEP_UP", &systVarWeight_BSEMILEP_UP);
                tree.SetBranchStatus("var_BSEMILEP_UP", 1);
                tree.SetBranchStatus ("var_BTAG_BSEMILEP_UP", 1);
                tree.SetBranchAddress("var_BTAG_BSEMILEP_UP", &systVarWeight_BTAG_BSEMILEP_UP);
            }
            if(sys.variation()==Systematic::Variation::down){
                tree.SetBranchStatus ("var_BSEMILEP_DOWN", 1);
                tree.SetBranchAddress("var_BSEMILEP_DOWN", &systVarWeight_BSEMILEP_DOWN);
                tree.SetBranchStatus ("var_BTAG_BSEMILEP_DOWN", 1);
                tree.SetBranchAddress("var_BTAG_BSEMILEP_DOWN", &systVarWeight_BTAG_BSEMILEP_DOWN);
            }
            break;
        case Systematic::pdf:
            break;
        case Systematic::psScaleWeight:
            // tree.SetBranchAddress("var_PS", &systVarWeight_PS);
            // tree.SetBranchStatus("var_PS", 1);
            if(sys.variationNumber()==4){
                if (sys.variation()==Systematic::Variation::up){
                    tree.SetBranchStatus ("var_PSSCALE_WEIGHT_4_UP", 1);
                    tree.SetBranchAddress("var_PSSCALE_WEIGHT_4_UP", &systVarWeight_PSSCALE_WEIGHT_4_UP);
                    tree.SetBranchStatus ("var_BTAG_PSSCALE_WEIGHT_4_UP", 1);
                    tree.SetBranchAddress("var_BTAG_PSSCALE_WEIGHT_4_UP", &systVarWeight_BTAG_PSSCALE_WEIGHT_4_UP);
                }
                if (sys.variation()==Systematic::Variation::down){
                    tree.SetBranchStatus ("var_PSSCALE_WEIGHT_4_DOWN", 1);
                    tree.SetBranchAddress("var_PSSCALE_WEIGHT_4_DOWN", &systVarWeight_PSSCALE_WEIGHT_4_DOWN);
                    tree.SetBranchStatus ("var_BTAG_PSSCALE_WEIGHT_4_DOWN", 1);
                    tree.SetBranchAddress("var_BTAG_PSSCALE_WEIGHT_4_DOWN", &systVarWeight_BTAG_PSSCALE_WEIGHT_4_DOWN);
                }
            }
            if(sys.variationNumber()==5){
                if (sys.variation()==Systematic::Variation::up){
                    tree.SetBranchStatus ("var_PSSCALE_WEIGHT_5_UP", 1);
                    tree.SetBranchAddress("var_PSSCALE_WEIGHT_5_UP", &systVarWeight_PSSCALE_WEIGHT_5_UP);
                    tree.SetBranchStatus ("var_BTAG_PSSCALE_WEIGHT_5_UP", 1);
                    tree.SetBranchAddress("var_BTAG_PSSCALE_WEIGHT_5_UP", &systVarWeight_BTAG_PSSCALE_WEIGHT_5_UP);
                }
                if (sys.variation()==Systematic::Variation::down){
                    tree.SetBranchStatus ("var_PSSCALE_WEIGHT_5_DOWN", 1);
                    tree.SetBranchAddress("var_PSSCALE_WEIGHT_5_DOWN", &systVarWeight_PSSCALE_WEIGHT_5_DOWN);
                    tree.SetBranchStatus ("var_BTAG_PSSCALE_WEIGHT_5_DOWN", 1);
                    tree.SetBranchAddress("var_BTAG_PSSCALE_WEIGHT_5_DOWN", &systVarWeight_BTAG_PSSCALE_WEIGHT_5_DOWN);
                }
            }
            break;
        case Systematic::alphasPdf:
            if(sys.variation()==Systematic::Variation::up){
                tree.SetBranchStatus ("var_PDF_ALPHAS_UP", 1);
                tree.SetBranchAddress("var_PDF_ALPHAS_UP", &systVarWeight_PDF_ALPHAS_UP);
                tree.SetBranchStatus ("var_BTAG_PDF_ALPHAS_UP", 1);
                tree.SetBranchAddress("var_BTAG_PDF_ALPHAS_UP", &systVarWeight_BTAG_PDF_ALPHAS_UP);
            }
            if(sys.variation()==Systematic::Variation::down){
                tree.SetBranchStatus ("var_PDF_ALPHAS_DOWN", 1);
                tree.SetBranchAddress("var_PDF_ALPHAS_DOWN", &systVarWeight_PDF_ALPHAS_DOWN);
                tree.SetBranchStatus ("var_BTAG_PDF_ALPHAS_DOWN", 1);
                tree.SetBranchAddress("var_BTAG_PDF_ALPHAS_DOWN", &systVarWeight_BTAG_PDF_ALPHAS_DOWN);
            }
            break;
        case Systematic::pu:
            if(sys.variation()==Systematic::Variation::up){
                tree.SetBranchStatus ("var_PU_UP", 1);
                tree.SetBranchAddress("var_PU_UP", &systVarWeight_PU_UP);
                tree.SetBranchStatus ("var_BTAG_PU_UP", 1);
                tree.SetBranchAddress("var_BTAG_PU_UP", &systVarWeight_BTAG_PU_UP);
            }
            if(sys.variation()==Systematic::Variation::down){
                tree.SetBranchStatus ("var_PU_DOWN", 1);
                tree.SetBranchAddress("var_PU_DOWN", &systVarWeight_PU_DOWN);
                tree.SetBranchStatus ("var_BTAG_PU_DOWN", 1);
                tree.SetBranchAddress("var_BTAG_PU_DOWN", &systVarWeight_BTAG_PU_DOWN);
            }
            break;
        case Systematic::trig:
            if(sys.variation()==Systematic::Variation::up){
                tree.SetBranchStatus ("var_TRIG_UP", 1);
                tree.SetBranchAddress("var_TRIG_UP", &systVarWeight_TRIG_UP);
                tree.SetBranchStatus ("var_BTAG_TRIG_UP", 1);
                tree.SetBranchAddress("var_BTAG_TRIG_UP", &systVarWeight_BTAG_TRIG_UP);
            }
            if(sys.variation()==Systematic::Variation::down){
                tree.SetBranchStatus ("var_TRIG_DOWN", 1);
                tree.SetBranchAddress("var_TRIG_DOWN", &systVarWeight_TRIG_DOWN);
                tree.SetBranchStatus ("var_BTAG_TRIG_DOWN", 1);
                tree.SetBranchAddress("var_BTAG_TRIG_DOWN", &systVarWeight_BTAG_TRIG_DOWN);
            }
            break;
        case Systematic::eleID:
            if(sys.variation()==Systematic::Variation::up){
                tree.SetBranchStatus ("var_ELE_ID_UP", 1);
                tree.SetBranchAddress("var_ELE_ID_UP", &systVarWeight_ELE_ID_UP);
                tree.SetBranchStatus ("var_BTAG_ELE_ID_UP", 1);
                tree.SetBranchAddress("var_BTAG_ELE_ID_UP", &systVarWeight_BTAG_ELE_ID_UP);
            }
            if(sys.variation()==Systematic::Variation::down){
                tree.SetBranchStatus ("var_ELE_ID_DOWN", 1);
                tree.SetBranchAddress("var_ELE_ID_DOWN", &systVarWeight_ELE_ID_DOWN);
                tree.SetBranchStatus ("var_BTAG_ELE_ID_DOWN", 1);
                tree.SetBranchAddress("var_BTAG_ELE_ID_DOWN", &systVarWeight_BTAG_ELE_ID_DOWN);
            }
            break;
        case Systematic::eleReco:
            if(sys.variation()==Systematic::Variation::up){
                tree.SetBranchStatus ("var_ELE_RECO_UP", 1);
                tree.SetBranchAddress("var_ELE_RECO_UP", &systVarWeight_ELE_RECO_UP);
                tree.SetBranchStatus ("var_BTAG_ELE_RECO_UP", 1);
                tree.SetBranchAddress("var_BTAG_ELE_RECO_UP", &systVarWeight_BTAG_ELE_RECO_UP);
        }
            if(sys.variation()==Systematic::Variation::down){
                tree.SetBranchStatus ("var_ELE_RECO_DOWN", 1);
                tree.SetBranchAddress("var_ELE_RECO_DOWN", &systVarWeight_ELE_RECO_DOWN);
                tree.SetBranchStatus ("var_BTAG_ELE_RECO_DOWN", 1);
                tree.SetBranchAddress("var_BTAG_ELE_RECO_DOWN", &systVarWeight_BTAG_ELE_RECO_DOWN);
            }
            break;
        case Systematic::muonID:
            if(sys.variation()==Systematic::Variation::up){
                tree.SetBranchStatus ("var_MUON_ID_UP", 1);
                tree.SetBranchAddress("var_MUON_ID_UP", &systVarWeight_MUON_ID_UP);
                tree.SetBranchStatus ("var_BTAG_MUON_ID_UP", 1);
                tree.SetBranchAddress("var_BTAG_MUON_ID_UP", &systVarWeight_BTAG_MUON_ID_UP);
            }
            if(sys.variation()==Systematic::Variation::down){
                tree.SetBranchStatus ("var_MUON_ID_DOWN", 1);
                tree.SetBranchAddress("var_MUON_ID_DOWN", &systVarWeight_MUON_ID_DOWN);
                tree.SetBranchStatus ("var_BTAG_MUON_ID_DOWN", 1);
                tree.SetBranchAddress("var_BTAG_MUON_ID_DOWN", &systVarWeight_BTAG_MUON_ID_DOWN);
            }
            break;
        case Systematic::muonIso:
            if(sys.variation()==Systematic::Variation::up){
                tree.SetBranchStatus ("var_MUON_ISO_UP", 1);
                tree.SetBranchAddress("var_MUON_ISO_UP", &systVarWeight_MUON_ISO_UP);
                tree.SetBranchStatus ("var_BTAG_MUON_ISO_UP", 1);
                tree.SetBranchAddress("var_BTAG_MUON_ISO_UP", &systVarWeight_BTAG_MUON_ISO_UP);
            }
            if(sys.variation()==Systematic::Variation::down){
                tree.SetBranchStatus ("var_MUON_ISO_DOWN", 1);
                tree.SetBranchAddress("var_MUON_ISO_DOWN", &systVarWeight_MUON_ISO_DOWN);
                tree.SetBranchStatus ("var_BTAG_MUON_ISO_DOWN", 1);
                tree.SetBranchAddress("var_BTAG_MUON_ISO_DOWN", &systVarWeight_BTAG_MUON_ISO_DOWN);
            }
            break;
        case Systematic::l1prefiring:
            if(sys.variation()==Systematic::Variation::up){
                tree.SetBranchStatus ("var_L1PREFIRING_UP", 1);
                tree.SetBranchAddress("var_L1PREFIRING_UP", &systVarWeight_L1PREFIRING_UP);
                tree.SetBranchStatus ("var_BTAG_L1PREFIRING_UP", 1);
                tree.SetBranchAddress("var_BTAG_L1PREFIRING_UP", &systVarWeight_BTAG_L1PREFIRING_UP);
            }
            if(sys.variation()==Systematic::Variation::down){
                tree.SetBranchStatus("var_L1PREFIRING_DOWN", 1);
                tree.SetBranchAddress("var_L1PREFIRING_DOWN", &systVarWeight_L1PREFIRING_DOWN);
                tree.SetBranchStatus ("var_BTAG_L1PREFIRING_DOWN", 1);
                tree.SetBranchAddress("var_BTAG_L1PREFIRING_DOWN", &systVarWeight_BTAG_L1PREFIRING_DOWN);
            }
            break;
        case Systematic::btag:
            if(sys.variation()==Systematic::Variation::up){
                tree.SetBranchAddress("var_BTAG_UP", &systVarWeight_BTAG_UP);
                tree.SetBranchStatus ("var_BTAG_UP", 1);
            }
            if(sys.variation()==Systematic::Variation::down){
                tree.SetBranchAddress("var_BTAG_DOWN", &systVarWeight_BTAG_DOWN);
                tree.SetBranchStatus ("var_BTAG_DOWN", 1);
            }
            break;
        case Systematic::btagLjet:
            if(sys.variation()==Systematic::Variation::up){
                tree.SetBranchAddress("var_BTAG_LJET_UP", &systVarWeight_BTAG_LJET_UP);
                tree.SetBranchStatus ("var_BTAG_LJET_UP", 1);
            }
            if(sys.variation()==Systematic::Variation::down){
                tree.SetBranchAddress("var_BTAG_LJET_DOWN", &systVarWeight_BTAG_LJET_DOWN);
                tree.SetBranchStatus ("var_BTAG_LJET_DOWN", 1);
            }
            break;
        default:
            tree.SetBranchStatus("var_*", 0);
            break;

    }
}
