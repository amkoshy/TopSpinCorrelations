#include "Math/GenVector/VectorUtil.h"

#include <MiniTreeReader.h>
#include <MiniTreeHelpers.h>
#include <sampleHelpers.h>
#include "../../common/include/utils.h"
#include "../../common/include/classes.h"
#include "../../common/include/classesFwd.h"
#include "PlotterConfigurationHelper.h"
#include "RootFileReader.h"




void Event::Clear(){
    weight = 1.;
    weightWithBtagSF=1.;
    genLevelWeight=1.;
    // hasLKRS=false;
    hasKRS=false;
    passStep3=false;
    // passStep4=false;
    // passStep5=false;
    // passStep6=false;
    // passStep7=false;
    // passStep8=false;
    // passStep8Loose=false;
    jets.clear();
    // jets.reserve(10);
    bjets.clear();
    // bjets.reserve(10);
    lepton1.SetCoordinates(0.,0.,0.,0.);
    lepton2.SetCoordinates(0.,0.,0.,0.);
    // leptons.clear();
    dilepton.SetCoordinates(0.,0.,0.,0.);
    // met.SetCoordinates(0.,0.,0.,0.);
    mlb_min = 0.;
    kinReco_top.SetCoordinates(0.,0.,0.,0.);
    kinReco_antitop.SetCoordinates(0.,0.,0.,0.);
    // kinReco_lepton.SetCoordinates(0.,0.,0.,0.);
    // kinReco_antilepton.SetCoordinates(0.,0.,0.,0.);
    // kinReco_b.SetCoordinates(0.,0.,0.,0.);
    // kinReco_antib.SetCoordinates(0.,0.,0.,0.);
    // kinReco_nonbjets.clear();
    // kinReco_ttbar.SetCoordinates(0.,0.,0.,0.);
    kinReco_rho=-999.;
    // kinReco_SF=1.;
    // looseKinReco_ttbar.SetCoordinates(0.,0.,0.,0.);
    // looseKinReco_lepton1.SetCoordinates(0.,0.,0.,0.);
    // looseKinReco_lepton2.SetCoordinates(0.,0.,0.,0.);
    // looseKinReco_b1.SetCoordinates(0.,0.,0.,0.);
    // looseKinReco_b2.SetCoordinates(0.,0.,0.,0.);
    // looseKinReco_nonbjets.clear();
    // looseKinReco_rho=-999.;
    // looseKinReco_SF=1.;
    // gen_additional_jets.clear();
    // gen_additional_jets.reserve(1);
    gen_additional_jet.SetCoordinates(0.,0.,0.,0.);

    gen_top.SetCoordinates(0.,0.,0.,0.);
    gen_antitop.SetCoordinates(0.,0.,0.,0.);
    // gen_lepton.SetCoordinates(0.,0.,0.,0.);
    // gen_antilepton.SetCoordinates(0.,0.,0.,0.);
    // gen_b.SetCoordinates(0.,0.,0.,0.);
    // gen_antib.SetCoordinates(0.,0.,0.,0.);
    // gen_ttbar.SetCoordinates(0.,0.,0.,0.);
    // gen_met.SetCoordinates(0.,0.,0.,0.);
    gen_rho=-999.;

    isTTJSignal=false;
    isTTSignalBackground=false;

    DNN_rho=-999.;
    final_rho=-999.;
    // DNN_rho_new=-999.;

    tempJet.SetCoordinates(0.,0.,0.,0.);
}

void Event::Fill(const MiniTreeReader &reader, const bool isMC, const bool isSignal, const Systematic::Systematic &sys, const bool hasPSWeights,
                const bool hasPDFWeights, const bool hasMEWeights, const bool hasBFragWeights, const bool hasBSemilepWeights){
    (void) hasPDFWeights;
    if(isMC){
        if(sys.type()==Systematic::nominal || sys.type()==Systematic::jes
            || sys.type()==Systematic::jer || sys.type()==Systematic::mass
            || sys.type()==Systematic::jerEta0 || sys.type()==Systematic::jerEta1
            || sys.type()==Systematic::ueTune || sys.type()==Systematic::unclustered
            || sys.type()==Systematic::match || sys.type()==Systematic::psFSRScale
            || sys.type()==Systematic::psISRScale){
            genLevelWeight = reader.weight;
            weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
            weightWithBtagSF=weight * reader.btagSF;
        }//TODO CHECK WEIGHTS FOR SYSTEMATICS
        // recoVariations
        if(sys.type()==Systematic::eleID){
            if (sys.variation()==Systematic::Variation::up){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.systVarWeight_ELE_ID_UP * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_ELE_ID_UP;
            }
            if (sys.variation()==Systematic::Variation::down){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.systVarWeight_ELE_ID_DOWN * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_ELE_ID_DOWN;
            }
        }
        if(sys.type()==Systematic::eleReco){
            if (sys.variation()==Systematic::Variation::up){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.systVarWeight_ELE_RECO_UP * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_ELE_RECO_UP;
            }
            if (sys.variation()==Systematic::Variation::down){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.systVarWeight_ELE_RECO_DOWN * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_ELE_RECO_DOWN;
            }
        }
        if(sys.type()==Systematic::muonID){
            if (sys.variation()==Systematic::Variation::up){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.systVarWeight_MUON_ID_UP * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_MUON_ID_UP;
            }
            if (sys.variation()==Systematic::Variation::down){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.systVarWeight_MUON_ID_DOWN * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_MUON_ID_DOWN;
            }
        }
        if(sys.type()==Systematic::muonIso){
            if (sys.variation()==Systematic::Variation::up){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.systVarWeight_MUON_ISO_UP * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_MUON_ISO_UP;
            }
            if (sys.variation()==Systematic::Variation::down){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.systVarWeight_MUON_ISO_DOWN * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_MUON_ISO_DOWN;
            }
        }
        if(sys.type()==Systematic::pu){
            if (sys.variation()==Systematic::Variation::up){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.leptonSF * reader.systVarWeight_PU_UP * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_PU_UP;
            }
            if (sys.variation()==Systematic::Variation::down){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.leptonSF * reader.systVarWeight_PU_DOWN * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_PU_DOWN;
            }
        }
        if(sys.type()==Systematic::trig){
            if (sys.variation()==Systematic::Variation::up){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.leptonSF * reader.pileupSF * reader.systVarWeight_TRIG_UP * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_TRIG_UP;
            }
            if (sys.variation()==Systematic::Variation::down){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.leptonSF * reader.pileupSF * reader.systVarWeight_TRIG_DOWN * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_TRIG_DOWN;
            }
        }
        if(sys.type()==Systematic::btag){
            if (sys.variation()==Systematic::Variation::up){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.leptonSF * reader.pileupSF * reader.triggerSF *  reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_UP;
            }
            if (sys.variation()==Systematic::Variation::down){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_DOWN ;
            }
        }
        if(sys.type()==Systematic::btagLjet){
            if (sys.variation()==Systematic::Variation::up){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight  * reader.systVarWeight_BTAG_LJET_UP ;
            }
            if (sys.variation()==Systematic::Variation::down){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_LJET_DOWN ;
            }
        }
        if(sys.type()==Systematic::l1prefiring){
            if (sys.variation()==Systematic::Variation::up){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.systVarWeight_L1PREFIRING_UP;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_L1PREFIRING_UP;
            }
            if (sys.variation()==Systematic::Variation::down){
                genLevelWeight = reader.weight;
                weight=reader.weight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.systVarWeight_L1PREFIRING_DOWN;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_L1PREFIRING_DOWN;
            }
        }
        //genVariations
        if(sys.type()==Systematic::meRenScale){
            if (sys.variation()==Systematic::Variation::up){
                genLevelWeight = reader.weight;
                // if(isSignal) genLevelWeight = genLevelWeight * reader.systVarWeight_MERENSCALE_UP);
                if(hasMEWeights) genLevelWeight = genLevelWeight * reader.systVarWeight_MERENSCALE_UP;
                weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_MERENSCALE_UP;
            }
            if (sys.variation()==Systematic::Variation::down){
                genLevelWeight =reader.weight;
                if(hasMEWeights) genLevelWeight = genLevelWeight * reader.systVarWeight_MERENSCALE_DOWN;
                weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_MERENSCALE_DOWN;
            }
        }
        if(sys.type()==Systematic::meFacScale){
            if (sys.variation()==Systematic::Variation::up){
                genLevelWeight = reader.weight;
                if(hasMEWeights) genLevelWeight = genLevelWeight * reader.systVarWeight_MEFACSCALE_UP;
                weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_MEFACSCALE_UP;
            }
            if (sys.variation()==Systematic::Variation::down){
                genLevelWeight =reader.weight;
                if(hasMEWeights) genLevelWeight = genLevelWeight * reader.systVarWeight_MEFACSCALE_DOWN;
                weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_MEFACSCALE_DOWN;
            }
        }
        if(sys.type()==Systematic::meScale){
            if (sys.variation()==Systematic::Variation::up){
                genLevelWeight = reader.weight;
                if(hasMEWeights) genLevelWeight = genLevelWeight * reader.systVarWeight_MESCALE_UP;
                weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_MESCALE_UP;
            }
            if (sys.variation()==Systematic::Variation::down){
                genLevelWeight =reader.weight;
                if(hasMEWeights) genLevelWeight = genLevelWeight * reader.systVarWeight_MESCALE_DOWN;
                weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_MESCALE_DOWN;
            }
        }
        if(sys.type()==Systematic::bSemilep){
            if (sys.variation()==Systematic::Variation::up){
                genLevelWeight = reader.weight;
                if(hasBSemilepWeights) genLevelWeight = genLevelWeight * reader.systVarWeight_BSEMILEP_UP;
                weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_BSEMILEP_UP;
            }
            if (sys.variation()==Systematic::Variation::down){
                genLevelWeight =reader.weight;
                if(hasBSemilepWeights) genLevelWeight = genLevelWeight * reader.systVarWeight_BSEMILEP_DOWN;
                weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_BSEMILEP_DOWN;
            }
        }
        if(sys.type()==Systematic::alphasPdf){
            if (sys.variation()==Systematic::Variation::up){
                genLevelWeight = reader.weight;
                if(hasMEWeights) genLevelWeight = genLevelWeight * reader.systVarWeight_PDF_ALPHAS_UP;
                weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_PDF_ALPHAS_UP;
            }
            if (sys.variation()==Systematic::Variation::down){
                genLevelWeight =reader.weight;
                if(hasMEWeights) genLevelWeight = genLevelWeight * reader.systVarWeight_PDF_ALPHAS_DOWN;
                weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_PDF_ALPHAS_DOWN;
            }
        }
        if(sys.type()==Systematic::bFrag){
            if (sys.variation()==Systematic::Variation::up){
                genLevelWeight = reader.weight;
                if(hasBFragWeights) genLevelWeight = genLevelWeight * reader.systVarWeight_BFRAG_UP;
                weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_BFRAG_UP;
            }
            if (sys.variation()==Systematic::Variation::down){
                genLevelWeight =reader.weight;
                if(hasBFragWeights) genLevelWeight = genLevelWeight * reader.systVarWeight_BFRAG_DOWN;
                weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_BFRAG_DOWN;
            }
        }
        if(sys.type()==Systematic::psScaleWeight){
                if(sys.variationNumber()==4){
                    if (sys.variation()==Systematic::Variation::up){
                        genLevelWeight = reader.weight;
                        // if(hasPSWeights) genLevelWeight = genLevelWeight * (reader.systVarWeight_PS->at(6));
                        if(hasPSWeights) genLevelWeight = genLevelWeight * (reader.systVarWeight_PSSCALE_WEIGHT_4_UP);
                        // if(hasPSWeights && reader.systVarWeight_PS->size()>5) genLevelWeight = genLevelWeight * (reader.systVarWeight_PS->at(6));
                        weight = genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                        weightWithBtagSF = weight * reader.systVarWeight_BTAG_PSSCALE_WEIGHT_4_UP;
                    }
                    if (sys.variation()==Systematic::Variation::down){
                        genLevelWeight =reader.weight;
                        // if(hasPSWeights) genLevelWeight = genLevelWeight * (reader.systVarWeight_PS->at(8));
                        if(hasPSWeights) genLevelWeight = genLevelWeight * (reader.systVarWeight_PSSCALE_WEIGHT_4_DOWN);
                        weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                        weightWithBtagSF=weight * reader.systVarWeight_BTAG_PSSCALE_WEIGHT_4_DOWN;
                    }
                }
                if(sys.variationNumber()==5){
                    if (sys.variation()==Systematic::Variation::up){
                        genLevelWeight = reader.weight;
                        if(hasPSWeights) genLevelWeight = genLevelWeight * (reader.systVarWeight_PSSCALE_WEIGHT_5_UP);
                        weight = genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                        weightWithBtagSF = weight * reader.systVarWeight_BTAG_PSSCALE_WEIGHT_5_UP;
                    }
                    if (sys.variation()==Systematic::Variation::down){
                        genLevelWeight =reader.weight;
                        if(hasPSWeights) genLevelWeight = genLevelWeight * (reader.systVarWeight_PSSCALE_WEIGHT_5_DOWN);
                        weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                        weightWithBtagSF=weight * reader.systVarWeight_BTAG_PSSCALE_WEIGHT_5_DOWN;
                    }
                }
        }
        if(sys.type()==Systematic::bFrag_Peterson){
                genLevelWeight = reader.weight;
                if(hasBFragWeights) genLevelWeight = genLevelWeight * reader.systVarWeight_BFRAG_PETERSON;
                weight=genLevelWeight * reader.leptonSF * reader.pileupSF * reader.triggerSF * reader.l1PrefiringWeight;
                weightWithBtagSF=weight * reader.systVarWeight_BTAG_BFRAG_PETERSON;
            }
    }else{
        genLevelWeight = 1.;
        weight = 1.;
    }

    passStep3 = reader.passStep3;
    // passStep4 = reader.passStep4;
    // passStep5 = reader.passStep5;
    // passStep6 = reader.passStep6;
    // passStep7 = reader.passStep7;
    // passStep8 = reader.passStep8;
    // passStep8Loose = reader.passStep8Loose;
    // hasLKRS=reader.hasLooseKinRecoSolution;
    hasKRS=reader.hasKinRecoSolution;

    // LV tempJet;
    for(unsigned int i = 0; i < reader.n_jets; i++){
        tempJet.SetCoordinates(reader.jets_pt[i],reader.jets_eta[i],reader.jets_phi[i],reader.jets_m[i]);
        jets.push_back(tempJet);
        if(reader.jets_btag[i]>0){
            bjets.push_back(tempJet);
        }
    }
    // for(unsigned int i = 0; i < reader.bjets_pt->size(); i++){
    //     tempJet.SetCoordinates(reader.bjets_pt->at(i),reader.bjets_eta->at(i),reader.bjets_phi->at(i),reader.bjets_m->at(i));
    //     bjets.push_back(tempJet);
    // }



    lepton1.SetCoordinates(reader.lepton1_pt ,reader.lepton1_eta,reader.lepton1_phi,reader.lepton1_m);
    lepton2.SetCoordinates(reader.lepton2_pt,reader.lepton2_eta,reader.lepton2_phi,reader.lepton2_m);
    dilepton=lepton1+lepton2;
    met.SetCoordinates(reader.met_pt,0.,reader.met_phi,0.);

    mlb_min = 9999.;
    float comb,comb1,comb2;
    for (unsigned int i = 0; i < bjets.size(); i++){
        comb1 = (lepton1+bjets.at(i)).M();
        comb2 = (lepton2+bjets.at(i)).M();
        comb = comb1<comb2? comb1 : comb2;
        if(comb<mlb_min){
            mlb_min=comb;
        }
    }

    if(hasKRS){
    //     kinReco_SF = reader.kinReco_weight;
    //
        kinReco_top.SetCoordinates(reader.kinReco_top_pt,reader.kinReco_top_eta,reader.kinReco_top_phi,reader.kinReco_top_m);
        kinReco_antitop.SetCoordinates(reader.kinReco_antitop_pt,reader.kinReco_antitop_eta,reader.kinReco_antitop_phi,reader.kinReco_antitop_m);
    //     kinReco_lepton.SetCoordinate(sreader.kinReco_lepton_pt,reader.kinReco_lepton_eta,reader.kinReco_lepton_phi,reader.kinReco_lepton_m);
    //     kinReco_antilepton.SetCoordinate(sreader.kinReco_antilepton_pt,reader.kinReco_antilepton_eta,reader.kinReco_antilepton_phi,reader.kinReco_antilepton_m);
    //     kinReco_b.SetCoordinates(reader.kinReco_b_pt,reader.kinReco_b_eta,reader.kinReco_b_phi,reader.kinReco_b_m);
    //     kinReco_antib.SetCoordinates(reader.kinReco_antib_pt,reader.kinReco_antib_eta,reader.kinReco_antib_phi,reader.kinReco_antib_m);
    //
    //     for(unsigned int i = 0; i < reader.kinReco_nonbjets_pt).size(); i++){
    //         tempJet.SetCoordinates(reader.kinReco_nonbjets_pt).at(i),reader.kinReco_nonbjets_eta).at(i),reader.kinReco_nonbjets_phi).at(i),reader.kinReco_nonbjets_m).at(i));
    //         kinReco_nonbjets.push_back(tempJet);
    //     }
        // kinReco_ttbar.SetCoordinates((kinReco_top+kinReco_antitop).Pt(),(kinReco_top+kinReco_antitop).Eta(),(kinReco_top+kinReco_antitop).Phi(),(kinReco_top+kinReco_antitop).M());
    //
    //     if(kinReco_nonbjets.size()>0){
    //         kinReco_rho = 340./((kinReco_nonbjets.at(0)+kinReco_ttbar).M());
    //     }else{
    //         kinReco_rho=-999.;
    //     }
        if(reader.kinReco_rho>0.){
            kinReco_nonbjet.SetCoordinates(reader.kinReco_nonbjet_pt,reader.kinReco_nonbjet_eta,reader.kinReco_nonbjet_phi,reader.kinReco_nonbjet_m);
            kinReco_rho = reader.kinReco_rho;
        }else{
            kinReco_nonbjet.SetCoordinates(0.,0.,0.,0.);
            kinReco_rho = -999.;
        }
    }
    // if(hasLKRS){
    //     looseKinReco_SF = reader.looseKinReco_weight;
    //     looseKinReco_ttbar.SetCoordinates(reader.looseKinReco_ttbar_pt,reader.looseKinReco_ttbar_eta,reader.looseKinReco_ttbar_phi,reader.looseKinReco_ttbar_m);
    //     looseKinReco_lepton1.SetCoordinates(reader.looseKinReco_l1_pt,reader.looseKinReco_l1_eta,reader.looseKinReco_l1_phi,reader.looseKinReco_l1_m);
    //     looseKinReco_lepton2.SetCoordinates(reader.looseKinReco_l2_pt,reader.looseKinReco_l2_eta,reader.looseKinReco_l2_phi,reader.looseKinReco_l2_m);
    //     looseKinReco_b1.SetCoordinates(reader.looseKinReco_b1_pt,reader.looseKinReco_b1_eta,reader.looseKinReco_b1_phi,reader.looseKinReco_b1_m);
    //     looseKinReco_b2.SetCoordinates(reader.looseKinReco_b2_pt,reader.looseKinReco_b2_eta,reader.looseKinReco_b2_phi,reader.looseKinReco_b2_m);
    //
    //     for(unsigned int i = 0; i < reader.looseKinReco_nonbjets_pt).size(); i++){
    //         tempJet.SetCoordinates(reader.looseKinReco_nonbjets_pt).at(i),reader.looseKinReco_nonbjets_eta).at(i),reader.looseKinReco_nonbjets_phi).at(i),reader.looseKinReco_nonbjets_m).at(i));
    //         looseKinReco_nonbjets.push_back(tempJet);
    //     }
    //     if(looseKinReco_nonbjets.size()>0){
    //         looseKinReco_rho = 340./((looseKinReco_nonbjets.at(0)+looseKinReco_ttbar).M());
    //     }else{
    //         looseKinReco_rho=-999.;
    //     }
    // }


    // DNN_rho = reader.DNN_rho;


    if(isSignal){

        gen_top.SetCoordinates(reader.gen_top_pt,reader.gen_top_eta,reader.gen_top_phi,reader.gen_top_m);
        gen_antitop.SetCoordinates(reader.gen_antitop_pt,reader.gen_antitop_eta,reader.gen_antitop_phi,reader.gen_antitop_m);
        // gen_lepton.SetCoordinates(reader.gen_lepton_pt,reader.gen_lepton_eta,reader.gen_lepton_phi,reader.gen_lepton_m);
        // gen_antilepton.SetCoordinates(reader.gen_antilepton_pt,reader.gen_antilepton_eta,reader.gen_antilepton_phi,reader.gen_antilepton_m);
        // gen_b.SetCoordinates(reader.gen_b_pt,reader.gen_b_eta,reader.gen_b_phi,reader.gen_b_m);
        // gen_antib.SetCoordinates(reader.gen_antib_pt,reader.gen_antib_eta,reader.gen_antib_phi,reader.gen_antib_m);
        // gen_ttbar.SetCoordinates((gen_top+gen_antitop).Pt(),(gen_top+gen_antitop).Eta(),(gen_top+gen_antitop).Phi(),(gen_top+gen_antitop).M());
        // gen_met.SetCoordinates(reader.gen_met_pt,reader.gen_met_eta,reader.gen_met_phi,reader.gen_met_m);

        // bool foundAddJet = false;
        // for(unsigned int i = 0; i < reader.gen_additional_jets_pt->size(); i++){
        //     if(!foundAddJet){
        //         gen_additional_jet.SetCoordinates(reader.gen_additional_jets_pt->at(i),reader.gen_additional_jets_eta->at(i),reader.gen_additional_jets_phi->at(i),reader.gen_additional_jets_m->at(i));
        //         foundAddJet=true;
        //     }
        // }

        // if(foundAddJet){
        //   if(gen_additional_jet.Pt()>50.){
        //     if(fabs(gen_additional_jet.Eta())<2.4){
        //           gen_rho =  reader.gen_rho;
        //
        //           isTTJSignal = true && isSignal;
        //         }else{
        //             gen_rho =  -999.;
        //             isTTSignalBackground = true;
        //         }
        //     }else{
        //         gen_rho =  -999.;
        //         isTTSignalBackground = true;
        //     }
        // }else{
        //     gen_rho =  -999.;
        //     isTTSignalBackground = true;
        // }

        if(reader.gen_partonLevel_rho>0.){
            gen_additional_jet.SetCoordinates(reader.gen_partonLevel_additional_jet_pt,reader.gen_partonLevel_additional_jet_eta,reader.gen_partonLevel_additional_jet_phi,reader.gen_partonLevel_additional_jet_m);
            if(gen_additional_jet.Pt()>50.){
                if(gen_additional_jet.Eta()<2.4){
                    gen_rho = reader.gen_partonLevel_rho;
                    isTTJSignal = true && isSignal;
                    isTTSignalBackground = false;
                }else{
                    gen_rho = -999.;
                    isTTJSignal = false;
                    isTTSignalBackground = true;
                }
            }else{
                gen_rho = -999.;
                isTTJSignal = false;
                isTTSignalBackground = true;
            }
            isTTJSignal = true;
        }else{
            gen_rho = -999.;
            isTTJSignal = false;
            isTTSignalBackground = true;
            gen_additional_jet.SetCoordinates(0.,0.,0.,0.);
        }
    }
}

BTagCategory::BTagCategory BTagCategory::convert(const TString& btag)
{
    if(btag=="ZeroBTag"){ return ZeroBTag;}
    else if(btag=="OneBTag"){ return OneBTag;}
    else if(btag=="TwoBTag"){ return TwoBTag;}
    else if(btag=="ThreeBTag"){ return ThreeBTag;}
    else if(btag=="FourBtag"){ return FourBtag;}
    else if(btag=="InclusiveBTag"){ return InclusiveBTag;}
    else if(btag=="GreaterZeroBTag"){ return GreaterZeroBTag;}
    else if(btag=="GreaterOneBTag"){ return GreaterOneBTag;}
    else if(btag=="GreaterTwoBTag"){ return GreaterTwoBTag;}
    else if(btag=="ZeroPlusGreaterTwoBTag"){ return ZeroPlusGreaterTwoBTag;}
    else{
      std::cerr<<"Error in BTagCategory::convert()! Following conversion is not implemented: "<<btag<<"\n...break\n"<<std::endl;
      exit(98);
    }
}

TString BTagCategory::convert(const BTagCategory& btag)
{
    switch(btag){
        case(ZeroBTag): return "ZeroBTag";
        case(OneBTag): return "OneBTag";
        case(TwoBTag): return "TwoBTag";
        case(ThreeBTag): return "ThreeBTag";
        case(FourBtag): return "FourBtag";
        case(InclusiveBTag): return "InclusiveBTag";
        case(GreaterZeroBTag): return "GreaterZeroBTag";
        case(GreaterOneBTag): return "GreaterOneBTag";
        case(GreaterTwoBTag): return "GreaterTwoBTag";
        case(ZeroPlusGreaterTwoBTag): return "ZeroPlusGreaterTwoBTag";
        default:
          std::cerr<<"Error in BTagCategory::convert()! Conversion is not implemented\n...break\n"<<std::endl;
          exit(98);
          break;
    }
}

std::vector<BTagCategory::BTagCategory> BTagCategory::convert(const std::vector<TString>& btags)
{
    std::vector<BTagCategory> v_btag;
    for(const auto& btag : btags) v_btag.push_back(convert(btag));
    return v_btag;
}

std::vector<BTagCategory::BTagCategory> BTagCategory::convert(const std::vector<std::string>& btags)
{
    std::vector<BTagCategory> v_btag;
    for(const auto& btag : btags) v_btag.push_back(convert(btag));
    return v_btag;
}

std::vector<TString> BTagCategory::convert(const std::vector<BTagCategory>& btags)
{
    std::vector<TString> v_btag;
    for(const auto& btag : btags) v_btag.push_back(convert(btag));
    return v_btag;
}

NJetCategory::NJetCategory NJetCategory::convert(const TString& njet)
{
        if(njet=="ZeroJet"){ return ZeroJet;}
        else if(njet=="OneJet"){ return OneJet;}
        else if(njet=="TwoJet"){ return TwoJet;}
        else if(njet=="ThreeJet"){ return ThreeJet;}
        else if(njet=="FourJet"){ return FourJet;}
        else if(njet=="InclusiveNJet"){ return InclusiveNJet;}
        else if(njet=="GreaterOneJet"){ return GreaterOneJet;}
        else if(njet=="GreaterTwoJet"){ return GreaterTwoJet;}
        else if(njet=="ZeroAndOneJet"){ return ZeroAndOneJet;}
        else{
          std::cerr<<"Error in NJetCategory::convert()! Following conversion is not implemented: "<<njet<<"\n...break\n"<<std::endl;
          exit(98);
    }
}

TString NJetCategory::convert(const NJetCategory& njet)
{
    switch(njet){
        case(ZeroJet): return "ZeroJet";
        case(OneJet): return "OneJet";
        case(TwoJet): return "TwoJet";
        case(ThreeJet): return "ThreeJet";
        case(FourJet): return "FourJet";
        case(InclusiveNJet): return "InclusiveNJet";
        case(GreaterOneJet): return "GreaterOneJet";
        case(GreaterTwoJet): return "GreaterTwoJet";
        case(ZeroAndOneJet): return "ZeroAndOneJet";
        default:
          std::cerr<<"Error in NJetCategory::convert()! Conversion is not implemented\n...break\n"<<std::endl;
          exit(98);
          break;
    }
}

std::vector<NJetCategory::NJetCategory> NJetCategory::convert(const std::vector<TString>& njets)
{
    std::vector<NJetCategory> v_njet;
    for(const auto& njet : njets) v_njet.push_back(convert(njet));
    return v_njet;
}

std::vector<NJetCategory::NJetCategory> NJetCategory::convert(const std::vector<std::string>& njets)
{
    std::vector<NJetCategory> v_njet;
    for(const auto& njet : njets) v_njet.push_back(convert(njet));
    return v_njet;
}

std::vector<TString> NJetCategory::convert(const std::vector<NJetCategory>& njets)
{
    std::vector<TString> v_njet;
    for(const auto& njet : njets) v_njet.push_back(convert(njet));
    return v_njet;
}






// NoRho, Rho0to01, Rho01to02, Rho02to03, Rho03to04, Rho04to05, Rho05to06, Rho06to07,
//  Rho07to08, Rho08to09, Rho09to1, undefined
RecoRhoBinCategory::RecoRhoBinCategory RecoRhoBinCategory::convert(const TString& rhobin)
{
    if(rhobin=="NoRho"){ return NoRho;}
    else if(rhobin=="InclusiveRho"){ return InclusiveRho;}
    else if(rhobin=="RhoStudies"){ return RhoStudies;}
    else if(rhobin=="Rho0to1"){ return Rho0to1;}
    // new binning
    else if(rhobin=="Rho0to02"){ return Rho0to02;}
    else if(rhobin=="Rho02to03"){ return Rho02to03;}
    else if(rhobin=="Rho03to04"){ return Rho03to04;}
    else if(rhobin=="Rho04to05"){ return Rho04to05;}
    else if(rhobin=="Rho05to06"){ return Rho05to06;}
    else if(rhobin=="Rho06to1"){ return Rho06to1;}
    else{
      std::cerr<<"Error in RecoRhoBinCategory::convert()! Following conversion is not implemented: "<<rhobin<<"\n...break\n"<<std::endl;
      exit(98);
    }
}
TString RecoRhoBinCategory::convert(const RecoRhoBinCategory& rhobin)
{
    switch(rhobin){
        case(NoRho): return "NoRho";
        case(InclusiveRho): return "InclusiveRho";
        case(RhoStudies): return "RhoStudies";
        case(Rho0to1): return "Rho0to1";
        // new binning
        case(Rho0to02): return "Rho0to02";
        case(Rho02to03): return "Rho02to03";
        case(Rho03to04): return "Rho03to04";
        case(Rho04to05): return "Rho04to05";
        case(Rho05to06): return "Rho05to06";
        case(Rho06to1): return "Rho06to1";
        default:
          std::cerr<<"Error in RecoRhoBinCategory::convert()! Conversion is not implemented\n...break\n"<<std::endl;
          exit(98);
          break;
    }
}

pair<float,float> RecoRhoBinCategory::getLowerUpperEdge(const RecoRhoBinCategory& rhobin)
{
    switch(rhobin){
        case(NoRho): return pair<float,float> (0.,99999999.);
        case(InclusiveRho): return pair<float,float> (-1.,999.);
        case(RhoStudies): return pair<float,float> (-1.,999.);
        case(Rho0to1): return pair<float,float> (0.,999.);
        // new binning
        case(Rho0to02): return pair<float,float> (0.,0.2);
        case(Rho02to03): return pair<float,float> (0.2,0.3);
        case(Rho03to04): return pair<float,float> (0.3,0.4);
        case(Rho04to05): return pair<float,float> (0.4,0.5);
        case(Rho05to06): return pair<float,float> (0.5,0.6);
        case(Rho06to1): return pair<float,float> (0.6,999.);
        default:
          std::cerr<<"Error in RecoRhoBinCategory::getLowerUpperEdge()! Conversion is not implemented\n...break\n"<<std::endl;
          exit(98);
          break;
    }
}

std::vector<RecoRhoBinCategory::RecoRhoBinCategory> RecoRhoBinCategory::convert(const std::vector<TString>& rhobins)
{
    std::vector<RecoRhoBinCategory> v_rhobin;
    for(const auto& rhobin : rhobins) v_rhobin.push_back(convert(rhobin));
    return v_rhobin;
}

std::vector<RecoRhoBinCategory::RecoRhoBinCategory> RecoRhoBinCategory::convert(const std::vector<std::string>& rhobins)
{
    std::vector<RecoRhoBinCategory> v_rhobin;
    for(const auto& rhobin : rhobins) v_rhobin.push_back(convert(rhobin));
    return v_rhobin;
}

std::vector<TString> RecoRhoBinCategory::convert(const std::vector<RecoRhoBinCategory>& rhobins)
{
    std::vector<TString> v_rhobin;
    for(const auto& rhobin : rhobins) v_rhobin.push_back(convert(rhobin));
    return v_rhobin;
}





Process::Process Process::convert(const TString& process)
{
    if(process=="ttbarsignal"){ return ttbarsignal;}
    else if(process=="ttbarsignal_genRho0to02"){ return ttbarsignal_genRho0to02;}
    else if(process=="ttbarsignal_genRho02to03"){ return ttbarsignal_genRho02to03;}
    else if(process=="ttbarsignal_genRho03to04"){ return ttbarsignal_genRho03to04;}
    else if(process=="ttbarsignal_genRho04to05"){ return ttbarsignal_genRho04to05;}
    else if(process=="ttbarsignal_genRho05to06"){ return ttbarsignal_genRho05to06;}
    else if(process=="ttbarsignal_genRho06to1"){ return ttbarsignal_genRho06to1;}
    else if(process=="ttbarsignal_genRho0to1"){ return ttbarsignal_genRho0to1;}
    else if(process=="ttbarsignal_BKGNoAddJet"){ return ttbarsignal_BKGNoAddJet;}
    else if(process=="ttbarbg"){ return ttbarbg;}
    else if(process=="ttbarX"){ return ttbarX;}
    else if(process=="dy"){ return dy;}
    else if(process=="singletop"){ return singletop;}
    else if(process=="diboson"){ return diboson;}
    else if(process=="wjets"){ return wjets;}
    else if(process=="data"){ return data;}
    else{
      std::cerr<<"Error in Process::convert()! Following conversion is not implemented: "<<process<<"\n...break\n"<<std::endl;
      exit(98);
    }
}

Process::Process Process::convertFromFilename(const TString& process)
{
    if(process.Contains("ttbarsignalplustau")){ return ttbarsignal;}
    else if(process.Contains("ttbarbg")){ return ttbarbg;}
    else if(process.Contains("ttbarW")||process.Contains("ttbarZ")||process.Contains("ttgjets")){ return ttbarX;}
    else if(process.Contains("dy")){ return dy;}
    else if(process.Contains("singletop")||process.Contains("singleantitop")){ return singletop;}
    else if(process.Contains("ww")||process.Contains("wz")||process.Contains("zz")){ return diboson;}
    else if(process.Contains("wtolnu")){ return wjets;}
    else if(process.Contains("run")){ return data;}
    else{
      std::cerr<<"Error in Process::convert()! Following conversion is not implemented: "<<process<<"\n...break\n"<<std::endl;
      exit(98);
    }
}

Process::Process Process::convertFromBinning(const float& bin){
    switch((int)(bin*1000.)){
    case(0):
        return ttbarsignal_genRho0to02;
    case(200):
        return ttbarsignal_genRho02to03;
    case(300):
        return ttbarsignal_genRho03to04;
    case(400):
        return ttbarsignal_genRho04to05;
    case(500):
        return ttbarsignal_genRho05to06;
    case(600):
        return ttbarsignal_genRho06to1;
    case(1000):
        return ttbarsignal_genRho06to1;
    default:
      std::cerr<<"Error in Process::convertFromBinning()! Conversion is not implemented "<<bin<<" "<<(int)(bin*1000.)<<" \n...break\n"<<std::endl;
      exit(98);
      break;
    }
    // if((int)(bin*1000.)==0){
    //     return ttbarsignal_genRho0to02;
    // }
    // else if((int)(bin*1000.)==200){
    //     return ttbarsignal_genRho02to03;
    // }
    // else if((int)(bin*1000.)==300){
    //     return ttbarsignal_genRho03to04;
    // }
    // else if((int)(bin*1000.)==400){
    //     return ttbarsignal_genRho04to05;
    // }
    // else if((int)(bin*1000.)==500){
    //     return ttbarsignal_genRho05to06;
    // }
    // else if((int)(bin*1000.)==600){
    //     return ttbarsignal_genRho06to1;
    // }
    // else if((int)(bin*1000.)==1000){
    //     return ttbarsignal_genRho06to1;
    // }else{
    //   std::cerr<<"Error in Process::convertFromBinning()! Conversion is not implemented "<<bin<<" "<<(int)(bin*1000.)<<" \n...break\n"<<std::endl;
    //   exit(98);
    // }
}

TString Process::convert(const Process& process)
{
    switch(process){
        case ttbarsignal: return "ttbarsignal";
        case ttbarsignal_genRho0to02: return "ttbarsignal_genRho0to02";
        case ttbarsignal_genRho02to03: return "ttbarsignal_genRho02to03";
        case ttbarsignal_genRho03to04: return "ttbarsignal_genRho03to04";
        case ttbarsignal_genRho04to05: return "ttbarsignal_genRho04to05";
        case ttbarsignal_genRho05to06: return "ttbarsignal_genRho05to06";
        case ttbarsignal_genRho06to1: return "ttbarsignal_genRho06to1";
        case ttbarsignal_genRho0to1: return "ttbarsignal_genRho0to1";
        case ttbarsignal_BKGNoAddJet: return "ttbarsignal_BKGNoAddJet";
        case ttbarbg: return "ttbarbg";
        case ttbarX: return "ttbarX";
        case dy: return "dy";
        case diboson: return "diboson";
        case singletop: return "singletop";
        case wjets: return "wjets";
        case data: return "data_obs";
        default:
          std::cerr<<"Error in Process::convert()! Conversion is not implemented\n...break\n"<<std::endl;
          exit(98);
          break;
    }
}

std::vector<Process::Process> Process::convert(const std::vector<TString>& processes)
{
    std::vector<Process> v_process;
    for(const auto& process : processes) v_process.push_back(convert(process));
    return v_process;
}

std::vector<Process::Process> Process::convert(const std::vector<std::string>& processes)
{
    std::vector<Process> v_process;
    for(const auto& process : processes) v_process.push_back(convert(process));
    return v_process;
}

std::vector<TString> Process::convert(const std::vector<Process>& processes)
{
    std::vector<TString> v_process;
    for(const auto& process : processes) v_process.push_back(convert(process));
    return v_process;
}



TString RenameSystematicString(const Systematic::Systematic& sys){
    TString out=sys.name();
    if(out.Contains("_UP")){
        out=out.ReplaceAll("_UP","Up");
    }
    if(out.Contains("_DOWN")){
        out=out.ReplaceAll("_DOWN","Down");
    }
    return out;
}


TString assignFolder(const char* baseDir, const Channel::Channel& channel, const TString& systematic, const char* subDir)
{
    TString path("");

    // Create all subdirectories contained in baseDir
    TObjArray* a_currentDir = TString(baseDir).Tokenize("/");
    for(Int_t iCurrentDir = 0; iCurrentDir < a_currentDir->GetEntriesFast(); ++iCurrentDir)
    {
        const TString& currentDir = a_currentDir->At(iCurrentDir)->GetName();
        path.Append(currentDir);
        path.Append("/");
        if(TString(baseDir).BeginsWith("/") && !path.BeginsWith("/")){ path.Prepend("/"); }
        gSystem->MakeDirectory(path);
    }

    // Create subdirectories for systematic and channel
    path.Append(systematic);
    path.Append("/");
    gSystem->MakeDirectory(path);
    path.Append(Channel::convert(channel));
    path.Append("/");
    gSystem->MakeDirectory(path);
    if(TString(subDir) != "")
    {
      path.Append(subDir);
      path.Append("/");
    }

    gSystem->MakeDirectory(path);

    return path;
}


const TString UserHisto::GetNameToKey(const UserHisto::HistoKey &key){
    TString tempString;
    if (key==nEvents)           {tempString="NEvents";}
    else if(key==genRho)        {tempString="genRho";}
    else if(key==genmTTbar)     {tempString="genTTbarM";}
    // else if(key==recoRhoKRS)    {tempString="kinRecoRho";}
    // else if(key==recoRhoLKRS)   {tempString="looseKinRecoRho";}
    else if(key==recoRhoDNN)   {tempString="recoRhoDNN";}
    // else if(key==mTTbarKRS)     {tempString="kinRecoTTbarM";}
    // else if(key==mTTbarLKRS)    {tempString="looseKinRecoTTbarM";}
    // else if(key==addJetPtKRS)   {tempString="kinRecoAddJetPt";}
    // else if(key==addJetPtLKRS)  {tempString="looseKinRecoAddJetPt";}
    else if(key==dileptonMass)  {tempString="dileptonMass";}
    else if(key==nJets)         {tempString="nJets";}
    else if(key==nBJets)        {tempString="nBJets";}
    else if(key==bJet1Pt)        {tempString="bJet1Pt";}
    else if(key==mlb_min)        {tempString="mlb_min";}

    // else if(key==genKRSRecoRhoDiff) {tempString="genKRSRecoRhoDiff";}
    // else if(key==genLKRSRecoRhoDiff) {tempString="genLKRSRecoRhoDiff";}
    else if(key==genRecoRhoDNNDiff) {tempString="genRecoRhoDNNDiff";}

    // MET
    else if(key==metPx)     {tempString="metPx";}
    else if(key==metPy)     {tempString="metPy";}
    else if(key==metPt)     {tempString="metPt";}
    else if(key==metM)      {tempString="metM";}
    else if(key==metPhi)    {tempString="metPhi";}
    else if(key==metE)      {tempString="metE";}
    // JETS
    else if(key==jet1Px)    {tempString="jet1Px";}
    else if(key==jet1Py)    {tempString="jet1Py";}
    else if(key==jet1Pz)    {tempString="jet1Pz";}
    else if(key==jet1M)     {tempString="jet1M";}
    else if(key==jet1Pt)    {tempString="jet1Pt";}
    else if(key==jet1Eta)   {tempString="jet1Eta";}
    else if(key==jet1Phi)   {tempString="jet1Phi";}
    else if(key==jet1E)     {tempString="jet1E";}

    else if(key==jet2Px)    {tempString="jet2Px";}
    else if(key==jet2Py)    {tempString="jet2Py";}
    else if(key==jet2Pz)    {tempString="jet2Pz";}
    else if(key==jet2M)     {tempString="jet2M";}
    else if(key==jet2Pt)    {tempString="jet2Pt";}
    else if(key==jet2Eta)   {tempString="jet2Eta";}
    else if(key==jet2Phi)   {tempString="jet2Phi";}
    else if(key==jet2E)     {tempString="jet2E";}

    else if(key==jet3Px)    {tempString="jet3Px";}
    else if(key==jet3Py)    {tempString="jet3Py";}
    else if(key==jet3Pz)    {tempString="jet3Pz";}
    else if(key==jet3M)     {tempString="jet3M";}
    else if(key==jet3Pt)    {tempString="jet3Pt";}
    else if(key==jet3Eta)   {tempString="jet3Eta";}
    else if(key==jet3Phi)   {tempString="jet3Phi";}
    else if(key==jet3E)     {tempString="jet3E";}
    else{
        std::cerr<<"Error in UserHisto::GetNameToKey()! Key unknown \n...break\n"<<std::endl;
    }
    return tempString;
}

const TString UserHisto::GetNameToKeys2D(const UserHisto::HistoKey &keyX, const UserHisto::HistoKey &keyY){
    TString tempStringX = GetNameToKey(keyX);
    TString tempStringY = GetNameToKey(keyY);

    return tempStringX+"_"+tempStringY;
}


float UserHisto::GetValueToKeyFromEvent(const Event &ev, const UserHisto::HistoKey &key){
    switch(key){
        case nEvents:       return (0.5);
        case genRho:        return (ev.gen_rho);
        // case genmTTbar:     return (ev.gen_ttbar.M());
        case genmTTbar:     return ((ev.gen_top+ev.gen_antitop).M());
        // case recoRhoDNN:    return (ev.DNN_rho);
        case recoRhoDNN:    return (ev.final_rho);
        case dileptonMass:  return (ev.dilepton.M());
        case nJets:         return (ev.jets.size());
        case nBJets:        return (ev.bjets.size());
        case bJet1Pt:       return (ev.bjets.at(0).Pt());
        case mlb_min:       return (ev.mlb_min);

        // case genRecoRhoDNNDiff:  return (ev.gen_rho-ev.DNN_rho);
        case genRecoRhoDNNDiff:  return (ev.gen_rho-ev.final_rho);
        // MET
        case metPx:     return (ev.met.Px());
        case metPy:     return (ev.met.Py());
        case metPt:     return (ev.met.Pt());
        case metM:      return (ev.met.M());
        case metPhi:    return (ev.met.Phi());
        case metE:      return (ev.met.E());
        // JETS
        case jet1Px:    return (ev.jets.at(0).Px());
        case jet1Py:    return (ev.jets.at(0).Py());
        case jet1Pz:    return (ev.jets.at(0).Pz());
        case jet1M:     return (ev.jets.at(0).M());
        case jet1Pt:    return (ev.jets.at(0).Pt());
        case jet1Eta:   return (ev.jets.at(0).Eta());
        case jet1Phi:   return (ev.jets.at(0).Phi());
        case jet1E:     return (ev.jets.at(0).E());

        case jet2Px:    return (ev.jets.at(1).Px());
        case jet2Py:    return (ev.jets.at(1).Py());
        case jet2Pz:    return (ev.jets.at(1).Pz());
        case jet2M:     return (ev.jets.at(1).M());
        case jet2Pt:    return (ev.jets.at(1).Pt());
        case jet2Eta:   return (ev.jets.at(1).Eta());
        case jet2Phi:   return (ev.jets.at(1).Phi());
        case jet2E:     return (ev.jets.at(1).E());

        case jet3Px:    return (ev.jets.at(2).Px());
        case jet3Py:    return (ev.jets.at(2).Py());
        case jet3Pz:    return (ev.jets.at(2).Pz());
        case jet3M:     return (ev.jets.at(2).M());
        case jet3Pt:    return (ev.jets.at(2).Pt());
        case jet3Eta:   return (ev.jets.at(2).Eta());
        case jet3Phi:   return (ev.jets.at(2).Phi());
        case jet3E:     return (ev.jets.at(2).E());
        default:
            std::cerr<<"Error in UserHisto::GetValueToKeyFromEvent()! Key unknown \n...break\n"<<std::endl;
            return -999.;
    }
}

bool UserHisto::CheckKeyConditionFromEvent(const Event &ev, const UserHisto::HistoKey &key){
    switch (key){
        case nEvents:       return true;
        case genRho:        return true;
        case genmTTbar:     return true;
        case recoRhoDNN:    return true;
        case dileptonMass:  return true;
        case nJets:         return true;
        case nBJets:        return true;
        case bJet1Pt:       return (ev.bjets.size()>0);
        case mlb_min:       return (ev.bjets.size()>0);

        case genRecoRhoDNNDiff: return (true);
    // MET
        case metPx:     return true;
        case metPy:     return true;
        case metPt:     return true;
        case metM:      return true;
        case metPhi:    return true;
        case metE:      return true;
    // JETS
        case jet1Px:    return (ev.jets.size()>0);
        case jet1Py:    return (ev.jets.size()>0);
        case jet1Pz:    return (ev.jets.size()>0);
        case jet1M:     return (ev.jets.size()>0);
        case jet1Pt:    return (ev.jets.size()>0);
        case jet1Eta:   return (ev.jets.size()>0);
        case jet1Phi:   return (ev.jets.size()>0);
        case jet1E:     return (ev.jets.size()>0);

        case jet2Px:    return (ev.jets.size()>1);
        case jet2Py:    return (ev.jets.size()>1);
        case jet2Pz:    return (ev.jets.size()>1);
        case jet2M:     return (ev.jets.size()>1);
        case jet2Pt:    return (ev.jets.size()>1);
        case jet2Eta:   return (ev.jets.size()>1);
        case jet2Phi:   return (ev.jets.size()>1);
        case jet2E:     return (ev.jets.size()>1);

        case jet3Px:    return (ev.jets.size()>2);
        case jet3Py:    return (ev.jets.size()>2);
        case jet3Pz:    return (ev.jets.size()>2);
        case jet3M:     return (ev.jets.size()>2);
        case jet3Pt:    return (ev.jets.size()>2);
        case jet3Eta:   return (ev.jets.size()>2);
        case jet3Phi:   return (ev.jets.size()>2);
        case jet3E:     return (ev.jets.size()>2);
        default:
            std::cout<<GetNameToKey(key)<<std::endl;
            std::cerr<<"Error in UserHisto::CheckKeyConditionFromEvent()! Key unknown \n...break\n"<<std::endl;
            return -999.;
    }
}
