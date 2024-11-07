#ifndef TopAnalysis_h
#define TopAnalysis_h

#include <map>

class TTree;
class TH1;
class TH2;

#include "../../common/include/AnalysisBase.h"
#include "../../common/include/classesFwd.h"
#include "analysisStructsFwd.h"
#include "utils.h"

#include "../../common/include/analysisObjectStructs.h"

class AnalysisConfig;
class AnalyzerBaseClass;
class TreeHandlerBase;
class EventMetadata;
class RecoObjects;
class CommonGenObjects;
class TopGenObjects;
class TopPseudoObjects;
class KinematicReconstructionSolutions;
class LooseKinRecoSolution;
namespace ttbar{
    class GenLevelWeights;
    class RecoLevelWeights;
    class GenObjectIndices;
    class RecoObjectIndices;
}




class TopAnalysis : public AnalysisBase
{

    /// Histograms

    // Pile-Up Reweighting (Step 3 and Step 8)

    TH1 *h_vertMulti_noPU_step3, *h_vertMulti_step3, *h_vertMulti_noPU_step8, *h_vertMulti_step8;
    
    // Full Dilepton System / pre-Zpeak cut (Step 3)

    TH1 *h_diLepMassFull;

    // Dilepton System / pre-jetMult cut (Step 4)

    TH1 *h_jetMulti_diLep, *h_jetMulti_Zwindow;
    TH1 *h_MET_diLep;
    TH1 *h_MuonpT_diLep, *h_MuonEta_diLep;
    TH1 *h_ElectronpT_diLep, *h_ElectronEta_diLep;
    TH1 *h_LeptonpT_diLep, *h_LeptonEta_diLep;
    TH1 *h_AntiLeptonpT_diLep, *h_AntiLeptonEta_diLep;

    // pre-MET cut (Step 5)
    TH1 *h_MET_preMETcut, *h_MET_Zwindow;

    // pre-btag cut (Step 6)
    TH1 *h_MuonpT_postMETcut, *h_MuonEta_postMETcut;
    TH1 *h_ElectronpT_postMETcut, *h_ElectronEta_postMETcut;
    TH1 *h_LeptonpT_postMETcut, *h_LeptonEta_postMETcut;
    TH1 *h_AntiLeptonpT_postMETcut, *h_AntiLeptonEta_postMETcut;
    TH1 *h_jetMulti_noBTag, *h_BjetMulti_noBTag, *h_BjetMulti_Zwindow;

    // pre-KinReco (Step 7)
    TH1 *h_jetMulti, *h_BjetMulti;
    TH1 *h_jetpT,*h_jetHT;
    TH1 *h_LeptonpT_BeforeKinReco, *h_BjetpT_BeforeKinReco;
    TH1 *h_MET_step7;
    TH1 *h_LeptonPt_BeforeKinReco, *h_LeptonEta_BeforeKinReco;
    TH1 *h_MuonpT_BeforeKinReco, *h_MuonEta_BeforeKinReco;
    TH1 *h_ElectronpT_BeforeKinReco, *h_ElectronEta_BeforeKinReco;
    TH1 *h_BjetEta_BeforeKinReco;

    // post-KinReco (Step 8)

    TH1 *h_diLepMassFull_fullSel;

    TH1 *h_MET_step8;

    TH1 *h_leptonPt_AfterKinReco, *h_leptonEta_AfterKinReco;
    TH1 *h_bjeteta_AfterKinReco;

    TH1 *h_LeptonpT, *h_LeptonEta;
    TH1 *h_AntiLeptonpT, *h_AntiLeptonEta;

    TH1 *h_jetMulti_kinReco, *h_bjetMulti_kinReco, *h_exjetMulti_kinReco;

    TH1 *h_HypAntiBJetpT, *h_HypAntiBJetEta, *h_HypAntiBJetRapidity,*h_HypAntiBJetCSVdiscriminator;
    TH1 *h_HypBJetpT, *h_HypBJetEta, *h_HypBJetRapidity, *h_HypBJetCSVdiscriminator;

    TH1 *h_HypAntiToppT, *h_HypAntiTopEta, *h_HypAntiTopMass,*h_HypAntiTopRapidity, *h_HypAntiTopRapidityAbs;
    TH1 *h_HypToppT, *h_HypTopEta,*h_HypTopMass, *h_HypTopRapidity, *h_HypTopRapidityAbs;

    TH1 *h_HypTTBarMass, *h_HypTTBarRapidity, *h_HypTTBarRapidityAbs, *h_HypTTBarpT;
    TH1 *h_HypTTBarMass_Loose, *h_HypTTBarRapidity_Loose, *h_HypTTBarRapidityAbs_Loose, *h_HypTTBarpT_Loose;
    TH1 *h_HypLLBarMass, *h_HypLLBarpT;
    TH1 *h_HypMet;

    TH1 *h_HypAntiLeptonpT, *h_HypAntiLeptonEta;
    TH1 *h_HypLeptonpT, *h_HypLeptonEta;

    TH1 *h_HypToppTLead,    *h_HypToppTNLead,    *h_HypTopRapidityLead, *h_HypTopRapidityNLead, *h_HypTopMassLead, *h_HypTopMassNLead;
    TH1 *h_HypLeptonpTLead, *h_HypLeptonpTNLead, *h_HypLeptonEtaLead,   *h_HypLeptonEtaNLead;
    TH1 *h_HypBJetpTLead,   *h_HypBJetpTNLead,   *h_HypBJetEtaLead,     *h_HypBJetEtaNLead;

    //    TH1 *h_exjetMulti_0, *h_exjetMulti_1, *h_exjetMulti_2;

    TH2 *h_GenRecoLeptonpT,*h_GenRecoAntiLeptonpT,*h_GenRecoLeptonEta,*h_GenRecoAntiLeptonEta, *h_GenRecoLLBarMass, *h_GenRecoLLBarpT;
    TH2 *h_GenRecoBJetpT,*h_GenRecoAntiBJetpT, *h_GenRecoBJetEta,*h_GenRecoAntiBJetEta, *h_GenRecoBJetRapidity, *h_GenRecoAntiBJetRapidity;
    TH2 *h_GenRecoToppT,*h_GenRecoAntiToppT,*h_GenRecoTopRapidity,*h_GenRecoAntiTopRapidity, *h_GenRecoTTBarMass, *h_GenRecoTTBarpT, *h_GenRecoTTBarRapidity;
    TH2 *h_GenRecoTopRapidityAbs,*h_GenRecoAntiTopRapidityAbs, *h_GenRecoTTBarRapidityAbs;
    TH2 *h_GenRecoMet;

    TH1 *h_NJetMatching;

    TH1 *h_jetMultiXSec,*h_jetMultiXSec_Loose,*h_jetMultiAll, *h_jetMultiNoPU, *h_jetMultiVisTop;
    TH1 *h_exjetMulti;
    TH1 *h_GenAll_RecoCuts, *h_GenAll_RecoCuts_noweight, *h_GenAll, *h_GenAll_noweight, *h_VisGenAll, *h_VisGenAll_noweight;

    TH1 *h_VisGenTTBarMass, *h_VisGenTTBarRapidity, *h_VisGenTTBarRapidityAbs, *h_VisGenTTBarpT;
    TH1 *h_VisGenTopRapidity,*h_VisGenAntiTopRapidity, *h_VisGenTopRapidityAbs, *h_VisGenAntiTopRapidityAbs;
    TH1 *h_VisGenLLBarMass,*h_VisGenLLBarpT;
    TH1 *h_VisGenMet;

    TH1 *h_RecoTTBarMass, *h_RecoTTBarRapidity, *h_RecoTTBarRapidityAbs, *h_RecoTTBarpT;
    TH1 *h_RecoToppT,*h_RecoAntiToppT,*h_RecoTopRapidity,*h_RecoAntiTopRapidity, *h_RecoTopRapidityAbs,*h_RecoAntiTopRapidityAbs;
    TH1 *h_RecoLLBarMass, *h_RecoLLBarpT;
    TH1 *h_RecoLeptonpT,*h_RecoAntiLeptonpT,*h_RecoLeptonEta,*h_RecoAntiLeptonEta;
    TH1 *h_RecoBJetpT,*h_RecoAntiBJetpT, *h_RecoBJetRapidity,*h_RecoAntiBJetRapidity,*h_RecoBJetEta,*h_RecoAntiBJetEta;
    TH1 *h_RecoMet;

    TH2 *h_GenRecoHT;
    TH1 *h_VisGenHT, *h_HypHT, *h_RecoHT;
    
    TH1 *h_HypNeutrinopT, *h_HypAntiNeutrinopT;
    TH1 *h_RecoNeutrinopT, *h_RecoAntiNeutrinopT;
    TH1 *h_VisGenNeutrinopT, *h_VisGenAntiNeutrinopT;
    TH2 *h_GenRecoNeutrinopT, *h_GenRecoAntiNeutrinopT;

    TH1 *h_VisGenAntiToppT, *h_VisGenAntiTopEta;
    TH1 *h_VisGenToppT, *h_VisGenTopEta;

    TH1 *h_VisGenAntiBJetpT, *h_VisGenAntiBJetEta, *h_VisGenAntiBJetRapidity;
    TH1 *h_VisGenBJetpT, *h_VisGenBJetEta, *h_VisGenBJetRapidity;

    TH1 *h_VisGenAntiBQuarkpT, *h_VisGenAntiBQuarkEta, *h_VisGenAntiBQuarkRapidity;
    TH1 *h_VisGenBQuarkpT, *h_VisGenBQuarkEta, *h_VisGenBQuarkRapidity;

    TH1 *h_VisGenAntiLeptonpT, *h_VisGenAntiLeptonEta;
    TH1 *h_VisGenLeptonpT, *h_VisGenLeptonEta;

    TH2 *h_GenRecoTTBarDeltaPhi, *h_GenRecoTTBarDeltaRapidity;
    TH1 *h_RecoTTBarDeltaPhi, *h_RecoTTBarDeltaRapidity;
    TH1 *h_HypTTBarDeltaPhi, *h_HypTTBarDeltaRapidity;
    TH1 *h_VisGenTTBarDeltaPhi, *h_VisGenTTBarDeltaRapidity;

    TH2 *h_GenRecoBBBarpT, *h_GenRecoBBBarMass, *h_GenRecoBBBarDPhi;
    TH1 *h_RecoBBBarpT, *h_RecoBBBarMass, *h_RecoBBBarDPhi;
    TH1 *h_HypBBBarpT, *h_HypBBBarMass, *h_HypBBBarDPhi;
    TH1 *h_VisGenBBBarpT, *h_VisGenBBBarMass, *h_VisGenBBBarDPhi;

    TH1 *h_HypToppTTTRestFrame, *h_HypAntiToppTTTRestFrame;
    TH1 *h_RecoToppTTTRestFrame, *h_RecoAntiToppTTTRestFrame;
    TH1 *h_VisGenToppTTTRestFrame, *h_VisGenAntiToppTTTRestFrame;
    TH2 *h_GenRecoToppTTTRestFrame, *h_GenRecoAntiToppTTTRestFrame;

    TH2 *h_GenRecoLLBarDPhi, *h_GenRecoLeptonantiBjetMass, *h_GenRecoAntiLeptonBjetMass, *h_GenRecoJetMult, *h_GenRecoLLBarDEta;
    TH1 *h_VisGenLLBarDPhi,  *h_VisGenLeptonantiBjetMass,  *h_VisGenAntiLeptonBjetMass,  *h_VisGenJetMult, *h_VisGenLLBarDEta;
    TH1 *h_HypLLBarDPhi,     *h_HypLeptonantiBjetMass,     *h_HypAntiLeptonBjetMass,     *h_HypJetMult, *h_HypLLBarDEta;
    TH1 *h_RecoLLBarDPhi,    *h_RecoLeptonantiBjetMass,    *h_RecoAntiLeptonBjetMass,    *h_RecoJetMult, *h_RecoLLBarDEta;

    TH1 *h_RecoToppTLead,    *h_RecoToppTNLead,    *h_RecoTopRapidityLead, *h_RecoTopRapidityNLead, *h_RecoTopMassLead, *h_RecoTopMassNLead;
    TH1 *h_RecoLeptonpTLead, *h_RecoLeptonpTNLead, *h_RecoLeptonEtaLead,   *h_RecoLeptonEtaNLead;
    TH1 *h_RecoBJetpTLead,   *h_RecoBJetpTNLead,   *h_RecoBJetEtaLead,     *h_RecoBJetEtaNLead;

    TH1 *h_VisGenToppTLead,    *h_VisGenToppTNLead,    *h_VisGenTopRapidityLead, *h_VisGenTopRapidityNLead, *h_VisGenTopMassLead, *h_VisGenTopMassNLead;
    TH1 *h_VisGenLeptonpTLead, *h_VisGenLeptonpTNLead, *h_VisGenLeptonEtaLead,   *h_VisGenLeptonEtaNLead;
    TH1 *h_VisGenBJetpTLead,   *h_VisGenBJetpTNLead,   *h_VisGenBJetEtaLead,     *h_VisGenBJetEtaNLead;

    TH2 *h_GenRecoToppTLead,    *h_GenRecoToppTNLead,    *h_GenRecoTopRapidityLead, *h_GenRecoTopRapidityNLead, *h_GenRecoTopMassLead, *h_GenRecoTopMassNLead;
    TH2 *h_GenRecoLeptonpTLead, *h_GenRecoLeptonpTNLead, *h_GenRecoLeptonEtaLead,   *h_GenRecoLeptonEtaNLead;
    TH2 *h_GenRecoBJetpTLead,   *h_GenRecoBJetpTNLead,   *h_GenRecoBJetEtaLead,     *h_GenRecoBJetEtaNLead;

    TH1 *h_RecoAbsDeltaPhiExtraJet12,*h_HypAbsDeltaPhiExtraJet12;
    TH1 *h_RecoAbsDeltaPhiExtraJet12_eta1,*h_HypAbsDeltaPhiExtraJet12_eta1;
    TH1 *h_RecoAbsDeltaPhiExtraJet12_eta2,*h_HypAbsDeltaPhiExtraJet12_eta2;
    TH1 *h_RecoAbsDeltaPhiExtraJet12_eta3,*h_HypAbsDeltaPhiExtraJet12_eta3;
    TH1 *h_RecoAbsDeltaPhiExtraJet12_eta4,*h_HypAbsDeltaPhiExtraJet12_eta4;
    TH1 *h_RecoAbsDeltaEtaExtraJet12,*h_HypAbsDeltaEtaExtraJet12;

    //Begin: Plots for Carmen
    TH1 *h_RecoLeadingJetpT,    *h_RecoNLeadingJetpT,    *h_RecoLeadingJetEta,    *h_RecoNLeadingJetEta;
    TH1 *h_HypLeadingJetpT,     *h_HypNLeadingJetpT,     *h_HypLeadingJetEta,     *h_HypNLeadingJetEta;
    TH2 *h_GenRecoLeadingJetpT, *h_GenRecoLeadingJetEta, *h_GenRecoNLeadingJetpT, *h_GenRecoNLeadingJetEta;
    TH1 *h_VisGenLeadingJetpT,  *h_VisGenLeadingJetEta,  *h_VisGenNLeadingJetpT,  *h_VisGenNLeadingJetEta;

    TH1 *h_RecoExtraJetpT,  *h_HypExtraJetpT, *h_VisGenExtraJetpT, *h_RecoExtraJetEta, *h_HypExtraJetEta, *h_VisGenExtraJetEta, *h_RecoExtraJetAbsEta, *h_HypExtraJetAbsEta, *h_VisGenExtraJetAbsEta;
    TH1 *h_RecoExtraJetpT2, *h_HypExtraJetpT2, *h_VisGenExtraJetpT2, *h_RecoExtraJetEta2, *h_HypExtraJetEta2, *h_VisGenExtraJetEta2, *h_RecoExtraJetAbsEta2, *h_HypExtraJetAbsEta2, *h_VisGenExtraJetAbsEta2;
    TH1 *h_RecoExtraJetpT3, *h_HypExtraJetpT3, *h_VisGenExtraJetpT3, *h_RecoExtraJetEta3, *h_HypExtraJetEta3, *h_VisGenExtraJetEta3, *h_RecoExtraJetAbsEta3, *h_HypExtraJetAbsEta3, *h_VisGenExtraJetAbsEta3;
    TH1 *h_RecoExtraJetpT4, *h_HypExtraJetpT4, *h_VisGenExtraJetpT4, *h_RecoExtraJetEta4, *h_HypExtraJetEta4, *h_VisGenExtraJetEta4, *h_RecoExtraJetAbsEta4, *h_HypExtraJetAbsEta4, *h_VisGenExtraJetAbsEta4;
    TH2 *h_GenRecoExtraJetpT, *h_GenRecoExtraJetEta, *h_GenRecoExtraJetpT2, *h_GenRecoExtraJetEta2, *h_GenRecoExtraJetpT3, *h_GenRecoExtraJetEta3, *h_GenRecoExtraJetpT4, *h_GenRecoExtraJetEta4, *h_GenRecoExtraJetAbsEta, *h_GenRecoExtraJetAbsEta2, *h_GenRecoExtraJetAbsEta3, *h_GenRecoExtraJetAbsEta4;

    TH1 *h_RecoExtraJetHT, *h_HypExtraJetHT, *h_VisGenExtraJetHT;
    TH2 *h_GenRecoExtraJetHT;

    TH1 *h_RecoJetMultpt30, *h_RecoJetMultpt40, *h_HypJetMultpt30, *h_VisGenJetMultpt30, *h_HypJetMultpt40, *h_VisGenJetMultpt40, *h_RecoJetMultpt60, *h_HypJetMultpt60, *h_VisGenJetMultpt60;
    TH1 *h_RecoJetMultpt100, *h_HypJetMultpt100, *h_VisGenJetMultpt100;
    TH2 *h_GenRecoJetMultpt30, *h_GenRecoJetMultpt40, *h_GenRecoJetMultpt60, *h_GenRecoJetMultpt100;

    TH1 *h_HypJetMultQ0, *h_RecoJetMultQ0, *h_VisGenJetMultQ0;
    TH1 *h_RecoJetMultTotal, *h_HypJetMultTotal, *h_VisGenJetMultTotal;
    TH1 *h_HypJetMultQsum, *h_RecoJetMultQsum, *h_VisGenJetMultQsum;
    TH1 *h_HypJetExtra2Q0, *h_RecoJetExtra2Q0, *h_VisGenJetExtra2Q0;
    TH2 *h_GenRecoJetMultQ0, *h_GenRecoJetExtra2Q0, *h_GenRecoJetMultQsum, *h_GenRecoJetMultTotal;

    TH2 *h_GenRecoDeltaRExtraJet12;
    TH1 *h_VisGenDeltaRExtraJet12, *h_RecoDeltaRExtraJet12, *h_HypDeltaRExtraJet12;
    TH2 *h_GenRecoDeltaPhiExtraJet12, *h_GenRecoPhiExtraJet12,*h_GenRecoTTBar1stJetMass, *h_GenRecoTTBar0Mass, *h_GenRecoMassExtraJet12;
    TH1 *h_VisGenDeltaPhiExtraJet12, *h_RecoDeltaPhiExtraJet12, *h_HypDeltaPhiExtraJet12, *h_VisGenPhiExtraJet12, *h_RecoPhiExtraJet12, *h_HypPhiExtraJet12;
    TH1 *h_VisGenMassExtraJet12, *h_RecoMassExtraJet12, *h_HypMassExtraJet12;

    TH1 *h_VisGenTTBar1stJetMass, *h_RecoTTBar1stJetMass, *h_HypTTBar1stJetMass;
    TH1 *h_VisGenTTBar0Mass, *h_RecoTTBar0Mass, *h_HypTTBar0Mass;

    TH1 *h_HypLLBarCosPhiInTopRestFrame, *h_VisGenLLBarCosPhiInTopRestFrame, *h_RecoLLBarCosPhiInTopRestFrame;
    TH2 *h_GenRecoLLBarCosPhiInTopRestFrame;

    //End: Plots for Carmen

    //=======================PSEUDO-TOP: BEGIN==========================//

    TH1 *h_VisPseudoAll, *h_VisPseudoAll_noweight, *h_VisPseudoNuNuBarpT;

    TH1 *h_VisPseudoJetMultpt30, *h_VisPseudoJetMultpt40, *h_VisPseudoJetMultpt60, *h_VisPseudoJetMultpt100;
    TH1 *h_VisPseudoLeptonpT, *h_VisPseudoAntiLeptonpT, *h_VisPseudoLeptonEta, *h_VisPseudoAntiLeptonEta;
    TH1 *h_VisPseudoLeptonpTLead, *h_VisPseudoLeptonpTNLead, *h_VisPseudoLeptonEtaLead, *h_VisPseudoLeptonEtaNLead;
    TH1 *h_VisPseudoBJetEta, *h_VisPseudoAntiBJetEta, *h_VisPseudoBJetRapidity, *h_VisPseudoAntiBJetRapidity, *h_VisPseudoBJetpT, *h_VisPseudoAntiBJetpT;
    TH1 *h_VisPseudoBJetpTLead, *h_VisPseudoBJetEtaLead, *h_VisPseudoBJetpTNLead, *h_VisPseudoBJetEtaNLead;
    TH1 *h_VisPseudoLLBarpT, *h_VisPseudoLLBarMass, *h_VisPseudoLLBarDPhi, *h_VisPseudoLLBarDEta, *h_VisPseudoLLBarCosPhiInTopRestFrame;
    TH1 *h_VisPseudoBBBarpT, *h_VisPseudoBBBarMass, *h_VisPseudoBBBarDPhi;

    TH1 *h_VisPseudoTTBarMass, *h_VisPseudoTTBarRapidity, *h_VisPseudoTTBarRapidityAbs, *h_VisPseudoTTBarpT, *h_VisPseudoTTBarDeltaPhi, *h_VisPseudoTTBarDeltaRapidity;
    TH1 *h_VisPseudoToppT, *h_VisPseudoAntiToppT, *h_VisPseudoTopRapidity, *h_VisPseudoAntiTopRapidity, *h_VisPseudoTopRapidityAbs, *h_VisPseudoAntiTopRapidityAbs;
    TH1 *h_VisPseudoToppTLead, *h_VisPseudoToppTNLead,    *h_VisPseudoTopRapidityLead, *h_VisPseudoTopRapidityNLead;

    TH2 *h_PseudoRecoJetMultpt30, *h_PseudoRecoJetMultpt40, *h_PseudoRecoJetMultpt60, *h_PseudoRecoJetMultpt100;
    TH2 *h_PseudoRecoLeptonpT, *h_PseudoRecoAntiLeptonpT, *h_PseudoRecoLeptonEta, *h_PseudoRecoAntiLeptonEta;
    TH2 *h_PseudoRecoLeptonpTLead, *h_PseudoRecoLeptonpTNLead, *h_PseudoRecoLeptonEtaLead, *h_PseudoRecoLeptonEtaNLead;
    TH2 *h_PseudoRecoBJetEta, *h_PseudoRecoAntiBJetEta, *h_PseudoRecoBJetRapidity, *h_PseudoRecoAntiBJetRapidity, *h_PseudoRecoBJetpT, *h_PseudoRecoAntiBJetpT;
    TH2 *h_PseudoRecoBJetpTLead, *h_PseudoRecoBJetEtaLead, *h_PseudoRecoBJetpTNLead, *h_PseudoRecoBJetEtaNLead;
    TH2 *h_PseudoRecoLLBarpT, *h_PseudoRecoLLBarMass, *h_PseudoRecoLLBarDPhi, *h_PseudoRecoLLBarDEta, *h_PseudoRecoLLBarCosPhiInTopRestFrame;
    TH2 *h_PseudoRecoBBBarpT, *h_PseudoRecoBBBarMass, *h_PseudoRecoBBBarDPhi;

    TH2 *h_PseudoRecoTTBarMass, *h_PseudoRecoTTBarRapidity, *h_PseudoRecoTTBarRapidityAbs, *h_PseudoRecoTTBarpT, *h_PseudoRecoTTBarDeltaPhi, *h_PseudoRecoTTBarDeltaRapidity;
    TH2 *h_PseudoRecoToppT, *h_PseudoRecoAntiToppT, *h_PseudoRecoTopRapidity, *h_PseudoRecoAntiTopRapidity, *h_PseudoRecoTopRapidityAbs, *h_PseudoRecoAntiTopRapidityAbs;
    TH2 *h_PseudoRecoToppTLead, *h_PseudoRecoToppTNLead, *h_PseudoRecoTopRapidityLead, *h_PseudoRecoTopRapidityNLead;

    TH2 *h_PseudoRecoLeptonantiBjetMass, *h_PseudoRecoAntiLeptonBjetMass;
    TH1 *h_VisPseudoLeptonantiBjetMass,  *h_VisPseudoAntiLeptonBjetMass;

    TH2 *h_PseudoRecoToppTTTRestFrame, *h_PseudoRecoAntiToppTTTRestFrame;
    TH1 *h_VisPseudoToppTTTRestFrame, *h_VisPseudoAntiToppTTTRestFrame;

    TH2 *h_GenPseudoToppT, *h_GenPseudoAntiToppT, *h_GenPseudoTopRapidity, *h_GenPseudoAntiTopRapidity, *h_GenPseudoTTBarRapidity, *h_GenPseudoTTBarpT, *h_GenPseudoTTBarMass;
    TH2 *h_GenPseudoLeptonpT, *h_GenPseudoAntiLeptonpT, *h_GenPseudoLeptonEta, *h_GenPseudoAntiLeptonEta;
    TH2 *h_GenPseudoBJetpT, *h_GenPseudoAntiBJetpT, *h_GenPseudoBJetEta, *h_GenPseudoAntiBJetEta, *h_GenPseudoJetMultpt30;

    TH1 *h_GenPseudoToppT_reso, *h_GenPseudoAntiToppT_reso, *h_GenPseudoTopRapidity_reso, *h_GenPseudoAntiTopRapidity_reso, *h_GenPseudoTTBarRapidity_reso, *h_GenPseudoTTBarpT_reso, *h_GenPseudoTTBarMass_reso;
    TH1 *h_GenPseudoLeptonpT_reso, *h_GenPseudoAntiLeptonpT_reso, *h_GenPseudoLeptonEta_reso, *h_GenPseudoAntiLeptonEta_reso;
    TH1 *h_GenPseudoBJetpT_reso, *h_GenPseudoAntiBJetpT_reso, *h_GenPseudoBJetEta_reso, *h_GenPseudoAntiBJetEta_reso, *h_GenPseudoJetMultpt30_reso;

    TH1 *h_HypToppT_OutOfSpace, *h_HypAntiToppT_OutOfSpace, *h_HypTopRapidity_OutOfSpace, *h_HypAntiTopRapidity_OutOfSpace, *h_HypTopRapidityAbs_OutOfSpace, *h_HypAntiTopRapidityAbs_OutOfSpace;
    TH1 *h_HypToppTLead_OutOfSpace, *h_HypToppTNLead_OutOfSpace, *h_HypTopRapidityLead_OutOfSpace, *h_HypTopRapidityNLead_OutOfSpace, *h_HypToppTTTRestFrame_OutOfSpace, *h_HypAntiToppTTTRestFrame_OutOfSpace;
    TH1 *h_HypTTBarMass_OutOfSpace, *h_HypTTBarpT_OutOfSpace, *h_HypTTBarRapidity_OutOfSpace, *h_HypTTBarRapidityAbs_OutOfSpace, *h_HypTTBarDeltaPhi_OutOfSpace, *h_HypTTBarDeltaRapidity_OutOfSpace;
    TH1 *h_HypJetMultpt30_OutOfSpace, *h_HypJetMultpt40_OutOfSpace, *h_HypJetMultpt60_OutOfSpace, *h_HypJetMultpt100_OutOfSpace;
    TH1 *h_HypLeptonpT_OutOfSpace, *h_HypAntiLeptonpT_OutOfSpace, *h_HypLeptonEta_OutOfSpace, *h_HypAntiLeptonEta_OutOfSpace;
    TH1 *h_HypLeptonpTLead_OutOfSpace, *h_HypLeptonpTNLead_OutOfSpace, *h_HypLeptonEtaLead_OutOfSpace, *h_HypLeptonEtaNLead_OutOfSpace;
    TH1 *h_HypLLBarMass_OutOfSpace, *h_HypLLBarpT_OutOfSpace, *h_HypLLBarDPhi_OutOfSpace, *h_HypLLBarDEta_OutOfSpace, *h_HypLLBarCosPhiInTopRestFrame_OutOfSpace;
    TH1 *h_HypBBBarpT_OutOfSpace, *h_HypBBBarMass_OutOfSpace, *h_HypBBBarDPhi_OutOfSpace;
    TH1 *h_HypBJetpT_OutOfSpace, *h_HypAntiBJetpT_OutOfSpace, *h_HypBJetEta_OutOfSpace, *h_HypAntiBJetEta_OutOfSpace, *h_HypBJetRapidity_OutOfSpace, *h_HypAntiBJetRapidity_OutOfSpace;
    TH1 *h_HypBJetpTLead_OutOfSpace, *h_HypBJetpTNLead_OutOfSpace, *h_HypBJetEtaLead_OutOfSpace, *h_HypBJetEtaNLead_OutOfSpace, *h_HypLeptonantiBjetMass_OutOfSpace, *h_HypAntiLeptonBjetMass_OutOfSpace;
    //=======================PSEUDO-TOP: END==========================//

    /// Plots for the parton momentum fraction defined by Olaf
    TH1 *h_HypTopPartonFraction, *h_HypAntiTopPartonFraction;
    TH1 *h_VisGenTopPartonFraction, *h_VisGenAntiTopPartonFraction;
    TH1 *h_RecoTopPartonFraction, *h_RecoAntiTopPartonFraction;
    TH2 *h_GenRecoTopPartonFraction, *h_GenRecoAntiTopPartonFraction;

    /// Histograms for event weights due to specific scale factor
    TH1 *h_PUSF, *h_TrigSF, *h_LepSF, *h_BTagSF, *h_KinRecoSF, *h_EventWeight, *h_PrefiringWeight, *h_ToppTWeight, *h_weights;


    /// Pointer to JetVetoMap instance
    const utils::JetVetoMap* jetVetoMap_;

    /// Histograms for resolution studies

    //TH1 *h_RMSvsGenToppT;
    //TH1 *h_RMSvsGenTopRapidity;
    //TH1 *h_RMSvsGenToppTLead;
    //TH1 *h_RMSvsGenTopRapidityLead;
    //TH1 *h_RMSvsGenToppTNLead;
    //TH1 *h_RMSvsGenTopRapidityNLead;

    //TH1 *h_RMSvsGenTTBarMass;
    //TH1 *h_RMSvsGenTTBarRapidity;
    //TH1 *h_RMSvsGenTTBarpT;
    //TH1 *h_RMSvsGenTTBarDPhi;
    //TH1 *h_RMSvsGenTTBarDeltaRapidity;

    TH2 *h_HypTTBarRapidityvsTTBarpT;

    TH2 *h_VisGenTTBarRapidityvsTTBarpT;
    /// ...

public:

    /// Constructor
    TopAnalysis(const AnalysisConfig& analysisConfig, TTree* =0);

    /// Inherited from AnalysisBase and overwritten for needs of TopAnalysis
    virtual void Begin(TTree*);
    virtual void SlaveBegin(TTree*);
    virtual void SlaveTerminate();
    virtual void Terminate();
    virtual Bool_t Process(Long64_t entry);

    /// Get a constant reference to nTuple branches for Top signal samples on generator level
    // Define own functionality and do not use AnalysisBase method
    const TopGenObjects& getTopGenObjects(const Long64_t& entry)const;

    /// Bool for separating direct dileptonic ttbar decays and decays via intermediate taus
    virtual void SetRunViaTau(const bool runViaTau);

    /// Set up closure test
    void SetClosureTest(TString closure, double slope);

    /// Class definition
    ClassDef(TopAnalysis, 0);

    /// Set up all analysers of type AnalyzerBaseClass
    void SetAllAnalyzers(std::vector<AnalyzerBaseClass*> v_analyzer);

    /// Set up all tree handlers of type TreeHandlerBase
    virtual void SetAllTreeHandlers(std::vector<TreeHandlerBase*> v_treeHandler);


    // ================================================================
    // virtual functions needed for PlainTree Analysis
    // void casts needed to suppress unused variable warnings
    // virtual void SetLeptonScaleFactors_Up(const LeptonScaleFactors* const scaleFactors){(void)scaleFactors;};
    // virtual void SetLeptonScaleFactors_Dn(const LeptonScaleFactors* const scaleFactors){(void)scaleFactors;};
    // virtual void SetLeptonScaleFactors_EleUp(const LeptonScaleFactors* const scaleFactors){(void)scaleFactors;};
    // virtual void SetLeptonScaleFactors_EleDn(const LeptonScaleFactors* const scaleFactors){(void)scaleFactors;};
    // virtual void SetLeptonScaleFactors_MuonUp(const LeptonScaleFactors* const scaleFactors){(void)scaleFactors;};
    // virtual void SetLeptonScaleFactors_MuonDn(const LeptonScaleFactors* const scaleFactors){(void)scaleFactors;};
    virtual void SetTriggerScaleFactors_Up(const TriggerScaleFactors* const scaleFactors){(void)scaleFactors;};
    virtual void SetTriggerScaleFactors_Dn(const TriggerScaleFactors* const scaleFactors){(void)scaleFactors;};
    virtual void SetTriggerScaleFactors_EtaUp(const TriggerScaleFactors* const scaleFactors){(void)scaleFactors;};
    virtual void SetTriggerScaleFactors_EtaDn(const TriggerScaleFactors* const scaleFactors){(void)scaleFactors;};
    virtual void SetBtagScaleFactors_Var_AllSyst(const std::map<std::tuple<Systematic::Type, Systematic::Variation, int>, BtagScaleFactors*>* mapAllSystBTagScaleFactors){(void)mapAllSystBTagScaleFactors;};

    void SetJetVetoMap(const utils::JetVetoMap* const ptr){ jetVetoMap_ = ptr;}
    // // ------------------------------

    virtual void SetKinematicReconstruction_Up(const KinematicReconstruction* const kinematicReconstruction,
                                    const KinematicReconstructionScaleFactors* const kinematicReconstructionScaleFactors){(void)kinematicReconstruction;(void)kinematicReconstructionScaleFactors;};
    virtual void SetKinematicReconstruction_Dn(const KinematicReconstruction* const kinematicReconstruction,
                                    const KinematicReconstructionScaleFactors* const kinematicReconstructionScaleFactors){(void)kinematicReconstruction;(void)kinematicReconstructionScaleFactors;};
    virtual void SetLooseKinReco_Up(const LooseKinReco* const looseKinReco,
                                    const LooseKinRecoScaleFactors* const looseKinRecoScaleFactors){(void)looseKinReco;(void)looseKinRecoScaleFactors;};

    virtual void SetLooseKinReco_Dn(const LooseKinReco* const looseKinReco,
                                    const LooseKinRecoScaleFactors* const looseKinRecoScaleFactors){(void)looseKinReco;(void)looseKinRecoScaleFactors;};

    virtual void SetPlainTreeOutputDirectory(const TString& plainTreeOutputDirectory) {(void)plainTreeOutputDirectory;};
    virtual void SetTreeOutputDirectory(const TString& treeOutputDirectory) {(void)treeOutputDirectory;};

    virtual void SetRecoTree(const bool& makeRecoTree){(void)makeRecoTree;};
    virtual void SetSimpleRecoTree(const bool& makeSimpleRecoTree){(void)makeSimpleRecoTree;};
    virtual void SetGenTree(const bool& makeGenTree){(void)makeGenTree;};

    virtual void SetTTBarSample(const bool& isTTBarMC){(void)isTTBarMC;};
    virtual void SetSystSample(const bool& isSystMC){(void)isSystMC;};
    virtual void SetSystVarSample(const bool& isSystVar){(void)isSystVar;};
// ================================================================

private:

    /// Create binned control plots
    // Create Nbins control plots for the differential distribution h_differential
    // Use h_control for the control plot name and binning
    void CreateBinnedControlPlots(TH1* h_differential, TH1* h_control, const bool fromHistoList =true);

    /// Fill binned control plots
    // h: differential distribution histogram
    // binvalue: the value of the quantity in the differential distribution histogram
    // the control plot histogram
    // the value for the control plot
    // weight: event weight
    void FillBinnedControlPlot(TH1* h_differential, double binvalue,
                               TH1 *h_control, double value, double weight);



    /// Get weight of closure test
    double calculateClosureTestWeight(const Long64_t& entry);

    /// Select events from Top signal samples which need to be removed due to generator selection
    virtual bool failsTopGeneratorSelection(const Long64_t& entry)const;

    void generatorTopEvent(LV& leadGenTop, LV& nLeadGenTop,
                           LV& leadGenLepton, LV& nLeadGenLepton,
                           LV& leadGenBJet, LV& nLeadGenBJet,
                           double& genHT,
                           const int bHadronIndex, const int antiBHadronIndex,
                           const double trueLevelWeightNoPileup, const double trueLevelWeight,
                           const TopGenObjects& topGenObjects);

    void generatorTTbarjetsEvent(double& jetHTGen,
                                 const int bHadronIndex, const int antiBHadronIndex,
                                 const double trueLevelWeight,
                                 std::vector<int>& genJetIndices,std::vector<int>& genJet40Indices, std::vector<int>& genJet60Indices,
                                 std::vector<int>& genJet100Indices, std::vector<int>& generatorExtraJetIndices,
                                 const TopGenObjects& topGenObjects);

    void pseudoTopEvent(LV& leadPseudoTop, LV& nLeadPseudoTop,
                        LV& leadPseudoLepton, LV& nLeadPseudoLepton,
                        LV& leadPseudoBJet, LV& nLeadPseudoBJet,
                        const double trueLevelWeightNoPileup,
                        const double trueLevelWeight,
                        std::vector<int>& pseudoJetIndices, std::vector<int>& pseudoJet40Indices, std::vector<int>& pseudoJet60Indices, std::vector<int>& pseudoJet100Indices,
                        const TopPseudoObjects& topPseudoObjects);

    void generatorVisJets(const TopGenObjects& topGenObjects,std::vector<int>& genVisJetIndices);

    void generatorExtraJets(const TopGenObjects& topGenObjects,std::vector<int>& genExtraJetIndices, const int bHadronIndex, const int antiBHadronIndex);
    /// Select from a vector of lepton indices those whose fulfill 2-dimensional isolation criterion
    void lepton2dIsolationIndices(std::vector<int>& v_index, const VLV& v_lvLep, const VLV& v_lvJet)const;

    void boostLeptonsToRespectiveTopQuarkRestFrames(LV& lv_lepton, LV& lv_antilepton, const LV& lv_top, const LV& lv_antitop);

    /// Fill all analysers and histograms in one method
    virtual void fillAll(const std::string& selectionStep,
                 const EventMetadata& eventMetadata,
                 const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                 const TopGenObjects& topGenObjects, const TopPseudoObjects& topPseudoObjects,
                 const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                 const LooseKinRecoSolution& looseKinRecoSolution,
                 const ttbar::GenObjectIndices& genObjectIndices, const ttbar::RecoObjectIndices& recoObjectIndices,
                 const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                 const double& defaultWeight, const double& defaultWeightLooseKinReco)const;

    /// Book all histograms of all analysers for all steps in one method
    void bookAll();

    /// Clear all analysers in one method
    void clearAll();

    /// Reference to the analysis config
    const AnalysisConfig& analysisConfig_;

    /// Do kinematic reconstruction on nTuple level
    bool kinRecoOnTheFly_;

    /// Histogram for total weight of closure test
    TH1 *h_ClosureTotalWeight;

    /// Histogram for total weight of PDF variations
    TH1 *h_PDFTotalWeight;

    /// Histogram for total weight of PS variations
    TH1 *h_PSTotalWeight;

    /// Whether to apply closure test
    bool doClosureTest_;

    /// Data for closure test
#ifndef __CINT__
    std::function<double(Long64_t)> closureFunction_;
#endif
    int closureMaxEvents_;

    /// Whether it is leptonic decays via tau in ttbar dilepton samples
    bool runViaTau_;

    /// Map holding binned control plots
    //binnedControlPlots contains:
    //map of name of differential distribution
    // -> pair( histogram with the binning of the differential distribution,
    //          vector(bin) -> map( control plot name -> TH1*))
    std::map<std::string, std::pair<TH1*, std::vector<std::map<std::string, TH1*> > > >* binnedControlPlots_;

    /// All analysers of type AnalyzerBaseClass
    std::vector<AnalyzerBaseClass*> v_analyzer_;

    /// All tree handlers of type TreeHandlerBase
    std::vector<TreeHandlerBase*> v_treeHandler_;

    // OZ 30.07.2017 to access everything about the event from there
    friend class TopAnalysisPlain;
    friend class MiniTreeTopAnalysis;

};


#endif
