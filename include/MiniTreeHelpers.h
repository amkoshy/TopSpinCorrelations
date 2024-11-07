#ifndef MINITREEHELPERS_H
#define MINITREEHELPERS_H

#include "Math/GenVector/VectorUtil.h"

#include <sampleHelpers.h>
#include "../../common/include/utils.h"
#include "../../common/include/classes.h"
#include "../../common/include/classesFwd.h"
#include <TSystem.h>
#include <TH1F.h>
#include <TH2F.h>

using namespace std;
using namespace ROOT;


TString assignFolder(const char* baseDir, const Channel::Channel& channel, const TString& systematic, const char* subDir ="");

namespace Process{
    enum Process{
        ttbarsignal, ttbarbg, ttbarX, dy, diboson, wjets, undefined, data, singletop,
        ttbarsignal_genRho0to02, ttbarsignal_genRho02to03, ttbarsignal_genRho03to04, ttbarsignal_genRho04to05, ttbarsignal_genRho05to06,
        ttbarsignal_genRho06to1,
        ttbarsignal_genRho0to1,
        ttbarsignal_BKGNoAddJet
    };
    /// Convert a BTags from string to enum
    Process convert(const TString& process);
    Process convertFromFilename(const TString& process);
    Process convertFromBinning(const float& bin);
    /// Convert a BTags from enum to string
    TString convert(const Process& process);
    /// Convert a vector of btags from string to enum
    std::vector<Process> convert(const std::vector<TString>& processes);
    /// Convert a vector of btags from string to enums
    std::vector<Process> convert(const std::vector<std::string>& processes);
    /// Convert a vector of btags from string to enum
    std::vector<TString> convert(const std::vector<Process>& processes);
}




namespace BTagCategory{
    /// All btag categories as needed in any part of the framework
    enum BTagCategory{
        ZeroBTag, OneBTag, TwoBTag, ThreeBTag, FourBtag, InclusiveBTag, undefined,
        GreaterZeroBTag, GreaterOneBTag, GreaterTwoBTag, ZeroPlusGreaterTwoBTag
    };
    /// Convert a BTags from string to enum
    BTagCategory convert(const TString& btag);
    /// Convert a BTags from enum to string
    TString convert(const BTagCategory& btag);
    /// Convert a vector of btags from string to enum
    std::vector<BTagCategory> convert(const std::vector<TString>& btags);
    /// Convert a vector of btags from string to enums
    std::vector<BTagCategory> convert(const std::vector<std::string>& btags);
    /// Convert a vector of btags from string to enum
    std::vector<TString> convert(const std::vector<BTagCategory>& btag);
}
namespace NJetCategory{
    /// All btag categories as needed in any part of the framework
    enum NJetCategory{
        ZeroJet, OneJet, TwoJet, ThreeJet, FourJet, InclusiveNJet, GreaterOneJet, undefined,GreaterTwoJet, ZeroAndOneJet
    };
    /// Convert a BTags from string to enum
    NJetCategory convert(const TString& njet);
    /// Convert a BTags from enum to string
    TString convert(const NJetCategory& njet);
    /// Convert a vector of btags from string to enum
    std::vector<NJetCategory> convert(const std::vector<TString>& njets);
    /// Convert a vector of btags from string to enums
    std::vector<NJetCategory> convert(const std::vector<std::string>& njets);
    /// Convert a vector of btags from string to enum
    std::vector<TString> convert(const std::vector<NJetCategory>& njet);
}
namespace RecoRhoBinCategory{
    /// All btag categories as needed in any part of the framework
    enum RecoRhoBinCategory{
        NoRho,
         undefined, InclusiveRho, RhoStudies,Rho0to1,
         Rho0to02, Rho02to03, Rho03to04,Rho04to05,Rho05to06,Rho06to1
    };
    /// Convert a BTags from string to enum
    RecoRhoBinCategory convert(const TString& rhobin);
    /// Convert a BTags from enum to string
    TString convert(const RecoRhoBinCategory& rhobin);
    /// Convert a vector of btags from string to enum
    std::vector<RecoRhoBinCategory> convert(const std::vector<TString>& rhobins);
    /// Convert a vector of btags from string to enums
    std::vector<RecoRhoBinCategory> convert(const std::vector<std::string>& rhobins);
    /// Convert a vector of btags from string to enum
    std::vector<TString> convert(const std::vector<RecoRhoBinCategory>& rhobins);

    pair<float,float> getLowerUpperEdge(const RecoRhoBinCategory& rhobin);
}

TString RenameSystematicString(const Systematic::Systematic& sys);

namespace Systematic{
    const std::vector<Type> allowedSystematics = {
        nominal,
        pu, trig,
        eleID, muonID,
        eleReco, muonIso,
        jer, jes,
        jerEta0, jerEta1,
        unclustered,
        btag,
        btagLjet,
        // lumi,
        cp5,
        topPt,
        match,
        meScale, meFacScale, meRenScale,
	    psISRScale, psFSRScale, psScale,
        alphasPdf, pdf, l1prefiring,
        bFrag, bFrag_Peterson,
        bSemilep,
        mass,
        erdon, erdonretune, gluonmovetune,
        ueTune,
        // independent variation of JES sources
        jesAbsoluteStat,
        jesAbsoluteScale,
        jesAbsoluteFlavMap,
        jesAbsoluteMPFBias,
        jesFragmentation,
        jesSinglePionECAL,
        jesSinglePionHCAL,
        jesFlavorQCD,
        jesFlavorRealistic,
        jesFlavorPureCharm,
        jesFlavorPureBottom,
        jesFlavorPureGluon,
        jesFlavorPureQuark,
        jesTimePtEta,
        jesRelativeBal,
        jesRelativeJEREC1,
        jesRelativePtBB,
        jesRelativePtEC1,
        jesRelativeFSR,
        jesRelativeStatFSR,
        jesRelativeStatEC,
        jesPileUpDataMC,
        jesPileUpPtRef,
        jesPileUpPtBB,
        jesPileUpPtEC1,
    };
}


struct Event{
    double weight;
    double weightWithBtagSF;
    double genLevelWeight;
    // bool hasLKRS;
    bool hasKRS;
    bool passStep3;
    // bool passStep4;
    // bool passStep5;
    // bool passStep6;
    // bool passStep7;
    // bool passStep8;
    // bool passStep8Loose;
    VLV jets;
    VLV bjets;
    LV lepton1;
    LV lepton2;
    // VLV leptons;
    LV dilepton;
    LV met;
    float mlb_min;
    LV kinReco_top;
    LV kinReco_antitop;
    // LV kinReco_lepton;
    // LV kinReco_antilepton;
    // LV kinReco_b;
    // LV kinReco_antib;
    // VLV kinReco_nonbjets;
    LV kinReco_nonbjet;
    // LV kinReco_ttbar;
    float kinReco_rho;
    // float kinReco_SF;
    // LV looseKinReco_ttbar;
    // LV looseKinReco_lepton1;
    // LV looseKinReco_lepton2;
    // LV looseKinReco_b1;
    // LV looseKinReco_b2;
    // VLV looseKinReco_nonbjets;
    // float looseKinReco_rho;
    // float looseKinReco_SF;
    LV gen_top;
    // VLV gen_additional_jets;
    // LV gen_additional_jet;
    LV gen_antitop;
    // LV gen_lepton;
    // LV gen_antilepton;
    // LV gen_b;
    // LV gen_antib;
    LV gen_ttbar;
    // LV gen_met;
    // float gen_rho;

    LV gen_additional_jet;
    float gen_rho;

    bool isTTJSignal;
    bool isTTSignalBackground;

    float DNN_rho;
    // float DNN_rho_new;

    float final_rho;

    LV tempJet;

    void Init(){
        jets.reserve(10);
        bjets.reserve(10);
    }
    void Clear();
    void Fill(const MiniTreeReader &reader_, const bool isMC, const bool isSignal, const Systematic::Systematic &sys, const bool hasPSWeights, const bool hasPDFWeights,
                const bool hasMEWeights, const bool hasBFragWeights, const bool hasBSemilepWeights);
};


namespace UserHisto{
    enum HistoKey{
        nEvents,
        // recoRhoKRS,recoRhoLKRS,
        recoRhoDNN,
        dileptonMass,
        jet1Px,jet1Py,jet1Pz,jet1M,jet1Pt,jet1Eta,jet1Phi,jet1E,
        jet2Px,jet2Py,jet2Pz,jet2M,jet2Pt,jet2Eta,jet2Phi,jet2E,
        jet3Px,jet3Py,jet3Pz,jet3M,jet3Pt,jet3Eta,jet3Phi,jet3E,
        nJets,nBJets,
        // mTTbarKRS,mTTbarLKRS,
        // addJetPtKRS,addJetPtLKRS,
        metPx,metPy,metM,metPt,metPhi,metE,
        bJet1Pt,mlb_min,

        genRho,genmTTbar,
        // genKRSRecoRhoDiff,genLKRSRecoRhoDiff,
        genRecoRhoDNNDiff
    };

    const TString GetNameToKey(const HistoKey &key);
    const TString GetNameToKeys2D(const HistoKey &keyX, const HistoKey &keyY);

    float GetValueToKeyFromEvent(const Event &ev, const HistoKey &key);
    bool CheckKeyConditionFromEvent(const Event &ev, const HistoKey &key);

    struct Histo{
        shared_ptr<TH1F> histogram;
        HistoKey variableKey;
        const TString GetName(){
            return UserHisto::GetNameToKey(variableKey);
        };
    };
    struct Histo2D{
        shared_ptr<TH2F> histogram;
        HistoKey variableKeyX;
        HistoKey variableKeyY;
        const TString GetName(){
            return UserHisto::GetNameToKeys2D(variableKeyX, variableKeyY);
        };
    };

}




#endif
