#ifndef AnalysisConfig_h
#define AnalysisConfig_h

#include <string>

#include "../../common/include/sampleHelpers.h"





class AnalysisConfig{

public:

    /// Struct for general info
    struct General{
        /// Constructor
        General();

        /// Print variables, return as string, and optionally print to screen
        std::string print(const bool screen =false)const;

        /// Analysis era
        Era::Era era_;

        /// Data luminosity in pb-1
        float luminosity_;

        /// Relative luminosity uncertainty
        float luminosityUncertainty_;
    };



    /// Struct for object/event corrections
    struct Corrections{
        /// Constructor
        Corrections();

        /// Print variables, return as string, and optionally print to screen
        std::string print(const bool screen =false)const;

        /// Pileup MC distribution file corresponding to MC sample in use
        std::string pileup_MC_File_;

        /// Pileup distribution histogram name corresponding to MC file given above
        std::string pileup_MC_Histogram_;

        /// Pileup distribution file corresponding to data sample in use
        std::string pileup_DATA_File_;

        /// Pileup distribution histogram name corresponding to data file given above
        std::string pileup_DATA_Histogram_;

        /// Pileup distribution files for syst. uncertainties corresponding to data sample in use
        std::string pileup_DATA_sysUP_File_;
        std::string pileup_DATA_sysDN_File_;


        /// File ending of dilepton trigger scale factors input file
        std::string triggerSFInputSuffix_;

        /// For now, there is no use, but the functionality for setting the mode is available
        /// As provided during readout of trigger eff. scale factors
        /// 0: no distortion, i.e. do nothing
        /// Updated: triggerSFReadoutMode is now used to set the case for a single input file (mode = 1)
        ///or for input files for ee, mumu, emu channels (mode = 2). Default is still 0 in which case
        ///the user will be reminded to specify the mode.
        int triggerSFReadoutMode_;
        std::string trigger_DoubleElectron_SF_File_;
        std::string trigger_DoubleElectron_SF_Histogram_;
        std::string trigger_DoubleElectron_SF_Histogram_Format_;
        std::string trigger_ElectronMuon_SF_File_;
        std::string trigger_ElectronMuon_SF_Histogram_;
        std::string trigger_ElectronMuon_SF_Histogram_Format_;
        std::string trigger_DoubleMuon_SF_File_;
        std::string trigger_DoubleMuon_SF_Histogram_;
        std::string trigger_DoubleMuon_SF_Histogram_Format_;

        /// Input file for electron ID scale factor
        std::string electronSFInputFile_;

        /// Input file for muon ID scale factor
        std::string muonSFInputFile_;

        /// Mirror lepton scale factors (SF) in eta
        /// Determine SF basing on absolute eta of the leptons
        /// To be configured in cases when one of the input files contains SF only with respect to absolute eta
        /// 0: Use SF as provided in the files (do nothing): sf for positive (negative) eta to be used for leptons with positive (negative) eta
        /// 1: Mirror SF for electrons only
        /// 2: Mirror SF for muons only
        /// 3: Mirror SF for electrons and muons
        int leptonSFReadoutMode_;

        // provide for both muon and electron SFs only one histogram (old workflow)
        bool combinedSF_;

        // electron ID SF File, Histogram name and readout format
        std::string SF_electron_ID_File_;
        std::string SF_electron_ID_Histo_;
        std::string SF_electron_ID_Format_;
        std::string SF_electron_Reco_File_;
        std::string SF_electron_Reco_Histo_;
        std::string SF_electron_Reco_Format_;
        std::string SF_electron_Reco_File_LowPt_;
        std::string SF_electron_Reco_Histo_LowPt_;
        std::string SF_electron_Reco_Format_LowPt_;
        // muon ID SF File, Histogram name and readout format
        std::string SF_muon_ID_File_;
        std::string SF_muon_ID_Histo_;
        std::string SF_muon_ID_Histo_Sys_;
        std::string SF_muon_ID_Histo_Stat_;
        std::string SF_muon_ID_Format_;
        std::string SF_muon_ISO_File_;
        std::string SF_muon_ISO_Histo_;
        std::string SF_muon_ISO_Histo_Sys_;
        std::string SF_muon_ISO_Histo_Stat_;
        std::string SF_muon_ISO_Format_;
        
        // extra DY extrapolation uncertainties applied on top of the lepton uncertainties
        float extra_muon_ISO_DY_extrapolation_unc_;
        float extra_electron_ID_DY_extrapolation_unc_;


        // Jet Veto Map File/Histograms
        std::string jetVetoMap_File_;
        std::string jetVetoMap_Histogram_;
        std::string jetVetoMap_Format_;

        // Apply Muon Energy Corrections
        bool applyMuonEnergyCorrections_;

        // Apply Electron Energy Corrections
        bool applyElectronEnergyCorrections_;

        // Apply Bjet Energy Regression
        bool applyBJetEnergyRegression_;

        // Propagate lepton and jet corrections to MET
        bool propagateCorrectionsToMET_;


        /// Types of jet corrections which should be applied
        /// JER applies only to MC, JES to data and MC
        /// Due to different precision of float and double, the JES re-calculation deviates by about 0.00001 GeV
        /// 0: Use corrected jets from ntuple, i.e. do nothing
        /// 1: Undo JER correction, i.e. jets as in ntuple without JER
        /// 2: Undo also JES correction, i.e. fully uncorrected jets
        /// 3: Apply JES correction, i.e. re-calculate JES correction for uncorrected jets
        /// 4: Apply JER correction, i.e. re-calculate JER correction for already JES re-corrected jets
        int jetCorrectionMode_;

        // Jet JER Corrections input(s)
        std::string jerUncertaintySourceName_;
        std::string jerSFBinning_;
        std::string jet_JERHybridMethod_PtResolutionFile_;
        float      jet_JERHybridMethod_jetConeRadius_;

        /// File containing the uncertainties associated to JES
        std::string jesUncertaintySourceFile_;

        /// File containing the reduced set of the uncertainties associated to JES
        std::string jesUncertaintyReducedSourceFile_;

        /// File containing the official JES scale factors in text (txt) format: L1 correction for Data
        std::string jesL1CorrectionDataFile_;

        /// File containing the official JES scale factors in text (txt) format: L2 correction for Data
        std::string jesL2CorrectionDataFile_;

        /// File containing the official JES scale factors in text (txt) format: L3 correction for Data
        std::string jesL3CorrectionDataFile_;

        /// File containing the official JES scale factors in text (txt) format: L2L3 residual correction for Data
        std::string jesL2L3CorrectionDataFile_;

        /// File containing the official JES scale factors in text (txt) format: L1 correction for MC
        std::string jesL1CorrectionMcFile_;

        /// File containing the official JES scale factors in text (txt) format: L2 correction for MC
        std::string jesL2CorrectionMcFile_;

        /// File containing the official JES scale factors in text (txt) format: L3 correction for MC
        std::string jesL3CorrectionMcFile_;


        /// Enable/disable the computation of JES uncertainty from the reduced set of sources
        bool useReducedJESSet_;

        /// Name of the signal ttbar file for b-tag efficienies
        std::string btagEfficiencyFilename_;

        /// Correction mode for the b-tagging
        Btag::CorrectionMode btagCorrectionMode_;

        /// File for the official heavy flavour scale factors for b-tag discriminator reweighting
        std::string btagHeavyFlavourFile_;

        /// File for the official light flavour scale factors for b-tag discriminator reweighting
        std::string btagLightFlavourFile_;

        /// File containing the official b-tag scale factors in comma separated value (csv) format
        std::string btagSFInputFile_;

        /// Perform correction of b-tagging scale factors from file via efficiencies for the selected systematic sources
        bool btagSFCorrectionViaEff_;

        /// Perform Top pt reweighting
        bool topPtReweighting_;

        /// Perform Top pt 2d reweighting
        bool applyTopSlope2DGenLevelReweighting_;

    };



    /// Struct for object selections
    struct Selections{
        /// Constructor
        Selections();

        /// Print variables, return as string, and optionally print to screen
        std::string print(const bool screen =false)const;

        /// Lepton eta selection (absolute value)
        float leptonEtaCut_;

        /// Lepton pt selection in GeV
        float leptonPtCut_;

        /// Leading lepton pt selection in GeV
        float leadingLeptonPtCut_;

        /// Muon isloation criteria in relative units of positive values above zero up to 1.0, as recommended for the moment of 20.02.2017 for PF-based combined relative, delta(Beta)-corrected algorithm:
        /// tight = 0.15
        /// loose = 0.25
        /// ... no application, if value is set to "-999.0", i.e. do nothing
        float muonIsoCut_;

        /// Jet eta selection (absolute value)
        float jetEtaCut_;

        /// Jet pt selection in GeV
        float jetPtCut_;

        /// Jet PF ID selection
        int jetPFIDWorkingPoint_;

        /// Jet pileup ID selection
        JetPileupId::WorkingPoint jetPileupIdWorkingPoint_;
        float jetPileupIdMinPt_;
        float jetPileupIdMaxPt_;

        /// Jet veto map selection
        float jetVetoMapCutValue_;

        /// Minimal deltaR for removal of jets that are close to leptons (if negative, no cut applied)
        float deltaRLeptonJetCut_;

        /// Leading 2 jet pt selection in GeV
        float lead2JetPtCut_;

        /// B-tag algorithm
        Btag::Algorithm btagAlgorithm_;

        /// B-tag working point
        Btag::WorkingPoint btagWorkingPoint_;

        /// MET algorithm
        Met::Algorithm metAlgorithm_;

        /// MET selection for same-flavour channels (ee, mumu)
        float metCut_;

        /// Generated jet eta selection (absolute value)
        float genJetEtaCut_;

        /// Generated jet pt selection in GeV
        float genJetPtCut_;

        /// Minimal deltaR for removal of generated jets that are close to leptons (if negative, no cut applied)
        float genDeltaRLeptonJetCut_;

        /// Whether veto events with more than two selected leptons
        /// 0: no veto, i.e. do nothing
        /// 1: veto events with more than two selected leptons
        int vetoEventsWithMoreThan2Leptons_;
    };



    /// Struct for sample compostion
    struct SampleComposition{
        /// Constructor
        SampleComposition();

        /// Print variables, return as string, and optionally print to screen
        std::string print(const bool screen =false)const;

        /// Whether to use single lepton trigger info alongside single lepton and dilepton data streams
        /// false or 0: use trigger menu corresponding for standalone usage of dilepton data streams
        /// true or 1: use trigger menu corresponding for usage of single lepton data streams alongside dilepton ones (supported starting from 2016 analysis, otherwise use "0")
        bool addSingleLeptonTriggerInfo_;

        /// Whether to use pseudodata and how
        /// 0: use real data
        /// 1: use pseudodata as stacksum
        int pseudodata_;

        /// Level of merging for different samples
        /// 0: no merging, use all defined processes individually
        /// 1-x: different merging schemes
        int mergeLevel_;

        /// Whether to use pseudo-top definiton (particle level definition within TOPLHCWG) and create/fill related collections. Mostly for use in 'diLeptonic' framework.
        /// 0: no usage of pseudo-top, running analysis solely with parton level top
        /// 1: in addition to parton level top, process and fill pseudo-top info in fiducial phase space; events outside of fiducial region are not considered (i.e. info is lost & forgotten)
        /// 2: in addition to parton level top, process and fill pseudo-top info in fiducial phase space; events outside of fiducial region (slightly looser definition wrt nominal fiducial)
        ///    are filled to underflow bins (in histograms needed for data unfolding) with negative weight of "-1000"
        /// 3: in addition to parton level top, process and fill pseudo-top info in fiducial phase space; events outside of fiducial region are filled to "OutOfSpace" histograms, as hypothesis reco. info,
        ///    for further treatment of such events as in case of the background source
        int pseudoTopMode_;
    };



    /// Struct for plot style
    struct PlotStyle{
        /// Constructor
        PlotStyle();

        /// Print variables, return as string, and optionally print to screen
        std::string print(const bool screen =false)const;

        /// CMS label
        /// -1: none
        ///  0: CMS
        ///  1: CMS preliminary
        ///  2: CMS private work
        int cmsLabel_;
    };



    /// Constructor reading in config file
    AnalysisConfig(const std::string& configfilename ="config.txt");

    /// Destructor
    ~AnalysisConfig(){}



    /// Return constant reference of struct general
    const General& general()const{return general_;}

    /// Return constant reference of struct object/event corrections
    const Corrections& corrections()const{return corrections_;}

    /// Return constant reference of struct object selections
    const Selections& selections()const{return selections_;}

    /// Return constant reference of struct sample compostion
    const SampleComposition& sampleComposition()const{return sampleComposition_;}

    /// Return constant reference of struct plot style
    const PlotStyle& plotStyle()const{return plotStyle_;}

    /// Print variables for all structs, return as string, and optionally print to screen
    std::string print(const bool screen =false)const;



private:

    /// Struct holding general info
    General general_;

    /// Struct holding object/event corrections
    Corrections corrections_;

    /// Struct holding object selections
    Selections selections_;

    /// Struct holding sample compostion
    SampleComposition sampleComposition_;

    /// Struct holding plot style
    PlotStyle plotStyle_;
};



#endif
