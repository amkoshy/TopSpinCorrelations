#ifndef PlotterConfigurationHelper_h
#define PlotterConfigurationHelper_h

//#include <vector>
#include <set>

#include <map>
#include <utility>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <vector>

#include <TMath.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TString.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TRandom3.h>
#include "TreePlain.h"
#include "utils.h"
#include <TChain.h>

class TString;
class TH1;
class TH2;

class RootFileReader;

class PlotterConfigurationHelper
{
  public:

    /// Full Constructor, use for all options
  PlotterConfigurationHelper(RootFileReader*  rootFileReader, const TString& year, const TString& yearCorr = "", bool doDYScale = false, bool setDiffXSec=false, bool flagProcSyst=false);
    /// Partial Constructor, use if only basic year options are needed: systematics vectors, luminosity, etc.
    PlotterConfigurationHelper(const TString& year, const TString& yearCorr = "", bool setDiffXSec = false, bool flagProcSyst = false);

    /// Default destructor
    ~PlotterConfigurationHelper(){}

    void setYearOptions(bool setDiffXSec, bool flagProcSyst);
    void setPlotterDiffXSecVars();
    static void fillVectorOfValidSystematics(std::vector<TString>& vect_all, const TString& year, const TString& yearCorr, const bool runningDiffXSec = false,
                                             const bool useAdditionalSystematics = false, const bool fill_as_Up_Down = true,
                                             const bool include_all_as_syst = true, const bool include_Nominal_as_syst = true,
                                             const bool include_vars_with_noUpDown = true, const bool only_include_vars_with_noUpDown = false,
                                             TString sys_type = "all", const bool printExclusionWarning = false);
    double SampleXSection(const TString& filename);
    static void fillLegendColorDataset(const TString& fileListName, std::vector<TString>& legends, std::vector<int>& colors, std::vector<TString>& dataset);
    double CalcLumiWeight(const TString& WhichSample, const TString in_year = "use initialized");
    void ApplyFlatWeights(TH1* varhists,   const double weight);
    void ApplyFlatWeights(TH2* varhists,   const double weight);
    void RescaleHistYearToYearCorr(TH1* varHist, TH1* nominalHist, const double year_corr_rescale_from_file);
    void RescaleHistYearToYearCorr(TH2* varHist, TH2* nominalHist, const double year_corr_rescale_from_file);
    std::vector<TString> InputNominalFileList(TString mode);
    void DYScaleFactor(TString SpetialComment,std::vector<double>& DYScale,TString name);
    void DYScaleFactorFromFit(TString SpetialComment,std::vector<double>& DYScale,TString name);
    double GetTrueLevelRenormalisationWeight(const TString& WhichSample, const TString &varName="Nominal");
    void DefaultSetup();

    // ======Porting of ttmd_settings========
     void fillPlainTreesFileListSampleVectors(TString in_year, std::vector<TString>& pattern_vector, std::vector<TString>& output_vector, TString sys = "", TString channel = "");
     void fillPlainTreesSampleIdVector(std::vector<TString>& fileListVector, std::vector<std::pair<TString, int>>& pattern_vector, std::vector<std::pair<TString, int>>& output_vector);
     static const TString ExtractFileNameFromPath(TString path, bool exclude_channel_label = false);
     const std::vector<TString> GetSysVars(const TString& sys, bool getEnvelopesTotVar=false);
     const std::vector<TString>& GetAllSys(const bool fill_Up_Down);
     const std::vector<TString>& GetJECSrc(const bool fill_Up_Down);
     const std::vector<TString>& GetExpSys(const bool fill_Up_Down);
     const std::vector<TString>& GetModPDFEig();
     bool IsModPDFEig(const TString& syst);
     const std::vector<TString>& GetModWSys(const bool fill_Up_Down);
     const std::vector<TString>& GetModISys(const bool fill_Up_Down);
     const std::vector<TString>& GetAltTheorSys(const bool fill_Up_Down);
     const std::vector<TString>& GetAllVar(const bool fill_Up_Down);
     bool IsExpSys(const TString& syst);
     bool IsModWSys(const TString& syst);
     bool IsModISys(const TString& syst);
     bool IsAltTheorSys(const TString& syst);
     bool IsNoUpDown(const TString& syst);
     static bool IsYearToYearCorr(const TString& year, const TString& yearCorr, const TString& sys);
     static const std::vector<TString> GetYearCorrList(const TString& year, const TString& yearCorr, const TString& sys, const bool isYearToYearCorr);
     static void AppendYearCorr(const TString& year, const TString& labelYearCorr, TString& sys, const bool isYearToYearCorr, const bool noUpDown = false);
     static void RemoveYearCorr(const TString& year, const TString& yearCorr, TString& sys, const bool noUpDown = false);
     const TString GetLabelYearCorr(const TString& var);

     /// Envelopes functions
     /**
      * @brief IsEnvelope: Checks if the input is an envelope name
      * @param sys: input systematic
      * @return : true if is an envelope
      */
     bool IsEnvelope(const TString& sys);

     /**
      * @brief GetEnvelopeVars: Get all the variations for a given envelope
      * @param sys: envelope name
      * @param get_UD_vars: if false only variations names will be obtained, excluding the _UP/_DOWN.
      * @return : list of systematics
      */
     const std::vector<TString> GetEnvelopeVars(const TString& sys, bool get_UD_vars=true);

     /**
      * @brief GetEnvelopeVars: Get the UP/DOWN variations for a given envelope, envelope_name_UP and envelope_name_DOWN
      * @param sys: envelope name
      * @return : envelope_name_UP and envelope_name_DOWN
      */
     const std::pair<TString,TString> GetTotEnvelopeVars(const TString& sys);

     /**
      * @brief GetEnvelopesList: Get the list of available envelopes
      * @param fill_tot_up_down_variations: if true -> UP/DOWN variation are included
      * @return : list of available envelopes
      */
     const std::vector<TString> GetEnvelopesList(const bool fill_tot_up_down_variations=false);

     /**
      * @brief GetVarEnvelope: Get the envelope for a given systematic
      * @param var: systematic variation
      * @return : the envelope name ("" if is not contained in any envelope)
      */
     const TString GetVarEnvelope(const TString& var);

     /**
      * @brief IsVarInEnvelope: Check if a variation is contained in an specific envelope
      * @param var: variation to check
      * @param envelope: envelope where the variation will be check
      * @return: true if is contained
      */
     bool IsVarInEnvelope(const TString& var, const TString& envelope);

     /**
      * @brief FoldSystematicVarsListIntoEnvelops: Takes a list of systematics and fold the systematics in their respective envelopes
      * @param list_of_systematics: list of systematics
      * @param list_with_variations_included: if true -> the input list should contain all _UP/DOWN variations.
      * @return list of systematics folded into the envelopes
      */
     const std::vector<TString> FoldSystematicVarsListIntoEnvelops(const std::vector<TString>& list_of_systematics, const bool list_with_variations_included = false);

     const TString VarToSys(const TString& var, int& sign);
     const TString VarToSys(const TString& var);
     const TString GetVarDirTree(const TString& var);
     const TString GetVarWeight(const TString& var);
     const TString GetVarDirTreePlain(const TString& var, const bool flagWeights = true);
     const TString GetVarDirPlot(const TString& var);
     TString GetShortVarName(const TString& var);
     TString GetShortVarNameForTex(const TString& var);
     static double GetYearToYearCorrRescale(const TString& year, const TString& yearCorr, const TString& complete_Var_name);
     double GetSysRescale(const TString& sys);
     double GetVarRescaleSample(const TString& var, TString mode = "default", const TString i_year = "");
     TString GetVarPlainTreeFileSuffix(const TString& var, const TString& complete_Var_name, const TString& use_year = "");
     TString GetSysBkgDir(const TString& sys);
     std::vector<TString> GetChannelsPlainTree(const TString& ch);
     double GetChannelBR(const TString& ch);
     TreePlain* GetTreeDataPlain(const TString& ch, const TString& dirSyst = "Nominal", const TString& wSyst = "");
     //long GetNEvents(const TString& nameFile, const TString& nameTree);
     TreePlain* GetTreeMCBgrPlain(const TString& complete_Var_name, const TString& ch, const TString& dirSyst = "Nominal", const TString& wSyst = "Nominal");
     TreePlain* GetTreeMCBttPlain(const TString& complete_Var_name, const TString& ch, const TString& dirSyst = "Nominal", const TString& wSyst = "Nominal", const bool flagMode = true);
     TreePlain* GetTreeMCRecPlain(const TString& complete_Var_name, const TString& ch, const TString& dirSyst = "Nominal", const TString& wSyst = "Nominal", const bool flagMode = true, const bool DoParticle = false);
     TreePlain* GetTreeMCGenPlain(const TString& complete_Var_name, const TString& ch, const TString& dirSyst = "Nominal", const TString& wSyst = "Nominal", const bool flagMode = true, const bool DoParticle = false);
     static std::vector<double> load_ttmd_DY_SFs(TString year);
     double get_DY_SF(TString year, TString channel);
     static std::pair<double, double> get_lumi_and_err(TString year);
     void get_year_from_file(const TString& filename);
     double get_year_corr_rescale_from_file(const TString& complete_Var_name, const TString& filepath);
     TString getCorrespondingNominalDataset(const TString& datasetPath, const TString& systematic, const TString& channel);
     void doYearToYearCorrRescaleSanityCheck(const TString& Var, const TString dirInput);

    std::vector<TString> _vectorOfSystematics;
    std::vector<TString> _vectorOfSystematics_noUpDown;
    std::vector<TString> _vectorOfSystematics_all;
    std::map<TString, std::vector<TString>> _envelopesMap;

    double energy, lumi, lumiError, topxsec, bRatio[4], bRatioError;
    TString lumi_year;
    std::ostringstream lumi_fb;
    std::vector<TString> years_list;
    TString energy_label;
    TString plotting_mode;
    double massUncReductionFactor;
    double fsrScaleUncReductionFactor;
    bool revokeAntiQuantityCombination;
    bool drawTotalMCErrorForCP;
    bool drawPlotRatio;
    bool makeControlPlotEventCountTable;
    bool runningDiffXSec;
    bool flagProcSyst;
    TString EventCountTableWhenHist;
    int cmsPrelimLabelMode;
    bool usePoissonSmearedPseudoData;
    bool estimateSysUncFromPseudoData;
    const TString baseFileList;
    const TString outputPlotsFolder;


    //Variables for TUnfold implementation
    TString gBaseDir;
    TString gPlainTreeVersion; // updated simple kin reco: configuration for the initial setup

    bool longTitleCrossSection; // Add extra information in the xsection title

    int gTTBARDebug;

    TString analysisSuffix; // optional analysis suffix for plots directory
    TString diLeptonic_dir;
    TString analysis_folder;
    TString fileList_prefix_plainTrees; // prefix in plainTrees FileList files
    TString fileList_format_plainTrees; // format extension in plainTrees FileList files
    //TString gTreeDir; // directory with plainTree: configuration for the initial setup
    TString gAnalDir; // directory with final plots
    TString gPseudoDir; // directory with toys
    TString outDir; // runtime active directory
    TString gCorrData;
    TString gForAN;
    TString gBTagEffData;
    TString gMGTabsDir;
    TString gNPCorrSubdir;
    int gFlagMG; // Mar18 -> Sep18

    TString gPlainTreeName;
    TString gPlainTreeName0;
    TString gPlainTreeNamePlain;
    TString gPlainTreeName0Plain;
    double gXsecTtbarTwiki;// / (0.02277);
    double gBrrTtbarToEMu;
    double gBrrTtbarToEMuUncRel;

    int gPlotPixelSmall;
    int gPlotPixelLarge;
    bool gSysPDF=0;
    int gSysJES; // 1 JES + MET (historical), 2 JES splitted + UNCLUSTERED
    int gOnlyDominantSyst; // 0 skip nothing, 1 skip non-dominant, 2 skip all

    // Systematic variation factor for Bkg and DY
    double gVarBkg, gVarDY;

    int gUnfUncRspTotal;

    // Plotting and output options
    bool ttmd_createEventSummaryTable;
    bool ttmd_plot_uncertainties_summary, ttmd_plot_only_Absolute_and_Regularized,
         ttmd_plot_only_Normalized_and_Regularized, ttmd_plot_only_Regularized,
         ttmd_plot_comparison_between_channels, ttmd_plot_comparison_between_years, ttmd_plot_systs_mig_matrix_comparison_between_years,
         ttmd_include_alternative_theory, ttmd_skip_extra_plots_for_systematics;
    // corrections options
    bool ttmd_apply_trueLevelWeights;
    bool doDYScaleFromFit;

    //Alternative Theory included in plots
    std::vector<TString> ttmd_alternative_MC_samples,ttmd_alternative_MC_samples_legend_names;

    TString GetCorrYearValue() {return yearCorr_; };

  private:

    /// File reader for accessing specific histogram from given file
    RootFileReader* fileReader;

    TString year_;
    TString yearCorr_;
    bool isSingleLeptonDatasetIncluded;
    bool doDYScale;
    bool flagUseZpeak;

    // *** Event summary ***
    std::vector<std::pair<TString, int> > _samples_event_summary_info;
    //Total events info variables
    unsigned int _data_events,_ttbar_gen_events, _ttbar_rec_events,_ttbar_bgviatau_events,_ttbar_bkg_events,_MCbkg_events;
    //Weighted events info variables
    double _ttbar_gen_events_w,_ttbar_rec_events_w,_ttbar_bgviatau_events_w,_ttbar_bkg_events_w,_MCbkg_events_w;

    //Sample pattern definition for ttmd
    std::vector<TString> _gVFilesData_pattern, _ttbar_bkg_pattern, _gVFilesMCbkg_pattern, _ttbar_signal_pattern, _ttbar_bgviatau_pattern, _DY_pattern;
    std::vector<std::pair<TString, int> > _gVFilesMCbkg_pattern_ID;
};

// random generator
class ZRandom
{
 public:
  static TRandom3& GetInstance()
  {
    printf("ZRandom::GetInstance()\n");
    static TRandom3 rnd;
    return rnd;
  }

 private:
  ZRandom() {printf("ZRandom::ZRandom()\n");}
  ZRandom(const ZRandom&);
  ZRandom& operator=(const ZRandom&);
};


// ======Porting of ttmd_settings========

namespace ZSet {

  enum Level
  {
    LevelRec = 0,
    LevelGen = 1,
    LevelRecSimpleKinRec = 2,
    LevelRecSimpleKinRecMlbCut = 3,
    LevelRecSimpleKinRecMlbCutCustom = 4,
    LevelRecSimpleKinRecMWXY = 5,
    LevelRecSimpleKinRecMWE = 6,
    LevelRecGeomKinRec = 7,
    LevelUndef = -1
  };

  enum AnalysisChannel
  {
    ChannelEM = -11 * 13,
    ChannelEE = -11 * 11,
    ChannelMM = -13 * 13,
    ChannelSE = 11,
    ChannelSM = 13,
    ChannelUndef = 0,
    ChannelDL = 1, // ee, mm, em
    ChannelAll = -1 // this is for true level MC all ttbar
  };

  enum SelectionStep
  {
    Step0 = 0,
    Step1 = 1,
    Step2 = 2,
    Step3 = 3,
    Step4 = 4,
    Step5 = 5,
    Step6 = 6,
    Step7 = 7,
    Step8 = 8,
    Step9 = 9,
    StepUndef = -1
  };
  //std::vector<SelectionStep> vSelectionStep;

  enum Variation
  {
    VarN = 0,
    VarU = 1,
    VarD = 2,
    VarUndef = -1
  };
}

#endif
