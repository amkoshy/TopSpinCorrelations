#ifndef MINITREEANALYZER_H
#define MINITREEANALYZER_H

#include "Math/GenVector/VectorUtil.h"

#include <MiniTreeReader.h>
#include <MiniTreeHelpers.h>
#include <sampleHelpers.h>
#include "../../common/include/utils.h"
#include "../../common/include/classes.h"
#include "../../common/include/classesFwd.h"
#include "PlotterConfigurationHelper.h"
#include "RootFileReader.h"

#include "TopAnalysis/ZTopUtils/interface/TFModel.h"

using namespace std;
using namespace ROOT;

// define methods nedded in unordered_map
namespace std {
    template<>
    struct hash<Channel::Channel>
    {
       typedef Channel::Channel argument_type;
       typedef size_t result_type;

       result_type operator () (const argument_type& x) const
       {
          using type = typename std::underlying_type<argument_type>::type;
          return std::hash<type>()(static_cast<type>(x));
       }
    };
    template<>
    struct hash<Process::Process>
    {
       typedef Process::Process argument_type;
       typedef size_t result_type;

       result_type operator () (const argument_type& x) const
       {
          using type = typename std::underlying_type<argument_type>::type;
          return std::hash<type>()(static_cast<type>(x));
       }
    };
    template<>
    struct hash<RecoRhoBinCategory::RecoRhoBinCategory>
    {
       typedef RecoRhoBinCategory::RecoRhoBinCategory argument_type;
       typedef size_t result_type;

       result_type operator () (const argument_type& x) const
       {
          using type = typename std::underlying_type<argument_type>::type;
          return std::hash<type>()(static_cast<type>(x));
       }
    };
    template<>
    struct hash<NJetCategory::NJetCategory>
    {
       typedef NJetCategory::NJetCategory argument_type;
       typedef size_t result_type;

       result_type operator () (const argument_type& x) const
       {
          using type = typename std::underlying_type<argument_type>::type;
          return std::hash<type>()(static_cast<type>(x));
       }
    };
    template<>
    struct hash<BTagCategory::BTagCategory>
    {
       typedef BTagCategory::BTagCategory argument_type;
       typedef size_t result_type;

       result_type operator () (const argument_type& x) const
       {
          using type = typename std::underlying_type<argument_type>::type;
          return std::hash<type>()(static_cast<type>(x));
       }
    };
}


class MiniTreeAnalyzer {
    public:
    vector<BTagCategory::BTagCategory> NBJetCategories = {
                        BTagCategory::ZeroPlusGreaterTwoBTag,
                        BTagCategory::OneBTag, BTagCategory::TwoBTag,
                        BTagCategory::InclusiveBTag};
    vector<NJetCategory::NJetCategory> NJetCategories = {
                        NJetCategory::ZeroAndOneJet, NJetCategory::TwoJet, NJetCategory::GreaterTwoJet,
                        NJetCategory::InclusiveNJet};
    vector<RecoRhoBinCategory::RecoRhoBinCategory> RecoRhoBinCategories = {
        RecoRhoBinCategory::NoRho, RecoRhoBinCategory::Rho0to02, RecoRhoBinCategory::Rho02to03,
        RecoRhoBinCategory::Rho03to04, RecoRhoBinCategory::Rho04to05,
        RecoRhoBinCategory::Rho05to06, RecoRhoBinCategory::Rho06to1,
        RecoRhoBinCategory::InclusiveRho,
        RecoRhoBinCategory::RhoStudies, RecoRhoBinCategory::Rho0to1};

    static const unsigned int nRecoRhoBins = 6;
    static const unsigned int nGenRhoBins = 6;
    float  recoRhoBinning[nRecoRhoBins+1]= {0., 0.2, 0.3, 0.4, 0.5,0.6, 1.0};
    float  genRhoBinning[nGenRhoBins+1] {0., 0.2, 0.3, 0.4, 0.5,0.6, 1.0};

    float scaleValue = 340.;

    float metCut = 40.;
    float mllCut = 20.;

    bool isData = true;
    bool isMC = false;
    bool isSignal = false;
    TString inputFileName = "";
    TString inputFilePath = "";
    TString outFilePath;

    bool hasPDFWeights = false;
    bool hasPSWeights = false;
    bool hasBFragWeights = false;
    bool hasMEWeights = false;
    bool hasBSemilepWeights = false;

    double lumiWeight = 1.;
    double renormWeight = 1.;

    uint maxSize = 3;

    ztop::TFModel *dnnModel;
    double dnn_inputs[21];
    vector<float> dnn_outputs;

    void InitDNNModel(TString& pathToModel);
    void EvaluateDNN(Event& event);

    // TH1D *h_eventsPassedStep3;
    // TH1D *h_eventsPassedStep4;
    // TH1D *h_eventsPassedStep5;
    // TH1D *h_eventsPassedStep6;
    // TH1D *h_eventsPassedStep7;
    // TH1D *h_eventsPassedStep8;
    // TH1D *h_eventsPassedStep8L;

    //Measure execution time
    chrono::steady_clock::time_point start;
    chrono::steady_clock::time_point end;

    int year;
    TString yearString;

    RootFileReader *fileReader;
    PlotterConfigurationHelper *configHelper;

    Event *event;

    Process::Process processType = Process::undefined;
    Systematic::Systematic systematic;
    Channel::Channel channel = Channel::undefined;

    public:
    MiniTreeAnalyzer(const TString &year_string, const TString &filename_string, const TString &systematic_string, const TString &channel_string);
    void SetOutpath(const TString &path_);

    bool GetZWindowCutDecisionStep4(const Event& ev);
    bool GetMETCutDecisionStep6(const Event& ev);

    void PushBackHistos(const BTagCategory::BTagCategory &btag,
            const RecoRhoBinCategory::RecoRhoBinCategory &rhobin, const NJetCategory::NJetCategory &njet, const TString histoname,
            const TString title, const int nbins, const float minX, const float maxX, const UserHisto::HistoKey &key, int &nHistos);
    void PushBackHistos(const BTagCategory::BTagCategory &btag,
            const RecoRhoBinCategory::RecoRhoBinCategory &rhobin, const NJetCategory::NJetCategory &njet, const TString histoname,
            const TString title, const int nbins, const float* xbins, const UserHisto::HistoKey &key, int &nHistos);
    void PushBackHistos2D(const BTagCategory::BTagCategory &btag,
            const RecoRhoBinCategory::RecoRhoBinCategory &rhobin, const NJetCategory::NJetCategory &njet, const TString histoname,
            const TString title, const int nbins, const float minX, const float maxX, const UserHisto::HistoKey &keyX,
            const int nbinsY, const float minY, const float maxY, const UserHisto::HistoKey &keyY, int &nHistos);

    bool GetCategoryFillDecision(const Event &ev, const BTagCategory::BTagCategory &btagcategory, const RecoRhoBinCategory::RecoRhoBinCategory &rhobin, const NJetCategory::NJetCategory &njetcategory);
    void FillHistograms(unordered_map<BTagCategory::BTagCategory,unordered_map<RecoRhoBinCategory::RecoRhoBinCategory,unordered_map<NJetCategory::NJetCategory,unordered_map<Process::Process,vector<shared_ptr<UserHisto::Histo>>>>>> &histMap);
    void FillHistograms2D(unordered_map<BTagCategory::BTagCategory,unordered_map<RecoRhoBinCategory::RecoRhoBinCategory,unordered_map<NJetCategory::NJetCategory,unordered_map<Process::Process,vector<shared_ptr<UserHisto::Histo2D>>>>>> &histMap);
    void InitHistograms(const vector<BTagCategory::BTagCategory> &btags,
                        const vector<RecoRhoBinCategory::RecoRhoBinCategory> &rhobins, const vector<NJetCategory::NJetCategory> &njets);

    unordered_map<BTagCategory::BTagCategory,unordered_map<RecoRhoBinCategory::RecoRhoBinCategory,unordered_map<NJetCategory::NJetCategory,unordered_map<Process::Process,vector<shared_ptr<UserHisto::Histo>>>>>> h1Map;
    unordered_map<BTagCategory::BTagCategory,unordered_map<RecoRhoBinCategory::RecoRhoBinCategory,unordered_map<NJetCategory::NJetCategory,unordered_map<Process::Process,vector<shared_ptr<UserHisto::Histo2D>>>>>> h2Map;

    bool Analyze();

    template<typename T>
    void save2File(unordered_map<BTagCategory::BTagCategory, unordered_map<RecoRhoBinCategory::RecoRhoBinCategory, unordered_map<NJetCategory::NJetCategory, unordered_map<Process::Process, vector<T>>>>> & hMaps, TFile& file);

    void FixHistogram(TH1& histo);
    void FixHistogram(TH2& histo);

    void ProgressBar(const int &progress);
    void WriteOutput();

};

#endif
