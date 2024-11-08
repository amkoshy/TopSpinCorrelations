#ifndef analysisHelpers_h
#define analysisHelpers_h

#include <vector>
#include <string>

class TString;




namespace AnalysisMode{

    /// Analysis modes for analysis
    enum AnalysisMode{
        cp,           // Basic control plots analyser
        gen,           // Basic gen level control plots analyser
        dda,          // Double differential analyser
        btop,         // Boosted top analyser
        ddaTree,      // Produce Tree for double differential analysis
        plainTree,    // Produce plain Tree with particle 4-momenta
        miniTree,    // Produce plain Tree with particle 4-momenta
        btopTree,     // Produce Tree for boosted top analysis
        kinReco,      // kin. Reco. efficiency analyser
      	looseKinReco, // Produce loose kin reco plots
      	kinRecoQualityStudies, // Produce plots for quality studies of full and loose kin recos
        dijet,        // Dijet analyser
        charge,       // Jet charge analyser
        match,        // Jet matching analyser
        playg,        // Playground analyser
        weight,       // Event weight analyser
        genEvent,     // Generator event analyser
        mvaTopP,      // Produce MVA input for top system jet assignment
        mvaTopA,      // Apply MVA weights for top system jet assignment
        mvaEventP,    // Produce MVA input for event classification
        mvaEventA,    // Apply MVA weights for event classification
        mvaChargeP,   // Produce MVA input for jet charge
        mvaChargeA,   // Apply MVA weights for jet charge
        spinCorr,     // Produce spin correlation histograms and branches
        undefined     // Undefined mode
    };

    /// All analysis modes allowed for analysis
    const std::vector<AnalysisMode> allowedAnalysisModes
    {cp, gen, dda, btop, ddaTree, plainTree, miniTree, btopTree, kinReco, looseKinReco, kinRecoQualityStudies, dijet, charge, match, playg, weight, genEvent, mvaTopP, mvaTopA, mvaEventP, mvaEventA, mvaChargeP, mvaChargeA, spinCorr};



    /// Convert an AnalysisMode from string to enum
    AnalysisMode convert(const TString& analysisMode);

    /// Convert an AnalysisMode from enum to string
    TString convert(const AnalysisMode& analysisMode);

    /// Convert a vector of AnalysisModes from string to enum
    std::vector<AnalysisMode> convert(const std::vector<TString>& analysisModes);

    /// Convert a vector of AnalysisModes from string to enum
    std::vector<AnalysisMode> convert(const std::vector<std::string>& analysisModes);

    /// Convert a vector of AnalysisModes from enum to string
    std::vector<TString> convert(const std::vector<AnalysisMode>& analysisModes);
}










#endif
