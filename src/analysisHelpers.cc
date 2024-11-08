#include <vector>
#include <iostream>
#include <cstdlib>

#include <TString.h>

#include "analysisHelpers.h"






AnalysisMode::AnalysisMode AnalysisMode::convert(const TString& analysisMode)
{
    if(analysisMode == "cp") return cp;
    if(analysisMode == "gen") return gen;
    if(analysisMode == "dda") return dda;
    if(analysisMode == "btop") return btop;
    if(analysisMode == "ddaTree") return ddaTree;
    if(analysisMode == "plainTree") return plainTree;
    if(analysisMode == "miniTree") return miniTree;
    if(analysisMode == "btopTree") return btopTree;
    if(analysisMode == "kinReco") return kinReco;
    if(analysisMode == "looseKinReco") return looseKinReco;
    if(analysisMode == "kinRecoQualityStudies") return kinRecoQualityStudies;
    if(analysisMode == "dijet") return dijet;
    if(analysisMode == "charge") return charge;
    if(analysisMode == "match") return match;
    if(analysisMode == "playg") return playg;
    if(analysisMode == "weight") return weight;
    if(analysisMode == "genEvent") return genEvent;
    if(analysisMode == "mvaTopP") return mvaTopP;
    if(analysisMode == "mvaTopA") return mvaTopA;
    if(analysisMode == "mvaEventP") return mvaEventP;
    if(analysisMode == "mvaEventA") return mvaEventA;
    if(analysisMode == "mvaChargeP") return mvaChargeP;
    if(analysisMode == "mvaChargeA") return mvaChargeA;
    if(analysisMode == "spinCorr") return spinCorr;
    std::cout<<"Warning! The following analysis mode conversion is not implemented: "<<analysisMode<<std::endl;
    return undefined;
}



TString AnalysisMode::convert(const AnalysisMode& analysisMode)
{
    if(analysisMode == cp) return "cp";
    if(analysisMode == gen) return "gen";
    if(analysisMode == dda) return "dda";
    if(analysisMode == btop) return "btop";
    if(analysisMode == ddaTree) return "ddaTree";
    if(analysisMode == plainTree) return "plainTree";
    if(analysisMode == miniTree) return "miniTree";
    if(analysisMode == btopTree) return "btopTree";
    if(analysisMode == kinReco) return "kinReco";
    if(analysisMode == looseKinReco) return "looseKinReco";
    if(analysisMode == kinRecoQualityStudies) return "kinRecoQualityStudies";
    if(analysisMode == dijet) return "dijet";
    if(analysisMode == charge) return "charge";
    if(analysisMode == match) return "match";
    if(analysisMode == playg) return "playg";
    if(analysisMode == weight) return "weight";
    if(analysisMode == genEvent) return "genEvent";
    if(analysisMode == mvaTopP) return "mvaTopP";
    if(analysisMode == mvaTopA) return "mvaTopA";
    if(analysisMode == mvaEventP) return "mvaEventP";
    if(analysisMode == mvaEventA) return "mvaEventA";
    if(analysisMode == mvaChargeP) return "mvaChargeP";
    if(analysisMode == mvaChargeA) return "mvaChargeA";
    if(analysisMode == spinCorr) return "spinCorr";
    if(analysisMode == undefined) return "";
    std::cout<<"Error! Analysis mode conversion is not implemented, break...\n"<<std::endl;
    exit(96);
}



std::vector<AnalysisMode::AnalysisMode> AnalysisMode::convert(const std::vector<TString>& analysisModes)
{
    std::vector<AnalysisMode> v_analysisMode;
    for(auto analysisMode : analysisModes) v_analysisMode.push_back(convert(analysisMode));
    return v_analysisMode;
}



std::vector<AnalysisMode::AnalysisMode> AnalysisMode::convert(const std::vector<std::string>& analysisModes)
{
    std::vector<AnalysisMode> v_analysisMode;
    for(auto analysisMode : analysisModes) v_analysisMode.push_back(convert(analysisMode));
    return v_analysisMode;
}



std::vector<TString> AnalysisMode::convert(const std::vector<AnalysisMode>& analysisModes)
{
    std::vector<TString> v_analysisMode;
    for(auto analysisMode : analysisModes) v_analysisMode.push_back(convert(analysisMode));
    return v_analysisMode;
}
