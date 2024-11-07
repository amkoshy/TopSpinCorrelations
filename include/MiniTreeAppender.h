#ifndef MINITREEAPPENDER_H
#define MINITREEAPPENDER_H

#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <map>
#include <cmath>

#include <thread>
#include <mutex>
#include <chrono>
#include <pthread.h>

#include <TROOT.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TRandom.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TMath.h>
#include "TSystem.h"

#include <MiniTreeReader.h>
#include <MiniTree.h>
#include <MiniTreeAnalyzer.h>

#include "TopAnalysis/ZTopUtils/interface/TFModel.h"


#include "Math/GenVector/VectorUtil.h"

#include <MiniTreeHelpers.h>
#include <sampleHelpers.h>
#include "../../common/include/utils.h"
#include "../../common/include/classes.h"
#include "../../common/include/classesFwd.h"


class MiniTreeAppender{
// public:

    TString inputFileName = "";
    TString inputFilePath = "";
    TString outFilePath;

    int year;
    TString yearString;

    Systematic::Systematic systematic;
    Channel::Channel channel = Channel::undefined;

    std::string oldFile_;
    std::string newFile_;
    std::string oldTree_;

    std::chrono::steady_clock::time_point start;
    std::chrono::steady_clock::time_point end;

    uint maxSize = 3;

    ztop::TFModel *dnnModel;
    double dnn_inputs[21];
    std::vector<float> dnn_outputs;
    // float dnn_output;

    void InitDNNModel(TString pathToModel);
    float EvaluateDNN(VLV &jets_, LV &lepton1_, LV &lepton2_, LV &met_);

    TChain* chain;
    TChain* chain_old;

    unsigned long long eventNumber;
    float b_lepton1_pt;
    float b_lepton1_phi;
    float b_lepton1_eta;
    float b_lepton1_m;
    float b_lepton2_pt;
    float b_lepton2_eta;
    float b_lepton2_phi;
    float b_lepton2_m;
    float b_met_pt;
    float b_met_phi;
    std::vector<float> *b_jets_pt=0;
    std::vector<float> *b_jets_eta=0;
    std::vector<float> *b_jets_phi=0;
    std::vector<float> *b_jets_m=0;
    // needed for DNN
    std::vector<float> jets_px;
    std::vector<float> jets_py;
    std::vector<float> jets_pz;
    std::vector<float> jets_e;
    float lepton1_px;
    float lepton1_py;
    float lepton1_pz;
    float lepton1_e;
    float lepton2_px;
    float lepton2_py;
    float lepton2_pz;
    float lepton2_e;
    float met_e;

public:
    MiniTreeAppender();
    MiniTreeAppender(const TString &year_string, const TString &filename_string, const TString &systematic_string, const TString &channel_string);
    void ProgressBar(const int &progress);

    void Append();
};




#endif
