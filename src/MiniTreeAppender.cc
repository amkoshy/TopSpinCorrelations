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

#include "TopAnalysis/ZTopUtils/interface/TFModel.h"
#include <MiniTreeAppender.h>

MiniTreeAppender::MiniTreeAppender() {}
MiniTreeAppender::MiniTreeAppender(const TString &year_string, const TString &filename_string, const TString &systematic_string, const TString &channel_string)

{
    systematic = Systematic::undefinedSystematic();
    systematic = Systematic::Systematic(systematic_string);

    Systematic::Systematic nominalSyst = Systematic::Systematic("Nominal");

    channel= Channel::undefined;
    channel = Channel::convert(channel_string);

    TString inpath;
    TString outpath = assignFolder("mergedMiniTree_"+year_string, channel, systematic.name());
    TString inputFilePath_ = outpath+Channel::convert(channel)+"_"+filename_string;
    TString outputFilePath_ = outpath+Channel::convert(channel)+"_appended_"+filename_string;
    InitDNNModel(TString("/nfs/dust/cms/user/sewuchte/TopRhoNetwork/rhoOutput/rhoRegModel_optim_2017.pb"));
    std::cout << "MiniTreeAppender(): Input file for analysis: " + inputFilePath_<< std::endl;
    oldFile_ = inputFilePath_;
    oldTree_ = "miniTree";
    newFile_ = outputFilePath_;
}

void MiniTreeAppender::ProgressBar(const int &progress){
    std::string progressBar = "[";
    for(int i = 0; i <= progress; i++){
        if(i%1 == 0) progressBar += "#";
    }
    for(int i = 0; i < 100 - progress; i++){
        if(i%1 == 0) progressBar += " ";
    }
    progressBar = progressBar + "] " + std::to_string(progress) + "% of Events processed";
    std::cout << "\r" << progressBar << std::flush;
    if(progress == 100) std::cout << std::endl;
}


void MiniTreeAppender::InitDNNModel(TString pathToModel){
    dnnModel = new ztop::TFModel( (std::string) pathToModel, 21, "input_1",1,  "dense_2/Sigmoid");
}


void MiniTreeAppender::Append(){

    TString branchName("DNN_rho");

    TFile* oldF = TFile::Open(oldFile_.c_str(), "update");
    TTree* oldT = (TTree*)oldF->Get(oldTree_.c_str());
    //Disables wished branches if they are already existing in old tree
    // if(oldT->GetListOfBranches()->Contains(branchName)){
        // oldT->SetBranchStatus(branchName, 0);
        // oldT->SetBranchStatus("*", 1);
    // }

    float dnn_output = -999.;

    TBranch *bDNN = oldT->Branch(branchName,&dnn_output);

    LV tempJet;
    LV lepton1;
    LV lepton2;
    VLV jets;
    LV met;

    oldT->SetBranchAddress("lepton1_pt",  &b_lepton1_pt);
    oldT->SetBranchAddress("lepton1_eta", &b_lepton1_eta);
    oldT->SetBranchAddress("lepton1_phi", &b_lepton1_phi);
    oldT->SetBranchAddress("lepton1_m",   &b_lepton1_m);
    oldT->SetBranchAddress("lepton2_pt",  &b_lepton2_pt);
    oldT->SetBranchAddress("lepton2_eta", &b_lepton2_eta);
    oldT->SetBranchAddress("lepton2_phi", &b_lepton2_phi);
    oldT->SetBranchAddress("lepton2_m",   &b_lepton2_m);
    oldT->SetBranchAddress("met_pt",   &b_met_pt);
    oldT->SetBranchAddress("met_phi",   &b_met_phi);

    oldT->SetBranchAddress("jets_pt",   &b_jets_pt);
    oldT->SetBranchAddress("jets_eta",   &b_jets_eta);
    oldT->SetBranchAddress("jets_phi",   &b_jets_phi);
    oldT->SetBranchAddress("jets_m",   &b_jets_m);

    Long64_t nentries = oldT->GetEntries();
    ProgressBar(0.);

    int processed = 0;
    for (Long64_t i=0;i<nentries;i++){
        lepton1.SetCoordinates(0., 0., 0., 0.);
        lepton2.SetCoordinates(0., 0., 0., 0.);
        met.SetCoordinates(0., 0., 0., 0.);
        jets.clear();

        oldT->GetEntry(i);

        for(unsigned int i = 0; i < b_jets_pt->size(); i++){
            tempJet.SetCoordinates(b_jets_pt->at(i),b_jets_eta->at(i),b_jets_phi->at(i),b_jets_m->at(i));
            jets.push_back(tempJet);
        }
        lepton1.SetCoordinates(b_lepton1_pt, b_lepton1_eta, b_lepton1_phi, b_lepton1_m);
        lepton2.SetCoordinates(b_lepton2_pt, b_lepton2_eta, b_lepton2_phi, b_lepton2_m);
        met.SetCoordinates(b_met_pt, 0., b_met_phi, 0.);

        dnn_output = EvaluateDNN(jets, lepton1, lepton2, met);

        bDNN->Fill();

        processed++;
        if(processed % 1000 == 0) {
            int progress = 100*(float)processed/oldT->GetEntries();
            ProgressBar(progress);
        }
    }

    ProgressBar(100);

    oldT->Write(0, TObject::kOverwrite, 0);
    delete oldF;
}

float MiniTreeAppender::EvaluateDNN(VLV &jets_, LV &lepton1_, LV &lepton2_, LV &met_)
{
    uint tempMaxsize = maxSize;
    dnn_inputs[0]= 0.; dnn_inputs[1]= 0.; dnn_inputs[2]= 0.; dnn_inputs[3]= 0.;
    dnn_inputs[4]= 0.; dnn_inputs[5]= 0.; dnn_inputs[6]= 0.; dnn_inputs[7]= 0.;
    dnn_inputs[8]= 0.; dnn_inputs[9]= 0.; dnn_inputs[10]= 0.; dnn_inputs[11]= 0.;
    dnn_inputs[12]= lepton1_.Px(); dnn_inputs[13]= lepton1_.Py(); dnn_inputs[14]= lepton1_.Pz(); dnn_inputs[15]= lepton1_.E();
    dnn_inputs[16]= lepton2_.Px(); dnn_inputs[17]= lepton2_.Py(); dnn_inputs[18]= lepton2_.Pz(); dnn_inputs[19]= lepton2_.E();
    dnn_inputs[20]= met_.E();

    if (jets_.size()<=tempMaxsize){
        tempMaxsize = jets_.size();
    }
    for(uint i=0; i < tempMaxsize; i++){
        dnn_inputs[i*4+0] = jets_.at(i).Px();
        dnn_inputs[i*4+1] = jets_.at(i).Py();
        dnn_inputs[i*4+2] = jets_.at(i).Pz();
        dnn_inputs[i*4+3] = jets_.at(i).E();
    }
    dnn_outputs = dnnModel->evaluate(dnn_inputs);
    return dnn_outputs.at(0);
}
