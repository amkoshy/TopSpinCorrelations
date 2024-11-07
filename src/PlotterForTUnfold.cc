// 
// This file taken from Ievgen Korol's code (/nfs/dust/cms/user/korol/CMSSW_5_3_27) on 14.01.2017 as used for 2D analysis TOP-14-013
// with purpose to run TUnfold. Modifed to work for 1D.
//
// Modified to work in 1D: "fake" first dimension with 1 bin has to be provided, see example in Plotter13TeV::CalcDiffXSec()
//
// Oleksandr Zenaiev (oleksandr.zenaiev@desy.de), 14.01.2017
//

#include <fstream>
#include <iostream>
// #include <cstdio>
#include <sstream>
// #include <cmath>
// #include <iomanip>

#include <TCanvas.h>
#include <TLegend.h>
// #include <TExec.h>
#include <TStyle.h>
#include <TSystem.h>
// #include <TMath.h>
// #include <TROOT.h>
#include <THStack.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMatrixD.h>
#include <TMatrixTSym.h>
#include <TDecompChol.h>

#include "utils.h"
#include "PlotterForTUnfold.h"
#include "../../common/include/utils.h"
#include "../../common/include/plotterUtils.h"
#include "UsefulTools.h"
#include "Output.h"

#include "TUnfold.h"
#include <TLine.h>
#include <TMath.h>
#include <TPaveText.h>
#include "TUnfoldBinning.h"
#include "TUnfoldDensity.h"





PlotterForTUnfold::PlotterForTUnfold(const Samples& samples, const double lumi, const double lumiUnc):

samples_(samples),
lumi_(lumi),
lumi0_(lumi),
lumiUnc_(lumiUnc),
topxsec_(-999),
nDataInc_(0),
nBgrInc_(0),
nRecoTtAllInc_(0),
nRecoSigInc_(0),
nRecoAllMC_(0),
nRecoTtOtherInc_(0),
nRecoSingleTop_(0),
nRecoWPlusJets_(0),
nRecoZee_(0),
nRecoZtautau_(0),
nRecoTtPlusZWG_(0),
nRecoDiboson_(0),
nRecoQCD_(0),

v_plotName_(std::vector<TString>(0)),
v_plotTitle_(std::vector<TString>(0)),
v_plotUnits_(std::vector<TString>(0)),
v_plotYunitscp_(std::vector<TString>(0)),
v_plotYlog_(0),

// OZ 17.01.2017
//plotsFolder_("Plots/"),
plotsFolder_("TUnfold/"),
tauFolder_(""),
ctFolder_(""),
systematicName_(""),

v_histOnReco_(0),
v_histOnGen_(0),

//Control plots//
v_cpNBins_(0),
v_R1_(std::vector<Double_t>(0)),
v_R2_(std::vector<Double_t>(0)),
v_ymax_(std::vector<Double_t>(0)),
v_leg2_(std::vector<Int_t>(0)),
v_eban_(std::vector<Int_t>(0)),
vv_SampleHist_(std::vector<std::vector<TH1D* > >(0)),
v_histRecoAllBins_(0),
vvv_SampleHist_(std::vector<std::vector<std::vector<TH1D* > > >(0)),
vv_UnderflowSampleHist_(std::vector<std::vector<TH1D* > >(0)),
vv_OverflowSampleHist_(std::vector<std::vector<TH1D* > >(0)),

//Unfolding Options//
regMethod_("rho"),
isTauTest_(false),
tau_(-999),
logTau_(-999),
isCT_(false),

var_(-999),
isVar_(false),

//Unfolding//
v_coarseBins_(std::vector<std::vector<Double_t> >(0)),
v_fineBins_(std::vector<std::vector<Double_t> >(0)),
v_uReco_(std::vector<int >(0)),
v_oReco_(std::vector<int >(0)),
v_uTrue_(std::vector<int >(0)),
v_oTrue_(std::vector<int >(0)),
detectorBinning_(0),
generatorBinning_(0),
detectorDistribution_(0),
generatorDistribution_(0),
histMigration_(0),
v_covX_(std::vector<double >(0)),
vX_(std::vector<double >(0)),
v_covXnorm_(std::vector<double >(0)),
vXnorm_(std::vector<double >(0)),
vXnormMC_(std::vector<double >(0)),
vXmc_(std::vector<double >(0)),

v_covY_(std::vector<double >(0)),
v_dY_(std::vector<double >(0)),
v_covYnorm_(std::vector<double >(0)),
v_dYnorm_(std::vector<double >(0)),
vYmc_(std::vector<double >(0)),
vYnormMC_(std::vector<double >(0)),

//.................................................................................................................................
histBgrUo_(0),
histBgrUoScaled_(0),
histBgr_(0),
histBgrScaled_(0),
histData_(0),
histDataScaled_(0),

histFsigReco_(0),
histTtAllReco_(0),
histTtSigReco_(0),
histAllReco_(0),
histAllRecoScaled_(0),
//.......................................................................................................................................
unfoldedData_(0),
unfoldedDataNorm_(0),
histUnfGenBin_(0),
histTauLog_(0),
histTauLogRhoScan_(0),
histTauLogLCurve_(0),

rhoIJtotal_(0),


//Closure test//
histSR_(0),
histSG_(0),
v_histGenRewUp_(0),
v_histGenRewDown_(0),
v_histRecoRewUp_(0),
v_histRecoRewDown_(0),

histGenPtRewUp_(0),
histGenPtRewDown_(0),
histGenRapidityRewUp_(0),
histGenRapidityRewDown_(0),
histRecoPtRewUp_(0),
histRecoPtRewDown_(0),
histRecoRapidityRewUp_(0),
histRecoRapidityRewDown_(0),



histGen_(0),
histGenNorm_(0),
histPurityAllBins_(0),
histStabilityAllBins_(0),
vectorPurityAllBins_(std::vector<double>(0)),
vectorStabilityAllBins_(std::vector<double>(0)),
histRecoGenAllBins_(0),
histEffAllBins_(0),

// Tree branches //
entry_(-999),
entry0_(-999),
eventWeight_(-999), 
trueLevelWeight_(-999),
trueLevelWeight0_(-999),
isSignal_(0),
branchVals_(std::vector<Float_t>(0)),
branchValsGen_(std::vector<Float_t>(0)),
branchValsGen0_(std::vector<Float_t>(0)),

isPDF_(0)

{
    r3_ = new TRandom3();
    r3_->SetSeed(3013);
    
    //FIXME: make it in better way
    v_BR_.push_back(0.01166);//ee
    v_BR_.push_back(0.02332);//emu
    v_BR_.push_back(0.01166);//mumu
    v_BR_.push_back(0.04666);//combined
     
     // switch on histogram errors
    TH1::SetDefaultSumw2();
     
}



void PlotterForTUnfold::setOptions(const std::vector<TString> v_plotName)
{
    
    //std::cout<<"[PlotterForTUnfold]:--- Beginning of plot options settings\n\n";
    
    //Clearing previous declarations  
    this->clearMemory();
    
    v_plotName_ = v_plotName;
    
    branchVals_ = std::vector<Float_t>(2,-999);
    branchValsGen_ = std::vector<Float_t>(2,-999);
    branchValsGen0_ = std::vector<Float_t>(2,-999);
    
    // reading Control Cards
    for(auto plotName : v_plotName_){
        ///check before start
        ///const std::string ccFileName(common::CMSSW_BASE() + "/src/TopAnalysis/Configuration/analysis/diLeptonic/ControlCards/" + plotName + ".txt");
        // OZ 15.01.2017 renamed input steering dir
        //const std::string ccFileName("./ControlCards/" + plotName + ".txt");
        const std::string ccFileName("./ControlCardsTUnfold/" + plotName + ".txt");
        std::ifstream ccFileStream(ccFileName.data(), std::ifstream::in);
        if (!ccFileStream.good()) {
          std::cerr<<"Error in PlotterForTUnfold! Cannot find file with name: "<< ccFileName <<"\n...break\n"<<std::endl;
          exit(12);
        }

       // Loop over all lines in ccFile
       bool isInBlock = false;
       while(ccFileStream.good()){
        // Read cc-File
            std::string line;
            getline(ccFileStream, line);
            line.erase(0, line.find_first_not_of(" \t"));
            if (line.size() == 0 || line[0] == '#') continue;

            // Loop over words in cc-File line and fill vWord
             std::vector<TString> vWord;
             std::string word;
             for (std::stringstream ss(line); ss >> word; ){
                vWord.push_back(word);
             }
             
             if(vWord.at(0) == "cpnbins")
             {
                     v_cpNBins_.push_back((vWord.at(1)).Atoi());
                     v_R1_.push_back((vWord.at(2)).Atof());
                     v_R2_.push_back((vWord.at(3)).Atof());
             }
             if(vWord.at(0) == "ymax")
             {
                     v_ymax_.push_back((vWord.at(1)).Atoi());
             }
             
             if(vWord.at(0) == "leg2")
             {
                     v_leg2_.push_back((vWord.at(1)).Atoi());
             }
             
             if(vWord.at(0) == "eban")
             {
                     v_eban_.push_back((vWord.at(1)).Atoi());
             }
             
             if(vWord.at(0) == "title")v_plotTitle_.push_back(vWord.at(1));
             
             if(vWord.at(0) == "units" ){
                if(vWord.size()>1) v_plotUnits_.push_back(vWord.at(1)); 
                if(vWord.size()<2) v_plotUnits_.push_back("");
             }
             
             if(vWord.at(0) == "yunitscp" ){
                if(vWord.size()>1){ 
                    vWord.at(1).ReplaceAll("~"," ");
                    v_plotYunitscp_.push_back(vWord.at(1));
                } 
                if(vWord.size()<2){ 
                    // OZ 15.01.2017 renamed input steering dir
                    //std::cout << "Warning: yunitscp should be define in ControlCards/" << plotName << ".txt" << std::endl;
                    std::cout << "Warning: yunitscp should be define in ControlCardsTUnfold/" << plotName << ".txt" << std::endl;
                    v_plotYunitscp_.push_back("");
                }
             }
             
             if(vWord.at(0) == "ylog")
             {
                     v_plotYlog_.push_back((vWord.at(1)).Atoi());
             }
             
             //Pass to dimensional block
             if( (int)vWord.size()==1 && vWord.at(0) == "2d" ){
                 if(isInBlock) isInBlock = false;
                 if(!isInBlock) isInBlock = true;
            }
             
             // Reading only from dimensional block we are interesting
             if(isInBlock){
                 std::vector<Double_t > binsVector;
                 
                if(vWord.at(0) == "coarse"){
                     for(auto word : vWord) binsVector.push_back(word.Atof());
                     binsVector.erase(binsVector.begin(), binsVector.begin() + 1);
                     v_coarseBins_.push_back(binsVector);
                     // OZ 17.01.2017
                     if(v_coarseBins_.size() == 2)
                     {
                       printf("PlotterForTUnfold::setOptions(): coarse [");
                       for(unsigned int b = 0; b < v_coarseBins_[0].size(); b++)
                         printf(" %.3f", v_coarseBins_[0][b]);
                       printf("] [");
                       for(unsigned int b = 0; b < v_coarseBins_[1].size(); b++)
                         printf(" %.3f", v_coarseBins_[1][b]);
                       printf("]\n");
                     }
                }
                if(vWord.at(0) == "fine"){
                     for(auto word : vWord) binsVector.push_back(word.Atof());
                     binsVector.erase(binsVector.begin(), binsVector.begin() + 1);
                     v_fineBins_.push_back(binsVector);
                     // OZ 17.01.2017
                     if(v_fineBins_.size() == 2)
                     {
                       printf("PlotterForTUnfold::setOptions(): fine [");
                       for(unsigned int b = 0; b < v_fineBins_[0].size(); b++)
                         printf(" %.3f", v_fineBins_[0][b]);
                       printf("] [");
                       for(unsigned int b = 0; b < v_fineBins_[1].size(); b++)
                         printf(" %.3f", v_fineBins_[1][b]);
                       printf("]\n");
                     }
                }
                if(vWord.at(0) == "uoTrue")
                {
                    v_uTrue_.push_back(vWord.at(1).Atoi());
                    v_oTrue_.push_back(vWord.at(2).Atoi());
                }
                if(vWord.at(0) == "uoReco")
                {
                    v_uReco_.push_back(vWord.at(1).Atoi());
                    v_oReco_.push_back(vWord.at(2).Atoi());
                }
                 
             }//read from dimensional block
             else continue;
             
       }//stream line
       
    }//dimension name

    //Set binning scheme
    //Definition of TUnfoldBinning no detector and generator level
    detectorBinning_ = new TUnfoldBinning("detector");
    detectorDistribution_ = detectorBinning_->AddBinning("detectordistribution");
    for(int i=0;i<2;++i)detectorDistribution_->AddAxis(v_plotName_.at(i),(int)v_fineBins_.at(i).size()-1,v_fineBins_.at(i).data(),v_uReco_.at(i), v_oReco_.at(i));
    
    generatorBinning_ = new TUnfoldBinning("generator");
    generatorDistribution_ = generatorBinning_->AddBinning("signal");
    for(int i=0;i<2;++i)generatorDistribution_->AddAxis("gen_"+v_plotName_.at(i),(int)v_coarseBins_.at(i).size()-1, v_coarseBins_.at(i).data() , v_uTrue_.at(i), v_oTrue_.at(i));
    
    //std::cout<<"[PlotterForTUnfold]:--- Finishing of plot options settings\n\n";
}



void PlotterForTUnfold::clearMemory()
{
    v_plotName_.clear();
    v_plotTitle_.clear();
    v_plotUnits_.clear();
    v_plotYunitscp_.clear();
    v_plotYlog_.clear();
    branchVals_.clear();
    branchValsGen_.clear();
    branchValsGen0_.clear();
    v_coarseBins_.clear();
    v_fineBins_.clear();
    v_uReco_.clear();
    v_oReco_.clear();
    v_uTrue_.clear();
    v_oTrue_.clear();
    v_cpNBins_.clear();
    v_R1_.clear();
    v_R2_.clear();
    v_ymax_.clear();
    v_leg2_.clear();
    
    if(detectorDistribution_)detectorDistribution_->Delete();
    if(detectorBinning_)detectorBinning_->Delete();
    if(generatorDistribution_)generatorDistribution_->Delete();
    if(generatorBinning_)generatorBinning_->Delete();
    
    
}

void PlotterForTUnfold::clearMemoryPerSystematicChannel()
{
    
    for(auto v_p : vv_SampleHist_)
     for(auto p : v_p)
         delete p;
    vv_SampleHist_.clear();

    for(auto v_p : vv_UnderflowSampleHist_)
     for(auto p : v_p)
         delete p;
    vv_UnderflowSampleHist_.clear();
    
    for(auto v_p : vv_OverflowSampleHist_)
     for(auto p : v_p)
         delete p;
    vv_OverflowSampleHist_.clear();
    
    for(auto vv_p : vvv_SampleHist_)
     for(auto v_p : vv_p)
      for(auto p : v_p)
         delete p;
    vvv_SampleHist_.clear();
    
    if(histMigration_)histMigration_->Delete();
    v_covX_.clear();
    vX_.clear();
    v_covXnorm_.clear();
    vXnorm_.clear();
    vXnormMC_.clear();
    vXmc_.clear();
    
    v_covY_.clear();
    v_dY_.clear();
    v_covYnorm_.clear();
    v_dYnorm_.clear();
    vYmc_.clear();
    vYnormMC_.clear();
    
    
    if(histBgrUo_)histBgrUo_->Delete();
    if(histBgrUoScaled_)histBgrUoScaled_->Delete();
    if(histBgr_)histBgr_->Delete();
    if(histBgrScaled_)histBgrScaled_->Delete();
    if(histData_)histData_->Delete();
    if(histDataScaled_)histDataScaled_->Delete();
    
    if(histFsigReco_)histFsigReco_->Delete();
    if(histTtAllReco_)histTtAllReco_->Delete();
    if(histTtSigReco_)histTtSigReco_->Delete();
    if(histAllReco_)histAllReco_->Delete();
    if(histAllRecoScaled_)histAllRecoScaled_->Delete();
    
    if(unfoldedData_)unfoldedData_->Delete();
    if(unfoldedDataNorm_)unfoldedDataNorm_->Delete();
    if(histUnfGenBin_)histUnfGenBin_->Delete();
    if(histTauLog_)histTauLog_->Delete();
    if(histTauLogRhoScan_)histTauLogRhoScan_->Delete();
    if(histTauLogLCurve_)histTauLogLCurve_->Delete();
    
    if(rhoIJtotal_)rhoIJtotal_->Delete();
    
    //Closurer test//
    if(histSR_)histSR_->Delete();
    if(histSG_)histSG_->Delete();
    
    for(auto p : v_histGenRewUp_)
        delete p;
    v_histGenRewUp_.clear();
    
    for(auto p : v_histGenRewDown_)
        delete p;
    v_histGenRewDown_.clear();
    
    for(auto p : v_histRecoRewUp_)
        delete p;
    v_histRecoRewUp_.clear();
    
    for(auto p : v_histRecoRewDown_)
        delete p;
    v_histRecoRewDown_.clear();
    
    if(histGenPtRewUp_)histGenPtRewUp_->Delete();
    if(histGenPtRewDown_)histGenPtRewDown_->Delete();
    if(histGenRapidityRewUp_)histGenRapidityRewUp_->Delete();
    if(histGenRapidityRewDown_)histGenRapidityRewDown_->Delete();
    if(histRecoPtRewUp_)histRecoPtRewUp_->Delete();
    if(histRecoPtRewDown_)histRecoPtRewDown_->Delete();
    if(histRecoRapidityRewUp_)histRecoRapidityRewUp_->Delete();
    if(histRecoRapidityRewDown_)histRecoRapidityRewDown_->Delete();
    
    for(auto p : v_histRecoAllBins_)
        delete p;
    v_histRecoAllBins_.clear();
    
    for(auto p : v_histOnReco_)
        delete p;
    v_histOnReco_.clear();
    
    for(auto p : v_histOnGen_)
        delete p;
    v_histOnGen_.clear();
    
    if(histGen_)histGen_->Delete();
    if(histGenNorm_)histGenNorm_->Delete();
    
    if(histEffAllBins_)histEffAllBins_->Delete();
    if(histPurityAllBins_)histPurityAllBins_->Delete();
    if(histStabilityAllBins_)histStabilityAllBins_->Delete();
    if(histRecoGenAllBins_)histRecoGenAllBins_->Delete();
    
    vectorPurityAllBins_.clear();
    vectorStabilityAllBins_.clear();
    
    plotsFolder_ = "";
    tauFolder_ = "";
    ctFolder_ = "";
    systematicName_ = "";
    
    r3_->SetSeed(3013);
}



void PlotterForTUnfold::prepareHistograms(const std::vector<Sample>& v_sample)
{
    for(int ind=0;ind<2;++ind){
        std::vector<TH1D* > v_SampleHist;
        std::vector<TH1D* > v_UnderflowSampleHist;
        std::vector<TH1D* > v_OverflowSampleHist;
        
        // OZ to eliminate unused var warning [preColor]
        //int preColor = -1;
        for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
            const auto& sample(v_sample.at(iSample));
            TH1D* sampleHist = new TH1D(v_plotName_.at(ind)+"_cp" + std::to_string(ind) + std::to_string((int)iSample)  ,v_plotTitle_.at(ind),v_cpNBins_.at(ind),v_R1_.at(ind),v_R2_.at(ind));//FIXME: remove this: + std::to_string(ind) + std::to_string((int)iSample)
            sampleHist->SetFillColor(sample.color());
            sampleHist->SetLineColor(sample.color());
            v_SampleHist.push_back(sampleHist);
            v_UnderflowSampleHist.push_back((TH1D*)(sampleHist->Clone(sampleHist->GetName() + TString("u"))));
            v_OverflowSampleHist.push_back((TH1D*)(sampleHist->Clone(sampleHist->GetName() + TString("o"))));
            
            // OZ to eliminate unused var warning [preColor]
            //if(iSample!=0)preColor=sample.color();
            
        }
        vv_SampleHist_.push_back(v_SampleHist);
        vv_UnderflowSampleHist_.push_back(v_UnderflowSampleHist);
        vv_OverflowSampleHist_.push_back(v_OverflowSampleHist);
    }

    for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
        const auto& sample(v_sample.at(iSample));
        TH1* sampleHist =generatorBinning_->CreateHistogram("histRecoAllBins"+(TString)std::to_string((int)iSample).data());
        sampleHist->SetFillColor(sample.color());
        sampleHist->SetLineColor(sample.color());
        v_histRecoAllBins_.push_back((TH1*)sampleHist->Clone());
        sampleHist->Delete();
    }
    
    //setinng vvv_SampleHist_
    for(int ind=0;ind<2;++ind){
        std::vector<std::vector<TH1D* > > vv_SampleHist;
        for(int jnd=0;jnd<2;++jnd){
            if(ind==jnd)continue;
            const auto& v_bins = v_coarseBins_.at(jnd);
            for(int ibin=0;ibin<(int)v_bins.size()-1;++ibin){
                std::vector<TH1D* > v_SampleHist;
                for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
                    const auto& sample(v_sample.at(iSample));
                    TH1D* sampleHist = new TH1D(v_plotName_.at(ind)+"_cpBin" + std::to_string(ind) + std::to_string(ibin) + std::to_string((int)iSample)  ,v_plotTitle_.at(ind),
                                                (int)v_coarseBins_.at(ind).size()-1,v_coarseBins_.at(ind).data());
                    sampleHist->SetFillColor(sample.color());
                    //sampleHist->SetLineColor(sample.color());
                    v_SampleHist.push_back(sampleHist);
                }
                vv_SampleHist.push_back(v_SampleHist);
            }
        }
        vvv_SampleHist_.push_back(vv_SampleHist);
    }
    
    ///check before start
// //     std::vector<int> colorNum={9,8,6,46};
// //     for(int i =0; i<(int)v_coarseBins_.at(0).size()-1 ;i++){
// //         TH1D* hist = new TH1D("name","tittle",(int)v_coarseBins_.at(0).size()-1,v_coarseBins_.at(0).data());
// //         TH1D* hist1 = new TH1D("name1","tittle1",20,0,500);
// //         //hist->SetFillColor(40+i*2);
// //         //hist1->SetFillColor(40+i*2);
// //         hist->SetFillColor(colorNum.at(i));
// //         hist1->SetFillColor(colorNum.at(i));
// //         hist->SetFillStyle(3001);
// //         hist1->SetFillStyle(3001);
// //         
// //         v_histOnReco_.push_back((TH1D*)hist1->Clone());
// //         v_histOnGen_.push_back((TH1D*)hist->Clone());
// //     }
    std::cout << "TEST:!!!!! 1" << std::endl;
    //Definition of unfolding histograms
    histMigration_ = TUnfoldBinning::CreateHistogramOfMigrations(generatorBinning_,detectorBinning_,"histMigration");
    
    histBgrUo_ = detectorBinning_->CreateHistogram("histBgrUo");
    histBgrUoScaled_ = detectorBinning_->CreateHistogram("histBgrUoScaled");
    histBgr_ = detectorBinning_->CreateHistogram("histBgr");
    histBgrScaled_ = detectorBinning_->CreateHistogram("histBgrScaled");
    histData_ = detectorBinning_->CreateHistogram("histData");
    histDataScaled_ = detectorBinning_->CreateHistogram("histDataScaled");
    
    histFsigReco_ = detectorBinning_->CreateHistogram("histFsigReco");
    histTtSigReco_ = detectorBinning_->CreateHistogram("histTtSigReco");
    histTtAllReco_ = detectorBinning_->CreateHistogram("histTtAllReco");
    histAllReco_ = detectorBinning_->CreateHistogram("histAllReco");
    histAllRecoScaled_ = detectorBinning_->CreateHistogram("histAllRecoScaled");
    
        histGen_ = generatorBinning_->CreateHistogram("histGen",kTRUE);
        histGenNorm_ = generatorBinning_->CreateHistogram("histGenNorm",kTRUE);
        
        histEffAllBins_ = generatorBinning_->CreateHistogram("histEffAllBins");
        histPurityAllBins_ = generatorBinning_->CreateHistogram("histPurityAllBins");
        histStabilityAllBins_ = generatorBinning_->CreateHistogram("histStabilityAllBins");
        histRecoGenAllBins_ = generatorBinning_->CreateHistogram("histRecoGenAllBins");
    
    //Closure test//
    histSR_ = detectorBinning_->CreateHistogram("histSR");
    histSG_ = generatorBinning_->CreateHistogram("histSG");
    for(int i=1;i<=histSG_->GetNbinsX();i++){
        v_histGenRewUp_.push_back((TH1*)histSG_->Clone());
        v_histGenRewDown_.push_back((TH1*)histSG_->Clone());
        v_histRecoRewUp_.push_back((TH1*)histSR_->Clone());
        v_histRecoRewDown_.push_back((TH1*)histSR_->Clone());
    }
    
    histGenPtRewUp_=(TH1*)histSG_->Clone();
    histGenPtRewDown_=(TH1*)histSG_->Clone();
    histGenRapidityRewUp_=(TH1*)histSG_->Clone();
    histGenRapidityRewDown_=(TH1*)histSG_->Clone();
    histRecoPtRewUp_=(TH1*)histSR_->Clone();
    histRecoPtRewDown_=(TH1*)histSR_->Clone();
    histRecoRapidityRewUp_=(TH1*)histSR_->Clone();
    histRecoRapidityRewDown_=(TH1*)histSR_->Clone();
    
    
    histTauLog_ = new TH1D("histTauLog_","Regularization strength;lg(#tau);N",200,-7,-2);
    histTauLogRhoScan_ = new TH1D("histTauLogRhoScan_","#rho scan;lg(#tau);N",200,-6,-2);
    histTauLogLCurve_ = new TH1D("histTauLogLCurve_","L curve;lg(#tau);N",200,-6,-2);
    
    std::cout << "TEST:!!!!! 2" << std::endl;
    
}



void PlotterForTUnfold::setRegMethod(TString method){
   regMethod_ = method;
}



void PlotterForTUnfold::setTauTest(bool isTauTest){
    isTauTest_=isTauTest;
}


void PlotterForTUnfold::setCT(bool isCT){
    isCT_=isCT;
}

void PlotterForTUnfold::setVar(double var){
    isVar_=true;
    var_=var;
}


void PlotterForTUnfold::producePlots()
{
    //std::cout<<"[PlotterForTUnfold]:--- Beginning of plot production\n\n";
    
    // Access correction factors
    const SystematicChannelFactors globalWeights = samples_.globalWeights("_step8").first;
    
    // Loop over all channels and systematics and produce plots
    const SystematicChannelSamples& m_systematicChannelSample(samples_.getSystematicChannelSamples());
    for(const auto& systematicChannelSamples : m_systematicChannelSample){                     //systematic loop
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        if((isTauTest_||isCT_)&&systematic.type()!=Systematic::nominal)continue;
        if(systematic.name().Contains("PDF"))isPDF_=1;
        else isPDF_=0;
        
        lumi_ = lumi0_;
        if(systematic.type() == Systematic::lumi){
            if(systematic.variation() == Systematic::up) lumi_ = lumi_*(1 + lumiUnc_);
            else if (systematic.variation() == Systematic::down)  lumi_ = lumi_*(1 - lumiUnc_);
        }
        for(const auto& channelSample : systematicChannelSamples.second){                     //channel loop
            this->clearMemoryPerSystematicChannel();
            const Channel::Channel& channel(channelSample.first);
            const std::vector<Sample>& v_sample(channelSample.second);
            const std::vector<double>& v_weight(globalWeights.at(systematic).at(channel));
            
            systematicName_ = systematic.name();
            
            //Set up folders:
            plotsFolder_ = common::assignFolder("Plots",channel,systematic,v_plotName_.at(0)+"-"+v_plotName_.at(1));
            tauFolder_ = "./Plots/preunfold/" + Channel::convert(channel) + "/" + v_plotName_.at(0)+"-"+v_plotName_.at(1) +"/";
            TString tauFile="";
            if(systematic.name().Contains("PDF")) tauFile = tauFolder_ + "optimalTau_PDF_0_CENTRAL.txt";
            else if(systematic.name().Contains("MCATNLO"))tauFile = tauFolder_ + "optimalTau_MCATNLO.txt";
            else if(systematic.name().Contains("POWHEG"))tauFile = tauFolder_ + "optimalTau_POWHEG.txt";
            else if(systematic.name().Contains("POWHEGHERWIG"))tauFile = tauFolder_ + "optimalTau_POWHEGHERWIG.txt"; 
            else tauFile = tauFolder_ + "optimalTau_Nominal.txt";
            gSystem->mkdir(tauFolder_,true);
            ctFolder_ = tauFolder_ + "/ct/";
            gSystem->mkdir(ctFolder_,true);
            
            
            //read nominal tau from file
            if(systematicName_!="Nominal"||isTauTest_||isCT_||isVar_){ /// check before start
                std::vector<double> tauTemp;
                tauFile = plotsFolder_;
                ///check before start
                tauFile = tauFile.ReplaceAll(systematicName_,"Nominal") + "/optimalTau.txt";
                ///tauFile = "./optimalTau.txt";
                std::cout << "try to read tau from: " << tauFile << std::endl;
                utils::readLineToVector(tauFile,"tau",tauTemp);
                tau_ = tauTemp.at(0);
                std::cout << "try to read tau status: successful !!! " << std::endl;
            }
                                    std::cout << "Test 41: !!!!!!!!!!" << std::endl ;                

            this->prepareHistograms(v_sample);
	                            std::cout << "Test 42: !!!!!!!!!!" << std::endl ;                

            
            TString sampleInfoFolder = common::assignFolder("Plots",channel,systematic) + "/sampleInfo.txt";
            Output sampleInfo("sample");
            
            nRecoSigInc_=0;
            nRecoSingleTop_=0;
            nRecoTtOtherInc_=0;
            nRecoTtAllInc_=0;
            nGenSigInc_= 0 ;
            nDataInc_ = 0;
            nBgrInc_=0;
            nRecoAllMC_=0;
            
            //TFile pdfWeightFile("./selectionRoot/PDF_0_CENTRAL/emu/emu_ttbarsignalplustau.root");
            //TH1D * pdfWeightHist = (TH1D*)pdfWeightFile.Get("pdfWeight");
            //double pdfWeightFromFile = 1.0/pdfWeightHist->GetMean();
            double pdfWeightFromFile = 1.0;
            //pdfWeightHist->Delete();
                                    std::cout << "Test 3: !!!!!!!!!!" << std::endl ;                

            for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){                  //sample loop
                const auto& sample(v_sample.at(iSample));
                
                isSignal_ = 0;
                if(sample.sampleType() == Sample::ttsignal)isSignal_=1;
                
                if(sample.legendEntry() == "QCD Multijet") continue; ///check before start
                //if(sample.legendEntry() == "W+Jets") continue;
                double pdfWeight = pdfWeightFromFile; ///check before start
                if(!((sample.inputFile().Contains("ttbarb")||sample.inputFile().Contains("ttbars"))&&systematicName_.Contains("PDF"))){
                    pdfWeight=1;
                }
                
                std::cout << "Reading file: " <<sample.inputFile() << std::endl;
                TFile* dataFile=new TFile(sample.inputFile().ReplaceAll("selectionRoot","ddaInput"));
                TTree *dataTree0=(TTree *) dataFile->Get("ttBar_treeVariables_step0"); 
                TTree *dataTree8=(TTree *) dataFile->Get("ttBar_treeVariables_step8");
                if(!dataTree0 || !dataTree8) {
                     std::cout<<"could not read 'ttBar_treeVariables_step' tree from " << dataFile->GetName() << "\n";
                 }
                // set branches
                dataTree0->ResetBranchAddresses();
                dataTree0->SetBranchAddress("entry",&entry0_);
                dataTree0->SetBranchAddress("trueLevelWeight",&trueLevelWeight0_);
                for(int i=0;i<2;++i)dataTree0->SetBranchAddress("gen_"+v_plotName_.at(i),&(branchValsGen0_.at(i)));
                
                dataTree8->ResetBranchAddresses();
                dataTree8->SetBranchAddress("entry",&entry_);
                dataTree8->SetBranchAddress("eventWeight",&eventWeight_);
                dataTree8->SetBranchAddress("trueLevelWeight",&trueLevelWeight_);
                dataTree8->SetBranchAddress("isKinReco",&isKinReco_);
                for(int i=0;i<2;++i)dataTree8->SetBranchAddress(v_plotName_.at(i),&(branchVals_.at(i)));
                for(int i=0;i<2;++i)dataTree8->SetBranchAddress("gen_"+v_plotName_.at(i),&(branchValsGen_.at(i)));
                float MttReco = -999;
                ///dataTree8->SetBranchAddress("ttbar_mass",&MttReco);
                
                sampleInfo.add("entries", utils::numToString(dataTree8->GetEntriesFast()));
                sampleInfo.add("weight",utils::numToString(v_weight.at(iSample)));
                sampleInfo.add("file",sample.inputFile());
                sampleInfo.add("name",sample.legendEntry());
                 std::cout << "Test 4: !!!!!!!!!!" << std::endl ;
                
                ///check before start
                double bgrWeight = 1;
                //if(isVar_&&(sample.sampleType() == Sample::dyee || sample.sampleType() == Sample::dymumu || sample.sampleType() == Sample::dytautau))bgrWeight = 1.*var_/100;
                //if(sample.sampleType() != Sample::data && sample.sampleType() != Sample::ttother && 
                       //sample.sampleType() != Sample::dyee && sample.sampleType() != Sample::dymumu && sample.sampleType() != Sample::dytautau && sample.sampleType() != Sample::ttsignal)bgrWeight = 1.*var_/100;
                ///if(sample.sampleType() != Sample::data&&isVar_)bgrWeight = 1.*var_/100;
                //if(isVar_)bgrWeight = 1.*var_/100;
                
                /// **************************************************  Filling Histograms ********************************** ///
                Int_t lastEvent=0;
                for(Int_t ievent=0;ievent<dataTree8->GetEntriesFast();ievent++){               //loop over tree8 events
                    if(dataTree8->GetEntry(ievent)<=0) break;
                    
                    eventWeight_ = eventWeight_*pdfWeight*bgrWeight;
                    trueLevelWeight_ = trueLevelWeight_*pdfWeight*bgrWeight;
                    
                    if(133909==ievent)continue;///std::cout << "Value: " << branchVals_.at(1) << std::endl;
                    if(MttReco>100000)continue;
//                     std::cout << "jenya: " << MttReco << std::endl;
//                     std::cout << "jenya: " << branchVals_.at(0) << " " << branchVals_.at(1) << std::endl;
                    
                    //Control plots// 
                    ///if(MttReco>600){
                    for(int ind=0;ind<2;++ind){
                        vv_SampleHist_.at(ind).at(iSample)->Fill(branchVals_.at(ind),eventWeight_*v_weight.at(iSample));
                        for(int jnd=0;jnd<2;++jnd){
                            if(ind==jnd)continue;
                            const auto& v_bins = v_coarseBins_.at(jnd);
                            if(branchVals_.at(jnd)<v_bins.at(0))vv_UnderflowSampleHist_.at(ind).at(iSample)->Fill(branchVals_.at(ind),eventWeight_*v_weight.at(iSample));
                            if(branchVals_.at(jnd)>v_bins.at((int)v_bins.size()-1))vv_OverflowSampleHist_.at(ind).at(iSample)->Fill(branchVals_.at(ind),eventWeight_*v_weight.at(iSample));
                            for(int ibin=0;ibin<(int)v_bins.size()-1;++ibin){
                                   if(branchVals_.at(jnd) >= v_bins.at(ibin) && branchVals_.at(jnd) < v_bins.at(ibin+1)){
                                        vvv_SampleHist_.at(ind).at(ibin).at(iSample)->Fill(branchVals_.at(ind),eventWeight_*v_weight.at(iSample));
                                    }
                            }
                        }
                    }
                    v_histRecoAllBins_.at(iSample)->Fill(genBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                    //}
                    // ... //
                                   std::cout << "Test 5: !!!!!!!!!!" << std::endl ;

                  if(isSignal_){
                      for(Int_t ievent0=lastEvent;ievent0<dataTree0->GetEntriesFast();ievent0++) {
                          if(dataTree0->GetEntry(ievent0)<=0) break;
                              trueLevelWeight0_ = trueLevelWeight0_*pdfWeight*bgrWeight;
                          
                              nGenSigInc_=nGenSigInc_+trueLevelWeight0_*v_weight.at(iSample);
                              ((TH2*)histGen_)->Fill(branchValsGen0_.at(0),branchValsGen0_.at(1),trueLevelWeight0_*v_weight.at(iSample));
                              histSG_->Fill(genBin_(branchValsGen0_),trueLevelWeight0_*v_weight.at(iSample));
                          
// 			      std::cout << "Test 51: !!!!!!!!!!" << std::endl ;

			  
                              ///check before start
//                               for(int i=0; i<(int)v_coarseBins_.at(0).size()-1;i++){
// 				std::cout << "Test 511: !!!!!!!!!!" << std::endl ;
//                                   double genVal = branchValsGen0_.at(0);
// 				std::cout << "Test 512: !!!!!!!!!!" << std::endl ;
//                                   if(genVal<v_coarseBins_.at(0).at(i+1)&&genVal>v_coarseBins_.at(0).at(i)){
// 				std::cout << "Test 513: !!!!!!!!!! " << v_histOnGen_.at(i) << std::endl ;
// 				  
// 				  v_histOnGen_.at(i)->Fill(genVal,trueLevelWeight0_*v_weight.at(iSample));
// 				std::cout << "Test 514: !!!!!!!!!!" << std::endl ;
//                                   }
//                               }
                              
                              
//                               std::cout << "Test 52: !!!!!!!!!!" << std::endl ;
                              
                              //CT filling gen//
                              int genBin = genBin_(branchValsGen0_);
                              for(int i=1;i<=histSG_->GetNbinsX();i++){
                                  if(genBin==i){
                                    v_histGenRewUp_.at(i-1)->Fill(genBin,(1.3)*trueLevelWeight0_*v_weight.at(iSample));
                                    v_histGenRewDown_.at(i-1)->Fill(genBin,(0.7)*trueLevelWeight0_*v_weight.at(iSample));
                                      
                                  }
                                  else{
                                      v_histGenRewUp_.at(i-1)->Fill(genBin,trueLevelWeight0_*v_weight.at(iSample));
                                      v_histGenRewDown_.at(i-1)->Fill(genBin,trueLevelWeight0_*v_weight.at(iSample));
                                  }
                              }
                              if(v_plotName_.at(0)=="top_arapidity" && v_plotName_.at(1)=="top_pt"){
                                  histGenPtRewUp_->Fill(genBin_(branchValsGen0_),rewTopPtUp(branchValsGen0_.at(1))*trueLevelWeight0_*v_weight.at(iSample));
                                  histGenPtRewDown_->Fill(genBin_(branchValsGen0_),rewTopPtDown(branchValsGen0_.at(1))*trueLevelWeight0_*v_weight.at(iSample));
                                  histGenRapidityRewUp_->Fill(genBin_(branchValsGen0_),rewTopRapidityUp(branchValsGen0_.at(0))*trueLevelWeight0_*v_weight.at(iSample));
                                  histGenRapidityRewDown_->Fill(genBin_(branchValsGen0_),rewTopRapidityDown(branchValsGen0_.at(0))*trueLevelWeight0_*v_weight.at(iSample));
                              }
                              // ... //
                              
                              std::cout << "Test 53: !!!!!!!!!!" << std::endl ;
                              
                              if(entry0_ != entry_)
                              {
                                  histMigration_->Fill(genBin_(branchValsGen0_),0.,trueLevelWeight0_*v_weight.at(iSample));
                              }
                              else if(entry0_ == entry_)
                              {
                                  lastEvent = ievent0+1;
                                  break;
                              }
                      }
                  }
                                                     std::cout << "Test 6: !!!!!!!!!!" << std::endl ;


                  if(isKinReco_){
                    if(v_sample.at(iSample).sampleType() == Sample::data){
                        histData_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                        nDataInc_=nDataInc_+eventWeight_*v_weight.at(iSample);
                    }
                    else
                    {
                        histAllReco_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                        nRecoAllMC_=nRecoAllMC_+eventWeight_*v_weight.at(iSample);
                        if(isSignal_){
                            
                                ///check before start
//                                 for(int i=0; i<(int)v_coarseBins_.at(0).size()-1;i++){
//                                 double genVal = branchValsGen_.at(0);
//                                 double recoVal = branchVals_.at(0);
//                                 if(genVal<v_coarseBins_.at(0).at(i+1)&&genVal>v_coarseBins_.at(0).at(i)){
//                                     v_histOnReco_.at(i)->Fill(recoVal,eventWeight_*v_weight.at(iSample));
//                                     }
//                                 }
                            
                            nRecoSigInc_=nRecoSigInc_+eventWeight_*v_weight.at(iSample);
                            
                            histTtSigReco_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                            histTtAllReco_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                            nRecoTtAllInc_=nRecoTtAllInc_+eventWeight_*v_weight.at(iSample);
                            
                            if( (branchValsGen_.at(0) <v_coarseBins_.at(0).at(0) && v_uTrue_.at(0)==0)
                                || (branchValsGen_.at(0)>v_coarseBins_.at(0).back() && v_oTrue_.at(0)==0)
                                || (branchValsGen_.at(1) <v_coarseBins_.at(1).at(0) && v_uTrue_.at(1)==0)
                                || (branchValsGen_.at(1)>v_coarseBins_.at(1).back() && v_oTrue_.at(1)==0) )
                            {
                                histBgrUo_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                            }
                            
                            histSR_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                            
                            
                            //CT filling reco//
                             int genBin = genBin_(branchValsGen_);
                              for(int i=1;i<=histSG_->GetNbinsX();i++){
                                  if(genBin==i){
                                    v_histRecoRewUp_.at(i-1)->Fill(recoBin_(branchVals_),(1.3)*eventWeight_*v_weight.at(iSample));
                                    v_histRecoRewDown_.at(i-1)->Fill(recoBin_(branchVals_),(0.7)*eventWeight_*v_weight.at(iSample));
                                  }
                                  else{
                                      v_histRecoRewUp_.at(i-1)->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                                      v_histRecoRewDown_.at(i-1)->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                                  }
                              }
                              if(v_plotName_.at(0)=="top_arapidity" && v_plotName_.at(1)=="top_pt"){
                                  histRecoPtRewUp_->Fill(recoBin_(branchVals_),rewTopPtUp(branchValsGen_.at(1))*eventWeight_*v_weight.at(iSample));
                                  histRecoPtRewDown_->Fill(recoBin_(branchVals_),rewTopPtDown(branchValsGen_.at(1))*eventWeight_*v_weight.at(iSample));
                                  histRecoRapidityRewUp_->Fill(recoBin_(branchVals_),rewTopRapidityUp(branchValsGen_.at(0))*eventWeight_*v_weight.at(iSample));
                                  histRecoRapidityRewDown_->Fill(recoBin_(branchVals_),rewTopRapidityDown(branchValsGen_.at(0))*eventWeight_*v_weight.at(iSample));
                              }
                              // ... // 
                            
                            histMigration_->Fill(genBin_(branchValsGen_),recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                            histMigration_->Fill(genBin_(branchValsGen_),0.,trueLevelWeight_*v_weight.at(iSample)-eventWeight_*v_weight.at(iSample));
                            
                            histEffAllBins_->Fill(genBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                            histPurityAllBins_->Fill(genBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                            histStabilityAllBins_->Fill(genBin_(branchValsGen_),eventWeight_*v_weight.at(iSample));
                            if(genBin_(branchVals_)==genBin_(branchValsGen_))histRecoGenAllBins_->Fill(genBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                            
                        }
                        else
                        {
                            if(sample.legendEntry()=="Single Top")nRecoSingleTop_ += eventWeight_*v_weight.at(iSample);
                            if(sample.legendEntry()=="W+Jets")nRecoWPlusJets_ += eventWeight_*v_weight.at(iSample);
                            if(sample.legendEntry()=="Z / #gamma* #rightarrow ee/#mu#mu")nRecoZee_ += eventWeight_*v_weight.at(iSample);
                            if(sample.legendEntry()=="Z / #gamma* #rightarrow #tau#tau")nRecoZtautau_ += eventWeight_*v_weight.at(iSample);
                            if(sample.legendEntry()=="t#bar{t}+Z/W/#gamma")nRecoTtPlusZWG_ += eventWeight_*v_weight.at(iSample);
                            if(sample.legendEntry()=="Diboson")nRecoDiboson_ += eventWeight_*v_weight.at(iSample);
                            if(sample.legendEntry()=="QCD Multijet")nRecoQCD_ += eventWeight_*v_weight.at(iSample);
                            
                            if(v_sample.at(iSample).sampleType() == Sample::ttother && ((isPDF_&&(! sample.inputFile().Contains("ttbarbg.root")))||(!isPDF_)) ){
                                histTtAllReco_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                                nRecoTtAllInc_=nRecoTtAllInc_+eventWeight_*v_weight.at(iSample);
                                nRecoTtOtherInc_=nRecoTtOtherInc_+eventWeight_*v_weight.at(iSample);
                            }
                            else
                            {
                                ///check before start
                                histBgr_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample)*bgrWeight);
                                nBgrInc_=nBgrInc_+eventWeight_*v_weight.at(iSample)*bgrWeight;
                                //histBgr_->Fill(recoBin_(branchVals_),eventWeight_*v_weight.at(iSample));
                                //nBgrInc_=nBgrInc_+eventWeight_*v_weight.at(iSample);
                            }
                        }
                    }
                  }
                    
                 /// ********************************************** Filling Histograms END  ******************************************** ///
                    
                }  // dataTree8 - loop
                
                dataTree0->Delete();
                dataTree8->Delete();
                dataFile->Close();
                dataFile->Delete();

            }//samples loop
            

            TCanvas* cRecoData = utils::setCanvas();
            TH1* hRecoData = plotNiceTitelReco(0,histData_->GetMaximum()*1.3);
            histData_->Draw("same e1");
            histAllReco_->SetLineColor(2);
            histAllReco_->Draw("same hist");
            writeCanvas(cRecoData,"recoData");
            
            histFsigReco_->Divide(histTtSigReco_, histTtAllReco_,1,1,"B");
            histFsigReco_->Draw();
            writeCanvas(cRecoData,"histFsigReco_");
            
            hRecoData->Delete();
            cRecoData->Delete();
            
            histAllRecoScaled_->Multiply(histAllReco_,histFsigReco_,1,1);
            histDataScaled_->Multiply(histData_,histFsigReco_,1,1);
            histBgrScaled_->Multiply(histBgr_,histFsigReco_,1,1);
            histBgrUoScaled_->Multiply(histBgrUo_,histFsigReco_,1,1);
            
            if(isTauTest_){
            /// ******  Tau scan START **** ///
                TH1D* histFluctuations = new TH1D("Fluctuations","Fluctuations;;Rms/meanErr",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax());
                TH1D* histDiff = new TH1D("Diff","Diff;;(mean-x0)",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax());
                TH2D* histUnfGen_Vs_Bin = new TH2D("UnfGen_Vs_Bin","UnfGen_Vs_Bin",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax(),1000,0,2*histSG_->GetMaximum());
                TH2D* histUnfGenErr_Vs_Bin = new TH2D("UnfGenErr_Vs_Bin","UnfGenErr_Vs_Bin",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax(),1000,0,0.2*histSG_->GetMaximum());
                
                TH2D* histUnfGen_Vs_Bin_Rho = new TH2D("UnfGen_Vs_Bin_Rho","UnfGen_Vs_Bin_Rho",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax(),1000,0,2*histSG_->GetMaximum());
                TH2D* histUnfGen_Vs_Bin_L = new TH2D("UnfGen_Vs_Bin_L","UnfGen_Vs_Bin_L",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax(),1000,0,2*histSG_->GetMaximum());
                TH2D* histUnfGen_Vs_Bin_0 = new TH2D("UnfGen_Vs_Bin_0","UnfGen_Vs_Bin_0",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax(),1000,0,2*histSG_->GetMaximum());
                
                TH1D* histFluctuations_Rho = new TH1D("Fluctuations_Pho","Fluctuations;;RMS(X)/<#DeltaX>",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax());
                TH1D* histFluctuations_L = new TH1D("Fluctuations_L","Fluctuations;;RMS(X)/<#DeltaX",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax());
                TH1D* histFluctuations_0 = new TH1D("Fluctuations_0","Fluctuations;;RMS(X)/<#DeltaX",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax());
                
                TProfile* histError_Rho = new TProfile("Error_Rho","mean of error;;<#DeltaX>",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax());
                TProfile* histError_L = new TProfile("Error_L","mean of error;;<#DeltaX>",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax());
                TProfile* histError_0 = new TProfile("Error_0","mean of error;;<#DeltaX>",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax());
                
                TProfile* histRelError_Rho = new TProfile("RelError_Rho","mean of relative error;;<#DeltaX/X>",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax());
                TProfile* histRelError_L = new TProfile("RelError_L","mean of relative error;;<#DeltaX/X>",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax());
                TProfile* histRelError_0 = new TProfile("RelError_0","mean of relative error;;<#DeltaX/X>",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax());
                
                TProfile* histRelDiff_Rho = new TProfile("Reldiff_Rho","mean of relative diff;;<(X-X0)/X0>",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax());
                TProfile* histRelDiff_L = new TProfile("Reldiff_L","mean of relative diff;;<(X-X0)/X0>",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax());
                TProfile* histRelDiff_0 = new TProfile("Reldiff_0","mean of relative diff;;<(X-X0)/X0>",histSG_->GetNbinsX(),histSG_->GetXaxis()->GetXmin(),histSG_->GetXaxis()->GetXmax());
                
                TH1D* histSigmaProb_Rho = new TH1D("histSigmaProb_Rho",";#pm#sigma range;%",5,0.5,5.5);
                TH1D* histSigmaProb_L = new TH1D("histSigmaProb_L",";#pm#sigma range;%",5,0.5,5.5);
                TH1D* histSigmaProb_0 = new TH1D("histSigmaProb_0",";#pm#sigma range;%",5,0.5,5.5);
                
                TH1D* histPDF_Rho = new TH1D("histPDF_Rho",";#chi2;",40,0,40);
                TH1D* histPDF_L = new TH1D("histPDF_L",";#chi2;",40,0,40);
                TH1D* histPDF_0 = new TH1D("histPDF_0",";#chi2;",40,0,40);
                histPDF_Rho->Sumw2();
                histPDF_L->Sumw2();
                
                
                int Nx = histSG_->GetNbinsX();
                
            for(int i=0;i<6000;i++){
            //for(int i=0;i<1000;i++){
                TH1* hist_Smeared=(TH1*)histAllRecoScaled_->Clone();
                performPoissonSmearing(histAllRecoScaled_,hist_Smeared);
                regMethod_="rho";
                runUnfolding(histMigration_,hist_Smeared,histBgrScaled_,histBgrUoScaled_);
                for(int j=1;j<=Nx;j++){
                    double x0 = histSG_->GetBinContent(j);
                    double x = histUnfGenBin_->GetBinContent(j);
                    double sigma = histUnfGenBin_->GetBinError(j);
                    double binX = histSG_->GetBinCenter(j);
                    for(int k=1;k<=5;k++){
                        if(x0<(x+sigma*k) && x0>(x-sigma*k))histSigmaProb_Rho->Fill(k);
                    }
                    histUnfGen_Vs_Bin->Fill(binX,x);
                    histUnfGenErr_Vs_Bin->Fill(binX,sigma);
                    histUnfGen_Vs_Bin_Rho->Fill(binX,x);
                    histRelError_Rho->Fill(binX,sigma/x);
                        histRelDiff_Rho->Fill(binX,(x-x0)/x0);
                    histError_Rho->Fill(binX,sigma);
                }
                if(1){
                    TMatrixD dV_X(1,(int)vX_.size());
                    dV_X.SetMatrixArray(utils::diffVect(vX_,vXmc_).data());
                    TMatrixD dVT_X = dV_X;
                    dVT_X.T();
                    TMatrixD covXtotalMatrix((int)vX_.size(),(int)vX_.size());
                    covXtotalMatrix.SetMatrixArray(v_covX_.data());
                    TDecompChol decompChol_X(covXtotalMatrix);
                    TMatrixD invertCovXtotalMatrix = ((TMatrixDSym)(decompChol_X.Invert()));
                    histPDF_Rho->Fill((dV_X*invertCovXtotalMatrix*dVT_X)(0,0));
                }
                
                
                
                regMethod_="l";
                runUnfolding(histMigration_,hist_Smeared,histBgrScaled_,histBgrUoScaled_);
                for(int j=1;j<=Nx;j++){
                    double x0 = histSG_->GetBinContent(j);
                    double x = histUnfGenBin_->GetBinContent(j);
                    double sigma = histUnfGenBin_->GetBinError(j);
                    double binX = histSG_->GetBinCenter(j);
                    for(int k=1;k<=5;k++){
                        if(x0<(x+sigma*k) && x0>(x-sigma*k))histSigmaProb_L->Fill(k);
                    }
                    histUnfGen_Vs_Bin_L->Fill(binX,x);
                    histRelError_L->Fill(binX,sigma/x);
                    histRelDiff_L->Fill(binX,(x-x0)/x0);
                    histError_L->Fill(binX,sigma);
                }
                if(1){
                    TMatrixD dV_X(1,(int)vX_.size());
                    dV_X.SetMatrixArray(utils::diffVect(vX_,vXmc_).data());
                    TMatrixD dVT_X = dV_X;
                    dVT_X.T();
                    TMatrixD covXtotalMatrix((int)vX_.size(),(int)vX_.size());
                    covXtotalMatrix.SetMatrixArray(v_covX_.data());
                    TDecompChol decompChol_X(covXtotalMatrix);
                    TMatrixD invertCovXtotalMatrix = ((TMatrixDSym)(decompChol_X.Invert()));
                    histPDF_L->Fill((dV_X*invertCovXtotalMatrix*dVT_X)(0,0));
                }
                
                regMethod_="0";
                runUnfolding(histMigration_,hist_Smeared,histBgrScaled_,histBgrUoScaled_);
                for(int j=1;j<=Nx;j++){
                    double x0 = histSG_->GetBinContent(j);
                    double x = histUnfGenBin_->GetBinContent(j);
                    double sigma = histUnfGenBin_->GetBinError(j);
                    // OZ to eliminate unused var warning
                    //double binX = histSG_->GetBinCenter(j);
                    for(int k=1;k<=5;k++){
                        if(x0<(x+sigma*k) && x0>(x-sigma*k))histSigmaProb_0->Fill(k);
                    }
                    histUnfGen_Vs_Bin_0->Fill(histSG_->GetBinCenter(j),x);
                    histRelError_0->Fill(histSG_->GetBinCenter(j),sigma/x);
                    histRelDiff_0->Fill(histSG_->GetBinCenter(j),(x-x0)/x0);
                    histError_0->Fill(histSG_->GetBinCenter(j),sigma);
                }
                if(1){
                    TMatrixD dV_X(1,(int)vX_.size());
                    dV_X.SetMatrixArray(utils::diffVect(vX_,vXmc_).data());
                    TMatrixD dVT_X = dV_X;
                    dVT_X.T();
                    TMatrixD covXtotalMatrix((int)vX_.size(),(int)vX_.size());
                    covXtotalMatrix.SetMatrixArray(v_covX_.data());
                    TDecompChol decompChol_X(covXtotalMatrix);
                    TMatrixD invertCovXtotalMatrix = ((TMatrixDSym)(decompChol_X.Invert()));
                    histPDF_0->Fill((dV_X*invertCovXtotalMatrix*dVT_X)(0,0));
                }
                

                if(hist_Smeared)hist_Smeared->Delete();
            }
            
            for(int ix=1;ix<=histSG_->GetNbinsX();ix++){
                TH1* histUnfGen_Vs_BinXi = ((TH1*)(histUnfGen_Vs_Bin->ProjectionY("UnfGen_Vs_BinXi",ix,ix,"e")));
                TH1* histUnfGen_Vs_BinXi_Rho = ((TH1*)(histUnfGen_Vs_Bin_Rho->ProjectionY("UnfGen_Vs_BinXi_Rho",ix,ix,"e")));
                TH1* histUnfGen_Vs_BinXi_L = ((TH1*)(histUnfGen_Vs_Bin_L->ProjectionY("UnfGen_Vs_BinXi_L",ix,ix,"e")));
                TH1* histUnfGen_Vs_BinXi_0 = ((TH1*)(histUnfGen_Vs_Bin_0->ProjectionY("UnfGen_Vs_BinXi_0",ix,ix,"e")));
                
                TH1* histUnfGenErr_Vs_BinXi = ((TH1*)(histUnfGenErr_Vs_Bin->ProjectionY("UnfGenErr_Vs_BinXi",ix,ix,"e")));
                double mean=histUnfGen_Vs_BinXi->GetMean();
                double meanE=histUnfGen_Vs_BinXi->GetMeanError();
                double rms  = histUnfGen_Vs_BinXi->GetRMS();
                double rmsE = histUnfGen_Vs_BinXi->GetRMSError();
                double meanErr=histUnfGenErr_Vs_BinXi->GetMean();
                double meanErrE=histUnfGenErr_Vs_BinXi->GetMeanError();
                histFluctuations->SetBinContent(ix,rms/meanErr);
                histFluctuations->SetBinError(ix,(rms/meanErr)*sqrt(rmsE*rmsE/rms/rms+meanErrE*meanErrE/meanErr/meanErr));
                
                histFluctuations_Rho->SetBinContent(ix,histUnfGen_Vs_BinXi_Rho->GetRMS()/histError_Rho->GetBinContent(ix));
                histFluctuations_L->SetBinContent(ix,histUnfGen_Vs_BinXi_L->GetRMS()/histError_L->GetBinContent(ix));
                histFluctuations_0->SetBinContent(ix,histUnfGen_Vs_BinXi_0->GetRMS()/histError_0->GetBinContent(ix));
                histFluctuations_Rho->SetBinError(ix,0.0001);
                histFluctuations_L->SetBinError(ix,0.0001);
                histFluctuations_0->SetBinError(ix,0.0001);
                
                
                histDiff->SetBinContent(ix,(mean-histSG_->GetBinContent(ix))/meanErr);
                histDiff->SetBinError(ix, sqrt(pow((meanE/meanErr),2)+(mean-histSG_->GetBinContent(ix))*pow(((meanErrE)/(meanErr*meanErr)),2)));
                histUnfGen_Vs_BinXi->Delete();
                histUnfGenErr_Vs_BinXi->Delete();
            }
            
            
            
            TCanvas* canvas = utils::setCanvas();
            TLegend* legend = utils::setLegend();
               histUnfGen_Vs_Bin->Draw("colz");
                canvas->Print(tauFolder_+ "UnfGen_Vs_Bin_"+systematic.name()+".png");
                canvas->Clear();
               histUnfGenErr_Vs_Bin->Draw("colz");
                canvas->Print(tauFolder_+ "UnfGenErr_Vs_Bin_"+systematic.name()+".png");
                canvas->Clear();
              histSG_->Draw();
                canvas->Print(tauFolder_+ "histSG_"+systematic.name()+".png");
                canvas->Clear();
                
               histTauLog_->DrawClone();
               // fit logTau //
//                double r1=-10,r2=10;
//                int maxBinTau = histTauLog_->GetMaximumBin();
//                double binContent=histTauLog_->GetBinContent(maxBinTau);
//                int binR1 = maxBinTau;
//                int binR2 = maxBinTau;
//                 while(binContent>0){
//                    binR1--;
//                    binR2++;
//                    binContent = histTauLog_->GetBinContent(binR1)<histTauLog_->GetBinContent(binR2) ? histTauLog_->GetBinContent(binR1) : histTauLog_->GetBinContent(binR2);
//                 }
//                 r1=histTauLog_->GetBinCenter(binR1+1);
//                 r2=histTauLog_->GetBinCenter(binR2-1);
//                TF1 tauFit("tauFit","gaus",r1,r2);
//                 tauFit.SetParameters(10,histTauLog_->GetBinCenter(maxBinTau),1);
//                 histTauLog_->Fit("tauFit","R");
//                 tau_ = pow(10,tauFit.GetParameter(1));
                runUnfolding(histMigration_,histDataScaled_,histBgrScaled_,histBgrUoScaled_);
                TLine logTauLine(logTau_,0,logTau_,histTauLog_->GetMaximum()*1.1);
                logTauLine.SetLineColor(2);
                //logTauLine.Draw("same");
                canvas->Print(tauFolder_+ "tauDist_"+systematic.name()+".png");
                canvas->Print(tauFolder_+ "tauDist_"+systematic.name()+".root");
                canvas->Clear();
                
                histTauLogRhoScan_->DrawClone();
                canvas->Print(tauFolder_+ "tauDistRhoScan_"+systematic.name()+".png");
                canvas->Print(tauFolder_+ "tauDistRhoScan_"+systematic.name()+".root");
                canvas->Clear();
                
                histTauLogLCurve_->DrawClone();
                canvas->Print(tauFolder_+ "tauDistLCurve_"+systematic.name()+".png");
                canvas->Print(tauFolder_+ "tauDistLCurve_"+systematic.name()+".root");
                canvas->Clear();
                
                histSigmaProb_Rho->SetMarkerColor(kBlue);
                histSigmaProb_L->SetMarkerColor(kRed);
                histSigmaProb_0->SetMarkerColor(kBlack);
                histSigmaProb_Rho->Scale(1.0/(10*histSG_->GetNbinsX()));
                histSigmaProb_L->Scale(1.0/(10*histSG_->GetNbinsX()));
                histSigmaProb_0->Scale(1.0/(10*histSG_->GetNbinsX()));
                histSigmaProb_Rho->SetStats(0);
                histSigmaProb_Rho->GetYaxis()->SetRangeUser(50,110);
                histSigmaProb_Rho->Draw("e1");
                histSigmaProb_L->Draw("same e1");
                histSigmaProb_0->Draw("same e1");
                legend->Clear();
                legend->AddEntry(histSigmaProb_0,"#tau = 0","pe");
                legend->AddEntry(histSigmaProb_L,"L-Scan","pe");
                legend->AddEntry(histSigmaProb_Rho,"#rho min","pe");
                legend->Draw("same");
                canvas->Print(tauFolder_+ "SigmaProb_"+systematic.name()+".png");
                canvas->Print(tauFolder_+ "SigmaProb_"+systematic.name()+".root");
                canvas->Clear();
                
                
                histPDF_Rho->SetMarkerColor(kBlue);
                histPDF_L->SetMarkerColor(kRed);
                histPDF_0->SetMarkerColor(kBlack);
                histPDF_Rho->Scale(1.0/histPDF_Rho->Integral());
                histPDF_L->Scale(1.0/histPDF_L->Integral());
                histPDF_0->Scale(1.0/histPDF_0->Integral());
                histPDF_Rho->SetStats(0);
                histPDF_Rho->Draw("e1");
                histPDF_L->Draw("same e1");
                histPDF_0->Draw("same e1");
                legend->Clear();
                legend->AddEntry(histPDF_0,"#tau = 0","pe");
                legend->AddEntry(histPDF_L,"L-Scan","pe");
                legend->AddEntry(histPDF_Rho,"#rho min","pe");
                legend->Draw("same");
                canvas->Print(tauFolder_+ "PDF_"+systematic.name()+".png");
                canvas->Print(tauFolder_+ "PDF_"+systematic.name()+".root");
                canvas->Clear();
                
                
                
               TH1D* hist = plotNiceTitel(-1,1,"(mean(X) - X_{0}) / mean(Stat.Err.)");
               histDiff->Draw("same e1");
                canvas->Print(tauFolder_+ "Diff_"+systematic.name()+".png");
                canvas->Clear();
                hist->Delete();
                
               hist = plotNiceTitel(0.5,1.5,"rms / mean(Stat.Err.)");
               histFluctuations->Draw("same e1");
               TLine horizontalLn(0,1,100,1);
               horizontalLn.DrawClone("same");
                canvas->Print(tauFolder_+ "fluctuation_"+systematic.name()+".png");
                canvas->Clear();
                hist->Delete();
                
                
               hist = plotNiceTitel(0.,0.4,"< #sigma(X) / X >");
               histRelError_Rho->SetMarkerColor(kRed);
               histRelError_L->SetMarkerColor(kBlue);
               histRelError_0->SetMarkerColor(kBlack);
               histRelError_0->Draw("same e1");
               histRelError_L->Draw("same e1");
               histRelError_Rho->Draw("same e1");
               //legend->SetX1(0.1);
               //legend->SetX2(0.4);
               legend->Clear();
               legend->AddEntry(histRelError_0,"#tau = 0","pe");
               legend->AddEntry(histRelError_L,"L-Scan","pe");
               legend->AddEntry(histRelError_Rho,"#rho min","pe");
               legend->Draw("same");
                canvas->Print(tauFolder_+ "relError_"+systematic.name()+".png");
                canvas->Print(tauFolder_+ "relError_"+systematic.name()+".pdf");
                canvas->Print(tauFolder_+ "relError_"+systematic.name()+".root");
                canvas->Clear();
                hist->Delete();
               
                
                
                
               hist = plotNiceTitel(-0.1,0.1,"<(X-X_{0})/X_{0}>");
               histRelDiff_Rho->SetMarkerColor(kBlue);
               histRelDiff_L->SetMarkerColor(kRed);
               histRelDiff_0->SetMarkerColor(kBlack); 
               //histRelDiff_0->Draw("same e1");
               histRelDiff_L->Draw("same e1");
               histRelDiff_Rho->Draw("same e1");
               //legend->SetX1(0.1);
               //legend->SetX2(0.4);
               legend->Clear();
               //legend->AddEntry(histRelDiff_0,"#tau = 0","pe");
               legend->AddEntry(histRelDiff_L,"L-Scan","pe");
               legend->AddEntry(histRelDiff_Rho,"#rho min","pe");
               legend->Draw("same");
                canvas->Print(tauFolder_+ "relDiff_"+systematic.name()+".png");
                canvas->Print(tauFolder_+ "relDiff_"+systematic.name()+".pdf");
                canvas->Print(tauFolder_+ "relDiff_"+systematic.name()+".root");
                canvas->Clear();
                hist->Delete();
                
                
               hist = plotNiceTitel(0.5,1.5,"RMS(X) / <#sigma(X)>");
               histFluctuations_Rho->SetMarkerColor(kBlue);
               histFluctuations_L->SetMarkerColor(kRed);
               histFluctuations_0->SetMarkerColor(kBlack);
               //histFluctuations_0->Draw("same e0");
               histFluctuations_L->Draw("same e0");
               histFluctuations_Rho->Draw("same e0");
               //TLine horizontalLn(0,1,100,1);
               horizontalLn.DrawClone("same");
               //legend->SetX1(0.1);
               //legend->SetX2(0.4);
               legend->Clear();
               //legend->AddEntry(histFluctuations_0,"#tau = 0","pe");
               legend->AddEntry(histFluctuations_L,"L-Scan","pe");
               legend->AddEntry(histFluctuations_Rho,"#rho min","pe");
               legend->Draw("same");
                canvas->Print(tauFolder_+ "rmsOverMeanError_"+systematic.name()+".png");
                canvas->Print(tauFolder_+ "rmsOverMeanError_"+systematic.name()+".pdf");
                canvas->Print(tauFolder_+ "rmsOverMeanError_"+systematic.name()+".root");
                canvas->Clear();
                hist->Delete();
                
                
                
            histFluctuations->Delete();
            histDiff->Delete();
            histUnfGen_Vs_Bin->Delete();
            histUnfGenErr_Vs_Bin->Delete();
            
            histRelError_Rho->Delete();
            histRelError_L->Delete();
            histRelError_0->Delete();
            
            canvas->Delete();
            
            /// ******  Tau scan END  ****** ///
            
            }
            else if(isCT_){
                /// *****  Closure Test START  ***** ///
//                 for(int i=1;i<=histSG_->GetNbinsX();i++){
//                     TCanvas c("","",600,600);
//                         TH1D* hist = plotNiceTitel(0,histSG_->GetMaximum()*1.2);
//                         v_histGenRewUp_.at(i-1)->SetStats(0);
//                         v_histGenRewUp_.at(i-1)->Draw("hist");
//                         v_histGenRewUp_.at(i-1)->SetLineColor(kMagenta);
//                         histSG_->SetLineColor(1);
//                         runUnfolding(histMigration_,v_histRecoRewUp_.at(i-1),0,0);
//                         histUnfGenBin_->SetMarkerColor(kMagenta);
//                         histUnfGenBin_->Draw("same e1");
//                         
//                         v_histGenRewDown_.at(i-1)->SetLineColor(kBlue);
//                         v_histGenRewDown_.at(i-1)->Draw("hist,same");
//                         runUnfolding(histMigration_,v_histRecoRewDown_.at(i-1),0,0);
//                         histUnfGenBin_->SetMarkerColor(kBlue);
//                         histUnfGenBin_->Draw("same e1");
//                         
//                         histSG_->Draw("hist,same");
//                     c.Print(ctFolder_+ "bin"+std::to_string(i)+"_"+systematic.name()+".png");
//                 }
                
                if(v_plotName_.at(0)=="top_arapidity" && v_plotName_.at(1)=="top_pt"){
                TCanvas* c = new TCanvas("","",600,600);
                TLegend* leg = utils::setLegend();
                
                TH1D* hist = plotNiceTitel(0,histSG_->GetMaximum()*1.2,v_plotYunitscp_.at(0));
                hist->SetTitleOffset(1.75,"Y");
                hist->SetTitleSize(hist->GetTitleSize("Y")*1.5,"Y");
                gPad->SetLeftMargin(0.2);
                histGenPtRewUp_->SetStats(0);
                histGenPtRewUp_->Draw("hist,same");
                histGenPtRewUp_->SetLineColor(kMagenta);
                histGenPtRewUp_->SetLineWidth(2);
                histGenPtRewUp_->SetLineStyle(2);
                histSG_->SetLineColor(2);
                histSG_->SetLineWidth(2);
                histSG_->SetLineColor(kRed-7);
                for(int i=1;i<=histRecoPtRewUp_->GetNbinsX();i++)histRecoPtRewUp_->SetBinError(i,sqrt(histRecoPtRewUp_->GetBinContent(i)));
                runUnfolding(histMigration_,histRecoPtRewUp_,0,0);
                histUnfGenBin_->Draw("same e1");
                histSG_->Draw("hist,same");
                
                leg->AddEntry(histUnfGenBin_,"Pseudo-Data","p");
                leg->AddEntry(histSG_,"Gen.","l");
                leg->AddEntry(histGenPtRewUp_,"Pseudo-Gen.","l");
                leg->Draw("same");
                c->Print(ctFolder_+ "PtRewUp_"+systematic.name()+".png");
                c->Print(ctFolder_+ "PtRewUp_"+systematic.name()+".pdf");
                c->Print(ctFolder_+ "PtRewUp_"+systematic.name()+".root");
                hist->Delete();
                
                hist = plotNiceTitel(0,histSG_->GetMaximum()*1.2,v_plotYunitscp_.at(0));
                hist->SetTitleOffset(1.75,"Y");
                hist->SetTitleSize(hist->GetTitleSize("Y")*1.5,"Y");
                gPad->SetLeftMargin(0.2);
                histGenPtRewDown_->SetStats(0);
                histGenPtRewDown_->Draw("hist,same");
                histGenPtRewDown_->SetLineColor(kMagenta);
                histGenPtRewDown_->SetLineWidth(2);
                histGenPtRewDown_->SetLineStyle(2);
                histSG_->SetLineColor(2);
                histSG_->SetLineWidth(2);
                histSG_->SetLineColor(kRed-7);
                for(int i=1;i<=histRecoPtRewDown_->GetNbinsX();i++)histRecoPtRewDown_->SetBinError(i,sqrt(histRecoPtRewDown_->GetBinContent(i)));
                runUnfolding(histMigration_,histRecoPtRewDown_,0,0);
                histUnfGenBin_->Draw("same e1");
                histSG_->Draw("hist,same");
                leg->Draw("same");
                c->Print(ctFolder_+ "PtRewDown_"+systematic.name()+".png");
                c->Print(ctFolder_+ "PtRewDown_"+systematic.name()+".pdf");
                c->Print(ctFolder_+ "PtRewDown_"+systematic.name()+".root");
                hist->Delete();
                
                
                hist = plotNiceTitel(0,histSG_->GetMaximum()*1.3,v_plotYunitscp_.at(0));
                hist->SetTitleOffset(1.75,"Y");
                hist->SetTitleSize(hist->GetTitleSize("Y")*1.5,"Y");
                gPad->SetLeftMargin(0.2);
                histGenRapidityRewUp_->SetStats(0);
                histGenRapidityRewUp_->Draw("hist,same");
                histGenRapidityRewUp_->SetLineColor(kMagenta);
                histGenRapidityRewUp_->SetLineWidth(2);
                histGenRapidityRewUp_->SetLineStyle(2);
                histSG_->SetLineColor(2);
                histSG_->SetLineWidth(2);
                histSG_->SetLineColor(kRed-7);
                for(int i=1;i<=histRecoRapidityRewUp_->GetNbinsX();i++)histRecoRapidityRewUp_->SetBinError(i,sqrt(histRecoRapidityRewUp_->GetBinContent(i)));
                runUnfolding(histMigration_,histRecoRapidityRewUp_,0,0);
                histUnfGenBin_->Draw("same e1");
                histSG_->Draw("hist,same");
                leg->Draw("same");
                c->Print(ctFolder_+ "RapidityRewUp_"+systematic.name()+".png");
                c->Print(ctFolder_+ "RapidityRewUp_"+systematic.name()+".pdf");
                c->Print(ctFolder_+ "RapidityRewUp_"+systematic.name()+".root");
                hist->Delete();
                
                hist = plotNiceTitel(0,histSG_->GetMaximum()*1.2,v_plotYunitscp_.at(0));
                hist->SetTitleOffset(1.75,"Y");
                hist->SetTitleSize(hist->GetTitleSize("Y")*1.5,"Y");
                gPad->SetLeftMargin(0.2);
                histGenRapidityRewDown_->SetStats(0);
                histGenRapidityRewDown_->Draw("hist,same");
                histGenRapidityRewDown_->SetLineColor(kMagenta);
                histGenRapidityRewDown_->SetLineWidth(2);
                histGenRapidityRewDown_->SetLineStyle(2);
                histSG_->SetLineColor(2);
                histSG_->SetLineWidth(2);
                histSG_->SetLineColor(kRed-7);
                for(int i=1;i<=histRecoRapidityRewDown_->GetNbinsX();i++)histRecoRapidityRewDown_->SetBinError(i,sqrt(histRecoRapidityRewDown_->GetBinContent(i)));
                runUnfolding(histMigration_,histRecoRapidityRewDown_,0,0);
                histUnfGenBin_->Draw("same e1");
                histSG_->Draw("hist,same");
                leg->Draw("same");
                c->Print(ctFolder_+ "RapidityRewDown_"+systematic.name()+".png");
                c->Print(ctFolder_+ "RapidityRewDown_"+systematic.name()+".pdf");
                c->Print(ctFolder_+ "RapidityRewDown_"+systematic.name()+".root");
                hist->Delete();
                
                c->Delete();
                }
                
                /// *****  Closure Test END ***** ///
            }
            else{
                
            /// ***** Cross Sections START ***** ///
            // Full //
                topxsec_ = ((nDataInc_ - nBgrInc_)*(nRecoSigInc_/nRecoTtAllInc_))/(nRecoSigInc_/nGenSigInc_)/lumi_/v_BR_.at(channel);
                sampleInfo.addLine("fullTopXSec = "+utils::numToString(topxsec_));
            sampleInfo.save(sampleInfoFolder);
            // ... //
            //Control plots//
            for(int ind=0;ind<2;++ind){
              writePlotCP( v_sample ,vv_SampleHist_.at(ind),ind);
                    //if(v_uTrue_.at(!ind))writePlotCP( v_sample ,vv_UnderflowSampleHist_.at(ind),ind,-1);
                    //if(v_oTrue_.at(!ind))writePlotCP( v_sample ,vv_OverflowSampleHist_.at(ind),ind,-2);
                for(int jnd=0;jnd<2;++jnd){
                  if(ind==jnd)continue;
                   const auto& v_bins = v_coarseBins_.at(jnd);
                   for(int ibin=0;ibin<(int)v_bins.size()-1;++ibin){
                      writePlotCP( v_sample ,vvv_SampleHist_.at(ind).at(ibin),ind,ibin);
                 }
                }
            }
            writePlotCPAllBins(v_sample,v_histRecoAllBins_);
            // ... //
            
    ///check before start
// //     TCanvas* can1 = utils::setCanvas();
// //     TLegend* leg1 = utils::setLegend(0.6,0.50,1.10,0.85);
// //     can1->cd();
// //     
// //     for(int i=0;i<(int)v_histOnGen_.size();i++){
// //         v_histOnGen_.at(i)->SetStats(0);
// //         v_histOnGen_.at(i)->SetTitle(";p_{T}(t), GeV;");
// //         if(i==0){
// //             v_histOnGen_.at(i)->Draw("hist");
// //             v_histOnGen_.at(i)->GetYaxis()->SetRangeUser(0,45000);
// //         }
// //         else v_histOnGen_.at(i)->Draw("hist,same");
// //         can1->Print(plotsFolder_+ "unfGen_" + v_plotName_.at(0) + std::to_string(i) + ".png");
// //     }
// //         
// //     THStack* stack1 = new THStack();
// //         stack1->SetTitle(";p_{T}(t), GeV;");
// //     for(int i=0;i<(int)v_histOnGen_.size();i++){
// //         //v_histOnReco_.at(i)->GetYaxis()->SetRangeUser(0,45000);
// //         stack1->Add(v_histOnReco_.at(i));
// //         stack1->Draw("HIST");
// //         can1->Print(plotsFolder_+ "unfReco_" + v_plotName_.at(0) + std::to_string(i) + ".png");
// //         v_histOnReco_.at(i)->SetLineColor(v_histOnReco_.at(i)->GetFillColor());
// //         
// //     }
// //     delete can1;
// //     delete leg1;
// //     // ...
            
            
            
            
            
            
            writePlotEPSAllBins();
            
            //Unfolding//
           cleanMigrationMatrix(histMigration_);
           
           ///check before start
           if(systematicName_!="Nominal"||isVar_)regMethod_="fix";
               //regMethod_="rho";
                //regMethod_="0";
           //tau_=0.00313052;
           //regMethod_="fix";
           
	   
	   std::cout << "Test 4: !!!!!!!!!!" << std::endl ;    
	   
           if(!isPDF_){runUnfolding(histMigration_,histDataScaled_,histBgrScaled_,histBgrUoScaled_);
           }
           else runUnfolding(histMigration_,histDataScaled_,histBgrScaled_,0);
           
            writePlotXSec(unfoldedDataNorm_,histGenNorm_,unfoldedData_,histGen_,channel);
            // ... //
           /// ***** Cross Sections END ***** ///
            } 
            
            
        }//channel loop
    }//systematics loop
    //std::cout<<"\n[PlotterForTUnfold]:=== Finishing of plot production\n\n";
}


int PlotterForTUnfold::genBin_(const std::vector<float>& val)
{
    int bin = -999;
        bin = generatorDistribution_->GetGlobalBinNumber(val.at(0),val.at(1));
    return bin;
}



int PlotterForTUnfold::recoBin_(const std::vector<float>& val)
{
    int bin = -999;
        bin = detectorDistribution_->GetGlobalBinNumber(val.at(0),val.at(1));
    return bin;
}



void PlotterForTUnfold::performPoissonSmearing(TH1* histIn,TH1*& histOut)
{
    for(int binIndex = 0; binIndex <= histIn->GetNbinsX()+1; ++binIndex) { // "0" and "+1" needed for under/overflow bins
        double binContent = histIn->GetBinContent(binIndex);
        double newContent = r3_->PoissonD(binContent);
        histOut->SetBinContent(binIndex, newContent);
        histOut->SetBinError(binIndex, sqrt(newContent));
    }
}


void PlotterForTUnfold::runUnfolding(const TH2* histMigration,const TH1* histInput,
                           const TH1* histBgr,const TH1* histBgrUo)
{

    printf("runUnfolding:\n \tUnfolding Options: \n");
    printf("\tisTauTest: %d\n",isTauTest_);
    printf("\tisCT: %d\n",isCT_);
    printf("\tregMethod: %s\n",regMethod_.Data());
    printf("\tsystematic: %s\n",systematicName_.Data()); 
    
    TUnfold::EConstraint constraintMode= TUnfold::kEConstraintArea;///check before start
    //TUnfold::EConstraint constraintMode= TUnfold::kEConstraintNone;
    
    TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
    //TUnfold::ERegMode regMode = TUnfold::kRegModeDerivative;
    //TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
    //TUnfold::ERegMode regMode = TUnfold::kRegModeNone;
    
    //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
    TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeNone;
    
    // detailed steering for regularisation
    const char *REGULARISATION_DISTRIBUTION=0;
    const char *REGULARISATION_AXISSTEERING=0;
    //detectorBinning_->PrintStream(std::cout);
    //generatorBinning_->PrintStream(std::cout);
    
    // check statistics in input detector distribution
    for(int b = 1; b <= histInput->GetNbinsX(); b++)
      if(histInput->GetBinContent(b) < 20)
        printf("Warning in PlotterForTUnfold::runUnfolding(): low statistics in bin[%d] = %f\n", b, histInput->GetBinContent(b));
    //printf("OZ Starting TUnfold\n");
    //printf("histMigration\n");
    //histMigration->Print("all");
    //printf("histInput\n");
    //histInput->Print("all");
    TUnfoldDensity unfold(histMigration,TUnfold::kHistMapOutputHoriz,
                        regMode,constraintMode,densityFlags,
                        generatorBinning_,detectorBinning_,
                        REGULARISATION_DISTRIBUTION,
                        REGULARISATION_AXISSTEERING);
    
    // special test to use input covariance matrix
    TH2D *inputEmatrix = detectorBinning_->CreateErrorMatrixHistogram("input_covar",true);
    for(int i=1;i<=inputEmatrix->GetNbinsX();i++) {
       Double_t e=histInput->GetBinError(i);
       inputEmatrix->SetBinContent(i,i,e*e);
       // test: non-zero covariance where variance is zero
       //   if(e<=0.) inputEmatrix->SetBinContent(i,i+1,1.0);
    }
    
    if(isPDF_)tau_=0;
    unfold.SetInput(histInput,1.0*(!isPDF_),0.0,inputEmatrix);
    
   // unfold.SetInput(histInput,0.0,0.0,inputEmatrix);
    
    if(histBgr){unfold.SubtractBackground(histBgr,"bgrSum");}
    if(histBgrUo){unfold.SubtractBackground(histBgrUo,"bgrUo");}
    
    
    Int_t nScan=50;
 TSpline *logTauX,*logTauY;
 TGraph *lCurve;
 Int_t iBest=0;
 if(isTauTest_){
 
    if(regMethod_ == "rho"){
        iBest=unfold.ScanTau(nScan,0.,0.,0,TUnfoldDensity::kEScanTauRhoAvg,0,0,&lCurve,&logTauX,&logTauY);
        histTauLogRhoScan_->Fill(TMath::Log10(unfold.GetTau()));
        }
    else if(regMethod_ == "l"){
        iBest=unfold.ScanLcurve(nScan,0.,0.,&lCurve,&logTauX,&logTauY);
        histTauLogLCurve_->Fill(TMath::Log10(unfold.GetTau()));
    }
    else if(regMethod_ == "0"){
        unfold.DoUnfold(0);
    }
    else if(regMethod_=="fix"){
        unfold.DoUnfold(tau_);
    }
 
    logTau_=TMath::Log10(unfold.GetTau());
    histTauLog_->Fill(TMath::Log10(unfold.GetTau()));
     
    
    v_covX_.clear();
    vX_.clear();
    vXmc_.clear();
    if(unfoldedData_)unfoldedData_->Delete();
    unfoldedData_  = unfold.GetOutput("unfoldedData_","unfolding result");
    for(int j=1;j<=unfoldedData_->GetNbinsY();j++){
        for(int i=1;i<=unfoldedData_->GetNbinsX();i++){
            vX_.push_back(unfoldedData_->GetBinContent(i,j));
            vXmc_.push_back(((TH2*)histGen_)->GetBinContent(i,j));
        }
    }
    TH2* t_Ematrix_2d = unfold.GetEmatrixTotal("t_Ematrix_2d", "error matrix (t)");
    for(int j=1;j<=t_Ematrix_2d->GetNbinsY();j++){
        for(int i=1;i<=t_Ematrix_2d->GetNbinsX();i++){
            v_covX_.push_back(t_Ematrix_2d->GetBinContent(i,j));
        }
    }
    t_Ematrix_2d->Delete();

 }
 else if(isCT_){
  iBest=unfold.ScanTau(nScan,0.,0.,0,TUnfoldDensity::kEScanTauRhoAvg,0,0,&lCurve,&logTauX,&logTauY);
 }
 else{
     
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = utils::setLegend();
     
    if(regMethod_ == "rho"){
        ///rho-scan
        TSpline  *logTauRho;
        iBest=unfold.ScanTau(nScan,0.0001,0.01,&logTauRho,TUnfoldDensity::kEScanTauRhoAvg,0,0,&lCurve,&logTauX,&logTauY); 
        logTauRho->SetTitle("Global correlation coef. vs lg(#tau);lg(#tau);Correlation coef.");
        logTauRho->Draw();
        double tempTitleofset = gStyle->GetTitleYOffset();
        gStyle->SetTitleYOffset(1.4);
        gPad->SetLeftMargin(0.2);
        Double_t t[1],rho[1];
        logTauRho->GetKnot(iBest,t[0],rho[0]);
        TGraph *bestRhoLogTau=new TGraph(1,t,rho);
        bestRhoLogTau->SetMarkerColor(kRed);
        bestRhoLogTau->Draw("*p,same");
        //logTauRho->Draw("same");
        writeCanvas(canvas,"logTauRho");
        logTauRho->Delete();
        bestRhoLogTau->Delete();
        gStyle->SetTitleYOffset(tempTitleofset);
    }
    else if(regMethod_ == "l"){
        ///L-curve
        iBest=unfold.ScanLcurve(nScan,0.,0.,&lCurve,&logTauX,&logTauY);
        Double_t x[1],y[1];
        lCurve->GetPoint(iBest,x[0],y[0]);
        TGraph *bestLCurve=new TGraph(1,x,y);
        //lCurve->SetTitle("L curve;log((y-Ax)^{T}V_{y}^{-1}(y-Ax));log((x-x0)L^{T}L(x-x0))");
        lCurve->SetTitle("L curve;log(L_{2}/#tau^{2});log(L_{1})");
        lCurve->GetYaxis()->SetTitleOffset(1.3);
        lCurve->GetYaxis()->SetTimeOffset(1.5);
        lCurve->GetXaxis()->SetTimeOffset(1.1);
        lCurve->Draw("AL");
        bestLCurve->SetMarkerColor(kRed);
        bestLCurve->Draw("*");
        writeCanvas(canvas,"Lcurve");
        lCurve->Delete();
        bestLCurve->Delete();
    }
    else if(regMethod_ == "fix"){
        unfold.DoUnfold(tau_);
    }
    else if(regMethod_ == "0"){
        unfold.DoUnfold(0);
    }
    
   // Probability matrix //TH2*
    TH2* probabilityMatrix = unfold.GetProbabilityMatrix("probabilityMatrix","probability matrix",kTRUE);
    probabilityMatrix->SetStats(0);
    probabilityMatrix->Draw("colz");
    //probabilityMatrix->Print("all");
    writeCanvas(canvas,probabilityMatrix->GetName());
    
    // rho ij total //TH2*
    gStyle->SetPaintTextFormat("0.2f");
    rhoIJtotal_ = unfold.GetRhoIJtotal("rhoIJtotal","rho ij total matrix",0,0,kTRUE);
    rhoIJtotal_->SetStats(0);
    rhoIJtotal_->Draw("colztext");
    writeCanvas(canvas,rhoIJtotal_->GetName());
    gStyle->SetPaintTextFormat("g");
    
    
    
    // Error matrices info
    TH1* t_Ematrix = generatorBinning_->CreateHistogram("t_Ematrix");//total error
    TH1* i_Ematrix = generatorBinning_->CreateHistogram("i_Ematrix");//input error
    TH1* bgr_Ematrix = generatorBinning_->CreateHistogram("bgr_Ematrix");//bgr error
    TH1* uo_Ematrix = generatorBinning_->CreateHistogram("uo_Ematrix");//uo error
    TH1* migrat_Ematrix = generatorBinning_->CreateHistogram("migrat_Ematrix");
    TH2* t_Ematrix_2d = unfold.GetEmatrixTotal("t_Ematrix_2d", "error matrix (t)");
    TH2* i_Ematrix_2d = unfold.GetEmatrixInput("i_Ematrix_2d", "error matrix (i)");
    TH2* bgr_Ematrix_2d = unfold.GetEmatrixSysBackgroundUncorr("bgrSum","bgr_Ematrix_2d", "error matrix (bgr)");
    TH2* uo_Ematrix_2d = unfold.GetEmatrixSysBackgroundUncorr("bgrUo","uo_Ematrix_2d", "error matrix (uo)");
    TH2* migrat_Ematrix_2d = unfold.GetEmatrixSysUncorr("migrat_Ematrix_2d", "error matrix (migrat)");
    
    TH1* unfoldedData = unfold.GetOutput("unfoldedData","unfolding result",0,0,kFALSE);
    //printf("OZ unfoldedData\n");
    //unfoldedData->Print("all");
    for(Int_t i=1;i<=i_Ematrix->GetNbinsX();i++) 
    {
        double data = unfoldedData->GetBinContent(i);
        if(data>=0){
            t_Ematrix->SetBinContent(i, TMath::Sqrt(t_Ematrix_2d->GetBinContent(i,i))/data );
            i_Ematrix->SetBinContent(i, TMath::Sqrt(i_Ematrix_2d->GetBinContent(i,i))/data );
            bgr_Ematrix->SetBinContent(i, TMath::Sqrt(bgr_Ematrix_2d->GetBinContent(i,i))/data );
            uo_Ematrix->SetBinContent(i, TMath::Sqrt(uo_Ematrix_2d->GetBinContent(i,i))/data );
            migrat_Ematrix->SetBinContent(i, TMath::Sqrt(migrat_Ematrix_2d->GetBinContent(i,i))/data );
        }
    }

    t_Ematrix->SetStats(0);
    t_Ematrix->GetYaxis()->SetRangeUser(-0.01,0.4);
    
    t_Ematrix->SetMarkerStyle(20);
    i_Ematrix->SetMarkerStyle(21);
    bgr_Ematrix->SetMarkerStyle(22);
    uo_Ematrix->SetMarkerStyle(23);
    migrat_Ematrix->SetMarkerStyle(33);
    
    t_Ematrix->SetMarkerColor(kBlack);
    i_Ematrix->SetMarkerColor(kRed);
    bgr_Ematrix->SetMarkerColor(kGreen);
    uo_Ematrix->SetMarkerColor(kBlue);
    migrat_Ematrix->SetMarkerColor(kMagenta-3);
    
    t_Ematrix->SetMarkerSize(0.6);
    i_Ematrix->SetMarkerSize(0.6);
    bgr_Ematrix->SetMarkerSize(0.6);
    uo_Ematrix->SetMarkerSize(0.6);
    migrat_Ematrix->SetMarkerSize(0.6);
    
    TH1* hTitleEmatrix =  plotNiceTitel(0,0.4,"rel. unc.");
    t_Ematrix->Draw("same,pe");
    i_Ematrix->Draw("same,pe");
    bgr_Ematrix->Draw("same,pe");
    uo_Ematrix->Draw("same,pe");
    migrat_Ematrix->Draw("same,pe");
    
    
    legend->SetX1(0.1);
    legend->SetX2(0.4);
    legend->Clear();
    legend->AddEntry(t_Ematrix,"Total","pe");
    legend->AddEntry(i_Ematrix,"Input","pe");
    legend->AddEntry(bgr_Ematrix,"Bgr.","pe");
    legend->AddEntry(uo_Ematrix,"Bgr. Uo.","pe");
    legend->AddEntry(migrat_Ematrix,"Migration Mat.","pe");
    
    legend->Draw("same");
    writeCanvas(canvas,"ematricesInfo");
    legend->Clear();
    
    hTitleEmatrix->Delete();
    t_Ematrix->Delete();
    i_Ematrix->Delete();
    bgr_Ematrix->Delete();
    uo_Ematrix->Delete();
    migrat_Ematrix->Delete();
    
    i_Ematrix_2d->Delete();
    bgr_Ematrix_2d->Delete();
    uo_Ematrix_2d->Delete();
    migrat_Ematrix_2d->Delete();
    // ...

    // OZ 16.01.2017 Hack to work with provided histograms, without calling prepareHistograms()
    {
      histMigration_ = TUnfoldBinning::CreateHistogramOfMigrations(generatorBinning_,detectorBinning_,"histMigration");
      
      histBgrUo_ = detectorBinning_->CreateHistogram("histBgrUo");
      histBgrUoScaled_ = detectorBinning_->CreateHistogram("histBgrUoScaled");
      histBgr_ = detectorBinning_->CreateHistogram("histBgr");
      histBgrScaled_ = detectorBinning_->CreateHistogram("histBgrScaled");
      histData_ = detectorBinning_->CreateHistogram("histData");
      histDataScaled_ = detectorBinning_->CreateHistogram("histDataScaled");
      
      histFsigReco_ = detectorBinning_->CreateHistogram("histFsigReco");
      histTtSigReco_ = detectorBinning_->CreateHistogram("histTtSigReco");
      histTtAllReco_ = detectorBinning_->CreateHistogram("histTtAllReco");
      histAllReco_ = detectorBinning_->CreateHistogram("histAllReco");
      histAllRecoScaled_ = detectorBinning_->CreateHistogram("histAllRecoScaled");
    
      // OZ 17.01.2017 these two histograms are provided separately via SetHistGen()
      //histGen_ = generatorBinning_->CreateHistogram("histGen",kTRUE);
      //histGenNorm_ = generatorBinning_->CreateHistogram("histGenNorm",kTRUE);
        
      histEffAllBins_ = generatorBinning_->CreateHistogram("histEffAllBins");
      histPurityAllBins_ = generatorBinning_->CreateHistogram("histPurityAllBins");
      histStabilityAllBins_ = generatorBinning_->CreateHistogram("histStabilityAllBins");
      histRecoGenAllBins_ = generatorBinning_->CreateHistogram("histRecoGenAllBins");
    }
        
    //Draw input correlation matrix
    // OZ 19.01.2017 it is empty
    /*if(histMigration_)
    {
      histMigration_->SetStats(0);
      histMigration_->Draw("colztext");
      writeCanvas(canvas,histMigration_->GetName());
    }*/
    
    //Delete objects
    unfoldedData->Delete();
    
    
    for(int j=1;j<=t_Ematrix_2d->GetNbinsY();j++){
        for(int i=1;i<=t_Ematrix_2d->GetNbinsX();i++){
            v_covX_.push_back(t_Ematrix_2d->GetBinContent(i,j));
        }
    }
    
    unfoldedData_  = unfold.GetOutput("unfoldedData_","unfolding result");
    for(int j=1;j<=unfoldedData_->GetNbinsY();j++){
        for(int i=1;i<=unfoldedData_->GetNbinsX();i++){
            float x = unfoldedData_->GetXaxis()->GetBinCenter(i);
            float y = unfoldedData_->GetYaxis()->GetBinCenter(j);
            std::vector<float> v;
            v.push_back(x);
            v.push_back(y);
            int bin = genBin_(v);
            double error = sqrt(t_Ematrix_2d->GetBinContent(bin,bin));
            unfoldedData_->SetBinError(i,j,error);
            vX_.push_back(unfoldedData_->GetBinContent(i,j));
        }
    }
    
        histGenNorm_ = (TH1*)histGen_->Clone();
        histGenNorm_->Scale(1./((TH2*)histGen_)->Integral());
    
    double integral = unfoldedData_->Integral();
    unfoldedDataNorm_ = unfold.GetOutput("unfoldedDataNorm_","unfolding result");
    unfoldedDataNorm_->Scale(1./integral);
    
    std::vector<double> vG;
    for(int j=0;j<(int)vX_.size();j++){
        for(int i=0;i<(int)vX_.size();i++){
            if(i==j)vG.push_back((integral-vX_.at(j))/(integral*integral));
            else vG.push_back(-vX_.at(j)/(integral*integral));
        }
    }
    TMatrixD matrixG((int)vX_.size(),(int)vX_.size());
    matrixG.SetMatrixArray(vG.data());
    TMatrixD matrixGT = matrixG;
    matrixGT.T();
    TMatrixD covMatrix((int)vX_.size(),(int)vX_.size());
    covMatrix.SetMatrixArray(v_covX_.data());
    TMatrixDSparse covMatrixNorm = matrixG*covMatrix*matrixGT;
    for(int j=1;j<=unfoldedDataNorm_->GetNbinsY();j++){
        for(int i=1;i<=unfoldedDataNorm_->GetNbinsX();i++){
            float x = unfoldedDataNorm_->GetXaxis()->GetBinCenter(i);
            float y = unfoldedDataNorm_->GetYaxis()->GetBinCenter(j);
            std::vector<float> v;
            v.push_back(x);
            v.push_back(y);
            int bin = genBin_(v);
            double error = sqrt(covMatrixNorm(bin-1,bin-1));
            unfoldedDataNorm_->SetBinError(i,j,error);
            vXnorm_.push_back(unfoldedDataNorm_->GetBinContent(i,j));
            vXnormMC_.push_back(((TH2*)histGenNorm_)->GetBinContent(i,j));
            vXmc_.push_back(((TH2*)histGen_)->GetBinContent(i,j));
        }
    }

    
    for(int j=0;j<(int)vX_.size();j++){
        for(int i=0;i<(int)vX_.size();i++){
            v_covXnorm_.push_back(covMatrixNorm(i,j));
        }
    }
    
    TMatrixD probMatrix(probabilityMatrix->GetNbinsY(),probabilityMatrix->GetNbinsX());
    TMatrixD probMatrixErr(probabilityMatrix->GetNbinsY(),probabilityMatrix->GetNbinsX());
      for(int j=1;j<=probabilityMatrix->GetNbinsY();j++){
        for(int i=1;i<=probabilityMatrix->GetNbinsX();i++){
            probMatrix(j-1,i-1)=probabilityMatrix->GetBinContent(i,j);
            probMatrixErr(j-1,i-1)=probabilityMatrix->GetBinError(i,j);
        }
    }
    TMatrixD vXmc((int)vXmc_.size(),1);
    vXmc.SetMatrixArray(vXmc_.data());
    TMatrixD vYmc = (probMatrix*vXmc);
    TMatrixD vYmcErr = (probMatrixErr*vXmc);
    TH1* histInputNorm = (TH1*)histInput->Clone();
    if(histBgr)histInputNorm->Add(histBgr,-1);
    if(histBgrUo)histInputNorm->Add(histBgrUo,-1);
    double histInputNormIntegral=histInputNorm->Integral();
    histInputNorm->Scale(1.0/histInputNormIntegral);
    TH1* histInputNormMC =(TH1*)histInput->Clone();
    for(int i=1;i<=histInputNormMC->GetNbinsX();i++){
        histInputNormMC->SetBinContent(i,vYmc(i-1,0));
        histInputNormMC->SetBinError(i,vYmcErr(i-1,0));
    }
    double histInputNormMCIntegral=histInputNormMC->Integral();
    histInputNormMC->Scale(1.0/histInputNormMCIntegral);
    for(int i=1;i<=histInput->GetNbinsX();i++){
        double histBgrUo_BinContent=0;
        if(histBgrUo)histBgrUo_BinContent=histBgrUo->GetBinContent(i);
        v_dY_.push_back(histInput->GetBinContent(i)-(histBgr ? histBgr->GetBinContent(i) : 0.0)-histBgrUo_BinContent-vYmc(i-1,0));
        double a=histInput->GetBinError(i);
        double b=(histBgr ? histBgr->GetBinError(i) : 0.0);
        double c=0;if(histBgrUo)c=histBgrUo->GetBinError(i);
        v_covY_.push_back(a*a+b*b+c*c+vYmcErr(i-1,0)*vYmcErr(i-1,0));
        v_dYnorm_.push_back(histInputNorm->GetBinContent(i)-histInputNormMC->GetBinContent(i));
        double d=histInputNorm->GetBinError(i);
        double e=histInputNormMC->GetBinError(i);
        v_covYnorm_.push_back(d*d+e*e);
        vYmc_.push_back(vYmc(i-1,0));
        vYnormMC_.push_back(histInputNormMC->GetBinContent(i));
    }
    
    histInputNorm->Delete();
    histInputNormMC->Delete();
    
    probabilityMatrix->Delete();
    t_Ematrix_2d->Delete();
    
    
    TH2D* histCovX = new TH2D("histCovX","Correlation Matrix; #bf{x} bins ; #bf{x} bins",(int)vX_.size(),0.5,(int)vX_.size()+0.5,(int)vX_.size(),0.5,(int)vX_.size()+0.5);
    histCovX->SetStats(0);
    for(int j=1;j<=histCovX->GetNbinsX();j++){
        for(int i=1;i<=histCovX->GetNbinsY();i++){
            int N = (int)vX_.size();
            histCovX->SetBinContent(i,j,(v_covX_.at((i-1)+(j-1)*N)/sqrt(v_covX_.at((i-1)+(i-1)*N)*v_covX_.at((j-1)+(j-1)*N))));
        }
    }
    gStyle->SetPaintTextFormat("0.2f");
    canvas->SetRightMargin(0.12);
    canvas->SetLeftMargin(0.08);
    histCovX->Draw("colztext");
    writeCanvas(canvas,histCovX->GetName());
    
    
    
    TH2D* histCovXnorm = new TH2D("histCovXnorm","Correlation Matrix; #bf{x} bins ; #bf{x} bins",(int)vX_.size(),0.5,(int)vX_.size()+0.5,(int)vX_.size(),0.5,(int)vX_.size()+0.5);
    histCovXnorm->SetContent(v_covXnorm_.data());
    histCovXnorm->SetStats(0);
    for(int i=1;i<=histCovXnorm->GetNbinsX();i++){
        for(int j=1;j<=histCovXnorm->GetNbinsY();j++){
            int N = (int)vX_.size();
            histCovXnorm->SetBinContent(i,j,(v_covXnorm_.at((i-1)+(j-1)*N)/sqrt(v_covXnorm_.at((i-1)+(i-1)*N)*v_covXnorm_.at((j-1)+(j-1)*N))));
        }
    }
    gStyle->SetPaintTextFormat("0.2f");
    histCovXnorm->Draw("colztext");
    writeCanvas(canvas,histCovXnorm->GetName());
    gStyle->SetPaintTextFormat("g");
    
    
    //Delete objects
    delete canvas;
    delete legend;
    
    
 }
 
   inputEmatrix->Delete();
   
  std::cout<<"chi**2="<<unfold.GetChi2A()<<"+"<<unfold.GetChi2L()
      <<" / "<<unfold.GetNdf()<<"\n";
    
    ///if(histUnfGenBin_)histUnfGenBin_->Delete(); ///check before start ///uncomment for toy experiments
    histUnfGenBin_ = unfold.GetOutput("histUnfGenBin_","unfolding result",0,0,kFALSE);
    
    if(systematicName_ == "Nominal" && !isCT_ && !isTauTest_ && !isVar_){ ///check before start
        Output saveTau("tau");
            saveTau.add("tau",std::vector<double>{unfold.GetTau()});
            saveTau.add("corr",std::vector<double>{unfold.GetScanVariable(TUnfoldDensity::kEScanTauRhoAvg,0,0)});
            saveTau.add("prt",vectorPurityAllBins_);
            saveTau.add("stb",vectorStabilityAllBins_);
            
            saveTau.save(plotsFolder_+ "/optimalTau.txt");
    }
    
}



void PlotterForTUnfold::writeCanvas( TCanvas* canvas, const TString name )
{
    canvas->Print(plotsFolder_+ name + "_" + v_plotName_.at(0) + "_vs_" + v_plotName_.at(1)  + ".pdf");
    canvas->Print(plotsFolder_+ name + "_" + v_plotName_.at(0) + "_vs_" + v_plotName_.at(1)  + ".png");
    canvas->Print(plotsFolder_+ name + "_" + v_plotName_.at(0) + "_vs_" + v_plotName_.at(1)  + ".root");
    canvas->Clear();
}


void PlotterForTUnfold::writePlotCP(const std::vector<Sample>& v_sample ,const std::vector<TH1D* >& v_SampleHist,const int& ind, const int binNum)
{
    int jnd = int((ind==0)?1:0);
    
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = NULL;
    TLegend* legend2 = NULL;
    std::cout << "test 1, " << std::endl;
    if(v_leg2_.at(ind)==1)legend = utils::setLegend(0.62,0.48,0.88,0.84);
    else if(v_leg2_.at(ind)==2){
        legend = utils::setLegend(0.51,0.55,0.77,0.91);
        legend2 = utils::setLegend(0.71,0.55,0.96,0.91);
        //legend2->AddEntry((TObject*)0, "", "");
        
    }
    canvas->cd();
    
    TH1D* histData = 0;

    //Fill data
    histData = (TH1D*)v_SampleHist.at(0)->Clone("name");
    for(size_t iSample = 1; iSample < v_sample.size(); ++iSample){
        if(v_sample.at(iSample).sampleType() == Sample::data){
            histData->Add((TH1D*)v_SampleHist.at(iSample)->Clone("name"));
        }
    }
    legend->AddEntry(histData,v_sample.at(0).legendEntry(),"pe");
    // ... //
    
    //Fill MC stack
    THStack* stack = new THStack();
    TH1D* histPerStackEntry = 0;
    for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
        if(v_sample.at(iSample).sampleType() == Sample::data)continue;
        TH1D* hist = (TH1D*)v_SampleHist.at(iSample)->Clone("name");
            hist->SetTitle(v_sample.at(iSample).legendEntry());
        if(!histPerStackEntry){
            histPerStackEntry = hist;
        }
        else if(v_sample.at(iSample-1).legendEntry()==v_sample.at(iSample).legendEntry())
        {
            histPerStackEntry->Add(hist);
        }
        else if(v_sample.at(iSample-1).legendEntry()!=v_sample.at(iSample).legendEntry()){
            histPerStackEntry->SetLineColor(1);
            stack->Add(histPerStackEntry);
            histPerStackEntry = hist;
        }
        if(iSample==v_sample.size()-1){
            histPerStackEntry->SetLineColor(1);
            stack->Add(histPerStackEntry);
        }
    }
    TList* listOfHistos = stack->GetHists();
    TIter next(listOfHistos,kIterBackward);
    TObject* object = 0;
    Int_t plotInSecondLeg = 0;
        while ((object = next()))
        {
            if(((TString)((TH1D*)object)->GetTitle())=="W+Jets")plotInSecondLeg = 1;
            if(plotInSecondLeg&&(v_leg2_.at(ind)==2)){
                legend2->AddEntry(((TH1D*)object),((TH1D*)object)->GetTitle(),"f");
                legend->AddEntry((TObject*)0, "", "");
            }
            else legend->AddEntry(((TH1D*)object),((TH1D*)object)->GetTitle(),"f");
            
        }
    plotInSecondLeg =  0;
    //Title
    TString plotTitle=systematicName_;
            plotTitle="";
    if(binNum == -1 ) plotTitle = v_plotTitle_.at(jnd) + " underflow bin" ;
    if(binNum == -2 ) plotTitle = v_plotTitle_.at(jnd) + " overflow bin" ;
    if(binNum >=0) plotTitle = utils::makeBinTitle(v_plotTitle_.at(jnd),v_coarseBins_.at(jnd).at(binNum),v_coarseBins_.at(jnd).at(binNum+1));
    if(binNum == -1 || binNum == -2)plotTitle = plotTitle + ";" + v_plotTitle_.at(ind) + ", " + v_plotUnits_.at(ind) + ";";
    else {
        if(v_plotUnits_.at(ind)!=""){
            plotTitle = plotTitle  + ";" + v_plotTitle_.at(ind) + ", " + v_plotUnits_.at(ind) + ";" + v_plotYunitscp_.at(ind);
        }
        else{
            plotTitle = plotTitle  + ";" + v_plotTitle_.at(ind) + " " + v_plotUnits_.at(ind) + ";" + v_plotYunitscp_.at(ind);
        }
    }
    
    //Draw
    TH1* stacksum = common::summedStackHisto(stack);
    if(v_ymax_.at(ind)<0)histData->SetAxisRange(0,( histData->GetMaximum() > stacksum->GetMaximum() ? histData->GetMaximum()*1.2 : stacksum->GetMaximum()*1.2 ),"Y");
    else histData->SetAxisRange(0,v_ymax_.at(ind),"Y");
    histData->SetStats(0);
    histData->SetTitle(plotTitle);
    histData->SetTitleOffset(1.65,"Y");
    histData->SetTitleSize(histData->GetTitleSize("Y")*1.4,"Y");
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.16);
    if(v_plotYlog_.at(ind)){
        
        histData->SetAxisRange(1,histData->GetBinContent(histData->GetMaximumBin())*1.3,"Y");
        gPad->SetLogy();
    }
    histData->Draw("e1");
    stack->Draw("same HIST");
    histData->Draw("e1 same");
    if(v_plotTitle_.at(ind)!="y^{t}"&&v_plotTitle_.at(ind)!="y^{tt}"&&v_plotYlog_.at(ind)==0)legend->Draw("same");
    legend->Draw("same");
    
    TPaveText *pt = new TPaveText(0.2, 0.945, 0.3, 1.0, "NDC");
        pt->SetFillColor( 0 );
        pt->SetFillStyle( 0 );
        pt->SetBorderSize( 0 );
        pt->SetTextSize(0.05);
        pt->AddText("CMS");
        pt->Draw("same");

        pt = new TPaveText(0.25, 0.95, 1., 1.0, "NDC");
        pt->SetFillColor( 0 );
        pt->SetFillStyle( 0 );
        pt->SetBorderSize( 0 );
        pt->SetTextSize(0.036);
        pt->AddText("                                     19.7 fb^{-1} (8 TeV)");
        pt->Draw("same");
    
    std::cout << "test 2, " << std::endl;
    
    //UsefulTools::DrawCMSLabels(lumi_,1);
    UsefulTools::DrawDecayChLabel("e#mu");
    //utils::drawRatio(histData, stacksum, NULL,0.5,1.5);
    
    if((!v_eban_.at(ind))&&(binNum==-999))utils::drawRatio2d(histData, stacksum, 0.5,1.5); 
    else if(binNum!=-999)utils::drawRatio2d(histData, stacksum, 0.5,1.5); 

    if(binNum==-999){
        Output cpInfo("xsec");
        std::vector<double> binContents;
        std::vector<double> binContentsErrors;
        for(int i=1;i<=stacksum->GetNbinsX();i++){
            binContents.push_back(stacksum->GetBinContent(i));
            binContentsErrors.push_back(stacksum->GetBinError(i));
        }
        cpInfo.add("cp",binContents);
        cpInfo.add("error",binContentsErrors);
        cpInfo.save(plotsFolder_+ "CP_" + v_plotName_.at(ind) + ".txt");
    }
    
        std::cout << "test 3, " << std::endl;

    
    if(v_eban_.at(ind)&&systematicName_=="Nominal")
    {
        
        std::vector<TString> v_systematic{
            
            "Nominal",
            
            "BG_UP","BTAG_ETA_UP","BTAG_LJET_ETA_UP","BTAG_LJET_PT_UP","BTAG_LJET_UP","BTAG_PT_UP","BTAG_UP","DY_UP",
            "JER_UP","JES_UP","KIN_UP","LEPT_UP","LUMI_UP","MASS_UP","MATCH_UP","PU_UP","SCALE_UP","TRIG_UP",
            "MCATNLO","POWHEG","POWHEGHERWIG",
            "PDF_0_CENTRAL","PDF_10_UP","PDF_11_UP","PDF_12_UP",
            "PDF_13_UP","PDF_14_UP","PDF_15_UP","PDF_16_UP",
            "PDF_17_UP","PDF_18_UP","PDF_19_UP","PDF_1_UP",
            "PDF_20_UP","PDF_21_UP","PDF_22_UP","PDF_23_UP","PDF_24_UP",
            "PDF_25_UP","PDF_26_UP","PDF_2_UP","PDF_3_UP","PDF_4_UP",
            "PDF_5_UP","PDF_6_UP","PDF_7_UP","PDF_8_UP","PDF_9_UP"
            };
            
                std::cout << "test 4, " << std::endl;

            std::vector<double> vXsecNominal;
            std::vector<double> vXsecNominalErrors;
            std::vector<double> vXsecNorm;
            std::vector<double> vXsecPdf;
            std::vector<double> vXsecPOWHEG;
            std::vector<double> vXsecPOWHEGHERWIG;
                            
            
            std::vector<double> vDiffUp;
            std::vector<double> vDiffDown;
            std::vector<double> vDiffHadrUp;
            std::vector<double> vDiffHadrDown;
            std::vector<double> vDiffHardUp;
            std::vector<double> vDiffHardDown;
            
            
            std::vector<double> vSystErrPos;
            std::vector<double> vSystErrNeg;
            std::vector<double> vSystErrAbsPos;
            std::vector<double> vSystErrAbsNeg;
            
            for(auto syst: v_systematic){ //systematic loop
                TString file = plotsFolder_+ "CP_" + v_plotName_.at(ind) + ".txt";
                file.ReplaceAll("Nominal",syst);
                std::ifstream tempStream(file.Data(), std::ifstream::in);
                if (!tempStream.good()) {
                    std::cerr<<"Error in PlotterForTUnfold::writePlotCP! Cannot find file with name: "<< file <<"\nrun: ./install/bin/Histo2 -c emu -p p -s " << syst <<  "\n"<<std::endl;
                }
                else{
                    
                    std::cout << "test 4.0, " << syst << std::endl;
                    
                    if(syst == "Nominal"){
                        utils::readLineToVector(file,"cp",vXsecNominal);
                        utils::readLineToVector(file,"error",vXsecNominalErrors);
                        vXsecNorm = vXsecNominal;
                        for(int i=0;i<(int)vXsecNominal.size();i++){
                            vSystErrPos.push_back(0);
                            vSystErrNeg.push_back(0);
                            vSystErrAbsPos.push_back(0);
                            vSystErrAbsNeg.push_back(0);
                        }
                        continue;
                    }
                    
                    std::vector<double> vSyst;
                    utils::readLineToVector(file,"cp",vSyst);
                    
                    std::cout << "test 4.01,"<<  std::endl;
                        if(syst == "POWHEG"){
                            utils::readLineToVector(file,"cp",vXsecPOWHEG);
                            continue;
                        }
                        if(syst == "POWHEGHERWIG"){
                            utils::readLineToVector(file,"cp",vXsecPOWHEGHERWIG);
                            continue;
                        }
                        
                        if(syst == "PDF_0_CENTRAL"){
                            utils::readLineToVector(file,"cp",vXsecPdf); 
                            continue;
                        }
                        double scale =1.0;
                        if(syst.Contains("PDF_")){
                            vXsecNorm = vXsecPdf;
                            scale = 1.0/1.645;
                        }
                    std::cout << "test 4.02,"<<  std::endl;
                    if(syst.Contains("MASS_"))scale=1.0/6;
                    
                    if(syst.Contains("_UP")){
                        std::cout << "test 4.0211,"<< vSyst.size() << " " << vXsecNorm.size() <<   std::endl;
                         vDiffUp = utils::divideVect(utils::diffVect(vSyst,vXsecNorm,scale),vXsecNorm,100);
                         std::cout << "test 4.021111111,"<<  std::endl;
                    }
                    std::cout << "test 4.0212,"<<  std::endl;
                    if(syst.Contains("_DOWN")){
                        std::cout << "test 4.021,"<<  std::endl;
                        vDiffDown = utils::divideVect(utils::diffVect(vSyst,vXsecNorm,scale),vXsecNorm,100);
                    }
                    std::cout << "test 4.03,"<<  std::endl;
                    if( vDiffUp.size()>0 && vDiffDown.size()>0 )
                    {
                        utils::fillSysUpDown(vDiffUp,vDiffDown,vSystErrPos,vSystErrNeg);
                        vDiffUp.clear();
                        vDiffDown.clear();
                    }
                    std::cout << "test 4.04,"<<  std::endl;
                }
                
            } //systematic loop
            
                    std::cout << "test 4.1, " << std::endl;

            
            vDiffHadrUp = utils::divideVect(utils::diffVect(vXsecPOWHEG,vXsecPOWHEGHERWIG),vXsecPOWHEG,100);
            vDiffHadrDown = utils::divideVect(utils::diffVect(vXsecPOWHEGHERWIG,vXsecPOWHEG),vXsecPOWHEG,100);
            utils::fillSysUpDown(vDiffHadrUp,vDiffHadrDown,vSystErrPos,vSystErrNeg);
            
            vDiffHardUp = utils::divideVect(utils::diffVect(vXsecNominal,vXsecPOWHEG),vXsecNominal,100);
            vDiffHardDown = utils::divideVect(utils::diffVect(vXsecPOWHEG,vXsecNominal),vXsecNominal,100);
            utils::fillSysUpDown(vDiffHardUp,vDiffHardDown,vSystErrPos,vSystErrNeg);
            
            std::vector<double> statErrors;
            statErrors = utils::divideVect(vXsecNominalErrors,vXsecNominal,100);
            for(int j=0;j<(int)vSystErrNeg.size();j++){
                vSystErrNeg.at(j) = sqrt( pow(vSystErrNeg.at(j),2) + pow(statErrors.at(j),2));
                vSystErrPos.at(j) = sqrt( pow(vSystErrPos.at(j),2) + pow(statErrors.at(j),2));
            }
            
            vSystErrAbsNeg=utils::multiplyVect(vSystErrNeg,vXsecNominal,0.01);
            vSystErrAbsPos=utils::multiplyVect(vSystErrPos,vXsecNominal,0.01);
            //utils::setPrecision(vSystErrNeg,1);
            //utils::setPrecision(vSystErrPos,1);
            
            std::vector<double> vx;
            std::vector<double> vx_error;
            std::vector<double> vy;
            std::vector<double> vy_norm;
            std::vector<double> v_zero;
            for(int i=1;i<=histData->GetNbinsX();i++){
                vx.push_back(histData->GetBinCenter(i));
                vy.push_back(stacksum->GetBinContent(i));
                vy_norm.push_back(1);
                v_zero.push_back(0);
                vx_error.push_back(histData->GetBinWidth(1)/2);
                vSystErrNeg.at(i-1)=vSystErrNeg.at(i-1)*0.01;
                vSystErrPos.at(i-1)=vSystErrPos.at(i-1)*0.01;
            }
            
            TGraphAsymmErrors* graph1 = new TGraphAsymmErrors(histData->GetNbinsX(),vx.data(),vy.data(),vx_error.data(),vx_error.data(),vSystErrAbsNeg.data(),vSystErrAbsPos.data());
            TGraphAsymmErrors* graph2 = new TGraphAsymmErrors(histData->GetNbinsX(),vx.data(),vy_norm.data(),vx_error.data(),vx_error.data(),vSystErrPos.data(),vSystErrNeg.data());
            gStyle->SetHatchesSpacing(1.5);
            graph1->SetFillColor(kBlack);
            graph1->SetFillStyle(3154);
            
            graph1->Draw("2");
            graph2->SetFillColor(kBlack);
            graph2->SetFillStyle(3154);
            
            if((v_leg2_.at(ind)==2))legend2->AddEntry(graph1,"Uncertainty","f");
            else legend->AddEntry(graph1,"Uncertainty","f");
            legend->Draw("same");
            if((v_leg2_.at(ind)==2)){
                legend2->AddEntry((TObject*)0, "", "");
                legend2->AddEntry((TObject*)0, "", "");
                legend2->AddEntry((TObject*)0, "", "");
                legend2->AddEntry((TObject*)0, "", "");
                legend2->AddEntry((TObject*)0, "", "");
                legend2->AddEntry((TObject*)0, "", "");
                legend2->AddEntry((TObject*)0, "", "");
                legend2->Draw("same");
            }
            
            utils::drawRatio2d_Eban(histData, stacksum, 0.5,1.5,"cp",graph1,graph2);
            
    }
    
    
    
    //Print
    if(binNum!=-999){
        if(binNum>=0)canvas->Print(plotsFolder_+ "CP_" + v_plotName_.at(ind) + "_IN_" + v_plotName_.at(jnd)  + "_" + std::to_string(binNum) + ".pdf");
        if(binNum==-1)canvas->Print(plotsFolder_+ "CP_" + v_plotName_.at(ind) + "_IN_" + v_plotName_.at(jnd)  + "_" + "u" + ".pdf");
        if(binNum==-2)canvas->Print(plotsFolder_+ "CP_" + v_plotName_.at(ind) + "_IN_" + v_plotName_.at(jnd)  + "_" + "o" + ".pdf");
    }
    if(binNum==-999)canvas->Print(plotsFolder_+ "CP_" + v_plotName_.at(ind) + ".png");
    if(binNum==-999)canvas->Print(plotsFolder_+ "CP_" + v_plotName_.at(ind) + ".pdf");
    if(binNum==-999)canvas->Print(plotsFolder_+ "CP_" + v_plotName_.at(ind) + ".eps");
    if(binNum==-999)canvas->Print(plotsFolder_+ "CP_" + v_plotName_.at(ind) + ".root");
    
        std::cout << "test 5, " << std::endl;

    
    //Delete
    stack->Delete();
    histData->Delete();
    delete canvas;
    delete legend;
    stacksum->Delete();
    
}

void PlotterForTUnfold::writePlotCPAllBins(const std::vector<Sample>& v_sample ,const std::vector<TH1* >& v_SampleHist)
{
    // Prepare canvas and legend
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = utils::setLegend();
    canvas->cd();
    
    //Fill data
    TH1D* histData = (TH1D*)v_SampleHist.at(0)->Clone("name");
    for(size_t iSample = 1; iSample < v_sample.size(); ++iSample)
    {
        if(v_sample.at(iSample).sampleType() == Sample::data){
            histData->Add((TH1D*)v_SampleHist.at(iSample)->Clone("name"));
        }
    }
    legend->AddEntry(histData,v_sample.at(0).legendEntry(),"pe");
    // ... //
    
    
    //Fill MC stack
    THStack* stack = new THStack();
    TH1D* histPerStackEntry = 0;
    for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){
        if(v_sample.at(iSample).sampleType() == Sample::data)continue;
        TH1D* hist = (TH1D*)v_SampleHist.at(iSample)->Clone("name");
            hist->SetTitle(v_sample.at(iSample).legendEntry());
        if(!histPerStackEntry){
            histPerStackEntry = hist;
        }
        else if(v_sample.at(iSample-1).legendEntry()==v_sample.at(iSample).legendEntry())
        {
            histPerStackEntry->Add(hist);
        }
        else if(v_sample.at(iSample-1).legendEntry()!=v_sample.at(iSample).legendEntry()){
            histPerStackEntry->SetLineColor(1);
            stack->Add(histPerStackEntry);
            histPerStackEntry = hist;
        }
        if(iSample==v_sample.size()-1){
            histPerStackEntry->SetLineColor(1);
            stack->Add(histPerStackEntry);
        }
    }
    TList* listOfHistos = stack->GetHists();
    TIter next(listOfHistos,kIterBackward);
    TObject* object = 0;
        while ((object = next()))
        {
            legend->AddEntry(((TH1D*)object),((TH1D*)object)->GetTitle(),"f");
        }
    // ... //
    
    //Draw
    TH1* stacksum = common::summedStackHisto(stack);
    double yMax = ( histData->GetMaximum() > stacksum->GetMaximum() ? histData->GetMaximum()*1.2 : stacksum->GetMaximum()*1.2 );
    
    TH1D* hist =  plotNiceTitel(0,yMax,v_plotYunitscp_.at(0));
    histData->SetTitleOffset(1.75,"Y");
    hist->SetTitleOffset(1.75,"Y");
    hist->SetTitleSize(hist->GetTitleSize("Y")*1.5,"Y");
    histData->SetAxisRange(0,( histData->GetMaximum() > stacksum->GetMaximum() ? histData->GetMaximum()*1.2 : stacksum->GetMaximum()*1.2 ),"Y");
    stack->Draw("same HIST");
    histData->Draw("e1 same");
    if(v_plotName_.at(0)=="top_arapidity"&&v_plotName_.at(1)=="top_pt")this->moveLegend(legend,0.13,-0.05,0.05,0.08);
    if(v_plotName_.at(0)=="top_pt"&&v_plotName_.at(1)=="ttbar_mass")this->moveLegend(legend,0.13,-0.05,0.05,0.08);
    if(v_plotName_.at(0)=="ttbar_arapidity"&&v_plotName_.at(1)=="ttbar_mass")this->moveLegend(legend,0.13,-0.02,0.05,0.06);
    if(v_plotName_.at(0)=="ttbar_arapidity"&&v_plotName_.at(1)=="ttbar_pt")this->moveLegend(legend,0.13,-0.05,0.05,0.08);
    if(v_plotName_.at(0)=="ttbar_delta_eta"&&v_plotName_.at(1)=="ttbar_mass")this->moveLegend(legend,0.13,-0.05,0.05,0.08);
    if(v_plotName_.at(0)=="ttbar_delta_phi"&&v_plotName_.at(1)=="ttbar_mass")this->moveLegend(legend,0.19,0.03,0.0,0.0);
    if(v_plotName_.at(0)=="ttbar_pt"&&v_plotName_.at(1)=="ttbar_mass")this->moveLegend(legend,0.13,-0.05,0.05,0.08);
    if(v_plotName_.at(0)=="ttbar_pt"&&v_plotName_.at(1)=="x1")this->moveLegend(legend,0.13,-0.05,0.05,0.08);
    if(v_plotName_.at(0)=="x1"&&v_plotName_.at(1)=="ttbar_mass")this->moveLegend(legend,0.13,-0.05,0.05,0.08);
    if(v_plotName_.at(0)=="top_arapidity"&&v_plotName_.at(1)=="ttbar_mass")this->moveLegend(legend,0.13,-0.05,0.05,0.08);
    legend->Draw("same");
    utils::drawRatio2d(stacksum, histData, 0.5,1.5,"cp",(TH1D*)hist->Clone());
    
    canvas->Print(plotsFolder_+ "CP_AllBins_" + v_plotName_.at(0) + "_vs_"+ v_plotName_.at(1) + ".pdf");
    canvas->Print(plotsFolder_+ "CP_AllBins_" + v_plotName_.at(0) + "_vs_"+ v_plotName_.at(1) + ".png");
    canvas->Print(plotsFolder_+ "CP_AllBins_" + v_plotName_.at(0) + "_vs_"+ v_plotName_.at(1) + ".root");
    
    //Delete objects
    hist->Delete();
    delete canvas;
    delete legend;
    stacksum->Delete();
}

void PlotterForTUnfold::moveLegend(TLegend*& legend,double moveX,double moveY,double plusWx,double plusWy){
    double x1=legend->GetX1();
    double x2=legend->GetX2();
    double y1=legend->GetY1();
    double y2=legend->GetY2();
    legend->SetX1(x1+moveX);
    legend->SetX2(x2+moveX+plusWx);
    legend->SetY1(y1+moveY);
    legend->SetY2(y2+moveY+plusWy);
}

void PlotterForTUnfold::writePlotEPSAllBins()
{
    // Prepare canvas and legend
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = utils::setLegend();
    legend->SetX1(0.1);
    legend->SetX2(0.4);
    legend->SetY1(0.7);
    legend->SetY2(0.85);
    
    canvas->cd();
    
    histEffAllBins_->Divide(histEffAllBins_,histSG_,1,1,"B");
    histEffAllBins_->SetStats(0);
    histEffAllBins_->SetMarkerStyle(20);
    histEffAllBins_->SetMarkerColor(3);
    legend->AddEntry(histEffAllBins_,"Efficiency","pe");
    
    histPurityAllBins_->Divide(histRecoGenAllBins_,histPurityAllBins_,1,1,"B");
    histPurityAllBins_->SetStats(0);
    histPurityAllBins_->SetMarkerStyle(23);
    histPurityAllBins_->SetMarkerColor(4);
    legend->AddEntry(histPurityAllBins_,"Purity","pe");
    
    histStabilityAllBins_->Divide(histRecoGenAllBins_,histStabilityAllBins_,1,1,"B");
    histStabilityAllBins_->SetStats(0);
    histStabilityAllBins_->SetMarkerStyle(22);
    histStabilityAllBins_->SetMarkerColor(2);
    legend->AddEntry(histStabilityAllBins_,"Stability","pe");
    
    for(int i=1;i<=histStabilityAllBins_->GetNbinsX();i++){
        vectorPurityAllBins_.push_back(histPurityAllBins_->GetBinContent(i));
        vectorStabilityAllBins_.push_back(histStabilityAllBins_->GetBinContent(i));
    }
    
    
    double yMax = histPurityAllBins_->GetMaximum() > histStabilityAllBins_->GetMaximum() ? histPurityAllBins_->GetMaximum()*1.15 : histStabilityAllBins_->GetMaximum()*1.15;
    histPurityAllBins_->SetAxisRange(0,yMax ,"Y");
    TH1D* hist = plotNiceTitel(0,1);
    
    histPurityAllBins_->Draw("e same");
    histStabilityAllBins_->Draw("e same");
    histEffAllBins_->Draw("e same");
    legend->Draw("same");
    canvas->Print(plotsFolder_+ "EPS_AllBins_" + v_plotName_.at(0) + "_vs_"+ v_plotName_.at(1) + ".pdf");
    canvas->Print(plotsFolder_+ "EPS_AllBins_" + v_plotName_.at(0) + "_vs_"+ v_plotName_.at(1) + ".root");
    canvas->Print(plotsFolder_+ "EPS_AllBins_" + v_plotName_.at(0) + "_vs_"+ v_plotName_.at(1) + ".png");
    
    //Delete objects
    hist->Delete();
    delete canvas;
    delete legend;
    
}



TH1D* PlotterForTUnfold::plotNiceTitel(double yAxisRangeMin,double yAxisRangeMax,TString yTitle){
    
    
    int Nbins1 = (int)(v_coarseBins_.at(1).size()-1)+(int)(v_oTrue_.at(1)+v_uTrue_.at(1));
    int Nbins0 = (int)(v_coarseBins_.at(0).size()-1)+(int)(v_oTrue_.at(0)+v_uTrue_.at(0));
    
    TH1D* hist = new TH1D("hist","",Nbins1,0.5,Nbins0*Nbins1+0.5);
    hist->SetTitle(";"+v_plotTitle_.at(1)+" "+v_plotUnits_.at(1)+";"+yTitle);
    hist->SetStats(0);
    hist->GetXaxis()->SetTickLength(0);
    hist->SetAxisRange(yAxisRangeMin,yAxisRangeMax ,"Y");
    
    
        for(int i=1;i<=Nbins1;++i)
        {
            TString binName = "";
            int isUnderflow = v_uTrue_.at(1);
            int isOverflow =  v_oTrue_.at(1);
            if(isUnderflow&&i==1) binName = "underflow";
            else if(isOverflow&&i==Nbins1){
                 binName = "overflow";
            }
            else{
                binName = utils::makeBinTitle("",v_coarseBins_.at(1).at(i-isUnderflow-1),v_coarseBins_.at(1).at(i-isUnderflow));
            }
            hist->GetXaxis()->SetBinLabel(i,binName.Data());
        }
        hist->SetTitle( utils::makeTitleBins(v_plotTitle_.at(0) +" "+ v_plotUnits_.at(0),v_coarseBins_.at(0),v_uTrue_.at(0),v_oTrue_.at(0) ) );
        hist->Draw();
        for(int i=1;i<Nbins1;++i)
        {
            double xLine = i*(Nbins0) + 0.5;
            TLine ln(xLine,yAxisRangeMin,xLine,yAxisRangeMax);
            ln.DrawClone("same");
        }
    
    return hist;
    
}


TH1D* PlotterForTUnfold::plotNiceTitelReco(double yAxisRangeMin,double yAxisRangeMax){
    
    
    int Nbins1 = (int)(v_fineBins_.at(1).size()-1)+(int)(v_oReco_.at(1)+v_uReco_.at(1));
    int Nbins0 = (int)(v_fineBins_.at(0).size()-1)+(int)(v_oReco_.at(0)+v_uReco_.at(0));
    
    TH1D* hist = new TH1D("hist","",Nbins1,0.5,Nbins0*Nbins1+0.5);
    hist->SetTitle(";"+v_plotTitle_.at(1)+" "+v_plotUnits_.at(1));
    hist->SetStats(0);
    hist->GetXaxis()->SetTickLength(0);
    hist->SetAxisRange(yAxisRangeMin,yAxisRangeMax ,"Y");
    
    
        for(int i=1;i<=Nbins1;++i)
        {
            TString binName = "";
            int isUnderflow = v_uReco_.at(1);
            int isOverflow =  v_oReco_.at(1);
            if(isUnderflow&&i==1) binName = "underflow";
            else if(isOverflow&&i==Nbins1){
                 binName = "overflow";
            }
            else{
                binName = utils::makeBinTitle("",v_fineBins_.at(1).at(i-isUnderflow-1),v_fineBins_.at(1).at(i-isUnderflow));
            }
            hist->GetXaxis()->SetBinLabel(i,binName.Data());
        }
        hist->SetTitle( utils::makeTitleBins(v_plotTitle_.at(0) +" "+ v_plotUnits_.at(0),v_fineBins_.at(0),v_uReco_.at(0),v_oReco_.at(0) ) );
        hist->DrawClone();
        for(int i=1;i<Nbins1;++i)
        {
            double xLine = i*(Nbins0) + 0.5;
            TLine ln(xLine,yAxisRangeMin,xLine,yAxisRangeMax);
            ln.DrawClone("same");
        }
    
    return hist;
    
}


void PlotterForTUnfold::writePlotXSec(const TH1* hData,const TH1* hMC,const TH1* hDataCount,const TH1* hMCCount,const int channel)
{
    // Prepare canvas and legend
    TCanvas* canvas = utils::setCanvas();
    TLegend* legend = utils::setLegend();
    canvas->SetBottomMargin(canvas->GetBottomMargin()+0.03);
    canvas->cd();
    
         TH1* histData = (TH1*)hData->Clone("data");
         //histData->Scale(1.0/lumi_);
         //histData->Scale(1.0/topxsec_);
         //histData->Scale(1.0/v_BR_.at(channel));
         histData->SetStats(0);
         
            TH1* histDataCount = (TH1*)hDataCount->Clone("dataCount");
            histDataCount->Scale(1.0/lumi_);
            histDataCount->Scale(1.0/v_BR_.at(channel));
         
         TH1* histMC   = (TH1*)hMC->Clone("mc");
         histMC->SetLineColor(2);
         
            TH1* histMCCount   = (TH1*)hMCCount->Clone("mcCount");
            histMCCount->Scale(1.0/lumi_);
            histMCCount->Scale(1.0/v_BR_.at(channel));
         
        TH2* h2dData = (TH2*)histData->Clone();
        TH2* h2dMC = (TH2*)histMC->Clone();
        
            TH2* h2dDataCount = (TH2*)histDataCount->Clone();
            TH2* h2dMCCount = (TH2*)histMCCount->Clone();

        for(int ix=-1;ix<=(int)v_coarseBins_.at(0).size()-1;ix++)
        {
            TString plotTitle = "";
            
            if(ix==-1)
            {
                if(v_uTrue_.at(0))plotTitle = v_plotTitle_.at(0) + " underflow bin";
                if(!v_uTrue_.at(0))continue;
            }
            else if(ix == (int)v_coarseBins_.at(0).size()-1)
            {
                if(v_oTrue_.at(0))plotTitle = v_plotTitle_.at(0) + " overflow bin";
                if(!v_oTrue_.at(0))continue;
            }
            else plotTitle = utils::makeBinTitle(v_plotTitle_.at(0),v_coarseBins_.at(0).at(ix),v_coarseBins_.at(0).at(ix+1));
            
            plotTitle = plotTitle + ";" + v_plotTitle_.at(1) + ", " + v_plotUnits_.at(1);
            
            TH1* h_temp = (TH1*)(h2dData->ProjectionY("tempData_ix",ix+1,ix+1,"e"));
            h_temp->Scale(1,"width");
            TH1* h_temp_mc =  (TH1*)(h2dMC->ProjectionY("tempMC_ix",ix+1,ix+1,"e"));
            h_temp_mc->Scale(1,"width");
            h_temp->SetStats(0);
            h_temp->SetTitle(plotTitle);
            
                TH1* h_tempCount = (TH1*)(h2dDataCount->ProjectionY("tempDataCount_ix",ix+1,ix+1,"e"));
                h_tempCount->Scale(1,"width");
                TH1* h_temp_mcCount =  (TH1*)(h2dMCCount->ProjectionY("tempMCCount_ix",ix+1,ix+1,"e"));
                h_temp_mcCount->Scale(1,"width");
            
            h_temp->SetAxisRange(0,( h_temp->GetMaximum() > h_temp_mc->GetMaximum() ? h_temp->GetMaximum()*1.15 : h_temp_mc->GetMaximum()*1.15 ),"Y");
            h_temp->GetYaxis()->SetLabelOffset(h_temp->GetYaxis()->GetLabelOffset()*1.2);
            h_temp->Draw("e");
            h_temp_mc->Draw("hist same");
            h_temp->Draw("e,same");
//          legend->AddEntry(h_temp,"Data","p");
//          legend->AddEntry(h_temp_mc,"MadGraph+Pythia","l");
//          legend->DrawClone("same");
            
            utils::drawRatio(h_temp, h_temp_mc, NULL,0.5,2);
            ///utils::drawRatioWithCorr(h_temp, h_temp_mc, rhoIJtotal_ ,0.49,1.51);
            
            if(ix==-1)canvas->Print(plotsFolder_+ "xSec_" + v_plotName_.at(1) + "_IN_"+ v_plotName_.at(0)+ "_" + "u" + ".pdf");
            else if(ix == (int)v_coarseBins_.at(0).size()-1)canvas->Print(plotsFolder_+ "xSec_" + v_plotName_.at(1) + "_IN_"+ v_plotName_.at(0)+ "_" + "o" + ".pdf");
            else{
                
                Output xsecInfo("xsec");
                    xsecInfo.add("lumi",std::vector<double>{lumi_});
                    xsecInfo.add("bin",std::vector<double>{v_coarseBins_.at(0).at(ix),v_coarseBins_.at(0).at(ix+1)});
                    xsecInfo.add(v_plotName_.at(0)+"-bins",v_coarseBins_.at(0));
                    xsecInfo.add(v_plotName_.at(1)+"-bins",v_coarseBins_.at(1));
                    std::vector<double> v_binWidth;
                    for(int i=1;i<=h_temp->GetNbinsX();i++)v_binWidth.push_back(h_temp->GetBinWidth(i));
                    xsecInfo.add("width",v_binWidth);
                    xsecInfo.add("inc.xsec",std::vector<double>{topxsec_});
                    xsecInfo.add("br",std::vector<double>{v_BR_.at(channel)});
                    xsecInfo.add("nDataInc_",std::vector<double>{nDataInc_});
                    xsecInfo.add("nRecoSigInc_",std::vector<double>{nRecoSigInc_});
                    xsecInfo.add("nRecoSingleTop_",std::vector<double>{nRecoSingleTop_});
                    xsecInfo.add("nRecoWPlusJets_",std::vector<double>{nRecoWPlusJets_});
                    xsecInfo.add("nRecoZee_",std::vector<double>{nRecoZee_});
                    xsecInfo.add("nRecoZtautau_",std::vector<double>{nRecoZtautau_});
                    xsecInfo.add("nRecoTtPlusZWG_",std::vector<double>{nRecoTtPlusZWG_});
                    xsecInfo.add("nRecoDiboson_",std::vector<double>{nRecoDiboson_});
                    xsecInfo.add("nRecoQCD_",std::vector<double>{nRecoQCD_});
                    xsecInfo.add("nRecoAllMC_",std::vector<double>{nRecoAllMC_});
                    xsecInfo.add("nRecoTtOtherInc_",std::vector<double>{nRecoTtOtherInc_});
                    xsecInfo.add("nGenSigInc_",std::vector<double>{nGenSigInc_});
                    xsecInfo.add("nBgrInc_",std::vector<double>{nBgrInc_});
                    xsecInfo.add("nRecoTtAllInc_",std::vector<double>{nRecoTtAllInc_});
                std::vector<double> v_xsec;
                std::vector<double> v_stat;
                std::vector<double> v_xsec_mc;
                    std::vector<double> v_xsecCount;
                    std::vector<double> v_statCount;
                    std::vector<double> v_xsec_mcCount;
                for(int i=1;i<=h_temp->GetNbinsX();i++){
                    v_xsec.push_back(h_temp->GetBinContent(i));
                    v_stat.push_back(h_temp->GetBinError(i));
                    v_xsec_mc.push_back(h_temp_mc->GetBinContent(i));
                        v_xsecCount.push_back(h_tempCount->GetBinContent(i));
                        v_statCount.push_back(h_tempCount->GetBinError(i));
                        v_xsec_mcCount.push_back(h_temp_mcCount->GetBinContent(i));
                }
                xsecInfo.add("stat",v_stat);
                xsecInfo.add("xsec",v_xsec);
                xsecInfo.add("mc",v_xsec_mc);
                xsecInfo.add("statCount",v_statCount);
                xsecInfo.add("xsecCount",v_xsecCount);
                xsecInfo.add("mcCount",v_xsec_mcCount);
                xsecInfo.add("vX",vX_);
                xsecInfo.add("vXmc",vXmc_);
                xsecInfo.add("covX",v_covX_);
                xsecInfo.add("vXnorm",vXnorm_);
                xsecInfo.add("vXnorm-mc",vXnormMC_);
                xsecInfo.add("covXnorm",v_covXnorm_);
                xsecInfo.add("vdY",v_dY_);
                xsecInfo.add("covY",v_covY_);
                xsecInfo.add("vdYnorm",v_dYnorm_);
                xsecInfo.add("covYnorm",v_covYnorm_);
                xsecInfo.add("vYmc",vYmc_);
                xsecInfo.add("vYnorm-mc",vYnormMC_);
                ///check before start
                    if(isVar_){xsecInfo.add("bgrWeight",std::vector<double>{var_});
                        xsecInfo.save(plotsFolder_+ "xSec_" + v_plotName_.at(1) + "_IN_"+ v_plotName_.at(0)+ "_" + std::to_string(ix) +"-"+ std::to_string((int)var_) + ".txt");}
                    else xsecInfo.save(plotsFolder_+ "xSec_" + v_plotName_.at(1) + "_IN_"+ v_plotName_.at(0)+ "_" + std::to_string(ix) + ".txt");
                    
                canvas->Print(plotsFolder_+ "xSec_" + v_plotName_.at(1) + "_IN_"+ v_plotName_.at(0)+ "_" + std::to_string(ix) + ".pdf");
                canvas->Print(plotsFolder_+ "xSec_" + v_plotName_.at(1) + "_IN_"+ v_plotName_.at(0)+ "_" + std::to_string(ix) + ".png");
            }
            canvas->Clear();
            h_temp->Delete();
            h_temp_mc->Delete();
            h_tempCount->Delete();
            h_temp_mcCount->Delete();
        }

      h2dData->Delete();
      h2dMC->Delete();
      h2dDataCount->Delete();
      h2dMCCount->Delete();

          histMC->Delete();
          histData->Delete();
          histMCCount->Delete();
          histDataCount->Delete();

    //Delete objects
    delete canvas;
    delete legend;
}



void PlotterForTUnfold::cleanMigrationMatrix(TH2*& histMigration)
{
    for(int i=1;i<=histMigration->GetYaxis()->GetNbins()+1;i++){
        histMigration->SetBinContent(0,i,0);
    }
}


double PlotterForTUnfold::rewTopPtUp(double pt)
{
    double w = 1;
    double b = 1.3;
    double a = -0.3/130.;
    w = a*pt + b;
    if(w<0.1)w=0.1;
    return w;
}


double PlotterForTUnfold::rewTopPtDown(double pt)
{
    double inW = rewTopPtUp(pt);
    double w = 2-inW;
    return w;
}


double PlotterForTUnfold::rewTopRapidityUp(double y)
{
    //func y = axx+b
    double w = 1;
    w = ((-0.2)/(1.1*1.1))*pow((y),2)+1.3;
    if(w<0.1)w=0.1;
    return w;
}


double PlotterForTUnfold::rewTopRapidityDown(double y)
{
    double inW = rewTopRapidityUp(y);
    double w = 2-inW;
    return w;
}

// OZ 15.01.2017 check equality of float numbers with tolerance
bool PlotterForTUnfold::IsEqual(const double val1, const double val2, const double eps/* = 1e-6*/) 
{   
  return (TMath::Abs(val1 - val2) < eps || TMath::Abs((val1 - val2) / val1) < eps); 
}

// OZ 15.01.2017 method to determine fine binning automatically
std::vector<double> PlotterForTUnfold::DetermineFineBinning(const TH1* hCoarse, const TH1* hInput, TString* message/* = NULL*/)
{
  // It takes hInput in very fine binning an find new binning to match hCoarse and to have at least minBinStat in each bin.
  // Bin edges of hInput must match hCoarse.
  const int minBinStat = 100;
  int debug = 1;
  
  // vector to store result binning
  std::vector<double> binning;
  
  printf("Input coarse binned histogram:\n");
  hCoarse->Print("all");
  printf("Input fine binned histogram:\n");
  hInput->Print("all");
  
  // check that coarse bin edges match fine bin edges
  int lastMatchedFineBin = 0;
  for(int coarseBin = 1; coarseBin <= hCoarse->GetNbinsX(); coarseBin++)
  {
    double edgeCoarseLow = hCoarse->GetBinLowEdge(coarseBin);
    double edgeCoarseUp = hCoarse->GetBinLowEdge(coarseBin + 1);
    if(debug)
      printf("coarse bin:  %f %f\n", edgeCoarseLow, edgeCoarseUp);
    bool matched = false;
    int events = 0;
    for(int fineBin = lastMatchedFineBin + 1; fineBin <= hInput->GetNbinsX() + 1; fineBin++)
    {
      double edgeFineLow = hInput->GetBinLowEdge(fineBin);
      double edgeFineUp = hInput->GetXaxis()->GetBinUpEdge(fineBin);
      if(debug)
        printf("fine bin:  %f %f\n", edgeFineLow, edgeFineUp);
      // skip fine bins below first coarse bin
      if(coarseBin == 1 && edgeFineLow < edgeCoarseLow)
        continue;
      // find first fine bin which matches first coarse bin
      if(coarseBin == 1 && IsEqual(edgeFineLow, edgeCoarseLow))
      {
        binning.push_back(edgeFineLow);
        if(debug)
          printf("pushing %f\n", edgeFineLow);
      }
      // now inside the range of the coarse binning: accumulate events
      events += hInput->GetBinContent(fineBin);
      if(debug)
        printf("events = %d\n", events);
      // next coarse bin edge is reached
      if(IsEqual(edgeFineUp, edgeCoarseUp))
      {
        // however if the accumulated event number is not enough needs to remove the previous bin edge (if it is not the previous coarse binning edge)
        if(coarseBin > 1 && events < minBinStat)
          if(!IsEqual(binning.back(), hCoarse->GetBinLowEdge(coarseBin)))
          {
            if(debug)
              printf("popping %f\n", binning.back());
            binning.pop_back();
          }
        matched = true;
        lastMatchedFineBin = fineBin;
        binning.push_back(edgeFineUp);
        if(debug)
          printf("pushing %f\n", edgeFineUp);
        break;
      }
      // minimum number of events is reached
      if(events >= minBinStat)
      {
        binning.push_back(edgeFineUp);
        if(debug)
          printf("pushing %f\n", edgeFineUp);
        events = 0;
      }
      // if already above the current coarse bin, then fine does not macth coarse binning
      if(edgeFineUp > edgeCoarseUp)
      {
        printf("Error in PlotterForTUnfold::DetermineFineBinning() bin %f %f from coarse binning does not match any bin edge in fine binning\n", edgeCoarseLow, edgeCoarseUp);
        exit(1);
      }
    }
    if(!matched)
    {
      printf("Error in PlotterForTUnfold::DetermineFineBinning() bin %f %f from coarse binning does not match any bin edge in fine binning\n", edgeCoarseLow, edgeCoarseUp);
      exit(1);
    }
  }

  if(debug)
  {
    printf("INPUT:  ");
    for(int b = 0; b <= hInput->GetNbinsX(); b++)
      printf(" %.3f", hInput->GetBinLowEdge(b + 1));
    printf("\n");
    printf("coarse");
    for(int b = 0; b <= hCoarse->GetNbinsX(); b++)
      printf(" %.4f", hCoarse->GetBinLowEdge(b + 1));
    printf("\n");
    printf("fine");
    for(unsigned int b = 0; b < binning.size(); b++)
      printf(" %.3f", binning[b]);
    printf("\n");
  }
  // output message
  if(message)
  {
    *message = "\n";
    *message += "title dummy\n";
    *message += "units dummy\n";
    *message += "yunits dummy\n";
    *message += "\n";
    *message += "textitle undef\n";
    *message += "texunits undef\n";
    *message += "\n";
    *message += "texxsecunitsnorm undef\n";
    *message += "texxsecunits undef\n";
    *message += "\n";
    *message += "yunitscp Top~quarks\n";
    *message += "ylog 0\n";
    *message += "\n";
    *message += "\n";
    *message += "2d\n";
    *message += "coarse";
    for(int b = 0; b <= hCoarse->GetNbinsX(); b++)
      *message += TString::Format(" %.4f", hCoarse->GetBinLowEdge(b + 1));
    *message += "\n";
    *message += "fine";
    for(unsigned int b = 0; b < binning.size(); b++)
      *message += TString::Format(" %.4f", binning[b]);
    *message += "\n";
    *message += "uoTrue 0 0\n";
    *message += "uoReco 0 0\n";
    *message += "2d\n";
  }
  
  return binning;
}

// OZ 15.01.2017 simplified unregularized unfolding, useful as cross check
// returns histogram with normalized differential cross section
TH1* PlotterForTUnfold::DoSimpleUnfolding(const TH1* hData, const TH2* hResp)
{
  int flagUFRec = 0;
  int flagOFRec = 0;
  int flagUFGen = 0;
  int flagOFGen = 0;
  
  TUnfoldBinning* detBinning = new TUnfoldBinning("det");
  detBinning->AddBinning("detdistribution")->AddAxis(*hResp->GetYaxis(), flagUFRec, flagOFRec);
  TUnfoldBinning* genBinning = new TUnfoldBinning("gen");
  genBinning->AddBinning("gendistribution")->AddAxis(*hResp->GetXaxis(), flagUFGen, flagOFGen);
  
  TUnfold::EConstraint constraintMode= TUnfold::kEConstraintArea;
  //TUnfold::EConstraint constraintMode= TUnfold::kEConstraintNone;
  
  TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
  //TUnfold::ERegMode regMode = TUnfold::kRegModeDerivative;
  //TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
  //TUnfold::ERegMode regMode = TUnfold::kRegModeNone;
  
  //TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeBinWidth;
  TUnfoldDensity::EDensityMode densityFlags = TUnfoldDensity::kDensityModeNone;
  
  // detailed steering for regularisation
  const char *REGDISTR=0;
  const char *REGAXISSTEER=0;
  //detectorBinning_->PrintStream(std::cout);
  //generatorBinning_->PrintStream(std::cout);

  TUnfoldDensity unfold(hResp, TUnfold::kHistMapOutputHoriz, regMode, constraintMode, densityFlags, genBinning, detBinning, REGDISTR, REGAXISSTEER);
  unfold.SetInput(hData);
  unfold.DoUnfold(0);  
  TH1* hOut  = unfold.GetOutput("unfoldedData_", "unfolding result");
  hOut->Scale(1.0 / hOut->Integral());
  return hOut;
}

TH2D* PlotterForTUnfold::H1dto2d(const TH1D* hIn)
{
  std::vector<double> bins;
  for(int b = 0; b <= hIn->GetNbinsX(); b++)
    bins.push_back(hIn->GetBinLowEdge(b + 1));
  std::vector<double> binsX;
  binsX.push_back(0.0);
  binsX.push_back(1.0);
  TH2D* hOut = new TH2D(hIn->GetName(), hIn->GetTitle(), 1, &binsX[0], hIn->GetNbinsX(), &bins[0]);
  for(int b = 0; b <= hIn->GetNbinsX() + 1; b++)
  {
    hOut->SetBinContent(1, b, hIn->GetBinContent(b));
    hOut->SetBinError(1, b, hIn->GetBinError(b));
  }
  return hOut;
}

// OZ 17.01.2017 2D (with 1 bin in 1st dim) to 1D histogram converter (overlow underflow in 1st dimension is lost)
TH1D* PlotterForTUnfold::H2dto1d(const TH2D* hIn)
{
  std::vector<double> bins;
  for(int b = 0; b <= hIn->GetNbinsY(); b++)
    bins.push_back(hIn->GetYaxis()->GetBinLowEdge(b + 1));
  TH1D* hOut = new TH1D(hIn->GetName(), hIn->GetTitle(), bins.size() - 1, &bins[0]);
  for(int b = 0; b <= hIn->GetNbinsY() + 1; b++)
  {
    hOut->SetBinContent(b, hIn->GetBinContent(1, b));
    hOut->SetBinError(b, hIn->GetBinError(1, b));
  }
  return hOut;
}

// OZ 17.01.2017 set generator level histogram
void PlotterForTUnfold::SetHistGen(TH2D* histGen)
{
  histGen_ = new TH2D(*histGen);
  /*for(int bx = 0; bx <= histGen->GetNbinsX() + 1; bx++)
  {
    for(int by = 0; by <= histGen->GetNbinsY() + 1; by++)
    {
      double bw = histGen->GetYaxis()->GetBinLowEdge(by + 1) - histGen->GetYaxis()->GetBinLowEdge(by);
      histGen_->SetBinContent(bx, by, histGen->GetBinContent(bx, by) / bw);
      histGen_->SetBinError(bx, by, histGen->GetBinError(bx, by) / bw);
    }
  }*/
  histGenNorm_ = new TH2D(*histGen);
  histGenNorm_->Scale(1.0 / histGenNorm_->Integral());
}
