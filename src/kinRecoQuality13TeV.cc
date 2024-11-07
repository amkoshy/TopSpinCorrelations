//#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>

#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1D.h>
#include <TFile.h>
#include <TSystem.h>

#include "utils.h"
#include "Sample.h"
#include "Samples.h"
#include "GlobalScaleFactors.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"
#include "../../common/include/utils.h"
//#include "EffHist.h"
#include "UsefulTools.h"

double sqr(double value) { return value*value; }




// int main(int argc, char *argv[])
int main()
{
    
    gSystem->Exec("mkdir Plots");
    
    RootFileReader *fileReader_ = RootFileReader::getInstance();

    UsefulTools* usefulTools = new UsefulTools(fileReader_,0,1);
    
    // Set up Lumi and xsec
    double Luminosity = usefulTools->lumi;
    
    // Set up channels
//     std::vector<Channel::Channel> v_channel(Channel::realChannels);
    std::vector<Channel::Channel> v_channel(1,Channel::emu);//emu
    std::cout << "Processing channels: ";
    for(auto channel : v_channel) std::cout << Channel::convert(channel) << " ";
    std::cout << "\n\n";
    
    // Set up systematics
    std::vector<Systematic::Systematic> v_systematic ={Systematic::nominalSystematic()}; 
    for(auto systematic : v_systematic) std::cout << systematic.name() << "\n ";
    std::cout << "\n\n";

    const bool dyCorrection = 1;
    const bool ttbbCorrection = 0;
    const GlobalScaleFactors* globalScaleFactors = new GlobalScaleFactors(v_channel, v_systematic, Luminosity, dyCorrection, ttbbCorrection);
    // Access all samples   
    const Samples samples("FileLists", v_channel, v_systematic, globalScaleFactors);
    
    // Access correction factors
    const TString fullStepname7 = ttbar::extractSelectionStepAndJetCategory("step7");
    const TString fullStepname8 = ttbar::extractSelectionStepAndJetCategory("step8");
    const SystematicChannelFactors globalWeights7 = samples.globalWeights(fullStepname7).first;
    const SystematicChannelFactors globalWeights8 = samples.globalWeights(fullStepname8).first;
    
//     std::vector<TString> histNames {"_Eff","_vs_JetMult","_vs_LepEta","_vs_JetEta",
//                                     "_vs_LeppT","_vs_MET","_vs_LepEta2","_vs_JetEta2",
//                                     "_vs_LeppT2","_vs_AntiLeptonpT","_vs_LeptonpT","_vs_LeptonEta",
//                                     "_vs_AntiLeptonEta","_vs_JetpT","_vs_JetpT2"};
//     TString histNamePrefix("KinReco_nRecoEvt");
    std::vector<TString> histNames {"HypToppT","HypTopRapidity","HypTTBarMass"};

    
    TString savePlotPath("");
    
    // Loop over all channels and systematics and produce plots
    const SystematicChannelSamples& m_systematicChannelSample(samples.getSystematicChannelSamples());
    for(const auto& systematicChannelSamples : m_systematicChannelSample){                     //systematic loop
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        for(const auto& channelSample : systematicChannelSamples.second){                     //channel loop
            const Channel::Channel& channel(channelSample.first);
            const std::vector<Sample>& v_sample(channelSample.second);
            const std::vector<double>& v_weight7(globalWeights7.at(systematic).at(channel));
            const std::vector<double>& v_weight8(globalWeights8.at(systematic).at(channel));
            
            
            for(TString histName : histNames){                                              //hist loop

                        savePlotPath = common::assignFolder("Plots",channel,systematic);
                        savePlotPath.Append("/kinRecoQuality_"+histName+".eps");

            TH1D* histSignal = 0;
            TH1D* histBg = 0;
                for(size_t iSample = 0; iSample < v_sample.size(); ++iSample){   //sample loop
                const Sample& sample(v_sample.at(iSample));
                    printf("%s  %s  %s  %s\n", systematic.name().Data(), Channel::convert(channel).Data(), histName.Data(), sample.inputFile().Data());
                if(sample.legendEntry() == "QCD multijet") continue;
                if(sample.legendEntry() == "W+jets") continue;
                if(sample.sampleType() == Sample::data) continue;

                TH1D* histTemp;
                TFile* dataFile=new TFile(sample.inputFile());
                    histTemp=(TH1D* ) dataFile->Get(histName);
                    histTemp->SetDirectory(0);
                dataFile->Close();
                    histTemp->Scale(v_weight8.at(iSample));
                    if(histName == "HypToppT")histTemp->Rebin(50);
                    if(histName == "HypTopRapidity")histTemp->Rebin(5);
                    if(histName == "HypTTBarMass")histTemp->Rebin(100);
                
                if(sample.inputFile().Contains("ttbarsignal")){ 
                    if(!histSignal)histSignal = (TH1D*)histTemp->Clone("histSignal");
                    else histSignal->Add(histTemp,1);
                }
                else{
                   if(!histBg) histBg  = (TH1D*)histTemp->Clone("histBg");
                   else histBg->Add(histTemp,1);
                }
                
                histTemp->Delete();
                
                }//samples loop
                
                TH1D* histS = (TH1D*)histSignal->Clone("histSignificance");
                for(int i = 1 ; i<=histS->GetNbinsX();i++)
                {
                    double S = histSignal->GetBinContent(i);
                    double dS = histSignal->GetBinError(i);
                    double B = histBg->GetBinContent(i);
                    double dB = histBg->GetBinError(i);
                    double Signif = S/sqrt(S+B);
                    double dSignif = sqrt(sqr(dS*(S+2*B))+sqr(dB*S))/(2*(S+B)*sqrt(S+B));
                    histS->SetBinContent(Signif,i);
                    histS->SetBinError(dSignif,i);
                    
                }




                TCanvas canvas("canvas","canvas");
                    histS->DrawCopy();
                canvas.Print(savePlotPath);
                
                if(histSignal)histSignal->Delete();
                if(histBg)histBg->Delete();
                if(histS)histS->Delete();
                
            }//hist loop
        }//channel loop
    }//systematics loop
    
//     TString rootFilePrePath(common::CMSSW_BASE()+"/src/TopAnalysis/Configuration/analysis/diLeptonic/");
//     

    
    return 0;
}
