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
#include <TRandom.h>
#include <TMath.h>
#include "TSystem.h"

#include <MiniTreeReader.h>
#include <MiniTree.h>
#include <MiniTreeAnalyzer.h>




using namespace std;


MiniTreeAnalyzer::MiniTreeAnalyzer(const TString &year_string, const TString &filename_string, const TString &systematic_string, const TString &channel_string)
{
    TH1::SetDefaultSumw2(kTRUE);

    systematic = Systematic::undefinedSystematic();
    systematic = Systematic::Systematic(systematic_string);

    Systematic::Systematic nominalSyst = Systematic::Systematic("Nominal");

    channel= Channel::undefined;
    channel = Channel::convert(channel_string);

    TString inpath;

    if( systematic.type()==Systematic::lept || systematic.type()==Systematic::trig ||
        systematic.type()==Systematic::ele || systematic.type()==Systematic::muon ||
        systematic.type()==Systematic::eleID || systematic.type()==Systematic::muonID ||
        systematic.type()==Systematic::eleReco || systematic.type()==Systematic::muonIso ||
        systematic.type()==Systematic::trigEta || systematic.type()==Systematic::pu ||
        systematic.type()==Systematic::kin || systematic.type()==Systematic::btag ||
        systematic.type()==Systematic::btagLjet || systematic.type()==Systematic::meScale ||
        systematic.type()==Systematic::meFacScale || systematic.type()==Systematic::meRenScale ||
        systematic.type()==Systematic::bFrag || systematic.type()==Systematic::bFrag_central ||
        systematic.type()==Systematic::bFrag_Peterson ||
        systematic.type()==Systematic::psScaleWeight ||
        systematic.type()==Systematic::pdf || systematic.type()==Systematic::alphasPdf ||
        systematic.type()==Systematic::bSemilep ||systematic.type()==Systematic::l1prefiring){
            inpath = common::accessFolder("mergedMiniTree_"+year_string, channel, nominalSyst);
        }else{
            inpath = common::accessFolder("mergedMiniTree_"+year_string, channel, systematic);
        }

    TString outpath = assignFolder("miniTreeHistograms_"+year_string, channel, RenameSystematicString(systematic.name()));
    TString inputFilePath_ = inpath+Channel::convert(channel)+"_"+filename_string;
    TString outputFilePath_ = outpath+Channel::convert(channel)+"_"+filename_string;

    TString DNNPath = "/nfs/dust/cms/user/sewuchte/TopRhoNetwork/rhoOutput/rhoRegModel_optim_2017.pb";
    InitDNNModel(DNNPath);

    SetOutpath(outputFilePath_);

    start = std::chrono::steady_clock::now();
    std::cout << "MiniTreeAnalyzer(): Input file for analysis: " + inputFilePath_<< std::endl;
    inputFilePath = inputFilePath_;
    inputFileName = filename_string;
    if (inputFileName.Contains("run")){
        isData = true;
        isMC = false;
    }else{
        isData = false;
        isMC = true;
    }
    if (inputFileName.Contains("ttbarsignalplustau")){
        if (!inputFileName.Contains("fromHadronic") && !inputFileName.Contains("fromLjets")){
            isSignal = true;
        }
    }

    bool isttbarsample = inputFileName.Contains("ttbarbg_fromDilepton")||inputFileName.Contains("ttbarbg_fromHadronic")||inputFileName.Contains("ttbarbg_fromLjets")
                        ||inputFileName.Contains("ttbarbgviatau")||inputFileName.Contains("ttbarsignalplustau");

    processType = Process::convertFromFilename(inputFileName);

    fileReader = RootFileReader::getInstance();
    configHelper = new PlotterConfigurationHelper(fileReader, year_string, "", false);

    lumiWeight = configHelper->CalcLumiWeight(inputFilePath_);

    TString nameTLWhist;
    TString nameTLNRWhist;
    hasPDFWeights=false;
    hasPSWeights=false;
    hasBFragWeights=false;
    hasMEWeights=false;
    hasBSemilepWeights=false;

    if((isMC)&&(
        systematic.type()==Systematic::pu ||
        systematic.type()==Systematic::meScale ||
        systematic.type()==Systematic::meFacScale || systematic.type()==Systematic::meRenScale ||
        systematic.type()==Systematic::bFrag ||
        systematic.type()==Systematic::psScaleWeight ||
        systematic.type()==Systematic::pdf || systematic.type()==Systematic::alphasPdf ||
        systematic.type()==Systematic::bSemilep))
    {
        nameTLWhist = "normalization_histograms/trueLevelWeightSum_hist_"+Systematic::convertType(systematic.type())+Systematic::convertVariation(systematic.variation());
        if(systematic.type()==Systematic::pdf){
            TString addition;
            addition.Form("_%i", systematic.variationNumber());
            nameTLWhist=nameTLWhist+addition;
        }
        if(systematic.type()==Systematic::psScaleWeight){
            TString addition;
            addition.Form("_%i", systematic.variationNumber());
            nameTLWhist=nameTLWhist+addition;
        }

        if(systematic.type()!=Systematic::pu){
            nameTLNRWhist = "trueLevelNoRenormalisationWeightSum_hist";
        }else{
            nameTLNRWhist = "normalization_histograms/trueLevelNoRenormalisationWeightSum_hist"+Systematic::convertType(systematic.type())+Systematic::convertVariation(systematic.variation());
        }

    }else{
        nameTLWhist = "trueLevelWeightSum_hist";
        nameTLNRWhist = "trueLevelNoRenormalisationWeightSum_hist";
    }
    if(systematic.type()==Systematic::psScaleWeight){
        const TH1D* h_trueLevelWeightSum_hist_temp = fileReader->Get<TH1D>(inputFilePath_, nameTLWhist, true);
        if(h_trueLevelWeightSum_hist_temp==NULL){
            hasPSWeights=false;
            nameTLWhist = "trueLevelWeightSum_hist";
            nameTLNRWhist = "trueLevelNoRenormalisationWeightSum_hist";
        }else{
            hasPSWeights=true;
            // if(inputFileName.Contains("ttbarZtollnunu")){
            if(!isttbarsample){
                hasPSWeights=false;
                nameTLWhist = "trueLevelWeightSum_hist";
                nameTLNRWhist = "trueLevelNoRenormalisationWeightSum_hist";
            }
        }
    }
    if(systematic.type()==Systematic::pdf){
        const TH1D* h_trueLevelWeightSum_hist_temp = fileReader->Get<TH1D>(inputFilePath_, nameTLWhist, true);
        if(h_trueLevelWeightSum_hist_temp==NULL){
            hasPDFWeights=false;
            nameTLWhist = "trueLevelWeightSum_hist";
            nameTLNRWhist = "trueLevelNoRenormalisationWeightSum_hist";
        }else{
            hasPDFWeights=true;
            if(!isttbarsample){
                hasPDFWeights=false;
                nameTLWhist = "trueLevelWeightSum_hist";
                nameTLNRWhist = "trueLevelNoRenormalisationWeightSum_hist";
            }
        }
    }
    if(systematic.type()==Systematic::bFrag){
        const TH1D* h_trueLevelWeightSum_hist_temp = fileReader->Get<TH1D>(inputFilePath_, nameTLWhist, true);
        if(h_trueLevelWeightSum_hist_temp==NULL){
            hasBFragWeights=false;
            nameTLWhist = "trueLevelWeightSum_hist";
            nameTLNRWhist = "trueLevelNoRenormalisationWeightSum_hist";
        }else{
            hasBFragWeights=true;
            if(!isttbarsample){
                hasBFragWeights=false;
                nameTLWhist = "trueLevelWeightSum_hist";
                nameTLNRWhist = "trueLevelNoRenormalisationWeightSum_hist";
            }
        }
    }
    if(systematic.type()==Systematic::bSemilep){
        const TH1D* h_trueLevelWeightSum_hist_temp = fileReader->Get<TH1D>(inputFilePath_, nameTLWhist, true);
        if(h_trueLevelWeightSum_hist_temp==NULL){
            hasBSemilepWeights=false;
            nameTLWhist = "trueLevelWeightSum_hist";
            nameTLNRWhist = "trueLevelNoRenormalisationWeightSum_hist";
        }else{
            hasBSemilepWeights=true;
            if(!isttbarsample){
                hasBSemilepWeights=false;
                nameTLWhist = "trueLevelWeightSum_hist";
                nameTLNRWhist = "trueLevelNoRenormalisationWeightSum_hist";
            }
        }
    }
    if(systematic.type()==Systematic::meScale ||
       systematic.type()==Systematic::meFacScale || systematic.type()==Systematic::meRenScale || systematic.type()==Systematic::alphasPdf){
        const TH1D* h_trueLevelWeightSum_hist_temp = fileReader->Get<TH1D>(inputFilePath_, nameTLWhist, true);
        if(h_trueLevelWeightSum_hist_temp==NULL){
            hasMEWeights=false;
            nameTLWhist = "trueLevelWeightSum_hist";
            nameTLNRWhist = "trueLevelNoRenormalisationWeightSum_hist";
        }else{
            hasMEWeights=true;
            if(!isttbarsample){
                hasMEWeights=false;
                nameTLWhist = "trueLevelWeightSum_hist";
                nameTLNRWhist = "trueLevelNoRenormalisationWeightSum_hist";
            }
        }
    }

    const TH1D* h_trueLevelWeightSum_hist = fileReader->Get<TH1D>(inputFilePath_, nameTLWhist);
    const TH1D* h_trueLevelNoRenormalisationWeightSum_hist = fileReader->Get<TH1D>(inputFilePath_, nameTLNRWhist);



    double trueLevelWeightSum_ = h_trueLevelWeightSum_hist->GetBinContent(1);
    double trueLevelNoRenormalisationWeightSum_ = h_trueLevelNoRenormalisationWeightSum_hist->GetBinContent(1);
    renormWeight = trueLevelWeightSum_ != 0. ? trueLevelNoRenormalisationWeightSum_/trueLevelWeightSum_ : -999.;
    renormWeight = trueLevelNoRenormalisationWeightSum_/trueLevelWeightSum_;
    std::cout<<"Read in h_trueLevelWeightSum_hist "<<nameTLWhist<<std::endl;
    std::cout<<"Read in h_trueLevelNoRenormalisationWeightSum_hist "<<nameTLNRWhist<<std::endl;

    if (renormWeight<0.){
        renormWeight = 1.;
        if(systematic.type()==Systematic::meScale ||
           systematic.type()==Systematic::meFacScale || systematic.type()==Systematic::meRenScale || systematic.type()==Systematic::alphasPdf){
               hasMEWeights=false;
           }
        if(systematic.type()==Systematic::psScaleWeight){
               hasPSWeights=false;
           }
        if(systematic.type()==Systematic::bFrag){
               hasBFragWeights=false;
           }
        if(systematic.type()==Systematic::bSemilep){
               hasBSemilepWeights=false;
           }
        if(systematic.type()==Systematic::pdf){
               hasPDFWeights=false;
           }
    }

    std::cout<<"luminosity weight: "<<lumiWeight<<" | renormalisation weight "<<renormWeight<<std::endl;

    event = new Event;
    event->Init();

    // h_eventsPassedStep3 = new TH1D("h_eventsPassedStep3","h_eventsPassedStep3",1,0.,1.);
    // h_eventsPassedStep4 = new TH1D("h_eventsPassedStep4","h_eventsPassedStep4",1,0.,1.);
    // h_eventsPassedStep5 = new TH1D("h_eventsPassedStep5","h_eventsPassedStep5",1,0.,1.);
    // h_eventsPassedStep6 = new TH1D("h_eventsPassedStep6","h_eventsPassedStep6",1,0.,1.);
    // h_eventsPassedStep7 = new TH1D("h_eventsPassedStep7","h_eventsPassedStep7",1,0.,1.);
    // h_eventsPassedStep8 = new TH1D("h_eventsPassedStep8","h_eventsPassedStep8",1,0.,1.);
    // h_eventsPassedStep8L = new TH1D("h_eventsPassedStep8L","h_eventsPassedStep8L",1,0.,1.);
}


void MiniTreeAnalyzer::ProgressBar(const int &progress){
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

void MiniTreeAnalyzer::SetOutpath(const TString &path_){
    outFilePath = path_;
}

void MiniTreeAnalyzer::InitDNNModel(TString& pathToModel){
    dnnModel = new ztop::TFModel( (std::string) pathToModel, 21, "input_1",1,  "dense_2/Sigmoid");
}

void MiniTreeAnalyzer::EvaluateDNN(Event& event)
{
    uint tempMaxsize = maxSize;
    dnn_inputs[0]= 0.; dnn_inputs[1]= 0.; dnn_inputs[2]= 0.; dnn_inputs[3]= 0.;
    dnn_inputs[4]= 0.; dnn_inputs[5]= 0.; dnn_inputs[6]= 0.; dnn_inputs[7]= 0.;
    dnn_inputs[8]= 0.; dnn_inputs[9]= 0.; dnn_inputs[10]= 0.; dnn_inputs[11]= 0.;
    dnn_inputs[12]= event.lepton1.Px(); dnn_inputs[13]= event.lepton1.Py(); dnn_inputs[14]= event.lepton1.Pz(); dnn_inputs[15]= event.lepton1.E();
    dnn_inputs[16]= event.lepton2.Px(); dnn_inputs[17]= event.lepton2.Py(); dnn_inputs[18]= event.lepton2.Pz(); dnn_inputs[19]= event.lepton2.E();
    dnn_inputs[20]= event.met.E();

    if (event.jets.size()<=tempMaxsize){
        tempMaxsize = event.jets.size();
    }
    for(uint i=0; i < tempMaxsize; i++){
        dnn_inputs[i*4+0] = event.jets.at(i).Px();
        dnn_inputs[i*4+1] = event.jets.at(i).Py();
        dnn_inputs[i*4+2] = event.jets.at(i).Pz();
        dnn_inputs[i*4+3] = event.jets.at(i).E();
    }
    dnn_outputs = dnnModel->evaluate(dnn_inputs);
    // event.DNN_rho_new = dnn_outputs.at(0);
    event.DNN_rho = dnn_outputs.at(0);
}


bool MiniTreeAnalyzer::GetZWindowCutDecisionStep4(const Event &ev){
    if(!(channel==Channel::emu)){
        bool isOnZ = (ev.dilepton.M()>76.) && (ev.dilepton.M()<106.);
        if(!isOnZ){
            return (ev.passStep3);
        }else{
            return false;
        }
    }else{
        return (ev.passStep3);
    }
}

bool MiniTreeAnalyzer::GetMETCutDecisionStep6(const Event &ev){
    if(!(channel==Channel::emu)){
        if(ev.met.Pt()>metCut){
            return (ev.passStep3);
        }else{
            return false;
        }
    }else{
        return (ev.passStep3);
    }
}



bool MiniTreeAnalyzer::GetCategoryFillDecision(const Event &ev, const BTagCategory::BTagCategory &btagcategory,
                                               const RecoRhoBinCategory::RecoRhoBinCategory &rhobin, const NJetCategory::NJetCategory &njetcategory){
    bool nJetPassed=false;
    bool nBJetPassed=false;
    bool rhoBinPassed=false;
    switch(njetcategory){
        case(NJetCategory::ZeroAndOneJet):
            if(ev.jets.size()==0 || ev.jets.size()==1){
                nJetPassed=true;
                break;
            }break;
        case(NJetCategory::TwoJet):
            if(ev.jets.size()==2){
                nJetPassed=true;
                break;
            }break;
        case(NJetCategory::GreaterTwoJet):
            if(ev.jets.size()>2){
                nJetPassed=true;
                break;
            }break;
        case(NJetCategory::InclusiveNJet):
                nJetPassed=true;
                break;
        default:
            std::cerr<<"Error in MiniTreeAnalyzer::GetFillDecision()! undefined NJetCategory\n...break\n"<<std::endl;
            exit(98);
            break;
    }

    switch(btagcategory){
        case(BTagCategory::ZeroPlusGreaterTwoBTag):
            if(ev.bjets.size()==0 || ev.bjets.size()>2){
                nBJetPassed=true;
                break;
            }break;
        case(BTagCategory::OneBTag):
            if(ev.bjets.size()==1){
                nBJetPassed=true;
                break;
            }break;
        case(BTagCategory::TwoBTag):
            if(ev.bjets.size()==2){
                nBJetPassed=true;
                break;
            }break;
        case(BTagCategory::InclusiveBTag):
                nBJetPassed=true;
                break;
        default:
            std::cerr<<"Error in MiniTreeAnalyzer::GetFillDecision()! undefined BTagCategory\n...break\n"<<std::endl;
            exit(98);
            break;
    }

    switch(rhobin){

        case(RecoRhoBinCategory::InclusiveRho):
            // rhoBinPassed=ev.passStep3;
            rhoBinPassed=ev.passStep3 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev);
            break;
        case(RecoRhoBinCategory::RhoStudies):
            // if(ev.jets.size()>2 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev) && ev.hasKRS){
            // if(ev.jets.size()>2 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev)){
            // if(ev.jets.size()>2 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev) && ev.DNN_rho>0.){
            if(ev.jets.size()>2 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev) && ev.final_rho>0.){
            // if(ev.jets.size()>2 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev) && ev.DNN_kinRecoSol.numberOfSolutions()){
            // if(ev.jets.size()>2 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev) && ev.hasLKRS){
                // if(ev.kinReco_nonbjets.at(0).Pt()>50.){
                // if(ev.kinReco_nonbjets.at(2).Pt()>50.){
                // if(ev.kinReco_nonbjets.at(0).Pt()>50. && fabs(ev.kinReco_nonbjets.at(0).Eta())<2.4){
                // if(ev.jets.at(0).Pt()>50.){
                if(ev.jets.at(2).Pt()>50.){
                // if(ev.looseKinReco_nonbjets.at(0).Pt()>50.){
                // if(ev.DNN_kinReco_nonbjets.at(0).Pt()>50.){
                    // float rhoValue = scaleValue/(ev.kinReco_ttbar+ev.kinReco_nonbjets.at(0)).M();
                    // rhoBinPassed = rhoValue > RecoRhoBinCategory::getLowerUpperEdge(rhobin).first && rhoValue < RecoRhoBinCategory::getLowerUpperEdge(rhobin).second;
                    rhoBinPassed = true;
                    break;
                }
            }
            break;
        case(RecoRhoBinCategory::NoRho):
            if(ev.jets.size()<3){
                if(ev.jets.size()<2){
                    // rhoBinPassed = ev.passStep3 && GetMETCutDecisionStep6(ev);
                    rhoBinPassed = ev.passStep3 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev);
                    break;
                }
                if(ev.jets.size()==2){
                    // rhoBinPassed = GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev) && ev.hasKRS;
                    rhoBinPassed = ev.passStep3 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev);
                    break;
                    // rhoBinPassed = ev.passStep3 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev) && ev.hasKRS;
                    // rhoBinPassed = GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev) && ev.hasLKRS;
                }
            }else{
              // if(ev.jets.size()>2 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev) && ev.hasKRS){
              // if(ev.jets.size()>2 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev) && ev.hasLKRS){
              if(ev.jets.size()>2 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev)){
                // if(ev.kinReco_nonbjets.at(0).Pt()<=50. && fabs(ev.kinReco_nonbjets.at(0).Eta())<2.4){
                // if(ev.looseKinReco_nonbjets.at(0).Pt()<=50. && fabs(ev.looseKinReco_nonbjets.at(0).Eta())<2.4){
                if(ev.jets.at(2).Pt()<=50. && fabs(ev.jets.at(2).Eta())<2.4){
                    rhoBinPassed = ev.passStep3;
                    break;
                }
              }
            }
            break;
        case(RecoRhoBinCategory::Rho0to02):
        case(RecoRhoBinCategory::Rho02to03):
        case(RecoRhoBinCategory::Rho03to04):
        case(RecoRhoBinCategory::Rho04to05):
        case(RecoRhoBinCategory::Rho05to06):
        case(RecoRhoBinCategory::Rho06to1):
        case(RecoRhoBinCategory::Rho0to1):
            // cout<<ev.jets.size()<<endl;
            // if(ev.jets.size()>2 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev) && ev.hasKRS){
            // if(ev.jets.size()>2 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev)){
            // if(ev.jets.size()>2 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev) && ev.DNN_rho>0.){
            if(ev.jets.size()>2 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev) && ev.final_rho>0.){
            // if(ev.jets.size()>2 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev) && ev.DNN_kinRecoSol.numberOfSolutions()){
            // if(ev.jets.size()>2 && GetZWindowCutDecisionStep4(ev) && GetMETCutDecisionStep6(ev) && ev.hasLKRS){
                // if(ev.kinReco_nonbjets.at(0).Pt()>50. && fabs(ev.kinReco_nonbjets.at(0).Eta())<2.4){
                // if(ev.jets.at(0).Pt()>50.){
                if(ev.jets.at(2).Pt()>50.){
                // if(ev.DNN_kinReco_nonbjets.at(0).Pt()>50.){
                // if(ev.looseKinReco_nonbjets.at(0).Pt()>50.){
                    // float rhoValue = scaleValue/(ev.kinReco_ttbar+ev.kinReco_nonbjets.at(0)).M();
                    // float rhoValue = scaleValue/(ev.looseKinReco_ttbar+ev.looseKinReco_nonbjets.at(0)).M();
                    // rhoBinPassed = rhoValue > RecoRhoBinCategory::getLowerUpperEdge(rhobin).first && rhoValue < RecoRhoBinCategory::getLowerUpperEdge(rhobin).second;
                    if(RecoRhoBinCategory::getLowerUpperEdge(rhobin).second>0.9999){
                        // rhoBinPassed = ev.looseKinReco_rho >= RecoRhoBinCategory::getLowerUpperEdge(rhobin).first;
                        // rhoBinPassed = ev.kinReco_rho >= RecoRhoBinCategory::getLowerUpperEdge(rhobin).first;
                        // rhoBinPassed = ev.DNN_rho > RecoRhoBinCategory::getLowerUpperEdge(rhobin).first;
                        rhoBinPassed = ev.final_rho > RecoRhoBinCategory::getLowerUpperEdge(rhobin).first;
                        break;
                    }else{
                        // rhoBinPassed = ev.looseKinReco_rho >= RecoRhoBinCategory::getLowerUpperEdge(rhobin).first && ev.looseKinReco_rho < RecoRhoBinCategory::getLowerUpperEdge(rhobin).second;
                        // rhoBinPassed = ev.kinReco_rho >= RecoRhoBinCategory::getLowerUpperEdge(rhobin).first && ev.kinReco_rho < RecoRhoBinCategory::getLowerUpperEdge(rhobin).second;
                        // rhoBinPassed = ev.DNN_rho > RecoRhoBinCategory::getLowerUpperEdge(rhobin).first && ev.DNN_rho < RecoRhoBinCategory::getLowerUpperEdge(rhobin).second;
                        rhoBinPassed = ev.final_rho > RecoRhoBinCategory::getLowerUpperEdge(rhobin).first && ev.final_rho < RecoRhoBinCategory::getLowerUpperEdge(rhobin).second;
                        break;
                    }
                }
            }
            break;
        default:
            std::cerr<<"Error in MiniTreeAnalyzer::GetFillDecision()! undefined RecoRhoBinCategory\n...break\n"<<std::endl;
            exit(98);
            break;
    }

    return (nJetPassed && nBJetPassed && rhoBinPassed);

}


void MiniTreeAnalyzer::PushBackHistos(const BTagCategory::BTagCategory &btag,
        const RecoRhoBinCategory::RecoRhoBinCategory &rhobin, const NJetCategory::NJetCategory &njet, const TString histoname,
        const TString title, const int nbins, const float minX, const float maxX, const UserHisto::HistoKey &key, int &nHistos){

    if(!isSignal){
        auto h1 = make_shared<UserHisto::Histo>();
        h1->histogram = make_shared<TH1F>(histoname+Form("_%i",nHistos), histoname+title, nbins, minX, maxX);
        h1->variableKey = key;
        h1Map[btag][rhobin][njet][processType].push_back(h1);
        nHistos=nHistos+1;
    }else{
        for(unsigned int i = 0; i < nGenRhoBins; i++){
            auto h1 = make_shared<UserHisto::Histo>();
            h1->histogram = make_shared<TH1F>(histoname+Form("_%i",nHistos), histoname+title, nbins, minX, maxX);
            h1->variableKey = key;
            h1Map[btag][rhobin][njet][Process::convertFromBinning(genRhoBinning[i])].push_back(h1);
            nHistos=nHistos+1;
        }
        auto h1 = make_shared<UserHisto::Histo>();
        h1->histogram = make_shared<TH1F>(histoname+Form("_%i",nHistos), histoname+title, nbins, minX,maxX);
        h1->variableKey = key;
        h1Map[btag][rhobin][njet][Process::ttbarsignal_BKGNoAddJet].push_back(h1);
        nHistos=nHistos+1;
    }
}

void MiniTreeAnalyzer::PushBackHistos(const BTagCategory::BTagCategory &btag,
        const RecoRhoBinCategory::RecoRhoBinCategory &rhobin, const NJetCategory::NJetCategory &njet, const TString histoname,
        const TString title, const int nbins, const float* xbins, const UserHisto::HistoKey &key, int &nHistos){
    if(!isSignal){
        auto h1 = make_shared<UserHisto::Histo>();
        h1->histogram = make_shared<TH1F>(histoname+Form("_%i",nHistos), histoname+title, nbins, xbins);
        h1->variableKey = key;
        h1Map[btag][rhobin][njet][processType].push_back(h1);
        nHistos=nHistos+1;
    }else{
        for(unsigned int i = 0; i < nGenRhoBins; i++){
            auto h1 = make_shared<UserHisto::Histo>();
            h1->histogram = make_shared<TH1F>(histoname+Form("_%i",nHistos), histoname+title, nbins, xbins);
            h1->variableKey = key;
            h1Map[btag][rhobin][njet][Process::convertFromBinning(genRhoBinning[i])].push_back(h1);
            nHistos=nHistos+1;
        }
        auto h1 = make_shared<UserHisto::Histo>();
        h1->histogram = make_shared<TH1F>(histoname+Form("_%i",nHistos), histoname+title, nbins, xbins);
        h1->variableKey = key;
        h1Map[btag][rhobin][njet][Process::ttbarsignal_BKGNoAddJet].push_back(h1);
        nHistos=nHistos+1;
    }
}
void MiniTreeAnalyzer::PushBackHistos2D(const BTagCategory::BTagCategory &btag,
        const RecoRhoBinCategory::RecoRhoBinCategory &rhobin, const NJetCategory::NJetCategory &njet, const TString histoname,
        const TString title, const int nbins, const float minX, const float maxX, const UserHisto::HistoKey &keyX, const int nbinsY, const float minY,
        const float maxY, const UserHisto::HistoKey &keyY, int &nHistos){
    if(!isSignal){
        auto h2 = make_shared<UserHisto::Histo2D>();
        h2->histogram = make_shared<TH2F>(histoname+Form("_%i",nHistos), histoname+title, nbins, minX, maxX, nbinsY, minY, maxY);
        h2->variableKeyX = keyX;
        h2->variableKeyY = keyY;
        h2Map[btag][rhobin][njet][processType].push_back(h2);
        nHistos=nHistos+1;
    }else{
        for(unsigned int i = 0; i < nGenRhoBins; i++){
            auto h2 = make_shared<UserHisto::Histo2D>();
            h2->histogram = make_shared<TH2F>(histoname+Form("_%i",nHistos), histoname+title, nbins, minX, maxX, nbinsY, minY, maxY);
            h2->variableKeyX = keyX;
            h2->variableKeyY = keyY;
            h2Map[btag][rhobin][njet][Process::convertFromBinning(genRhoBinning[i])].push_back(h2);
            nHistos=nHistos+1;
        }
        auto h2 = make_shared<UserHisto::Histo2D>();
        h2->histogram = make_shared<TH2F>(histoname+Form("_%i",nHistos), histoname+title, nbins, minX, maxX, nbinsY, minY, maxY);
        h2->variableKeyX = keyX;
        h2->variableKeyY = keyY;
        h2Map[btag][rhobin][njet][Process::ttbarsignal_BKGNoAddJet].push_back(h2);
        nHistos=nHistos+1;
    }
}




void MiniTreeAnalyzer::InitHistograms(const vector<BTagCategory::BTagCategory> &btags,
                                    const vector<RecoRhoBinCategory::RecoRhoBinCategory> &rhobins, const vector<NJetCategory::NJetCategory> &njets){
    int nHistos(0);
    for(auto btag : btags){
        for(auto rhobin : rhobins){
            for(auto njet : njets){
                if(rhobin==RecoRhoBinCategory::NoRho){
                    if(njet==NJetCategory::ZeroAndOneJet){
                        if(btag==BTagCategory::OneBTag){
                            PushBackHistos(btag, rhobin, njet,"final", ";p_{T}^{b};Events",17, 30., 200., UserHisto::bJet1Pt, nHistos);
                        }
                        PushBackHistos(btag, rhobin, njet,"final", ";p_{t}^{first Jet};N_{Events}",12, 30, 150, UserHisto::jet1Pt, nHistos);
                        PushBackHistos(btag, rhobin, njet, "final", ";;N_{Events}",1, 0., 10000., UserHisto::nEvents, nHistos);
                        PushBackHistos(btag, rhobin, njet, "final", ";nJets;N_{Events}",7, 0, 7, UserHisto::nJets, nHistos);
                        PushBackHistos(btag, rhobin, njet, "final", ";nBJets;N_{Events}",5, 0, 5, UserHisto::nBJets, nHistos);
                        PushBackHistos(btag, rhobin, njet,"final", ";m_{ll};Events", 28, 20., 300., UserHisto::dileptonMass, nHistos);
                        PushBackHistos(btag, rhobin, njet,"final", ";MET;Events", 25, 0., 250., UserHisto::metPt, nHistos);

                    }
                    if(njet==NJetCategory::InclusiveNJet){
                        PushBackHistos(btag, rhobin, njet, "final", ";;N_{Events}",1, 0., 10000., UserHisto::nEvents, nHistos);
                        PushBackHistos(btag, rhobin, njet, "final", ";nJets;N_{Events}",7, 0, 7, UserHisto::nJets, nHistos);
                        PushBackHistos(btag, rhobin, njet, "final", ";nBJets;N_{Events}",5, 0, 5, UserHisto::nBJets, nHistos);
                        PushBackHistos(btag, rhobin, njet,"final", ";m_{ll};Events", 18, 20., 200., UserHisto::dileptonMass, nHistos);
                        PushBackHistos(btag, rhobin, njet,"final", ";p_{t}^{first Jet};N_{Events}",27, 30, 300, UserHisto::jet1Pt, nHistos);
                        PushBackHistos(btag, rhobin, njet, "final", ";p_{t}^{MET}; N_{Events}", 60, 0, 300, UserHisto::metPt, nHistos);
                        PushBackHistos(btag, rhobin, njet, "final", ";p_{t}^{second Jet};N_{Events}", 60, 0, 300, UserHisto::jet2Pt, nHistos);
                        PushBackHistos(btag, rhobin, njet, "final", ";p_{t}^{third Jet};N_{Events}", 60, 0, 300, UserHisto::jet3Pt, nHistos);
                        PushBackHistos(btag, rhobin, njet,"final", ";p_{T}^{b};Events",17, 30., 200., UserHisto::bJet1Pt, nHistos);

                    }
                    if(njet==NJetCategory::TwoJet){
                        if(btag==BTagCategory::ZeroPlusGreaterTwoBTag){
                            PushBackHistos(btag, rhobin, njet,"final", ";;N_{Events}",1, 0., 10000., UserHisto::nEvents, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nJets;N_{Events}",7, 0, 7, UserHisto::nJets, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nBJets;N_{Events}",5, 0, 5, UserHisto::nBJets, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";m_{ll};N_{Events}", 18, 20., 200., UserHisto::dileptonMass, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";p_{t}^{second Jet};N_{Events}",12, 30, 150, UserHisto::jet2Pt, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";p_{t}^{first Jet};N_{Events}",12, 30, 150, UserHisto::jet1Pt, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";MET;Events", 25, 0., 250., UserHisto::metPt, nHistos);
                        }
                        if(btag==BTagCategory::OneBTag){
                            PushBackHistos(btag, rhobin, njet,"final", ";;N_{Events}",1, 0., 10000., UserHisto::nEvents, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nJets;N_{Events}",7, 0, 7, UserHisto::nJets, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nBJets;N_{Events}",5, 0, 5, UserHisto::nBJets, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";p_{T}^{b};Events",17, 30., 200., UserHisto::bJet1Pt, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";MET;Events", 25, 0., 250., UserHisto::metPt, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";m_{lb}^{min};Events", 6, 0., 120., UserHisto::mlb_min, nHistos);
                        }
                        if(btag==BTagCategory::TwoBTag){
                            PushBackHistos(btag, rhobin, njet,"final", ";;N_{Events}",1, 0., 10000., UserHisto::nEvents, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nJets;N_{Events}",7, 0, 7, UserHisto::nJets, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nBJets;N_{Events}",5, 0, 5, UserHisto::nBJets, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";m_{lb}^{min};Events", 6, 0., 120., UserHisto::mlb_min, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";MET;Events", 25, 0., 250., UserHisto::metPt, nHistos);
                        }
                        if(btag==BTagCategory::InclusiveBTag){
                            PushBackHistos(btag, rhobin, njet,"final", ";;N_{Events}",1, 0., 10000., UserHisto::nEvents, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nJets;N_{Events}",7, 0, 7, UserHisto::nJets, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nBJets;N_{Events}",5, 0, 5, UserHisto::nBJets, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";MET;Events", 25, 0., 250., UserHisto::metPt, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";m_{lb}^{min};Events", 6, 0., 120., UserHisto::mlb_min, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";p_{T}^{b};Events",17, 30., 200., UserHisto::bJet1Pt, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";p_{t}^{second Jet};N_{Events}",12, 30, 150, UserHisto::jet2Pt, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";p_{t}^{first Jet};N_{Events}",12, 30, 150, UserHisto::jet1Pt, nHistos);
                        }
                    }
                }
                if(rhobin==RecoRhoBinCategory::Rho0to02 || rhobin==RecoRhoBinCategory::Rho02to03 ||
                rhobin==RecoRhoBinCategory::Rho03to04 || rhobin==RecoRhoBinCategory::Rho04to05 ||
                rhobin==RecoRhoBinCategory::Rho05to06 ||
                rhobin==RecoRhoBinCategory::Rho06to1 ||
                rhobin==RecoRhoBinCategory::Rho0to1){
                    if(njet==NJetCategory::GreaterTwoJet){
                        if(btag==BTagCategory::ZeroPlusGreaterTwoBTag){
                            PushBackHistos(btag, rhobin, njet,"final", ";;N_{Events}",1, 0., 10000., UserHisto::nEvents, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nJets;N_{Events}",7, 0, 7, UserHisto::nJets, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nBJets;N_{Events}",5, 0, 5, UserHisto::nBJets, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";m_{ll};Events", 4, 0., 160., UserHisto::dileptonMass, nHistos);
                            // PushBackHistos(btag, rhobin, njet, "RhoStudies_RecoRhoDNN", ";#rho", 1000., 0., 1., UserHisto::recoRhoDNN, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";p_{T}^{third Jet};Events", 2, 50., 100., UserHisto::jet3Pt, nHistos);
                            PushBackHistos(btag, rhobin, njet, "RhoStudies_RecoRhoDNN", ";#rho", nRecoRhoBins, recoRhoBinning, UserHisto::recoRhoDNN, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";MET;Events", 25, 0., 250., UserHisto::metPt, nHistos);
                        }
                        if(btag==BTagCategory::InclusiveBTag){
                            PushBackHistos(btag, rhobin, njet,"final", ";;N_{Events}",1, 0., 10000., UserHisto::nEvents, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nJets;N_{Events}",7, 0, 7, UserHisto::nJets, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nBJets;N_{Events}",5, 0, 5, UserHisto::nBJets, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";m_{ll};Events", 9, 20., 200., UserHisto::dileptonMass, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";p_{T}^{third Jet};Events", 5, 50., 150., UserHisto::jet3Pt, nHistos);
                            PushBackHistos(btag, rhobin, njet, "RhoStudies_RecoRhoDNN", ";#rho", nRecoRhoBins, recoRhoBinning, UserHisto::recoRhoDNN, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";MET;Events", 25, 0., 250., UserHisto::metPt, nHistos);
                        }
                        if(btag==BTagCategory::OneBTag){
                            PushBackHistos(btag, rhobin, njet,"final", ";;N_{Events}",1, 0., 10000., UserHisto::nEvents, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nJets;N_{Events}",7, 0, 7, UserHisto::nJets, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nBJets;N_{Events}",5, 0, 5, UserHisto::nBJets, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";p_{T}^{third Jet};Events", 2, 50., 100., UserHisto::jet3Pt, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";m_{ll};Events", 4, 0., 160., UserHisto::dileptonMass, nHistos);
                            PushBackHistos(btag, rhobin, njet, "RhoStudies_RecoRhoDNN", ";#rho", nRecoRhoBins, recoRhoBinning, UserHisto::recoRhoDNN, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";MET;Events", 25, 0., 250., UserHisto::metPt, nHistos);
                        }
                        if(btag==BTagCategory::TwoBTag){
                            PushBackHistos(btag, rhobin, njet,"final", ";;N_{Events}", 1, 0., 10000., UserHisto::nEvents, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nJets;N_{Events}", 7, 0, 7, UserHisto::nJets, nHistos);
                            PushBackHistos(btag, rhobin, njet, "final", ";nBJets;N_{Events}", 5, 0, 5, UserHisto::nBJets, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";m_{lb};Events", 2, 0., 140., UserHisto::mlb_min, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";p_{T}^{third Jet};Events", 2, 50., 100., UserHisto::jet3Pt, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";m_{ll};Events", 4, 0., 160., UserHisto::dileptonMass, nHistos);
                            PushBackHistos(btag, rhobin, njet, "RhoStudies_RecoRhoDNN", ";#rho", nRecoRhoBins, recoRhoBinning, UserHisto::recoRhoDNN, nHistos);
                            PushBackHistos(btag, rhobin, njet,"final", ";MET;Events", 25, 0., 250., UserHisto::metPt, nHistos);
                        }
                    }

                }
                if(rhobin==RecoRhoBinCategory::InclusiveRho){
                    PushBackHistos(btag, rhobin, njet, "final", ";;N_{Events}",1, 0., 10000., UserHisto::nEvents, nHistos);
                    PushBackHistos(btag, rhobin, njet, "final", ";nJets;N_{Events}",7, 0, 7, UserHisto::nJets, nHistos);
                    PushBackHistos(btag, rhobin, njet, "final", ";nBJets;N_{Events}",5, 0, 5, UserHisto::nBJets, nHistos);
                    // CPs Inputs
                    PushBackHistos(btag, rhobin, njet, "final", ";p_{t}^{MET}; N_{Events}", 60, 0, 300, UserHisto::metPt, nHistos);
                    PushBackHistos(btag, rhobin, njet, "final", ";p_{t}^{first Jet};N_{Events}", 60, 0, 300, UserHisto::jet1Pt, nHistos);
                    PushBackHistos(btag, rhobin, njet, "final", ";p_{t}^{second Jet};N_{Events}", 60, 0, 300, UserHisto::jet2Pt, nHistos);
                    PushBackHistos(btag, rhobin, njet, "final", ";p_{t}^{third Jet};N_{Events}", 60, 0, 300, UserHisto::jet3Pt, nHistos);



                }
                if(rhobin==RecoRhoBinCategory::RhoStudies){
                    // PushBackHistos2D(btag, rhobin, njet, "RhoStudies_GenRecoDiff2D", ";#rho_{True};#rho_{True}-#rho_{Reco}",100., 0., 1., UserHisto::genRho, 200., -1., 1., UserHisto::genLKRSRecoRhoDiff, nHistos);
                    // PushBackHistos2D(btag, rhobin, njet, "RhoStudies_GenReco2D", ";#rho_{True};#rho_{Reco}", 100., 0., 1., UserHisto::genRho, 100., 0., 1., UserHisto::recoRhoLKRS, nHistos);
                    // PushBackHistos(btag, rhobin, njet, "RhoStudies_RecoRho", ";#rho", 1000., 0., 1., UserHisto::recoRhoLKRS, nHistos);
                    // PushBackHistos2D(btag, rhobin, njet, "RhoStudies_GenRecoDiff2D", ";#rho_{True};#rho_{True}-#rho_{Reco}",100., 0., 1., UserHisto::genRho, 200., -1., 1., UserHisto::genKRSRecoRhoDiff, nHistos);
                    // PushBackHistos2D(btag, rhobin, njet, "RhoStudies_GenReco2D", ";#rho_{True};#rho_{Reco}", 100., 0., 1., UserHisto::genRho, 100., 0., 1., UserHisto::recoRhoKRS, nHistos);
                    // PushBackHistos(btag, rhobin, njet, "RhoStudies_RecoRho", ";#rho", 1000., 0., 1., UserHisto::recoRhoKRS, nHistos);
                    PushBackHistos2D(btag, rhobin, njet, "RhoStudies_GenRecoDiff2D", ";#rho_{True};#rho_{True}-#rho_{Reco}",100., 0., 1., UserHisto::genRho, 200., -1., 1., UserHisto::genRecoRhoDNNDiff, nHistos);
                    PushBackHistos2D(btag, rhobin, njet, "RhoStudies_GenReco2D", ";#rho_{True};#rho_{Reco}", 100., 0., 1., UserHisto::genRho, 100., 0., 1., UserHisto::recoRhoDNN, nHistos);
                    PushBackHistos(btag, rhobin, njet, "RhoStudies_RecoRhoDNN", ";#rho", 1000., 0., 1., UserHisto::recoRhoDNN, nHistos);
                    // gen plots
                    PushBackHistos2D(btag, rhobin, njet, "GenMttbar_vs_genRho", ";#rho_{True};#m_{ttbar}^{true}", 100., 0., 1., UserHisto::genRho, 1000., 0., 1000., UserHisto::genmTTbar, nHistos);

                }
            }
        }
    }
}


void MiniTreeAnalyzer::FillHistograms(unordered_map<BTagCategory::BTagCategory,unordered_map<RecoRhoBinCategory::RecoRhoBinCategory,unordered_map<NJetCategory::NJetCategory,unordered_map<Process::Process,vector<shared_ptr<UserHisto::Histo>>>>>> &histMap){
        for(auto &hMapIt : histMap){
            for(auto &hMapIt2 : hMapIt.second){
                    for(auto &hMapIt3: hMapIt2.second){
                        for(auto &foldername: hMapIt3.second){
                            for(auto &histo: foldername.second){
                                if(GetCategoryFillDecision(*event, hMapIt.first, hMapIt2.first, hMapIt3.first)){
                                    if(event->isTTJSignal){
                                        for(unsigned int i = 0; i < nGenRhoBins; i++){
                                        if(foldername.first==Process::convertFromBinning(genRhoBinning[i])){
                                            float cutvalueGenRhoLow =genRhoBinning[i];
                                            float cutvalueGenRhoHigh =genRhoBinning[i+1];
                                            if((event->gen_rho < cutvalueGenRhoHigh && event->gen_rho >= cutvalueGenRhoLow)||(((int) 10.*cutvalueGenRhoHigh)==10 && event->gen_rho>=cutvalueGenRhoHigh)){
                                                if(hMapIt2.first==RecoRhoBinCategory::NoRho){
                                                    // if(hMapIt3.first == NJetCategory::NJetCategory::ZeroAndOneJet){
                                                    if(UserHisto::CheckKeyConditionFromEvent(*event,histo->variableKey)){
                                                        histo->histogram->Fill(UserHisto::GetValueToKeyFromEvent(*event,histo->variableKey),event->weightWithBtagSF);
                                                    }
                                                        // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                    // }
                                                    // if(hMapIt3.first == NJetCategory::NJetCategory::TwoJet){
                                                        // if(hMapIt.first==BTagCategory::ZeroPlusGreaterTwoBTag){
                                                            // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                        // }
                                                        // if(hMapIt.first==BTagCategory::OneBTag){
                                                            // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                        // }
                                                        // if(hMapIt.first==BTagCategory::TwoBTag){
                                                            // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                        // }
                                                    // }
                                                }
                                                if(
                                                    hMapIt2.first==RecoRhoBinCategory::Rho0to02 || hMapIt2.first==RecoRhoBinCategory::Rho02to03 ||
                                                    hMapIt2.first==RecoRhoBinCategory::Rho03to04 || hMapIt2.first==RecoRhoBinCategory::Rho04to05 ||
                                                    hMapIt2.first==RecoRhoBinCategory::Rho05to06 ||
                                                    hMapIt2.first==RecoRhoBinCategory::Rho06to1 ||
                                                    hMapIt2.first==RecoRhoBinCategory::Rho0to1){
                                                    // if(hMapIt.first==BTagCategory::ZeroPlusGreaterTwoBTag){
                                                    if(UserHisto::CheckKeyConditionFromEvent(*event,histo->variableKey)){
                                                        histo->histogram->Fill(UserHisto::GetValueToKeyFromEvent(*event,histo->variableKey),event->weightWithBtagSF);
                                                    }
                                                    // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                    // }
                                                    // if(hMapIt.first==BTagCategory::OneBTag){
                                                        // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                    // }
                                                    // if(hMapIt.first==BTagCategory::TwoBTag){
                                                        // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                    // }
                                                }
                                                if(hMapIt2.first==RecoRhoBinCategory::InclusiveRho){
                                                    // if(CheckKeyConditionFromEvent(*event,histo.variableKey)){
                                                        if(UserHisto::CheckKeyConditionFromEvent(*event,histo->variableKey)){
                                                            histo->histogram->Fill(UserHisto::GetValueToKeyFromEvent(*event,histo->variableKey),event->weightWithBtagSF);
                                                        }
                                                        // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                    // }
                                                }
                                                if(hMapIt2.first==RecoRhoBinCategory::RhoStudies){
                                                    if(UserHisto::CheckKeyConditionFromEvent(*event,histo->variableKey)){
                                                        histo->histogram->Fill(UserHisto::GetValueToKeyFromEvent(*event,histo->variableKey),event->weightWithBtagSF);
                                                    }
                                                    // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                }
                                            }
                                        }
                                    }
                                    }else{
                                        if(event->isTTSignalBackground){
                                        if(foldername.first==Process::ttbarsignal_BKGNoAddJet){
                                            if(hMapIt2.first==RecoRhoBinCategory::NoRho){
                                                // if(hMapIt3.first == NJetCategory::NJetCategory::ZeroAndOneJet){
                                                if(UserHisto::CheckKeyConditionFromEvent(*event,histo->variableKey)){
                                                    histo->histogram->Fill(UserHisto::GetValueToKeyFromEvent(*event,histo->variableKey),event->weightWithBtagSF);
                                                }
                                                    // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                // }
                                                // if(hMapIt3.first == NJetCategory::NJetCategory::TwoJet){
                                                    // if(hMapIt.first==BTagCategory::ZeroPlusGreaterTwoBTag){
                                                        // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                    // }
                                                    // if(hMapIt.first==BTagCategory::OneBTag){
                                                        // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                    // }
                                                    // if(hMapIt.first==BTagCategory::TwoBTag){
                                                        // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                    // }
                                                // }
                                            }if(hMapIt2.first==RecoRhoBinCategory::Rho0to02 || hMapIt2.first==RecoRhoBinCategory::Rho02to03 ||
                                            hMapIt2.first==RecoRhoBinCategory::Rho03to04 || hMapIt2.first==RecoRhoBinCategory::Rho04to05 ||
                                            hMapIt2.first==RecoRhoBinCategory::Rho05to06 ||
                                            hMapIt2.first==RecoRhoBinCategory::Rho06to1 ||
                                            hMapIt2.first==RecoRhoBinCategory::Rho0to1){
                                                // if(hMapIt.first==BTagCategory::ZeroPlusGreaterTwoBTag){
                                                if(UserHisto::CheckKeyConditionFromEvent(*event,histo->variableKey)){
                                                    histo->histogram->Fill(UserHisto::GetValueToKeyFromEvent(*event,histo->variableKey),event->weightWithBtagSF);
                                                }
                                                    // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                // }
                                                // if(hMapIt.first==BTagCategory::OneBTag){
                                                //     histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                // }
                                                // if(hMapIt.first==BTagCategory::TwoBTag){
                                                //     histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                // }
                                            }
                                            if(hMapIt2.first==RecoRhoBinCategory::InclusiveRho){
                                                // if(CheckKeyConditionFromEvent(*event,histo.variableKey)){
                                                    if(UserHisto::CheckKeyConditionFromEvent(*event,histo->variableKey)){
                                                        histo->histogram->Fill(UserHisto::GetValueToKeyFromEvent(*event,histo->variableKey),event->weightWithBtagSF);
                                                    }
                                                    // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                // }
                                                    // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                            }
                                            if(hMapIt2.first==RecoRhoBinCategory::RhoStudies){
                                                if(UserHisto::CheckKeyConditionFromEvent(*event,histo->variableKey)){
                                                    histo->histogram->Fill(UserHisto::GetValueToKeyFromEvent(*event,histo->variableKey),event->weightWithBtagSF);
                                                }
                                                // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                            }
                                        }}else{
                                            if(hMapIt2.first==RecoRhoBinCategory::NoRho){
                                                // if(hMapIt3.first == NJetCategory::NJetCategory::ZeroAndOneJet){
                                                if(UserHisto::CheckKeyConditionFromEvent(*event,histo->variableKey)){
                                                    histo->histogram->Fill(UserHisto::GetValueToKeyFromEvent(*event,histo->variableKey),event->weightWithBtagSF);
                                                }
                                                    // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                // }
                                                // if(hMapIt3.first == NJetCategory::NJetCategory::TwoJet){
                                                //     if(hMapIt.first==BTagCategory::ZeroPlusGreaterTwoBTag){
                                                //         histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                //     }
                                                //     if(hMapIt.first==BTagCategory::OneBTag){
                                                //         histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                //     }
                                                //     if(hMapIt.first==BTagCategory::TwoBTag){
                                                //         histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                //     }
                                                // }
                                            }if(hMapIt2.first==RecoRhoBinCategory::Rho0to02 || hMapIt2.first==RecoRhoBinCategory::Rho02to03 ||
                                            hMapIt2.first==RecoRhoBinCategory::Rho03to04 || hMapIt2.first==RecoRhoBinCategory::Rho04to05 ||
                                            hMapIt2.first==RecoRhoBinCategory::Rho05to06 ||
                                            hMapIt2.first==RecoRhoBinCategory::Rho06to1 ||
                                            hMapIt2.first==RecoRhoBinCategory::Rho0to1){
                                                // if(hMapIt.first==BTagCategory::ZeroPlusGreaterTwoBTag){
                                                // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                if(UserHisto::CheckKeyConditionFromEvent(*event,histo->variableKey)){
                                                    histo->histogram->Fill(UserHisto::GetValueToKeyFromEvent(*event,histo->variableKey),event->weightWithBtagSF);
                                                }
                                                // }
                                                // if(hMapIt.first==BTagCategory::OneBTag){
                                                //     histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                // }
                                                // if(hMapIt.first==BTagCategory::TwoBTag){
                                                //     histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                                // }
                                            }
                                            if(hMapIt2.first==RecoRhoBinCategory::InclusiveRho){
                                                if(UserHisto::CheckKeyConditionFromEvent(*event,histo->variableKey)){
                                                    histo->histogram->Fill(UserHisto::GetValueToKeyFromEvent(*event,histo->variableKey),event->weightWithBtagSF);
                                                }
                                                        // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                            }
                                            if(hMapIt2.first==RecoRhoBinCategory::RhoStudies){
                                                if(UserHisto::CheckKeyConditionFromEvent(*event,histo->variableKey)){
                                                    histo->histogram->Fill(UserHisto::GetValueToKeyFromEvent(*event,histo->variableKey),event->weightWithBtagSF);
                                                }
                                                // histo.histogram.Fill(UserHisto::GetValueToKeyFromEvent(*event,histo.variableKey),event->weightWithBtagSF);
                                            }
                                        }
                                }
                        }
                    }
                }
                }
            }
        }
}

void MiniTreeAnalyzer::FillHistograms2D(unordered_map<BTagCategory::BTagCategory,unordered_map<RecoRhoBinCategory::RecoRhoBinCategory,unordered_map<NJetCategory::NJetCategory,unordered_map<Process::Process,vector<shared_ptr<UserHisto::Histo2D>>>>>> &histMap){
    for(auto &hMapIt : histMap){
        for(auto &hMapIt2 : hMapIt.second){
            for(auto &hMapIt3: hMapIt2.second){
                for(auto &foldername: hMapIt3.second){
                    for(auto &histo: foldername.second){
                        if(GetCategoryFillDecision(*event, hMapIt.first, hMapIt2.first, hMapIt3.first)){
                            if(event->isTTJSignal){
                                for(unsigned int i = 0; i < nGenRhoBins; i++){
                                    if(foldername.first==Process::convertFromBinning(genRhoBinning[i])){
                                        float cutvalueGenRhoLow =genRhoBinning[i];
                                        float cutvalueGenRhoHigh =genRhoBinning[i+1];
                                        if((event->gen_rho < cutvalueGenRhoHigh && event->gen_rho >= cutvalueGenRhoLow)||(cutvalueGenRhoHigh==1. && event->gen_rho>=cutvalueGenRhoHigh)){
                                            if(hMapIt2.first==RecoRhoBinCategory::Rho0to02 || hMapIt2.first==RecoRhoBinCategory::Rho02to03 ||
                                            hMapIt2.first==RecoRhoBinCategory::Rho03to04 || hMapIt2.first==RecoRhoBinCategory::Rho04to05 ||
                                            hMapIt2.first==RecoRhoBinCategory::Rho05to06 ||
                                            hMapIt2.first==RecoRhoBinCategory::Rho06to1 ||
                                            hMapIt2.first==RecoRhoBinCategory::Rho0to1){
                                            }
                                            if(hMapIt2.first==RecoRhoBinCategory::RhoStudies){
                                                histo->histogram->Fill(UserHisto::GetValueToKeyFromEvent(*event,histo->variableKeyX),UserHisto::GetValueToKeyFromEvent(*event,histo->variableKeyY),event->weightWithBtagSF);
                                                histo->histogram->Fill(UserHisto::GetValueToKeyFromEvent(*event,histo->variableKeyX),UserHisto::GetValueToKeyFromEvent(*event,histo->variableKeyY),event->weightWithBtagSF);
                                            }
                                        }
                                    }
                                }
                            }else{
                            }
                        }
                    }
                }
            }
        }
    }
}



bool MiniTreeAnalyzer::Analyze(){

    TFile* file = TFile::Open(inputFilePath, "READ");
    TTree* eventTree = (TTree*)file->Get("miniTree");
    MiniTreeReader myReader(*eventTree, systematic);

    InitHistograms(NBJetCategories, RecoRhoBinCategories, NJetCategories);

    std::cout<<"All histograms initialized, begin with event loop:"<<std::endl;
    ProgressBar(0.);
    int processed = 0;
    for(Int_t i=0; i<eventTree->GetEntries(); ++i) {
        eventTree->GetEntry(i);
        event->Clear();
        event->Fill(myReader, isMC, isSignal, systematic, hasPSWeights, hasPDFWeights, hasMEWeights, hasBFragWeights,hasBSemilepWeights);
        if(event->hasKRS && event->kinReco_rho>0.){
            event->final_rho = event->kinReco_rho;
        }else{
            EvaluateDNN(*event);
            event->final_rho = event->DNN_rho;
        }
        // if(event->passStep3){
            // h_eventsPassedStep3->Fill(0.5,event->weight);
            // if(event->passStep4){
            //     h_eventsPassedStep4->Fill(0.5,event->weight);
            //     if(event->passStep5){
            //         h_eventsPassedStep5->Fill(0.5,event->weight);
            //         if(event->passStep6){
            //             h_eventsPassedStep6->Fill(0.5,event->weight);
            //             if(event->passStep7){
            //                 h_eventsPassedStep7->Fill(0.5,event->weightWithBtagSF);
            //                 if(event->passStep8){
            //                     h_eventsPassedStep8->Fill(0.5,event->weightWithBtagSF * event->kinReco_SF);
            //                 }
            //                 if(event->passStep8Loose){
            //                     h_eventsPassedStep8L->Fill(0.5,event->weightWithBtagSF * event->looseKinReco_SF);
            //                 }
            //             }
            //         }
            //     }
            // }
        // }

        FillHistograms(h1Map);
        FillHistograms2D(h2Map);

        processed++;
        if(processed % 1000 == 0) {
            int progress = 100*(float)processed/eventTree->GetEntries();
            ProgressBar(progress);
        }
        continue;
    }
    ProgressBar(100);
    std::chrono::steady_clock::time_point endLoop = std::chrono::steady_clock::now();
    std::cout<<"Event loop finished in "<<std::chrono::duration_cast<std::chrono::seconds>(endLoop - start).count()<<" seconds for "<<eventTree->GetEntries()<<" Events - (ave "<<(float)eventTree->GetEntries()/((float)std::chrono::duration_cast<std::chrono::seconds>(endLoop - start).count())<<" Events/s)"<<std::endl;

    eventTree->SetBranchStatus("*", 1);

    file->Close();
    return true;

}
template<typename T>
void MiniTreeAnalyzer::save2File(unordered_map<BTagCategory::BTagCategory,unordered_map<RecoRhoBinCategory::RecoRhoBinCategory,unordered_map<NJetCategory::NJetCategory,unordered_map<Process::Process,vector<T>>>>> & hMaps, TFile& file)
{
    if (!file.Get(RenameSystematicString(systematic))) {
        file.mkdir(RenameSystematicString(systematic));
    }
    file.cd(RenameSystematicString(systematic));
    if (!file.Get(RenameSystematicString(systematic)+"/"+Channel::convert(channel))) {
            file.mkdir(RenameSystematicString(systematic)+"/"+Channel::convert(channel));
    }
    file.cd(RenameSystematicString(systematic)+"/"+Channel::convert(channel));
    for(auto& hMapIt3 : hMaps) {
        if (!file.Get(RenameSystematicString(systematic)+"/"+Channel::convert(channel)+"/"+BTagCategory::convert(hMapIt3.first))) {
                file.mkdir(RenameSystematicString(systematic)+"/"+Channel::convert(channel)+"/"+BTagCategory::convert(hMapIt3.first));
        }
        file.cd(RenameSystematicString(systematic)+"/"+Channel::convert(channel)+"/"+BTagCategory::convert(hMapIt3.first));
        for(auto& hMapIt4 : hMapIt3.second) {
            if (!file.Get(RenameSystematicString(systematic)+"/"+Channel::convert(channel)+"/"+BTagCategory::convert(hMapIt3.first)+"/"+RecoRhoBinCategory::convert(hMapIt4.first))){
                    file.mkdir(RenameSystematicString(systematic)+"/"+Channel::convert(channel)+"/"+BTagCategory::convert(hMapIt3.first)+"/"+RecoRhoBinCategory::convert(hMapIt4.first));
            }
            file.cd(RenameSystematicString(systematic)+"/"+Channel::convert(channel)+"/"+BTagCategory::convert(hMapIt3.first)+"/"+RecoRhoBinCategory::convert(hMapIt4.first));
            for(auto& hMapIt5 : hMapIt4.second) {
                if (!file.Get(RenameSystematicString(systematic)+"/"+Channel::convert(channel)+"/"+BTagCategory::convert(hMapIt3.first)+"/"+RecoRhoBinCategory::convert(hMapIt4.first)+"/"+NJetCategory::convert(hMapIt5.first))){
                    file.mkdir(RenameSystematicString(systematic)+"/"+Channel::convert(channel)+"/"+BTagCategory::convert(hMapIt3.first)+"/"+RecoRhoBinCategory::convert(hMapIt4.first)+"/"+NJetCategory::convert(hMapIt5.first));
                }
                file.cd(RenameSystematicString(systematic)+"/"+Channel::convert(channel)+"/"+BTagCategory::convert(hMapIt3.first)+"/"+RecoRhoBinCategory::convert(hMapIt4.first)+"/"+NJetCategory::convert(hMapIt5.first));
                for (auto& hMapIt6: hMapIt5.second) {
                    if (!file.Get(RenameSystematicString(systematic)+"/"+Channel::convert(channel)+"/"+BTagCategory::convert(hMapIt3.first)+"/"+RecoRhoBinCategory::convert(hMapIt4.first)+"/"+NJetCategory::convert(hMapIt5.first)+"/"+Process::convert(hMapIt6.first))){
                        file.mkdir(RenameSystematicString(systematic)+"/"+Channel::convert(channel)+"/"+BTagCategory::convert(hMapIt3.first)+"/"+RecoRhoBinCategory::convert(hMapIt4.first)+"/"+NJetCategory::convert(hMapIt5.first)+"/"+Process::convert(hMapIt6.first));
                    }
                    file.cd(RenameSystematicString(systematic)+"/"+Channel::convert(channel)+"/"+BTagCategory::convert(hMapIt3.first)+"/"+RecoRhoBinCategory::convert(hMapIt4.first)+"/"+NJetCategory::convert(hMapIt5.first)+"/"+Process::convert(hMapIt6.first));
                    for (auto& histo: hMapIt6.second) {
                        histo->histogram->Scale(lumiWeight);
                        histo->histogram->Scale(renormWeight);
                        FixHistogram(*histo->histogram);
                        histo->histogram->Write(histo->GetName(), TObject::kWriteDelete);

                    }
                    file.cd();
                }
            }
        }
    }
}

void MiniTreeAnalyzer::FixHistogram(TH1& histo){
    if(histo.InheritsFrom(TH1::Class())){
        for(int x=0; x<histo.GetNbinsX()+2; x++){
            if(systematic.type()==Systematic::nominal){
            if(histo.GetBinContent(x) < 0.){
                histo.SetBinContent(x, 0.);
            }
        }else{
            if(histo.GetBinContent(x) < 0.01){
                histo.SetBinContent(x, 0.01);
            }
        }
        }
    }else{
        const std::string log("No TH1, should not be here!");
        throw std::runtime_error("MiniTreeAnalyzer::FixHistogram -- "+log);
    }
}

void MiniTreeAnalyzer::FixHistogram(TH2& histo){
    if(histo.InheritsFrom(TH2::Class())){
        for(int x=0; x<histo.GetNbinsX()+2; x++){
            for(int y=0; y<histo.GetNbinsY()+2; y++){
                if(histo.GetBinContent(x, y) < 0.){
                    histo.SetBinContent(x, y, 0.);
                }
            }
        }
    }else{
        const std::string log("No TH2, should not be here!");
        throw std::runtime_error("MiniTreeAnalyzer::FixHistogram -- "+log);
    }
}



void MiniTreeAnalyzer::WriteOutput(){
    std::cout<<"Writing Output..."<<std::endl;
    TFile* file = TFile::Open(outFilePath, "RECREATE");
    save2File(h1Map, *file);
    save2File(h2Map, *file);

    // h_eventsPassedStep3->Scale(lumiWeight*renormWeight);
    // h_eventsPassedStep4->Scale(lumiWeight*renormWeight);
    // h_eventsPassedStep5->Scale(lumiWeight*renormWeight);
    // h_eventsPassedStep6->Scale(lumiWeight*renormWeight);
    // h_eventsPassedStep7->Scale(lumiWeight*renormWeight);
    // h_eventsPassedStep8->Scale(lumiWeight*renormWeight);
    // h_eventsPassedStep8L->Scale(lumiWeight*renormWeight);
    //
    // h_eventsPassedStep3->Write();
    // h_eventsPassedStep4->Write();
    // h_eventsPassedStep5->Write();
    // h_eventsPassedStep6->Write();
    // h_eventsPassedStep7->Write();
    // h_eventsPassedStep8->Write();
    // h_eventsPassedStep8L->Write();

    file->Write();
    file->Close();
    end = std::chrono::steady_clock::now();
    std::cout << "Finished analyzer in " << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() <<" seconds!"<< std::endl;
    std::cout << "Output file created: " + outFilePath << std::endl;
    std::cout<<std::endl;
}
