#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>

#include <TString.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1D.h>
#include <TFile.h>
#include <TSystem.h>

#include "../../common/include/RootFileReader.h"
#include "../../common/include/utils.h"
#include "EffHist.h"
#include "PlotterConfigurationHelper.h"
#include "../../common/include/CommandLineParameters.h"


using namespace std;

double sqr(double value) { return value*value; }



int main(int argc, char** argv)
{

  CLParameter<std::string> opt_c("c", "Specify a certain channel (ee, emu, mumu). No channel specified = run on all channels", false, 1, 1, common::makeStringCheck(Channel::convert(Channel::allowedChannelsAnalysis)));
  CLParameter<std::string> opt_year("y", "Specify year, valid: 2016preVFP, 2016postVFP, 2017, 2018 ", true, 1, 1, common::makeStringCheck({"2016preVFP", "2016postVFP", "2017", "2018"}));

  CLAnalyser::interpretGlobal(argc, argv);

  TString year;
  if (opt_year.isSet()) year = opt_year.getArguments()[0];
  std::cout << "Processing setup for year: " << year << std::endl;

  std::vector<std::string> channels {"emu", "ee", "mumu"};
  if (opt_c.isSet()) channels = opt_c.getArguments();
  std::cout << "Processing channels: ";
  for (auto ch: channels) std::cout << ch << " ";
  std::cout << "\n";

  TString kinRecoSF("const std::map<Channel::Channel, double> m_sfNominal { ");
  TString kinRecoSFerr("const std::map<Channel::Channel, double> m_sfUnc { ");

  gSystem->Exec("mkdir -p KinRecoEffSFPlots_"+year);

  RootFileReader *fileReader_ = RootFileReader::getInstance();

  vector<TString> systematics {"Nominal"};

  TString rootFilePrePath("");

  vector<TString> effHistNames {"KinReco_nRecoEvt","KinReco_nRecoEvt","_Eff","_vs_JetMult","_vs_LepEta","_vs_JetEta","_vs_LeppT","_vs_MET","_vs_LepEta2","_vs_JetEta2","_vs_LeppT2","_vs_AntiLeptonpT","_vs_LeptonpT","_vs_LeptonEta","_vs_AntiLeptonEta","_vs_JetpT","_vs_JetpT2"};

  double NdataAfter = 0;
  double NdataBefore = 0;
  double dataEff = 0;
  double dataEffUnc = 0;

  double NmcAfter = 0;
  double NmcBefore = 0;
  double allmcEff =  0;
  double allmcEffUnc = 0;

  double sfUnc = 0;

  char dataEffString[100];
  char allmcEffString[100];
  char sfString[100];

  PlotterConfigurationHelper* plotterConfigHelper = new PlotterConfigurationHelper(fileReader_, year, year, 1, 1, 0);

  for ( auto systematic : systematics ){
              gSystem->Exec("mkdir -p KinRecoEffSFPlots_"+year+"/"+systematic);
    for ( TString channel : channels ){
        gSystem->Exec("mkdir -p KinRecoEffSFPlots_"+year+"/"+systematic+"/"+channel);

        std::vector<TString> dataset;
        std::vector<TString> legends;
        int channelType(-1);
        TString channelLabel;
        if(channel.Contains("ee")){channelLabel="ee";channelType=0;}
        if(channel.Contains("mumu")){channelLabel="#mu#mu";channelType=1;}
        if(channel.Contains("emu")){channelLabel="e#mu";channelType=2;}
        if(channel.Contains("combined")){channelLabel="Dilepton Combined";channelType=3;}

        if (channelType == -1){
          printf("ERROR: No channel defined!\n");
          break;
        }

        std::vector<double> DYScale={1,1,1,1};
        plotterConfigHelper->DYScaleFactor("Standard",DYScale,"");


        for(Size_t i=0;i<DYScale.size();++i) std::cout << "DYScale " << i << " " << DYScale.at(i) << std::endl;

        TString fileListName(common::CMSSW_BASE()+"/src/TopAnalysis/Configuration/analysis/diLeptonic/FileLists_"+year+"/HistoFileList_");
        fileListName.Append(systematic+"_"+channel+".txt");
        ifstream FileList(fileListName);
        if (FileList.fail()) std::cerr << "Error reading " << fileListName << std::endl;

        for ( int j=2; j<(int)effHistNames.size();j++ ){

            TString filename;
            while(!FileList.eof()){
                FileList>>filename;
                if(filename==""){continue;}//Skip empty lines
                filename.Prepend(rootFilePrePath);
                dataset.push_back(filename);
                if(filename.Contains("run")){legends.push_back("Data"); }
                else if(filename.Contains("ttbarsignal")){legends.push_back("t#bar{t} signal"); }
                else if(filename.Contains("ttbarbg")){legends.push_back("t#bar{t} other"); }
                else if(filename.Contains("single")){legends.push_back("Single t"); }
                else if(filename.Contains("ww") ||filename.Contains("wz")||filename.Contains("zz")){legends.push_back("Diboson"); }
                else if(filename.Contains("dytautau")){legends.push_back("Z+jets"); }
                else if(filename.Contains("dymumu")||filename.Contains("dyee")){legends.push_back("Z+jets"); }
                else if(filename.Contains("wtolnu")){legends.push_back("W+jets"); }
                else if(filename.Contains("qcd")){legends.push_back("QCD multijet"); }
                else if(filename.Contains("ttbarZ") || filename.Contains("ttgjets")|| filename.Contains("ttbarW")){legends.push_back("t#bar{t}+Z/W");}
            }
            TString DYEntry = "Z+jets";

            TString h_NumName=effHistNames.at(0);
            TString h_DeNumName=effHistNames.at(1);
            h_NumName+=effHistNames.at(j);
            h_DeNumName+=effHistNames.at(j);
            h_NumName.Append("_step8");
            h_DeNumName.Append("_step7");
            std::vector<TH1D> histsNum;
            std::vector<TH1D> histsDeNum;

        for(unsigned int i=0; i<dataset.size(); i++){

            TFile *f = new TFile(dataset.at(i));

            TH1D * histNum = (TH1D*)f->Get(h_NumName);
            TH1D * histDeNum = (TH1D*)f->Get(h_DeNumName);

            if (!histNum) std::cout<< "Error: Can not read histo: "<< h_NumName << " from file: " << dataset.at(i) << "\n";
            if (!histDeNum) std::cout<< "Error: Can not read histo: "<< h_DeNumName << " from file: " << dataset.at(i) << "\n";

            //Rescaling to the data luminosity
            double LumiWeight = plotterConfigHelper->CalcLumiWeight(dataset.at(i));
            plotterConfigHelper->ApplyFlatWeights(histNum, LumiWeight);
            plotterConfigHelper->ApplyFlatWeights(histDeNum, LumiWeight);

            if((legends.at(i) == DYEntry)&&(channelType!=2)){
                histNum->Scale(DYScale[channelType]);
                histDeNum->Scale(DYScale[channelType]);
            }

            for(int ibin = 0; ibin <= histNum->GetNbinsX(); ++ibin) {
            //Check visible bins + underflow
                if(histNum->GetBinContent(ibin) < 0.0) {
                histNum->SetBinContent(ibin, 0.0);// set prediction to 0.0
                histNum->SetBinError(ibin, 1.0);// set error to 1 in order to match expectation of +/-1 event
                }
            }
            if(histNum->GetBinContent(histNum->GetNbinsX()+1) < 0.0) {
            histNum->SetBinContent(histNum->GetNbinsX()+1, 0.0);// set prediction to 0.0
            histNum->SetBinError(histNum->GetNbinsX()+1, 1.0);// set error to 1 in order to match expectation of +/-1 event
            }

            for(int ibin = 0; ibin <= histDeNum->GetNbinsX(); ++ibin) {
            //Check visible bins + underflow
                if(histDeNum->GetBinContent(ibin) < 0.0) {
                histDeNum->SetBinContent(ibin, 0.0);// set prediction to 0.0
                histDeNum->SetBinError(ibin, 1.0);// set error to 1 in order to match expectation of +/-1 event
                }
            }
            if(histDeNum->GetBinContent(histDeNum->GetNbinsX()+1) < 0.0) {
            histDeNum->SetBinContent(histDeNum->GetNbinsX()+1, 0.0);// set prediction to 0.0
            histDeNum->SetBinError(histDeNum->GetNbinsX()+1, 1.0);// set error to 1 in order to match expectation of +/-1 event
            }

            histsNum.push_back(*histNum);
            histsDeNum.push_back(*histDeNum);

            f->Close();

        }

        if (histsNum.size() == 0)std::cerr << "***ERROR! No num-histograms available! "<< std::endl;
        if (histsDeNum.size() == 0)std::cerr << "***ERROR! No denum-histograms available! "<< std::endl;

        TH1 *mcNumHist = NULL;
        TH1 *mcDeNumHist = NULL;
        TH1 *dataNumHist = NULL;
        TH1 *dataDeNumHist = NULL;

        if(legends.at(0) != "Data") std::cout << "***ERROR! legends_.at(0) != Data " << std::endl;
        for(unsigned int i=0; i<histsNum.size(); i++)
        {
            if(legends.at(i) == "Data"){
                if(dataNumHist==NULL){
                    dataNumHist=&histsNum.at(i);
                    dataDeNumHist=&histsDeNum.at(i);
                }
                else{
                    dataNumHist->Add(&histsNum.at(i));
                    dataDeNumHist->Add(&histsDeNum.at(i));
                }
            }
            else{
           //else if(legends.at(i) == "t#bar{t} signal"){
                if(mcNumHist==NULL){
                    mcNumHist=&histsNum.at(i);
                    mcDeNumHist=&histsDeNum.at(i);
                } else {
                    mcNumHist->Add(&histsNum.at(i));
                    mcDeNumHist->Add(&histsDeNum.at(i));
                }

            //}
            }
        }

// In case you want to rebin your histograms
//         mcNumHist -> Rebin(1);
//         mcDeNumHist -> Rebin(1) ;
//         dataNumHist-> Rebin(1) ;
//         dataDeNumHist-> Rebin(1) ;

          EffHist effHist(effHistNames.at(j));

          TString savepath;
          savepath.Append(+"./KinRecoEffSFPlots_"+year+"/"+systematic+"/"+channel+"/"+"KinRecoEff_"+effHistNames.at(j)+".pdf");
          savepath.ReplaceAll("_vs_",4,"",0);

            //TString title(channel);
            TString title("");
            if(effHistNames.at(j)=="_vs_JetMult")title.Append("Kinematic reconstruction efficiency (jet multiplicity); N_{jets}; a.u.");
            if(effHistNames.at(j)=="_vs_LepEta")title.Append("Kinematic reconstruction efficiency (leading lepton); #eta_{l}; a.u.");
            if(effHistNames.at(j)=="_vs_LepEta2")title.Append("Kinematic reconstruction efficiency (trailing lepton); #eta_{l}; a.u.");
            if(effHistNames.at(j)=="_vs_JetEta")title.Append("Kinematic reconstruction efficiency (leading jet); #eta_{jet}; a.u.");
            if(effHistNames.at(j)=="_vs_JetEta2")title.Append("Kinematic reconstruction efficiency (subleading jet); #eta_{jet}; a.u.");
            if(effHistNames.at(j)=="_vs_JetpT")title.Append("Kinematic reconstruction efficiency (leading jet); p_{T}^{jet} [GeV]; a.u.");
            if(effHistNames.at(j)=="_vs_JetpT2")title.Append("Kinematic reconstruction efficiency (subleading jet); p_{T}^{jet} [GeV]; a.u.");
            if(effHistNames.at(j)=="_vs_LeppT")title.Append("Kinematic reconstruction efficiency (leading lepton); p_{T}^{l} [GeV]; a.u.");
            if(effHistNames.at(j)=="_vs_LeppT2")title.Append("Kinematic reconstruction efficiency (trailing lepton); p_{T}^{l} [GeV]; a.u.");
            if(effHistNames.at(j)=="_vs_AntiLeptonpT")title.Append("Kinematic reconstruction efficiency (lepton(+)); p_{T}^{#bar{l}} [GeV]; a.u.");
            if(effHistNames.at(j)=="_vs_LeptonpT")title.Append("Kinematic reconstruction efficiency (lepton(-));p_{T}^{l} [GeV]; a.u.");
            if(effHistNames.at(j)=="_vs_LeptonEta")title.Append("Kinematic reconstruction efficiency (lepton(-)); #eta_{l}; a.u.");
            if(effHistNames.at(j)=="_vs_AntiLeptonEta")title.Append("Kinematic reconstruction efficiency (lepton(+)); #eta_{#bar{l}}; a.u.");
            if(effHistNames.at(j)=="_vs_MET")title.Append("Kinematic reconstruction efficiency (missing transverse energy); E_{T}^{miss} [GeV]; a.u.");
            effHist.setTitle(title);
//             std::cout << title << std::endl;


            if(j==2){
                NdataAfter=dataNumHist->GetBinContent(1);
                NdataBefore=dataDeNumHist->GetBinContent(1);
                NmcAfter=mcNumHist->GetBinContent(1);
                NmcBefore=mcDeNumHist->GetBinContent(1);
                dataEff = NdataAfter / NdataBefore;
                dataEffUnc = sqrt(dataEff * (1-dataEff) / NdataBefore);
                allmcEff =  NmcAfter / NmcBefore;
                allmcEffUnc = sqrt(allmcEff * (1-allmcEff) / NmcBefore);
                sfUnc = sqrt(sqr(dataEffUnc/dataEff) + sqr(allmcEffUnc/allmcEff));

                sprintf(dataEffString, "%.2f%%", 100*dataEff);
                sprintf(allmcEffString, "%.2f%%", 100*allmcEff);
                sprintf(sfString, "%.4f #pm %.4f", dataEff/allmcEff, sfUnc);

                char tempCharVal[100];
                char tempCharErr[100];
                sprintf(tempCharVal,"{Channel::%s, %.4f},",channel.Data(),dataEff/allmcEff);
                sprintf(tempCharErr,"{Channel::%s, %.4f},",channel.Data(),sfUnc);
                printf("{Channel::%s, %s},\n",channel.Data(),sfString);
                kinRecoSF.Append(tempCharVal);
                kinRecoSFerr.Append(tempCharErr);
            }

            if(j>2){
                effHist.addData((TH1D*)dataNumHist,(TH1D*)dataDeNumHist);
                effHist.addMc((TH1D*)mcNumHist,(TH1D*)mcDeNumHist);

                effHist.savePlotEffSF(savepath,dataEffString,allmcEffString,sfString);
            }

      }

    }
  }



  printf("INFO: The next two lines must be copied into KinematicReconstruction.cc :\n");
  kinRecoSF.Append("};");
  kinRecoSFerr.Append("};");
  printf("%s\n%s\n",kinRecoSF.Data(),kinRecoSFerr.Data());


  return 0;
}
