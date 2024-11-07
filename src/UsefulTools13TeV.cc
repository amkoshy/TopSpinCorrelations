#include <map>
#include <utility>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <TH1D.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TPaveText.h>


#include "UsefulTools13TeV.h"
#include "../../common/include/RootFileReader.h"

  ///Please put the variation of each systematics one after each other, starting from the UP variation.
  ///    NOT valid example: MATCH_UP, MASS_DOWN, MASS_UP
  ///    It is important to keep "BTAG_PT" and "BTAG_ETA" in consecutive order for BCJets and LJets, in example: "BTAG_PT_UP", "BTAG_PT_DOWN", "BTAG_ETA_UP", "BTAG_ETA_DOWN"
   const std::vector<const char*> VectorOfValidSystematics 
    {"Nominal",
    "JER_UP", "JER_DOWN",

    //"JES_UP", "JES_DOWN", // cumulative jes uncertainty
    "JESAbsoluteStat_UP", "JESAbsoluteStat_DOWN", // this JES source should be always the first, the treatment of the cumulative uncertainty due to jes sources is associated with this entry    
    "JESAbsoluteScale_UP", "JESAbsoluteScale_DOWN", 
    // // "JESAbsoluteFlavMap_UP", "JESAbsoluteFlavMap_DOWN",
    "JESAbsoluteMPFBias_UP", "JESAbsoluteMPFBias_DOWN", 
    "JESFragmentation_UP", "JESFragmentation_DOWN", 
    "JESSinglePionECAL_UP", "JESSinglePionECAL_DOWN", 
    "JESSinglePionHCAL_UP", "JESSinglePionHCAL_DOWN", 
    "JESFlavorQCD_UP", "JESFlavorQCD_DOWN", 
    "JESTimePtEta_UP", "JESTimePtEta_DOWN",
    "JESRelativeBal_UP", "JESRelativeBal_DOWN", 
    "JESRelativeJEREC1_UP", "JESRelativeJEREC1_DOWN", 
    // // "JESRelativeJEREC2_UP", "JESRelativeJEREC2_DOWN", 
    // // "JESRelativeJERHF_UP", "JESRelativeJERHF_DOWN", 
    "JESRelativePtBB_UP", "JESRelativePtBB_DOWN", 
    "JESRelativePtEC1_UP", "JESRelativePtEC1_DOWN", 
    // // "JESRelativePtEC2_UP", "JESRelativePtEC2_DOWN", 
    // // "JESRelativePtHF_UP", "JESRelativePtHF_DOWN", 
    "JESRelativeFSR_UP", "JESRelativeFSR_DOWN", 
    "JESRelativeStatFSR_UP", "JESRelativeStatFSR_DOWN", 
    "JESRelativeStatEC_UP", "JESRelativeStatEC_DOWN", 
    // // "JESRelativeStatHF_UP", "JESRelativeStatHF_DOWN", 
    "JESPileUpDataMC_UP", "JESPileUpDataMC_DOWN", 
    "JESPileUpPtRef_UP", "JESPileUpPtRef_DOWN", 
    "JESPileUpPtEC1_UP", "JESPileUpPtEC1_DOWN", 
    // // "JESPileUpPtEC2_UP", "JESPileUpPtEC2_DOWN", 
    // // "JESPileUpPtHF_UP", "JESPileUpPtHF_DOWN", 
    "JESPileUpPtBB_UP", "JESPileUpPtBB_DOWN",
    
    "UNCLUSTERED_UP", "UNCLUSTERED_DOWN", "PU_UP", "PU_DOWN",
    "TRIG_UP", "TRIG_DOWN", "TRIG_ETA_UP", "TRIG_ETA_DOWN", "LEPT_UP", "LEPT_DOWN",
    "DY_UP", "DY_DOWN", "BG_UP", "BG_DOWN", "KIN_UP", "KIN_DOWN",
    "BTAG_UP", "BTAG_DOWN", "BTAG_LJET_UP", "BTAG_LJET_DOWN",
    "BTAG_PT_UP", "BTAG_PT_DOWN", "BTAG_ETA_UP", "BTAG_ETA_DOWN",
    "BTAG_LJET_PT_UP", "BTAG_LJET_PT_DOWN", "BTAG_LJET_ETA_UP", "BTAG_LJET_ETA_DOWN",
    //"BTAG_BEFF_UP", "BTAG_BEFF_DOWN", "BTAG_CEFF_UP", "BTAG_CEFF_DOWN", "BTAG_LEFF_UP", "BTAG_LEFF_DOWN",
    "MESCALE_UP", "MESCALE_DOWN", "MEFACSCALE_UP", "MEFACSCALE_DOWN", "MERENSCALE_UP", "MERENSCALE_DOWN",
    "PSISRSCALE_UP", "PSISRSCALE_DOWN", "PSFSRSCALE_UP", "PSFSRSCALE_DOWN",
    "BFRAG_UP", "BFRAG_DOWN", "BFRAG_PETERSON", //"BFRAG_CENTRAL",
    "BSEMILEP_UP", "BSEMILEP_DOWN", 
    "ERDON", "ERDONRETUNE", "GLUONMOVETUNE",
    "UETUNE_UP", "UETUNE_DOWN",
    //"SCALE_UP", "SCALE_DOWN",
    "MATCH_UP", "MATCH_DOWN",
    "MASS_UP", "MASS_DOWN",
    //"MCATNLO", "POWHEG", //"POWHEGHERWIG", //"POWHEGV2", //"POWHEGV2HERWIG", //"AMCATNLOFXFX", //"MADGRAPHMLM",// "PERUGIA11", // "SPINCORR",
    "PDF_ALPHAS_UP", "PDF_ALPHAS_DOWN",
    "PDF_UP", "PDF_DOWN",

    "all"};

//
// N.B. This is only for top mass extraction analysis. Systematics for differential cross-section analysis should be added to VectorOfValidSystematics
//
const std::vector<const char*> VectorOfSystematicsForTopMassAnalysis
{   "ELE_UP","ELE_DOWN",
    "MUON_UP","MUON_DOWN",
    "TOP_PT", "POWHEG",
     "TOP_PT_UP","TOP_PT_DOWN",

    "BTAGDISCR_BPURITY_UP", "BTAGDISCR_BPURITY_DOWN", 
    "BTAGDISCR_LPURITY_UP", "BTAGDISCR_LPURITY_DOWN", 
    "BTAGDISCR_BSTAT1_UP", "BTAGDISCR_BSTAT1_DOWN",   
    "BTAGDISCR_BSTAT2_UP", "BTAGDISCR_BSTAT2_DOWN",   
    "BTAGDISCR_LSTAT1_UP", "BTAGDISCR_LSTAT1_DOWN",   
    "BTAGDISCR_LSTAT2_UP", "BTAGDISCR_LSTAT2_DOWN",   
    "BTAGDISCR_CERR1_UP", "BTAGDISCR_CERR1_DOWN",     
    "BTAGDISCR_CERR2_UP", "BTAGDISCR_CERR2_DOWN",     

    "JESAbsoluteStat_UP", "JESAbsoluteStat_DOWN",     
    "JESAbsoluteScale_UP", "JESAbsoluteScale_DOWN", 
    "JESAbsoluteFlavMap_UP", "JESAbsoluteFlavMap_DOWN", 
    "JESAbsoluteMPFBias_UP", "JESAbsoluteMPFBias_DOWN", 
    "JESFragmentation_UP", "JESFragmentation_DOWN", 
    "JESSinglePionECAL_UP", "JESSinglePionECAL_DOWN", 
    "JESSinglePionHCAL_UP", "JESSinglePionHCAL_DOWN", 
    "JESFlavorQCD_UP", "JESFlavorQCD_DOWN", 
    "JESTimePtEta_UP", "JESTimePtEta_DOWN", 
    "JESRelativeJEREC1_UP", "JESRelativeJEREC1_DOWN", 
    "JESRelativeJEREC2_UP", "JESRelativeJEREC2_DOWN", 
    "JESRelativeJERHF_UP", "JESRelativeJERHF_DOWN", 
    "JESRelativePtBB_UP", "JESRelativePtBB_DOWN", 
    "JESRelativePtEC1_UP", "JESRelativePtEC1_DOWN", 
    "JESRelativePtEC2_UP", "JESRelativePtEC2_DOWN", 
    "JESRelativePtHF_UP", "JESRelativePtHF_DOWN", 
    "JESRelativeFSR_UP", "JESRelativeFSR_DOWN", 
    "JESRelativeStatFSR_UP", "JESRelativeStatFSR_DOWN", 
    "JESRelativeStatEC_UP", "JESRelativeStatEC_DOWN", 
    "JESRelativeStatHF_UP", "JESRelativeStatHF_DOWN", 
    "JESPileUpDataMC_UP", "JESPileUpDataMC_DOWN", 
    "JESPileUpPtRef_UP", "JESPileUpPtRef_DOWN", 
    "JESPileUpPtEC1_UP", "JESPileUpPtEC1_DOWN", 
    "JESPileUpPtEC2_UP", "JESPileUpPtEC2_DOWN", 
    "JESPileUpPtHF_UP", "JESPileUpPtHF_DOWN", 
    "JESPileUpPtBB_UP", "JESPileUpPtBB_DOWN",

     "PDF_1_UP", "PDF_1_DOWN",
     "PDF_2_UP", "PDF_2_DOWN",
     "PDF_3_UP", "PDF_3_DOWN",
     "PDF_4_UP", "PDF_4_DOWN",
     "PDF_5_UP", "PDF_5_DOWN",
     "PDF_6_UP", "PDF_6_DOWN",
     "PDF_7_UP", "PDF_7_DOWN",
     "PDF_8_UP", "PDF_8_DOWN",
     "PDF_9_UP", "PDF_9_DOWN",
     "PDF_10_UP", "PDF_10_DOWN",
     "PDF_11_UP", "PDF_11_DOWN",
     "PDF_12_UP", "PDF_12_DOWN",
     "PDF_13_UP", "PDF_13_DOWN",
     "PDF_14_UP", "PDF_14_DOWN",
     "PDF_15_UP", "PDF_15_DOWN",
     "PDF_16_UP", "PDF_16_DOWN",
     "PDF_17_UP", "PDF_17_DOWN",
     "PDF_18_UP", "PDF_18_DOWN",
     "PDF_19_UP", "PDF_19_DOWN",
     "PDF_20_UP", "PDF_20_DOWN",
     "PDF_21_UP", "PDF_21_DOWN",
     "PDF_22_UP", "PDF_22_DOWN",
     "PDF_23_UP", "PDF_23_DOWN",
     "PDF_24_UP", "PDF_24_DOWN",
     "PDF_25_UP", "PDF_25_DOWN",
     "PDF_26_UP", "PDF_26_DOWN",
     "PDF_27_UP", "PDF_27_DOWN",
     "PDF_28_UP", "PDF_28_DOWN",
     "PDF_29_UP", "PDF_29_DOWN"

};

UsefulTools13TeV::UsefulTools13TeV( RootFileReader* rootFileReader,bool isClosureTest,bool isDYScale):

energy(13.0), // centre-of-mass energy in TeV
lumi(35922.), // data luminosity in pb-1:  Full 2016: 35687. ;  2015 C,D (non-EA dataset): 2225.51

topxsec(831.76), //https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO

fileReader(rootFileReader),
is2015(false), // make sure also to update lumi-values and its error accordingly
is2016(true), //  make sure also to update lumi-values and its error accordingly
isSingleLeptonDatasetIncluded(true),// supported only for 2016 data
doClosureTest(isClosureTest),
doDYScale(isDYScale),
isPoissonSmearedPseudoData(false) // in case of real data available use "false"
{
    for (auto s: VectorOfValidSystematics) ListOfSyst.insert(s);
}

void UsefulTools13TeV::fillSetListOfSystematics(std::set<TString>& set)
{

    for (auto s: VectorOfValidSystematics) set.insert(s);
}

void UsefulTools13TeV::fillVectorOfValidSystematics(std::vector<const char*>& vect, const bool useTopMassSetup)
{
    vect.clear();
    
    for (auto s: VectorOfValidSystematics) vect.push_back(s);
    
    if (useTopMassSetup){
        for (auto s: VectorOfSystematicsForTopMassAnalysis){
            if (std::find(begin(VectorOfValidSystematics), end(VectorOfValidSystematics), s) == end(VectorOfValidSystematics))
                vect.push_back(s);
        }
    }
}

double UsefulTools13TeV::SampleXSection(const TString& filename){
    
    //13 TeV MC cross sections taken from:
    //  https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
    //  Ones without "updated" mark still from 8TeV 
    
    if(filename.Contains("run"))              {if(isPoissonSmearedPseudoData) return topxsec;
                                               else return 1.;}
    else if(filename.Contains("fromDilepton"))    {return topxsec * 0.10706;}// updated PDG 2017
    else if(filename.Contains("FullLept"))    {return topxsec * 0.1049;}
    else if(filename.Contains("SemiLept"))    {return topxsec * 0.4380;}
    else if(filename.Contains("Hadronic"))    {return topxsec * 0.4570;}
    else if(filename.Contains("Perugia11") &&
        filename.Contains("signal"))          {return topxsec * 0.1049;}
    else if(filename.Contains("ttbar") && !filename.Contains("ttbarW") && 
        !filename.Contains("ttbarZ"))         {return topxsec;}// updated
    else if(filename.Contains("single"))      {return 35.85;}// updated
    else if(filename.Contains("ww"))          {return 118.7;}// updated
    else if(filename.Contains("wz"))          {return 47.13;}// updated 04.09.17: NLO from MCFM
    else if(filename.Contains("zz"))          {return 16.523;}// updated 04.09.17: NLO from MCFM
    else if(filename.Contains("1050"))        {return 3.*7545.03;}// updated
    else if(filename.Contains("50inf"))       {return 3.*1921.8;}// updated 04.09.17: NNLO FEWZ 3.1.b2
    else if(filename.Contains("wtolnu"))      {return 3.*20508.9;}// updated
    //else if(filename.Contains("qcdmu15"))     {return 866600000*0.00044;}// update needed
    //else if(filename.Contains("qcd2em30inf"))     {return 54120000*0.002;}// update needed
    //else if(filename.Contains("qcd2em40inf"))     {return 866600000*0.00044;}
    //else if(filename.Contains("qcd2em3040"))     {return 866600000*0.00044;}
    //else if(filename.Contains("qcdmu2030"))   {return 2.870E8*6.500E-3;}// update needed
    //else if(filename.Contains("qcdmu3050"))   {return 6.609E7*12.20E-3;}// update needed
    //else if(filename.Contains("qcdmu5080"))   {return 8.802E6*21.80E-3;}// update needed
    //else if(filename.Contains("qcdmu80120"))  {return 1.024E6*39.50E-3;}// update needed
    //else if(filename.Contains("qcdmu120170")) {return 1.578E5*47.30E-3;}// update needed
    //else if(filename.Contains("qcdem2030"))   {return 2.886E8*10.10E-3;}// update needed
    //else if(filename.Contains("qcdem3080"))   {return 7.433E7*62.10E-3;}// update needed
    //else if(filename.Contains("qcdem80170"))  {return 1.191E6*153.9E-3;}// update needed
    //else if(filename.Contains("qcdbcem2030")) {return 2.886E8*5.800E-4;}// update needed
    //else if(filename.Contains("qcdbcem3080")) {return 7.424E7*2.250E-3;}// update needed
    //else if(filename.Contains("qcdbcem80170")){return 1.191E6*10.90E-3;}// update needed
    //else if(filename.Contains("ttbarW"))      {return 1.152;}// update needed
    //else if(filename.Contains("ttbarZ"))      {return 2.232;}// update needed
    //else if(filename.Contains("ttgjets"))     {return 1.8;}// update needed
    else if(filename.Contains("ttbarWjetstolnu"))      {return 0.2043;}// updated: nlo
    else if(filename.Contains("ttbarWjetstoqq"))      {return 0.4062;}// updated: nlo
    else if(filename.Contains("ttbarZtollnunu"))      {return 0.2529;}// updated: nlo
    else if(filename.Contains("ttbarZtoqq"))     {return 0.5297;}// updated: nlo 
    
    return -1;
}


void UsefulTools13TeV::fillLegendColorDataset(const TString& fileListName, std::vector<TString>& legends, std::vector<int>& colors, std::vector<TString>& dataset){

        std::cout << "reading " << fileListName << std::endl;
    
        std::ifstream FileList(fileListName);
        if (FileList.fail()){
            std::cerr << "Error reading " << fileListName << std::endl;
            exit(1);
        }
        
        TString filename;
        
        dataset.clear();
        legends.clear();
        colors.clear();
    
        while(!FileList.eof()){ 
        FileList>>filename;
        if(filename==""){continue;}//Skip empty lines
        dataset.push_back(filename);
        if(filename.Contains("run")){legends.push_back("Data"); colors.push_back(kBlack);}
        else if(filename.Contains("ttbarsignal")){legends.push_back("t#bar{t} signal"); colors.push_back(kRed+1);}
        else if(filename.Contains("ttbarbg")){legends.push_back("t#bar{t} other"); colors.push_back(kRed-7);}
        else if(filename.Contains("single")){legends.push_back("Single t"); colors.push_back(kMagenta);}
        else if(filename.Contains("ww") ||filename.Contains("wz")||filename.Contains("zz")){legends.push_back("Diboson"); colors.push_back(10);}
        //else if(filename.Contains("dytautau")){legends.push_back("Z / #gamma* #rightarrow #tau#tau"); colors.push_back(kAzure+8);}
        else if(filename.Contains("dymumu")||filename.Contains("dyee")||filename.Contains("dytautau")){legends.push_back("Z+jets"); colors.push_back(kAzure-2);}
        else if(filename.Contains("wtolnu")){legends.push_back("W+jets"); colors.push_back(kGreen-3);}
        else if(filename.Contains("qcd")){legends.push_back("QCD multijet"); colors.push_back(kYellow);}
        else if(filename.Contains("ttbarZ") || filename.Contains("ttbarW")){legends.push_back("t#bar{t}+Z/W"); colors.push_back(kOrange-2);}
    }
    FileList.close();
}


double UsefulTools13TeV::CalcLumiWeight(const TString& WhichSample){
    if (WhichSample.Contains("run") && !isPoissonSmearedPseudoData) return 1;
    double lumiWeight=0;
    if(WhichSample!=""){
        double XSection = SampleXSection(WhichSample);
        if(XSection <= 0.){
            std::cout<<"Sample XSection is <0. Can't calculate luminosity weight!! returning"<<std::endl;
            return 0;
        }
        //From 'filename' get the number of weighted (MC weights) event processed.
        const TH1 *h_NrOfEvts = fileReader->Get<TH1>(WhichSample, "weightedEvents");
        double NrOfEvts = h_NrOfEvts->GetBinContent(1);
        if (!WhichSample.Contains("run")) lumiWeight = lumi*XSection/NrOfEvts;
        if (WhichSample.Contains("run") && isPoissonSmearedPseudoData) lumiWeight = 2*lumi*XSection/NrOfEvts; // if ttbarsignlaplustau is merged with ttbarbgviatau --> 2
        //std::cout<<WhichSample<<" "<<lumiWeight<<" xsec:"<<XSection<<std::endl;
    }

    if (lumiWeight == 0) {
    std::cout << WhichSample << " has lumi weight 0\n";
    }
    return lumiWeight;
}

std::vector<TString> UsefulTools13TeV::InputNominalFileList(TString mode)
{
    //Hard code the input nominal file list
    
    std::vector<TString> FileVector;
    FileVector.clear();
    
    if( mode.CompareTo("combined") && mode.CompareTo("ee") && mode.CompareTo("emu") && mode.CompareTo("mumu")){
        std::cout<<"The decay channel you provided is not supported."<<std::endl;
        std::cout<<"Please use: ee, emu, mumu, combined"<<std::endl;
        return FileVector;
    }
    
    if(!mode.CompareTo("combined")){
        std::vector<TString> eemode   = UsefulTools13TeV::InputNominalFileList(TString("ee"));
        std::vector<TString> emumode  = UsefulTools13TeV::InputNominalFileList(TString("emu"));
        std::vector<TString> mumumode = UsefulTools13TeV::InputNominalFileList(TString("mumu"));
        FileVector.insert(FileVector.end(), eemode.begin(), eemode.end());
        FileVector.insert(FileVector.end(), emumode.begin(), emumode.end());
        FileVector.insert(FileVector.end(), mumumode.begin(), mumumode.end());
        return FileVector;
    }

    //data is only stored in the Nominal directory
    TString nominalPath = TString("selectionRoot/Nominal/") + mode + "/" + mode;
    if (!doClosureTest || !isPoissonSmearedPseudoData) {      
        
        if (is2015) {
            FileVector.push_back(nominalPath + "_run2015D.root");
            FileVector.push_back(nominalPath + "_run2015C.root");
        }

        if (is2016) {
            FileVector.push_back(nominalPath + "_run2016B.root");
            FileVector.push_back(nominalPath + "_run2016C.root");
            FileVector.push_back(nominalPath + "_run2016D.root");
            FileVector.push_back(nominalPath + "_run2016E.root");
            //FileVector.push_back(nominalPath + "_run2016F.root");
            FileVector.push_back(nominalPath + "_run2016F1.root");
            FileVector.push_back(nominalPath + "_run2016F2.root");
            FileVector.push_back(nominalPath + "_run2016G.root");
            FileVector.push_back(nominalPath + "_run2016H.root");
            if (isSingleLeptonDatasetIncluded) {
                FileVector.push_back(nominalPath + "_se_run2016B.root");
                FileVector.push_back(nominalPath + "_se_run2016C.root");
                FileVector.push_back(nominalPath + "_se_run2016D.root");
                FileVector.push_back(nominalPath + "_se_run2016E.root");
                //FileVector.push_back(nominalPath + "_se_run2016F.root");
                FileVector.push_back(nominalPath + "_se_run2016F1.root");
                FileVector.push_back(nominalPath + "_se_run2016F2.root");
                FileVector.push_back(nominalPath + "_se_run2016G.root");
                FileVector.push_back(nominalPath + "_se_run2016H.root");

                FileVector.push_back(nominalPath + "_smu_run2016B.root");
                FileVector.push_back(nominalPath + "_smu_run2016C.root");
                FileVector.push_back(nominalPath + "_smu_run2016D.root");
                FileVector.push_back(nominalPath + "_smu_run2016E.root");
                //FileVector.push_back(nominalPath + "_smu_run2016F.root");
                FileVector.push_back(nominalPath + "_smu_run2016F1.root");
                FileVector.push_back(nominalPath + "_smu_run2016F2.root");
                FileVector.push_back(nominalPath + "_smu_run2016G.root");
                FileVector.push_back(nominalPath + "_smu_run2016H.root");
            }
        }     
    } else {
        FileVector.push_back(nominalPath + "_ttbarsignalplustau_fakerun_nominal.root");
    }
      
    /*FileVector.push_back(nominalPath + "_dyee1050.root");
    FileVector.push_back(nominalPath + "_dyee50inf.root");
    FileVector.push_back(nominalPath + "_dymumu1050.root");
    FileVector.push_back(nominalPath + "_dymumu50inf.root");
    FileVector.push_back(nominalPath + "_dytautau1050.root");
    FileVector.push_back(nominalPath + "_dytautau50inf.root");*/
    FileVector.push_back(nominalPath + "_dyee1050_amcatnlofxfx.root");
    FileVector.push_back(nominalPath + "_dyee50inf_amcatnlofxfx.root");
    FileVector.push_back(nominalPath + "_dymumu1050_amcatnlofxfx.root");
    FileVector.push_back(nominalPath + "_dymumu50inf_amcatnlofxfx.root");
    FileVector.push_back(nominalPath + "_dytautau1050_amcatnlofxfx.root");
    FileVector.push_back(nominalPath + "_dytautau50inf_amcatnlofxfx.root");
    FileVector.push_back(nominalPath + "_singleantitop_tw.root");
    FileVector.push_back(nominalPath + "_singletop_tw.root");
    FileVector.push_back(nominalPath + "_wtolnu.root");
    FileVector.push_back(nominalPath + "_wwtoall.root");
    FileVector.push_back(nominalPath + "_wztoall.root");
    FileVector.push_back(nominalPath + "_zztoall.root");
    FileVector.push_back(nominalPath + "_ttbarWjetstolnu.root");
    FileVector.push_back(nominalPath + "_ttbarWjetstoqq.root");
    FileVector.push_back(nominalPath + "_ttbarZtollnunu.root");
    FileVector.push_back(nominalPath + "_ttbarZtoqq.root");
    FileVector.push_back(nominalPath + "_ttbarsignalplustau.root");
    FileVector.push_back(nominalPath + "_ttbarbgviatau.root");
    FileVector.push_back(nominalPath + "_ttbarbg.root");
    
    return FileVector;
    
}

void UsefulTools13TeV::ApplyFlatWeights(TH1* varhists, const double weight){

    if(weight == 0) {std::cout<<"Warning: the weight your applying is 0. This will remove your distribution."<<std::endl;}
    //if(weight >=1e3){std::cout<<"Warning: the weight your applying is >= 1e3. This will enlarge too much your distribution."<<std::endl;}
    varhists->Scale(weight);
}

void UsefulTools13TeV::ApplyFlatWeights(TH2* varhists, const double weight){

    if(weight == 0) {std::cout<<"Warning: the weight your applying is 0. This will remove your distribution."<<std::endl;}
    //if(weight >=1e3){std::cout<<"Warning: the weight your applying is >= 1e3. This will enlarge too much your distribution."<<std::endl;}
    varhists->Scale(weight);
}

void UsefulTools13TeV::DYScaleFactor(TString SpecialComment,std::vector<double>& DYScale,TString name){

    DYScale = {1,1,1,1};

    if(!doDYScale || doClosureTest || isPoissonSmearedPseudoData) return; //need to make a switch for control plots that don't want DYScale

    TString nameAppendix = "";
    if ( !SpecialComment.BeginsWith("_post") &&  SpecialComment != "Standard" ){
        std::cout<<"\n\n*******************************************************************"<<std::endl;
        std::cout<<"ERROR: When calculating the DY Scale factor you must specify in which step you want to calculate the DY SF:"<<std::endl;
        std::cout<<" '_postZcut', '_post2jets', '_postMET', '_post1btag', '_postKinReco' or 'Standard' = _postKinReco"<<std::endl;
        std::cout<<"*******************************************************************\n\n"<<std::endl;
        exit(444);
    }
    if (SpecialComment.BeginsWith("_post")){
        if(SpecialComment.EqualTo("_postZcut"))nameAppendix = "_step4";
        if(SpecialComment.EqualTo("_post2jets"))nameAppendix = "_step5";
        if(SpecialComment.EqualTo("_postMET"))nameAppendix = "_step6";
        if(SpecialComment.EqualTo("_post1btag"))nameAppendix = "_step7";
        if(SpecialComment.EqualTo("_postKinReco"))nameAppendix = "_step8";
        
    } else if ( SpecialComment == "Standard") {
        nameAppendix = "_step5";
    }

    std::cout<<"\n\nBegin DYSCALE FACTOR calculation at selection step "<<nameAppendix<<std::endl;
    
    std::vector<TString> Vec_Files = InputNominalFileList("combined");//Read the hardcoded list of files
    if(Vec_Files.size()<1) {std::cout<<"WARNING(in DYScaleFactor)!!! No datasets available to calculate DY SF. EXITING!!"<<std::endl; return;}
    
    double NoutEEDYMC=0, NinEEDYMC=0, NoutMuMuDYMC=0, NinMuMuDYMC=0;//Number of events in/out of z-veto region for the DY MC
    double NoutEMuDYMC=0;// same as above, but for EMU channel, needed for correct combination of SFs for the "combined" channel
    double NinEE=0, NinMuMu=0, NinEMu=0;//Number of events in z-veto region for data
    double NinEEloose=0, NinMuMuloose=0;//Number of data events in Z-Veto region with MET cut
    double NinEEMC=0, NinMuMuMC=0;//All other MC events

    for(size_t i=0; i < Vec_Files.size(); i++){
        double LumiWeight = CalcLumiWeight(Vec_Files.at(i));
        double allWeights=LumiWeight;//calculate here all the flat-weights we apply: Lumi*others*...
        if(Vec_Files.at(i).Contains("ee_") || Vec_Files.at(i).Contains("mumu_")){
            if(Vec_Files.at(i).Contains("run")){
                TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Zh1").Append(nameAppendix)).Append("zWindow"));
                TH1D *htemp1 = fileReader->GetClone<TH1D>(Vec_Files.at(i), "dyScaling_Looseh1");
                ApplyFlatWeights(htemp, allWeights);
                ApplyFlatWeights(htemp1, allWeights);
                if(Vec_Files.at(i).Contains("ee_")){
                    NinEE+=htemp->Integral();
                    NinEEloose+=htemp1->Integral();
                }
                if(Vec_Files.at(i).Contains("mumu_")){
                    NinMuMu+=htemp->Integral();
                    NinMuMuloose+=htemp1->Integral();
                }
                delete htemp; delete htemp1;
            }
            else if(Vec_Files.at(i).Contains("dy")){
                if(Vec_Files.at(i).Contains("50inf")){
                    TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Zh1").Append(nameAppendix)).Append("zWindow"));
                    TH1D *htemp1 = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("dyScaling_TTh1").Append(nameAppendix));
                    ApplyFlatWeights(htemp, LumiWeight);
                    ApplyFlatWeights(htemp1, LumiWeight);
                    if(Vec_Files.at(i).Contains("ee_")){
                        NinEEDYMC+=htemp->Integral();
                        NoutEEDYMC+=htemp1->Integral();
                    }
                    if(Vec_Files.at(i).Contains("mumu_")){
                        NinMuMuDYMC+=htemp->Integral();
                        NoutMuMuDYMC+=htemp1->Integral();
                    }
                    delete htemp; delete htemp1;
                }
                else{
                    TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("dyScaling_TTh1").Append(nameAppendix));
                    ApplyFlatWeights(htemp, LumiWeight);
                    if(Vec_Files.at(i).Contains("ee_")){   NoutEEDYMC+=htemp->Integral();}
                    if(Vec_Files.at(i).Contains("mumu_")){ NoutMuMuDYMC+=htemp->Integral();}
                    delete htemp;
                }
            }
            else{
                TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Zh1").Append(nameAppendix)).Append("zWindow"));
                ApplyFlatWeights(htemp, LumiWeight);
                if(Vec_Files.at(i).Contains("ee_")){   NinEEMC+=htemp->Integral();   }
                if(Vec_Files.at(i).Contains("mumu_")){ NinMuMuMC+=htemp->Integral(); }
                delete htemp;
            }
        }
        
        if(Vec_Files.at(i).Contains("emu_")){
            if(Vec_Files.at(i).Contains("run")) {
                TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString(TString("dyScaling_Zh1").Append(nameAppendix)).Append("zWindow"));
                ApplyFlatWeights(htemp, LumiWeight);
                NinEMu+=htemp->Integral();
                delete htemp;
            }
            else if(Vec_Files.at(i).Contains("dy")){
                if(Vec_Files.at(i).Contains("50inf")){
                    TH1D *htemp1 = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("dyScaling_TTh1").Append(nameAppendix));
                    ApplyFlatWeights(htemp1, LumiWeight);
                    NoutEMuDYMC+=htemp1->Integral();
                    delete htemp1;
                }
                else{
                    TH1D *htemp = fileReader->GetClone<TH1D>(Vec_Files.at(i), TString("dyScaling_TTh1").Append(nameAppendix));
                    ApplyFlatWeights(htemp, LumiWeight);
                    NoutEMuDYMC+=htemp->Integral();
                    delete htemp;
                }
            }
        }

    }
    
    
    double kee = sqrt(NinEEloose/NinMuMuloose);
//     printf("kee = sqrt(%.2f/%.2f) = %.2f\n", NinEEloose, NinMuMuloose, kee);
    
    double kmumu = sqrt(NinMuMuloose/NinEEloose);
//     printf("kmumu = sqrt(%.2f/%.2f) = %.2f\n", NinMuMuloose, NinEEloose, kmumu);
    
    double RoutinEE = NoutEEDYMC/NinEEDYMC;
//     printf("RoutinEE = %.2f/%.2f = %.2f\n", NoutEEDYMC, NinEEDYMC, RoutinEE);
    
    double RoutinMuMu = NoutMuMuDYMC/NinMuMuDYMC;
//     printf("RoutinMuMu = %.2f/%.2f = %.2f\n", NoutMuMuDYMC, NinMuMuDYMC, RoutinMuMu);
    
    double NoutMCEE = RoutinEE*(NinEE - 0.5*NinEMu*kee);
    double NoutMCMuMu = RoutinMuMu*(NinMuMu - 0.5*NinEMu*kmumu);

    double DYSFEE = NoutMCEE/NoutEEDYMC;
    double DYSFMuMu = NoutMCMuMu/NoutMuMuDYMC;
    double DYSFEMu = std::sqrt(DYSFEE * DYSFMuMu);

    std::cout << std::endl;
    std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;
    std::cout << "Calculation of DY Scale Factors for '" << name << "'  at selection step "<<nameAppendix << std::endl;

    std::cout<<"DYSFEE:                 "<<DYSFEE<<std::endl;
    std::cout<<"DYSFMuMu:               "<<DYSFMuMu<<std::endl;
    std::cout<<"DYSFEMu:                "<<DYSFEMu<<std::endl;

    std::cout<<"NinEEloose:             "<<NinEEloose<<std::endl;
    std::cout<<"NinMMloose:             "<<NinMuMuloose<<std::endl;

    std::cout<<"kee:                    "<<kee<<" +- "<<kee*0.5*TMath::Sqrt(1./NinMuMuloose + 1./NinEEloose)<<std::endl;
    std::cout<<"kmumu:                  "<<kmumu<<" +- "<<kmumu*0.5*TMath::Sqrt(1./NinMuMuloose + 1./NinEEloose)<<std::endl;

    std::cout<<"Rout/Rin ee:            "<<RoutinEE<<std::endl;
    std::cout<<"Rout/Rin Mumu:          "<<RoutinMuMu<<std::endl;

    std::cout<<"Est. From Data(ee):     "<<NoutMCEE<<std::endl;
    std::cout<<"Est. From Data(mumu):   "<<NoutMCMuMu<<std::endl;

    std::cout<<"Est. From MC(ee):       "<<NoutEEDYMC<<std::endl;
    std::cout<<"Est. From MC(mumu):     "<<NoutMuMuDYMC<<std::endl;
    std::cout<<"Est. From MC(emu):      "<<NoutEMuDYMC<<std::endl;

    std::cout << "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << std::endl;
    std::cout << std::endl;

    if (DYSFEE < 0.2 || DYSFMuMu < 0.2) {
        std::cout << "The DY SF is too low (below 0.2). Something is probably wrong.\n";
        std::exit(1);
    }
    
    DYScale.at(0)=DYSFEE;
    DYScale.at(1)=DYSFMuMu;
    DYScale.at(2)=DYSFEMu;
    DYScale.at(3)=(DYSFEE*NoutEEDYMC+DYSFMuMu*NoutMuMuDYMC+DYSFEMu*NoutEMuDYMC)/(NoutEEDYMC+NoutMuMuDYMC+NoutEMuDYMC);

    std::cout<<"End DYSCALE FACTOR calculation\n"<<std::endl;

}

// Draw official labels (CMS Preliminary, luminosity and CM energy) above plot
void UsefulTools13TeV::DrawCMSLabels(double lumi,int cmsprelim, double energy, double textSize) {

    const char *text;
    if(cmsprelim ==2 ) {//Private work for PhDs students
        text = "Private Work, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV";
    } else if (cmsprelim==1) {//CMS preliminary label
        text = "CMS Preliminary, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV";
    } else {//CMS label
        text = "CMS, %2.1f fb^{-1} at #sqrt{s} = %2.f TeV";
    }
    
    TPaveText *label = new TPaveText();
    label->SetX1NDC(gStyle->GetPadLeftMargin());
    //label->SetY1NDC(1.0-gStyle->GetPadTopMargin());
    label->SetX2NDC(1.0-gStyle->GetPadRightMargin());
    //label->SetY2NDC(1.0);
    label->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.02 );
    label->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() + 0.03 );
    
    label->SetTextFont(42);
    label->AddText(Form(text, lumi/1000, energy));
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    if (textSize!=0) label->SetTextSize(textSize);
    label->SetTextAlign(32);
    label->Draw("same");
}

// Draw label for Decay Channel in upper left corner of plot
void UsefulTools13TeV::DrawDecayChLabel(TString decaychannel, double textSize) {

    TPaveText *decch = new TPaveText();

    decch->AddText(decaychannel);

    decch->SetX1NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength()        );
    decch->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 );
    decch->SetX2NDC(      gStyle->GetPadLeftMargin() + gStyle->GetTickLength() + 0.15 );
    decch->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength()        );

    decch->SetFillStyle(0);
    decch->SetBorderSize(0);
    if (textSize!=0) decch->SetTextSize(textSize);
    decch->SetTextAlign(12);
    decch->Draw("same");
}

void UsefulTools13TeV::setStyle(TH1 *hist, TString Axis)
{
//     hist->SetLineColor(2);
//     hist->SetLineWidth(1);
//     hist->GetXaxis()->SetLabelFont(42);
//     hist->GetYaxis()->SetLabelFont(42);
//     hist->GetXaxis()->SetTitleFont(42);
//     hist->GetYaxis()->SetTitleFont(42);
//     hist->GetXaxis()->SetTitleSize(0.05);
//     hist->GetYaxis()->SetTitleSize(0.05);
//     hist->GetXaxis()->SetTitleOffset(1.08);
//     hist->GetYaxis()->SetTitleOffset(1.7);
//     hist->GetXaxis()->SetLabelOffset(0.007);
//     hist->GetYaxis()->SetLabelOffset(0.007);
// 
// 
//         hist->SetFillColor(0);
//     hist->SetMarkerStyle(20);
        if ((Axis.Contains("p_{T}", TString::kIgnoreCase) || Axis.Contains("m(t#bar{t})", TString::kIgnoreCase)) && 
            (!name.Contains("1st") && !name.Contains("Rapidity", TString::kIgnoreCase) && !name.Contains("Eta", TString::kIgnoreCase) && !name.Contains("Phi", TString::kIgnoreCase) && !name.Contains("JetMult"))) {
            hist->GetXaxis()->SetTitle(Axis+" #left[GeV#right]");
            hist->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d"+Axis+"}"+" #left[GeV^{-1}#right]"); 
            
        } else {
            hist->GetXaxis()->SetTitle(Axis);
            hist->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d"+Axis+"}");
        }

}

