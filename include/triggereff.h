#include <TH1F.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TChain.h>
#include <TCanvas.h>
#include <iostream>
#include "TLorentzVector.h"
#include "TMath.h"
#include <vector>
#include "../../common/include/classesFwd.h"
#include "../../common/include/classes.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../../../ZTopUtils/interface/BitsManager.h"
#include "../../../../ZTopUtils/interface/PUReweighter.h"
#include "../../common/include/ScaleFactors.h"
#include "../../common/include/sampleHelpers.h"

using namespace std;

namespace ztop{
    class PUReweighter;
}
/// Pointer to the pileup reweighter instance
ztop::PUReweighter* const puReweighter_(new ztop::PUReweighter());


string era_;
string MCorData_;
string sample_;

bool isMC;

BitsManager<string> HLTBitsManager_; 
TLorentzVector LeptonA;
TLorentzVector LeptonB;

TString filename;
string suffix = "";
string DATA_PATH_COMMON("../common/data");

string SF_electron_ID_File; 
string SF_electron_ID_Histo = "EGamma_SF2D";
string SF_electron_ID_Format = "pt_vs_etaSc";

string SF_electron_Reco_File;
string SF_electron_Reco_Histo  = "EGamma_SF2D";
string SF_electron_Reco_Format = "pt_vs_etaSc";

string SF_electron_Reco_File_LowPt;
string SF_electron_Reco_Histo_LowPt  = "EGamma_SF2D";
string SF_electron_Reco_Format_LowPt = "pt_vs_etaSc";

string SF_muon_ID_File;
string SF_muon_ID_Histo;
string SF_muon_ID_Format;
string SF_muon_ID_Histo_Stat;
string SF_muon_ID_Histo_Sys;

string SF_muon_ISO_File;
string SF_muon_ISO_Histo;
string SF_muon_ISO_Format;
string SF_muon_ISO_Histo_Stat;
string SF_muon_ISO_Histo_Sys;

string pileup_MC_File = "@{INPUT_FILE}";
string pileup_MC_Histogram = "PileupMCTemplateMaker/MC_TrueNIntBX0"; 

string pileup_DATA_File;
string pileup_DATA_Histogram = "pileup";

vector<string> HLTPathsOR_ee;

vector<string> HLTPathsOR_emu;

vector<string> HLTPathsOR_mumu;

vector<string> HLTPathsOR_met;


TH2* h2_electron_ID_histo_(0);
TH2* h2_electron_Reco_histo_(0);
TH2* h2_electron_Reco_histo_LowPt_(0);
TH2* h2_muon_ID_histo_(0);
TH2* h2_muon_ISO_histo_(0);


enum readoutFormat{undefined, pt_vs_eta, eta_vs_pt, pt_vs_absEta, absEta_vs_pt,
                      pt_vs_etaSc, etaSc_vs_pt, pt_vs_absEtaSc, absEtaSc_vs_pt};

/// format of the provided histogram in case of seperate histograms
readoutFormat ele_format_(undefined);
readoutFormat muon_format_(undefined);

int leadinglepton = -999;
int subleadinglepton = -999;

unsigned long long timestamp = 0;
VLV* v_leptons = 0;  
VLV* v_jets = 0;
LV* met = 0;  
vector<int>* lepPdgId = 0;
vector<float>* lepPfIso = 0;
vector<int>*   lepID_MuonTight = 0;
vector<int>*   lepID_ElecCutBased = 0;
vector<float>* lepSCEta = 0;
vector<int>*   lepID_ElecMVAIso = 0;
vector<int>*   lepID_MVATTH = 0;
vector<unsigned long long>* HLTBits = 0 ;
vector<int>*   jetPFID = 0;
Bool_t metFilters;
Int_t vertMultiTrue_ = -999;
vector<float>* jetBTagDeepCSV = 0;
vector<float>* jetBTagDeepJet = 0;
vector<float>* lepEnergyCorr = 0;
Int_t NPV_all = -999;  // number of reconstructed PVs
Int_t NPV_good = -999; // number of reconstructed PVs passing reconstruction-quality cuts
float GenWeight = 0;

double Omega_UNIX = 2*TMath::Pi()/ 86164.09053083288;
double Omega_Sidereal = 2*TMath::Pi()/86400.;
double phi_UNIX = 1.7493;
double phi_longitude = 1.5336;
double t0;
double t_sidereal;
double hour_sidereal;


TBranch* b_leptons = 0;   
TBranch* b_jets = 0;
TBranch* b_met = 0;   
TBranch* b_lepPdgId = 0;
TBranch* b_lepPfIso = 0;
TBranch* b_lepID_MuonTight = 0;
TBranch* b_lepID_ElecCutBased = 0; 
TBranch* b_lepSCEta = 0;
TBranch* b_lepID_ElecMVAIso = 0;
TBranch* b_lepID_MVATTH = 0;
TBranch* b_HLTBits = 0;
TBranch* b_jetPFID = 0;
TBranch* b_metFilters = 0;
TBranch* b_vertMultiTrue = 0;
TBranch* b_NPV_all = 0;
TBranch* b_NPV_good = 0;
TBranch* b_jetBTagDeepCSV = 0;
TBranch* b_jetBTagDeepJet = 0;
TBranch* b_lepEnergyCorr = 0;
TBranch* b_timestamp = 0;

// nTuple branch for generator event weight
TBranch* b_weightGenerator = 0;



/// Functions ///

int sign(int x) { return (x > 0) - (x < 0); }


Double_t deltaPhi(Double_t phi1, Double_t phi2) {
  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
  return dPhi;
}


Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {
  Double_t dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);
  return sqrt(dEta*dEta+dPhi*dPhi);
}


// To get HLTPaths 
void SetHLTBitsManager(const BitsManager<string>& hbm){ HLTBitsManager_ = hbm; }


readoutFormat convertReadoutFormat(string readoutFormat_)
{

    if(readoutFormat_=="pt_vs_eta") return pt_vs_eta;
    if(readoutFormat_=="eta_vs_pt") return eta_vs_pt;
    if(readoutFormat_=="pt_vs_absEta") return pt_vs_absEta;
    if(readoutFormat_=="absEta_vs_pt") return absEta_vs_pt;
    if(readoutFormat_=="pt_vs_etaSc") return pt_vs_etaSc;
    if(readoutFormat_=="etaSc_vs_pt") return etaSc_vs_pt;
    if(readoutFormat_=="pt_vs_absEtaSc") return pt_vs_absEtaSc;
    if(readoutFormat_=="absEtaSc_vs_pt") return absEtaSc_vs_pt;
    cerr<<"ERROR in constructor of LeptonScaleFactors! No readout format given!"
             <<"\n...break\n"<<endl;
    exit(764);
    return undefined;
}


/// Get 2-dimensional scale factor from histogram
float get2DSF(const TH2* const histo, const float& x, const float& y)
{
    int xbin, ybin, dummy;
    histo->GetBinXYZ(histo->FindFixBin(x, y), xbin, ybin, dummy);
    //overflow to last bin
    xbin = min(xbin, histo->GetNbinsX());
    ybin = min(ybin, histo->GetNbinsY());
    return histo->GetBinContent(xbin, ybin);
}


/// Return the scale factor histogram
TH2* prepareSF(string sfInputFileName, string histogramName) //const SystematicInternal& systematic)const
{
    if(sfInputFileName.empty()){
        cout<<"Input filename is empty - Set scale factors = 1.\n";
        return 0;
    }

    // Check if file containing scale factors exists, and open in case
    string inputFileName(DATA_PATH_COMMON);
    inputFileName.append("/");
    inputFileName.append(sfInputFileName);
    ifstream inputfile(inputFileName);
    if(!inputfile.is_open()){
        cerr<<"Error in LeptonScaleFactors::prepareSF()! File containing lepton ID/iso scale factors not found: "<<inputFileName
                 <<"\n...break\n"<<endl;
        exit(39);
    }
    inputfile.close();
    TFile scaleFactorFile(inputFileName.c_str());

    // Access histogram containing scale factors
    TH2* h_scaleFactorPtEta(0);
    h_scaleFactorPtEta = dynamic_cast<TH2*>(scaleFactorFile.Get(histogramName.c_str()));
    if(!h_scaleFactorPtEta){
        cerr<<"Error in LeptonScaleFactors::prepareSF()! TH2 for lepton Id/Iso scale factors not found: "<<histogramName
                 <<"\n...break\n"<<endl;
        exit(39);
    }

    // Store histogram in memory and close file
    h_scaleFactorPtEta->SetDirectory(0);
    scaleFactorFile.Close();

    return h_scaleFactorPtEta;
}


TH2* combineUncHistos(const TH2* nomHisto, const TH2* statHisto, const TH2* sysHisto)
{
  TH2* outputHist= (TH2*) nomHisto->Clone();
  int nX = nomHisto->GetNbinsX();
  int nY = nomHisto->GetNbinsY();
  for(int x=0; x<nX+2; x++){
    for(int y=0; y<nY+2; y++){
      if(nomHisto->GetBinContent(x, y)&&statHisto->GetBinContent(x, y)&&sysHisto->GetBinContent(x, y)){
        float sysErr=sysHisto->GetBinError(x, y);
        float statErr=statHisto->GetBinError(x, y);
        float totalErr=TMath::Sqrt((sysErr*sysErr) + (statErr*statErr));
        outputHist->SetBinError(x,y,totalErr);
      }
    }
  }
  //Fill Under/Overflow in case they didn't exist before
  for (int x=0; x<nX+2; x++) {
    if (!outputHist->GetBinContent(x, 0)){
      outputHist->SetBinContent(x, 0, outputHist->GetBinContent(x, 1));
      outputHist->SetBinError(x, 0, outputHist->GetBinError(x, 1));
    }
    if (!outputHist->GetBinContent(x, nY+1)){
      outputHist->SetBinContent(x, nY+1, outputHist->GetBinContent(x, nY));
      outputHist->SetBinError(x, nY+1, outputHist->GetBinError(x, nY));
    }
  }
  for (int y=0; y<nY+2; y++) {
          if (!outputHist->GetBinContent(0, y)){
            outputHist->SetBinContent(0, y, outputHist->GetBinContent(1, y));
            outputHist->SetBinError(0, y, outputHist->GetBinError(1, y));
          }
          if (!outputHist->GetBinContent(nX+1, y)){
            outputHist->SetBinContent(nX+1, y, outputHist->GetBinContent(nX, y));
            outputHist->SetBinError(nX+1, y, outputHist->GetBinError(nX, y));
          }
  }
  outputHist->SetDirectory(0);
  return outputHist;
}


TH2* appendFlowBins(const TH2* histo)
{
  TH2* outputHist= (TH2*) histo->Clone();
  int nX = histo->GetNbinsX();
  int nY = histo->GetNbinsY();
  for(int x=0; x<nX+2; x++){
    for(int y=0; y<nY+2; y++){
      if(histo->GetBinContent(x, y)){
        float val=histo->GetBinContent(x, y);
        float err=histo->GetBinError(x, y);
        outputHist->SetBinContent(x,y,val);
        outputHist->SetBinError(x,y,err);
      }
    }
  }
  //Fill Under/Overflow in case they didn't exist before
  for (int x=0; x<nX+2; x++) {
    if (!outputHist->GetBinContent(x, 0)){
      outputHist->SetBinContent(x, 0, outputHist->GetBinContent(x, 1));
      outputHist->SetBinError(x, 0, outputHist->GetBinError(x, 1));
    }
    if (!outputHist->GetBinContent(x, nY+1)){
      outputHist->SetBinContent(x, nY+1, outputHist->GetBinContent(x, nY));
      outputHist->SetBinError(x, nY+1, outputHist->GetBinError(x, nY));
    }
  }
  for (int y=0; y<nY+2; y++) {
          if (!outputHist->GetBinContent(0, y)){
            outputHist->SetBinContent(0, y, outputHist->GetBinContent(1, y));
            outputHist->SetBinError(0, y, outputHist->GetBinError(1, y));
          }
          if (!outputHist->GetBinContent(nX+1, y)){
            outputHist->SetBinContent(nX+1, y, outputHist->GetBinContent(nX, y));
            outputHist->SetBinError(nX+1, y, outputHist->GetBinError(nX, y));
          }
  }
  outputHist->SetDirectory(0);
  return outputHist;
}



void prepareLeptonSF() {
    cout<<"--- Beginning preparation of lepton scale factors\n";

    // Access electron scale factor
    cout<<"Electron:\n";
    h2_electron_ID_histo_ = prepareSF(SF_electron_ID_File, SF_electron_ID_Histo);
    if(h2_electron_ID_histo_) cout<<"ID scale factors found - will be used\n";
    // cout<<"histo(5,1) val = " << h2_electron_ID_histo_->GetBinContent(5, 1) << endl; 

    h2_electron_Reco_histo_ = prepareSF(SF_electron_Reco_File, SF_electron_Reco_Histo);
    if(h2_electron_Reco_histo_) cout<<"Reconstruction scale factors found - will be used\n";
    
    h2_electron_Reco_histo_LowPt_ = prepareSF(SF_electron_Reco_File_LowPt, SF_electron_Reco_Histo_LowPt);
    if(h2_electron_Reco_histo_LowPt_) cout<<"Reconstruction scale factors (lowPt) found - will be used\n";

    h2_electron_ID_histo_ = appendFlowBins(h2_electron_ID_histo_);
    h2_electron_Reco_histo_ = appendFlowBins(h2_electron_Reco_histo_);
    h2_electron_Reco_histo_LowPt_ = appendFlowBins(h2_electron_Reco_histo_LowPt_);


    // Access muon scale factor
    cout<<"Muon:\n";
    h2_muon_ID_histo_ = prepareSF(SF_muon_ID_File, SF_muon_ID_Histo);
    cout<<"ID nominal factors found\n";

    TH2* h2_muon_ID_histo_Stat = prepareSF(SF_muon_ID_File, SF_muon_ID_Histo_Stat);
    cout<<"ID stat. Unc. scale factors found\n";

    TH2* h2_muon_ID_histo_Sys = prepareSF(SF_muon_ID_File, SF_muon_ID_Histo_Sys);
    cout<<"ID syst. Unc. factors found\n";

    TH2* h2_muon_ID_histo_Err = combineUncHistos(h2_muon_ID_histo_, h2_muon_ID_histo_Stat, h2_muon_ID_histo_Sys);
    if(h2_muon_ID_histo_Err) cout<<"ID scale factors found - will be used\n";

    h2_muon_ISO_histo_ = prepareSF(SF_muon_ISO_File, SF_muon_ISO_Histo);
    cout<<"ISO nominal factors found\n";

    TH2* h2_muon_ISO_histo_Stat = prepareSF(SF_muon_ISO_File, SF_muon_ISO_Histo_Stat);
    cout<<"ISO stat. Unc. scale factors found\n";

    TH2* h2_muon_ISO_histo_Sys = prepareSF(SF_muon_ISO_File, SF_muon_ISO_Histo_Sys);
    cout<<"ISO syst. Unc. factors found\n";

    TH2* h2_muon_ISO_histo_Err = combineUncHistos(h2_muon_ISO_histo_, h2_muon_ISO_histo_Stat, h2_muon_ISO_histo_Sys);
    if(h2_muon_ISO_histo_Err) cout<<"iso scale factors found - will be used\n";

    h2_muon_ID_histo_ = appendFlowBins(h2_muon_ID_histo_Err);
    h2_muon_ISO_histo_ = appendFlowBins(h2_muon_ISO_histo_Err);

    cout<<"--- Finishing preparation of lepton scale factors\n\n";

    return;
}



/// Get per-lepton scale factor
float getLeptonSF(const LV lepton, int pdgId, float SCEta = -999.) {

    float result(1.);
    int absPdgId(abs(pdgId));
    const TH2* histo(0);
    const TH2* histo2(0);


    ele_format_ = convertReadoutFormat(SF_electron_ID_Format); 
    // cout<< "ele_format = " << ele_format_ << endl;
    muon_format_ = convertReadoutFormat(SF_muon_ID_Format);


    if(absPdgId == 11){
      histo = (lepton.Pt()<20.) ? h2_electron_Reco_histo_LowPt_ : h2_electron_Reco_histo_;
      histo2 = h2_electron_ID_histo_;
    }
    else if(absPdgId == 13){
      histo = h2_muon_ISO_histo_;
      histo2 = h2_muon_ID_histo_;
    }
    else{
        cout<<"WARNING in method scaleFactorAllLeptons! LeptonPdgId is not as expected: "
                <<absPdgId<<"\n...will use scale factor = 1.\n";
        return 1.;
    }
    if(!histo) return 1.; 
    if(absPdgId==11) {
      switch(ele_format_)
      {
        case pt_vs_eta:
          result *= get2DSF(histo, lepton.eta(), lepton.pt());
          result *= get2DSF(histo2, lepton.eta(), lepton.pt());
          break;
        case eta_vs_pt:
          result *= get2DSF(histo, lepton.pt(), lepton.eta());
          result *= get2DSF(histo2, lepton.pt(), lepton.eta());
          break;
        case pt_vs_absEta:
          result *= get2DSF(histo, fabs(lepton.eta()), lepton.pt());
          result *= get2DSF(histo2, fabs(lepton.eta()), lepton.pt());
          break;
        case absEta_vs_pt:
          result *= get2DSF(histo, lepton.pt(), fabs(lepton.eta()));
          result *= get2DSF(histo2, lepton.pt(), fabs(lepton.eta()));
          break;
        case pt_vs_etaSc:
          result *= get2DSF(histo, SCEta, lepton.pt());
        //   cout<<"SF pt_vs_etaSc reco = " << get2DSF(histo, SCEta, lepton.pt()) << endl; 
        //   cout<<"result pt_vs_etaSc reco = " << result << endl; 
          result *= get2DSF(histo2, SCEta, lepton.pt());
        //   cout<<"SF pt_vs_etaSc id = " << get2DSF(histo2, SCEta, lepton.pt()) << endl; 
        //   cout<<"result pt_vs_etaSc id = " << result << endl; 
          break;
        case etaSc_vs_pt:
          result *= get2DSF(histo, lepton.pt(), SCEta);
          result *= get2DSF(histo2, lepton.pt(), SCEta);
          // cout<<"result etaSc_vs_pt = " << result << endl; 
          break;
        case pt_vs_absEtaSc:
          result *= get2DSF(histo, fabs(SCEta), lepton.pt());
          result *= get2DSF(histo2, fabs(SCEta), lepton.pt());
          break;
        case absEtaSc_vs_pt:
          result *= get2DSF(histo, lepton.pt(), fabs(SCEta));
          result *= get2DSF(histo2, lepton.pt(), fabs(SCEta));
          break;
        default:
          break;
      }
    }
      if(absPdgId==13) {
        switch(muon_format_)
        {
          case pt_vs_eta:
            result *= get2DSF(histo, lepton.eta(), lepton.pt());
            result *= get2DSF(histo2, lepton.eta(), lepton.pt());
            break;
          case eta_vs_pt:
            result *= get2DSF(histo, lepton.pt(), lepton.eta());
            result *= get2DSF(histo2, lepton.pt(), lepton.eta());
            break;
          case pt_vs_absEta:
            result *= get2DSF(histo, fabs(lepton.eta()), lepton.pt());
            result *= get2DSF(histo2, fabs(lepton.eta()), lepton.pt());
            break;
          case absEta_vs_pt:
            result *= get2DSF(histo, lepton.pt(), fabs(lepton.eta()));
            // cout<<"SF absEta_vs_pt iso = " << get2DSF(histo, lepton.pt(), fabs(lepton.eta())) << endl; 
            // cout<<"result absEta_vs_pt iso = " << result << endl; 
            result *= get2DSF(histo2, lepton.pt(), fabs(lepton.eta()));
            // cout<<"SF absEta_vs_pt id = " << get2DSF(histo2, lepton.pt(), fabs(lepton.eta())) << endl; 
            // cout<<"result absEta_vs_pt id = " << result << endl; 
            break;
          case pt_vs_etaSc:
            result *= get2DSF(histo, SCEta, lepton.pt());
            result *= get2DSF(histo2, SCEta, lepton.pt());
            break;
          case etaSc_vs_pt:
            result *= get2DSF(histo, lepton.pt(), SCEta);
            result *= get2DSF(histo2, lepton.pt(), SCEta);
            break;
          case pt_vs_absEtaSc:
            result *= get2DSF(histo, fabs(SCEta), lepton.pt());
            result *= get2DSF(histo2, fabs(SCEta), lepton.pt());
            break;
          case absEtaSc_vs_pt:
            result *= get2DSF(histo, lepton.pt(), fabs(SCEta));
            result *= get2DSF(histo2, lepton.pt(), fabs(SCEta));
            break;
          default:
            break;
        }
      }

    return result;
}


float getWeightLeptonSF(int leadingLeptonIndex, int nLeadingLeptonIndex, VLV* leptons, vector<int>* lepPdgIds, vector<float>* lepSCEta){
    float result(1.);
    if(!isMC) return 1.;  
    if(leadingLeptonIndex<0 || nLeadingLeptonIndex<0) return 1.; 

    const LV leadingLepton(leptons->at(leadingLeptonIndex));
    const LV nLeadingLepton(leptons->at(nLeadingLeptonIndex));
    int leadingPdgId(lepPdgIds->at(leadingLeptonIndex));
    int nLeadingPdgId(lepPdgIds->at(nLeadingLeptonIndex));
    float leadingSCEta = -999.;
    float nLeadingSCEta = -999.;
    if(   ( leadingLeptonIndex < int(lepSCEta->size()))
       && (nLeadingLeptonIndex < int(lepSCEta->size())))
    {
       leadingSCEta  = lepSCEta->at(leadingLeptonIndex);
       nLeadingSCEta = lepSCEta->at(nLeadingLeptonIndex);
    }
    result *= getLeptonSF(leadingLepton, leadingPdgId, leadingSCEta);
    result *= getLeptonSF(nLeadingLepton, nLeadingPdgId, nLeadingSCEta);
    return result;
}


void preparePUSF(TString filename){

  std::cout<<"--- Beginning preparation of pileup reweighter\n";
  
  // Set data pileup distribution
  TString pileupInput(DATA_PATH_COMMON);
  TString pileupInputHisto("");
  pileupInput.Append("/").Append(pileup_DATA_File);
  pileupInputHisto.Append(pileup_DATA_Histogram);
  cout<<"Using data PU input file:\n"<<pileupInput<<endl;
  cout<<"Using data PU input histo:\n"<<pileupInputHisto<<endl;
  puReweighter_->setDataTruePUInput(pileupInput.Data(), pileupInputHisto.Data()); 

  // Set MC pileup distribution
  TString pileupInputMC("");
  TString pileupInputMCHisto("");
  string PU_MC_File = (pileup_MC_File == "@{INPUT_FILE}") ? string(filename) : DATA_PATH_COMMON+"/"+pileup_MC_File; 
  pileupInputMC.Append(PU_MC_File);
  pileupInputMCHisto.Append(pileup_MC_Histogram);
  cout<<"Using MC PU input file:\n"<<pileupInputMC<<endl;
  cout<<"Using MC PU input histo:\n"<<pileupInputMCHisto<<endl;
  puReweighter_->setMCTruePUInput(pileupInputMC.Data(), pileupInputMCHisto.Data()); 
  
  cout<<"--- Finishing preparation of pileup reweighter\n\n";

  return;
}

/// Get weight due to pileup reweighting
double getWeightPileup(const Long64_t& entry){
    double result(1.);
    if(!isMC) return 1.; 
    b_vertMultiTrue->GetEntry(entry); 
    result *= puReweighter_->getPUweight(vertMultiTrue_);
    return result;
}


double GetSiderealTime(double Omega_UNIX_, double t_CMS_, double t0_, double phi_UNIX_, double phi_longitude_, double Omega_Sidereal_){

  double t_sidereal_ = (Omega_UNIX_*(t_CMS_ - t0_) + phi_UNIX_ + phi_longitude_)/Omega_Sidereal_;
  return t_sidereal_;
}
