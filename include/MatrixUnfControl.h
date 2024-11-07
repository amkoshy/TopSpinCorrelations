#ifndef MatrixUnfControl_h
#define MatrixUnfControl_h

/***************************************
 * file: MatrixUnfControl.h
 * 
 * Header file for MatrixUnfControl
 * Author: Ajeeta Khatiwada, Jacob Linacre
 */
#include "stdInc.h"
#include "TRandom3.h"
#include "THStack.h"
#include <vector>
#include <string.h>
// for the unfolding!
#include "TUnfoldDensity.h"
#include "MatrixUnf.h"

const Int_t kMaxNumDataFiles = 200;
const Int_t kMaxNumBgFiles   = 200;
const Int_t kMaxNumSysFiles  = 200; 

class MatrixUnfControl{
 protected:
  //TString *fOutputName;
  //TFile *fMCFile; 
  TObjArray *fSysFiles; 
  TObjArray *fBgFiles; 
  TObjArray *fDataFiles;
  TObjArray *fMCFiles;
  
  //Double_t fXSectionSig;
  std::vector<Double_t>fXSectionData;
  std::vector<Double_t>fXSectionMC;
  std::vector<Double_t>fXSectionBg;
  std::vector<Double_t>fXSectionSys;

  std::vector<Double_t>fLumiData;
  std::vector<Double_t>fLumiMC;
  std::vector<Double_t>fLumiBg;
  std::vector<Double_t>fLumiSys;

  //std::vector<TString>fVariableNames;
  std::vector<TH1D*>fPseudoUnfResults;
  std::vector<TH1D*>fPseudoUnfResultsRebinnedA;
  std::vector<TH1D*>fPseudoUnfResultsRebinnedB;
  std::vector<TH1D*>fPseudoUnfAfbResultsRebinnedA;
  std::vector<TH1D*>fPseudoUnfAfbResultsRebinnedB;
  std::vector<TH1D*>fPseudoUnfOtherResultsRebinnedA;
  std::vector<TH1D*>fPseudoUnfOtherResultsRebinnedB;

  TObjArray *fDataArray;
  TObjArray *fRecArray;
  TObjArray *fGenArray;
  TObjArray *fVisGenArray;
  TObjArray *fMigMatArray;
  TObjArray *fBgArray;
  TObjArray *fBgSamplesArray;
  TObjArray *fMigMatSysArray;
  TObjArray *fTrueDataArray; //for closure tests this is gen info of the MC used as Data
  
  std::string fEra; 
  double fLumi;
  double fTauMin;//change the values from 1e-20(small regularization)
  double fTauMax;// to 1(completely regularized), 0 searches for itself
  float fTauSys;
  bool fMinRhoAVG;
  bool fClosureTest;
  bool fSelfConsistency;
  bool fEnsembleTest;
  bool fLinearityTest;
  bool fAllowRebinOption;
  bool fCombinedSumptt;
  bool fConstrainTau;
  bool fSubtrBgMCFromData;
  bool fDoSubtrBg;
  Int_t fUseAreaConstrain;
  bool fUseBiasInUnfold;
  Int_t fRegMode;//1=size,2=derivative,3=curvature,4=mixed,5=densityCurvature
  bool fsetSyst; 


  //histograms to store linearity summary results (showing all variables in one plot)
  TH1D* hLinearityMaxBiasFrac;
  TH1D* hLinearityMaxBiasFracStat;
  TH1D* hLinearityMaxBiasSlope;
  TH1D* hLinearityMaxBiasFracMax;
  TH1D* hLinearityMaxBiasFracStatMax;
  TH1D* hLinearityMaxBiasSlopeMax;
  TH1D* hLinearityAvgBiasFrac;
  TH1D* hLinearityAvgBiasFracStat;
  TH1D* hLinearityAvgBiasSlope;

  TH1D* hLinearityMaxBiasFrac_rebinnedA;
  TH1D* hLinearityMaxBiasFracStat_rebinnedA;
  TH1D* hLinearityMaxBiasSlope_rebinnedA;
  TH1D* hLinearityMaxBiasFracMax_rebinnedA;
  TH1D* hLinearityMaxBiasFracStatMax_rebinnedA;
  TH1D* hLinearityMaxBiasSlopeMax_rebinnedA;
  TH1D* hLinearityAvgBiasFrac_rebinnedA;
  TH1D* hLinearityAvgBiasFracStat_rebinnedA;
  TH1D* hLinearityAvgBiasSlope_rebinnedA;

  TH1D* hLinearityMaxBiasFrac_rebinnedB;
  TH1D* hLinearityMaxBiasFracStat_rebinnedB;
  TH1D* hLinearityMaxBiasSlope_rebinnedB;
  TH1D* hLinearityMaxBiasFracMax_rebinnedB;
  TH1D* hLinearityMaxBiasFracStatMax_rebinnedB;
  TH1D* hLinearityMaxBiasSlopeMax_rebinnedB;
  TH1D* hLinearityAvgBiasFrac_rebinnedB;
  TH1D* hLinearityAvgBiasFracStat_rebinnedB;
  TH1D* hLinearityAvgBiasSlope_rebinnedB;


  TH1D* hAfbLinearityMaxBiasFrac_rebinnedA;
  TH1D* hAfbLinearityMaxBiasFracStat_rebinnedA;
  TH1D* hAfbLinearityMaxBiasSlope_rebinnedA;
  TH1D* hAfbLinearityMaxBiasFracMax_rebinnedA;
  TH1D* hAfbLinearityMaxBiasFracStatMax_rebinnedA;
  TH1D* hAfbLinearityMaxBiasSlopeMax_rebinnedA;
  TH1D* hAfbLinearityAvgBiasFrac_rebinnedA;
  TH1D* hAfbLinearityAvgBiasFracStat_rebinnedA;
  TH1D* hAfbLinearityAvgBiasSlope_rebinnedA;

  TH1D* hAfbLinearityMaxBiasFrac_rebinnedB;
  TH1D* hAfbLinearityMaxBiasFracStat_rebinnedB;
  TH1D* hAfbLinearityMaxBiasSlope_rebinnedB;
  TH1D* hAfbLinearityMaxBiasFracMax_rebinnedB;
  TH1D* hAfbLinearityMaxBiasFracStatMax_rebinnedB;
  TH1D* hAfbLinearityMaxBiasSlopeMax_rebinnedB;
  TH1D* hAfbLinearityAvgBiasFrac_rebinnedB;
  TH1D* hAfbLinearityAvgBiasFracStat_rebinnedB;
  TH1D* hAfbLinearityAvgBiasSlope_rebinnedB;


  //unused TH1D* hAfbConstantMaxBiasFrac_rebinnedA;
  //unused TH1D* hAfbConstantMaxBiasFracStat_rebinnedA;
  TH1D* hAfbConstantMaxBiasSlope_rebinnedA;
  //unused TH1D* hAfbConstantMaxBiasFracMax_rebinnedA;
  //unused TH1D* hAfbConstantMaxBiasFracStatMax_rebinnedA;
  //unused TH1D* hAfbConstantMaxBiasSlopeMax_rebinnedA;
  //unused TH1D* hAfbConstantAvgBiasFrac_rebinnedA;
  //unused TH1D* hAfbConstantAvgBiasFracStat_rebinnedA;
  TH1D* hAfbConstantAvgBiasSlope_rebinnedA;

  //unused TH1D* hAfbConstantMaxBiasFrac_rebinnedB;
  //unused TH1D* hAfbConstantMaxBiasFracStat_rebinnedB;
  TH1D* hAfbConstantMaxBiasSlope_rebinnedB;
  //unused TH1D* hAfbConstantMaxBiasFracMax_rebinnedB;
  //unused TH1D* hAfbConstantMaxBiasFracStatMax_rebinnedB;
  //unused TH1D* hAfbConstantMaxBiasSlopeMax_rebinnedB;
  //unused TH1D* hAfbConstantAvgBiasFrac_rebinnedB;
  //unused TH1D* hAfbConstantAvgBiasFracStat_rebinnedB;
  TH1D* hAfbConstantAvgBiasSlope_rebinnedB;


  bool AddFile(TString filename, double xsec_norm, TObjArray*&filelist, std::vector<Double_t>&norm_list, std::vector<Double_t>&lumi_list);  
  //bool AddFile(TString filename, double xsec_norm, double nevt, TObjArray *filelist,std::vector<Double_t>norm_list,std::vector<Double_t>evts_list);  
  //void SetOutputFileName(TString outputname){fOutputName = outputname;}

 public:
  
  //MatrixUnfControl(bool constrainTau=kTRUE,float tausys=1.,double taumin=1e-3, double taumax=7e-1,bool allowrebinopt=kTRUE,bool combinedsumptt=kTRUE,bool usesigfracbgsubt=kFALSE,bool subtbgmcfromdata=kTRUE,bool dosubtrbg=kFALSE, bool usebiasintunfold=kTRUE,int useareaconstrain = 0,int regmode=2){}//default constructor
  MatrixUnfControl();
  void DoAllUnfold(TString preunfoldedfile, TString outputfile, TString channel);
  void MakeAllPlots(TString unfoldedfile, TString plotfile, TString channel, bool doOptTotalUnc);
  void doPlotPulls(TString channel, TString plotfile);
  bool SetTauOptions(bool constrainTau=kTRUE,float tausys=1.,double taumin=1e-3, double taumax=7e-1, bool minrhoavg=kFALSE);
  bool SetOptions(std::string era, double lumi, bool closuretest=kFALSE, bool selfconsistency=kFALSE, bool ensembletest=kFALSE, bool linearitytest=kFALSE, bool allowrebinopt=kTRUE,bool combinedsumptt=kTRUE,bool subtbgmcfromdata=kTRUE,bool dosubtrbg=kFALSE, bool usebiasintunfold = kTRUE,int useareaconstrain = 0,int regmode=3, bool setsys = kFALSE);
  bool AddUnfoldVariable(std::string varname);
  bool AddUnfoldVariables(TString varlist);
  bool AddSystematics(std::vector<std::string> vecinputsyst, std::vector<std::string> vecsyst);
  bool AddBgFile(TString filename, double xsec_norm){AddFile(filename,xsec_norm,fBgFiles,fXSectionBg,fLumiBg); return kTRUE;}
  bool AddDataFile(TString filename, double xsec_norm){AddFile(filename,xsec_norm,fDataFiles,fXSectionData,fLumiData); return kTRUE;}
  bool AddMCFile(TString filename, double xsec_norm){AddFile(filename,xsec_norm,fMCFiles,fXSectionMC,fLumiMC); return kTRUE;}
  //bool OpenMCFile(TString filename, Double_t xsec_norm);
  bool AddSysFile(TString filename, Double_t xsec_norm){AddFile(filename, xsec_norm, fSysFiles, fXSectionSys, fLumiSys); return kTRUE;}
  bool MakePseudoHistos(TH1D* inputhisto, TH1D*& pseudohisto, TString varname);
  void CalculateHistRange(std::vector<TH1D*>PseudoUnfResults, int i_bin, double &mean, double &rms, bool iserror, bool wholerange=true);
  void CreateLinearitySummaryHistos(int numCoefficients);
  bool MakeLinearityHistos(TH2D *inputmighisto, TH1D*& linearityhisto, TH1D*& linearitygenhisto, TString varname, float slope, const TUnfoldBinning* rebinnedBinning, const TUnfoldBinning* originalBinning, int numVariables, int rebin=1, bool reweightusingwidebincentres=false, int binToReweight=-1);
  void simpleFit(TGraph* graph, TH1D* hist, TH1D* histMax, int i_bin, double DeltaCtotal, double &param, double &parammax, bool isslope = false);
  void fillLinearityGraphs(std::vector<TH1D*>PseudoUnfResults, int rebin, bool isasym, TString nametag, int nLinearity, int ith_variable, int unf_nbins, TH1D* hLinearityGen[], double delta_coefficient[], TH1D* hLinearityDBiasFractional, TH1D* hLinearityDBiasFractionalMax, TH1D* hLinearityDBias, TH1D* hLinearityDBiasMax, TH1D* hLinearityDInjectedShape, TH1D* hLinearityDInjectedShapeMax, TH1D* hLinearityDMeasuredShape, TH1D* hLinearityDMeasuredShapeMax, TH1D* hLinearitySlope, TH1D* hLinearitySlopeMax, TH1D* hLinearityDBiasOverStatErr, TH1D* hLinearityDBiasOverStatErrMax, TH1D* hLinearityCBias, TH1D* hLinearityCBiasFractional, TH1D* hLinearityCBiasOverStatErr,    TH1D* hLinearityMaxBiasFrac, TH1D* hLinearityMaxBiasFracStat, TH1D* hLinearityMaxBiasSlope, TH1D* hLinearityMaxBiasFracMax, TH1D* hLinearityMaxBiasFracStatMax, TH1D* hLinearityMaxBiasSlopeMax, TH1D* hLinearityAvgBiasFrac, TH1D* hLinearityAvgBiasFracStat, TH1D* hLinearityAvgBiasSlope, TH1D* hConstantMaxBiasSlope, TH1D* hConstantAvgBiasSlope, int numCoefficients);
  void calculateCoefficientUncertainty(TH2D* hsystcovmat, TH1D* tunfresult_temp, TH1D* genhistogram, std::vector<double> &a_optimal_reco, TString VariableName, double *binsum_temp, int numCoefficients, double Cerr[]);
  double calculateBinSumUncertaintySnake(TH2D* hsystcovmat, TH1D* tunfresult_temp);
  TH2D* calculateCoefficientSystCorrmAllVars(TString unfoldedfile, TH2D* hsystcovmat, TH1D* hist_acombfactor[], TString NameTag, int inorm);
  TH2D* calculateBinSumSystCorrmAllVars(TString unfoldedfile, TH2D* hsystcovmat, TH1D* hist_acombfactor[], TString NameTag, int inorm);
  void FillCorrelationMatrixAllVariables(int nbins, TString namestring, std::vector<std::vector<std::vector<TH1D*>>>hPseudoUnfDeltapull, TH2D* &h_RhoVV, TMatrixD &m_RhoVV);
  void calculateCorrelationMatrices(int nVar, int nbinsrhoi, int nPE, bool docorrectpullwidth, TString rebinnedstring, TFile* BootstrapCorrelationMatricesFile, TH1D* hAllVarUnf, TH1D* hist_acombfactor[], std::vector<std::vector<std::vector<float>>> &BinPE, std::vector<std::vector<std::vector<float>>> &CoefPE, std::vector<double> &pullwidthBins, std::vector<double> &pullwidthVBins);
  void CalculateChannelFractions(double bkgsubdata_yields[], double sig_yields[], double channel_SF[]);
  void CorrectChannelFractions(TH2D*& htemprecoVsgen, TString MCFileName, double channel_SF[]);
  void sigfracBgSubtMethod(TH1D*& data_subt, TH1D* sig, TH1D* bg, TH1D* data);
  void regularBgSubtMethod(TH1D*& data_subt, TH1D* bg, TH1D* data);
  Double_t getLumiWeight(TFile*f,Double_t XSection,Double_t Lumi);
  //Int_t unfoldData(TString NameTUnfObjX,TString xAxisTitle,TString NameTUnfObjY,TString yAxisTitle,TString yAxisTitle2,TString binning,TH1D* input, TH2D* MigMat, TH1D* truth, TH1D* visgen, TH1D* reco, TH1D* bg, Double_t tauMin, Double_t tauMax, bool ensTest, TString channel);
  //Int_t unfoldData(TString NameTUnfObjX,TString xAxisTitle,TString NameTUnfObjY,TString yAxisTitle,TString yAxisTitle2,TString binning,int rebinfine,TH1D* input, TH2D* MigMat, TH1D* truth, TH1D* visgen, TH1D* reco, TH1D* bg, Double_t tauMin, Double_t tauMax, bool ensTest, TString channel, std::vector<TH2D*>sys);
  Int_t unfoldData(MatrixUnf* pttMatUnf, TString NameTUnfObjX, TH1D* input, TH2D* MigMat, TH1D* truth, TH1D* visgen, TH1D* reco, std::vector<TH1D*> vecbg, Double_t tauMin, Double_t tauMax, bool ensTest, TString channel, std::vector<TH2D*>sys);
  bool FillHistoArrays(TString CHAN);
  std::vector<TString>fVariableNames;
  std::vector<std::string>fVectorOfValidSystematics;  
  std::vector<std::string>fVectorOfInputSystematics;  
  TH1D* GetEnvelopeForUnfoldedResults( TString Channel, TString SystName, TString Variable, TFile* unffile, Bool_t normXsec, TString NameTag);
  TH1D* CalcUpDownDifference ( TString Channel, TString SystName, TString Variable, TFile* unffile, Bool_t normXsec, TString NameTag);
  void DrawCMSLabels(double textSize, double lumi = 35922.0, int cmsprelim = 1, double energy = 13);
  void DrawDecayChLabel(TString decaychannel, double textSize);
  TH2D* GetCombinedCovarianceMatrix( TString Channel, TString SystName, TString Variable, TFile* unffile, TH1D* htot_reldelta, Bool_t normXsec, TString NameTag );

  void GetDeltasAllVars(TH1D*& hnom_AllVar, TH1D*& htot_reldelta_AllVar, TString Channel, TString SystName, TString unfoldedfile, Bool_t normXsec, TString NameTag );
  TH2D* GetCombinedCovarianceMatrixAllVars(TH1D *hnom_AllVar, TH1D *htot_reldelta_AllVar, TString Channel, TString SystName, TString unfoldedfile, Bool_t normXsec, TString NameTag );

  void setTheoryStyleAndFillLegend(TH1 *histo, TString theoryName, TLegend *leg = 0);
  void setResultLegendStyle(TLegend *leg, TString name);
  void setResultLegendStyleN(TLegend *leg, TString name);

  // Amandeep : For N-D bin sizes
  Double_t GetPhysicalBinWidth(TString VariableName, int BinNumber);
  void SetAxis(TString VariableName, TAxis *Axis);
  void Draw2DLabels(TString VariableName, double y1_NDC = 0.71, double text_size = 0.025);
  void Draw2DLabelsVertical(TString VariableName, double x1_NDC = 0.82, double text_size = 0.025);
  void DrawCov2DAxisLabels(TString name);
  void DrawResMat2DAxisLabels(TString name);
  void DrawTotEff2DAxisLabels(TString name);
  void DrawPurStab2DAxisLabels(TString name);

  // End
};


#endif
