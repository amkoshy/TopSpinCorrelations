#ifndef TopAnalysisSpinVars_h
#define TopAnalysisSpinVars_h
#include <TH1.h>
#include <TH2.h>
#include "TopAnalysis.h"

// Ajeeta 2017.10.10
// Class for adding Spin correlation related variables in the TopAnalysis. 
// Overrides a few methods of the base class to store more variables.


class TopAnalysisSpinVars : public TopAnalysis
{
  using TopAnalysis::TopAnalysis;
  
    //Histograms
    TH1 *h_HypAntiLeptonBk, *h_HypAntiLeptonBj, *h_HypAntiLeptonBr, *h_HypAntiLeptonBq, *h_HypAntiLeptonBn;
    TH1 *h_HypLeptonBk, *h_HypLeptonBj, *h_HypLeptonBr, *h_HypLeptonBq, *h_HypLeptonBn;
    TH1 *h_HypLLBarBPnn,  *h_HypLLBarBMnn,  *h_HypLLBarBPrr,  *h_HypLLBarBMrr,  *h_HypLLBarBPkk,  *h_HypLLBarBMkk;
    TH1 *h_HypLLBarBPjj,  *h_HypLLBarBMjj,  *h_HypLLBarBPqq,  *h_HypLLBarBMqq;
    TH1 *h_HypLLBarCkk,  *h_HypLLBarCrr,  *h_HypLLBarCnn;
    TH1 *h_HypLLBarCrk,  *h_HypLLBarCkr,  *h_HypLLBarCnr,  *h_HypLLBarCrn,  *h_HypLLBarCnk,  *h_HypLLBarCkn;
    TH1 *h_HypLLBarCPrk,  *h_HypLLBarCMrk,  *h_HypLLBarCPnr,  *h_HypLLBarCMnr,  *h_HypLLBarCPnk,  *h_HypLLBarCMnk;
    TH1 *h_HypLLBarcHel,  *h_HypLLBarcLab,  *h_HypLLBarkNorm,  *h_HypLLBarrNorm,  *h_HypLLBarnNorm;

    TH2 *h_GenRecoAntiLeptonBk, *h_GenRecoAntiLeptonBj, *h_GenRecoAntiLeptonBr, *h_GenRecoAntiLeptonBq, *h_GenRecoAntiLeptonBn;
    TH2 *h_GenRecoLeptonBk, *h_GenRecoLeptonBj, *h_GenRecoLeptonBr, *h_GenRecoLeptonBq, *h_GenRecoLeptonBn;
    TH2 *h_GenRecoLLBarBPnn,  *h_GenRecoLLBarBMnn,  *h_GenRecoLLBarBPrr,  *h_GenRecoLLBarBMrr,  *h_GenRecoLLBarBPkk,  *h_GenRecoLLBarBMkk;
    TH2 *h_GenRecoLLBarBPjj,  *h_GenRecoLLBarBMjj,  *h_GenRecoLLBarBPqq,  *h_GenRecoLLBarBMqq;
    TH2 *h_GenRecoLLBarCkk,  *h_GenRecoLLBarCrr,  *h_GenRecoLLBarCnn;
    TH2 *h_GenRecoLLBarCrk,  *h_GenRecoLLBarCkr,  *h_GenRecoLLBarCnr,  *h_GenRecoLLBarCrn,  *h_GenRecoLLBarCnk,  *h_GenRecoLLBarCkn;
    TH2 *h_GenRecoLLBarCPrk,  *h_GenRecoLLBarCMrk,  *h_GenRecoLLBarCPnr,  *h_GenRecoLLBarCMnr,  *h_GenRecoLLBarCPnk,  *h_GenRecoLLBarCMnk;
    TH2 *h_GenRecoLLBarcHel,  *h_GenRecoLLBarcLab,  *h_GenRecoLLBarkNorm,  *h_GenRecoLLBarrNorm,  *h_GenRecoLLBarnNorm;

    TH1 *h_VisGenAntiLeptonBk, *h_VisGenAntiLeptonBj, *h_VisGenAntiLeptonBr, *h_VisGenAntiLeptonBq, *h_VisGenAntiLeptonBn;
    TH1 *h_VisGenLeptonBk, *h_VisGenLeptonBj, *h_VisGenLeptonBr, *h_VisGenLeptonBq, *h_VisGenLeptonBn;
    TH1 *h_VisGenLLBarBPnn,  *h_VisGenLLBarBMnn,  *h_VisGenLLBarBPrr,  *h_VisGenLLBarBMrr,  *h_VisGenLLBarBPkk,  *h_VisGenLLBarBMkk;
    TH1 *h_VisGenLLBarBPjj,  *h_VisGenLLBarBMjj,  *h_VisGenLLBarBPqq,  *h_VisGenLLBarBMqq;
    TH1 *h_VisGenLLBarCkk,  *h_VisGenLLBarCrr,  *h_VisGenLLBarCnn;
    TH1 *h_VisGenLLBarCrk,  *h_VisGenLLBarCkr,  *h_VisGenLLBarCnr,  *h_VisGenLLBarCrn,  *h_VisGenLLBarCnk,  *h_VisGenLLBarCkn;
    TH1 *h_VisGenLLBarCPrk,  *h_VisGenLLBarCMrk,  *h_VisGenLLBarCPnr,  *h_VisGenLLBarCMnr,  *h_VisGenLLBarCPnk,  *h_VisGenLLBarCMnk;
    TH1 *h_VisGenLLBarcHel,  *h_VisGenLLBarcLab,  *h_VisGenLLBarkNorm,  *h_VisGenLLBarrNorm,  *h_VisGenLLBarnNorm;

    TH1 *h_RecoAntiLeptonBk, *h_RecoAntiLeptonBj, *h_RecoAntiLeptonBr, *h_RecoAntiLeptonBq, *h_RecoAntiLeptonBn;
    TH1 *h_RecoLeptonBk, *h_RecoLeptonBj, *h_RecoLeptonBr, *h_RecoLeptonBq, *h_RecoLeptonBn;
    TH1 *h_RecoLLBarBPnn,  *h_RecoLLBarBMnn,  *h_RecoLLBarBPrr,  *h_RecoLLBarBMrr,  *h_RecoLLBarBPkk,  *h_RecoLLBarBMkk;
    TH1 *h_RecoLLBarBPjj,  *h_RecoLLBarBMjj,  *h_RecoLLBarBPqq,  *h_RecoLLBarBMqq;
    TH1 *h_RecoLLBarCkk,  *h_RecoLLBarCrr,  *h_RecoLLBarCnn;
    TH1 *h_RecoLLBarCrk,  *h_RecoLLBarCkr,  *h_RecoLLBarCnr,  *h_RecoLLBarCrn,  *h_RecoLLBarCnk,  *h_RecoLLBarCkn;
    TH1 *h_RecoLLBarCPrk,  *h_RecoLLBarCMrk,  *h_RecoLLBarCPnr,  *h_RecoLLBarCMnr,  *h_RecoLLBarCPnk,  *h_RecoLLBarCMnk;
    TH1 *h_RecoLLBarcHel,  *h_RecoLLBarcLab,  *h_RecoLLBarkNorm,  *h_RecoLLBarrNorm,  *h_RecoLLBarnNorm;
  
  private:
    virtual void SlaveBegin(TTree*);
    virtual Bool_t Process(Long64_t entry);

};

#endif
