#include "analysisStructs.h"





ttbar::GenLevelWeights::GenLevelWeights(const double& weightMadgraphCorrection, const double& weightPileup, const double& weightGenerator,
					const double& weightPDF, const double& weightPS, const double& weightAlphasPDF,
					const double& weightMEScale, const double& weightMEFacScale, const double& weightMERenScale,
					const double& weightBFrag, const double& weightBSemilep, const double& weightTopPt,
					const double& trueLevelWeightNoPileup, const double& trueLevelWeight):
weightMadgraphCorrection_(weightMadgraphCorrection),
weightPileup_(weightPileup),
weightGenerator_(weightGenerator),
weightPDF_(weightPDF),
weightPS_(weightPS),
weightAlphasPDF_(weightAlphasPDF),
weightMEScale_(weightMEScale),
weightMEFacScale_(weightMEFacScale),
weightMERenScale_(weightMERenScale),
weightBFrag_(weightBFrag),
weightBSemilep_(weightBSemilep),
weightTopPt_(weightTopPt),
trueLevelWeightNoPileup_(trueLevelWeightNoPileup),
trueLevelWeight_(trueLevelWeight)
{}



ttbar::RecoLevelWeights::RecoLevelWeights(const double& weightLeptonSF, const double& weightTriggerSF, const double& weightBtagSF,
                                        const double& weightNoPileup, const double& weightPU,const double& weight):
weightLeptonSF_(weightLeptonSF),
weightTriggerSF_(weightTriggerSF),
weightBtagSF_(weightBtagSF),
weightNoPileup_(weightNoPileup),
weightPU_(weightPU),
weight_(weight)
{}



ttbar::GenObjectIndices::GenObjectIndices(const int& genBjetFromTopIndex, const int& genAntiBjetFromTopIndex,
                                        const int& recoBjetFromTopIndex, const int& recoAntiBjetFromTopIndex,
                                        const int& genBjetFromHiggsIndex, const int& genAntiBjetFromHiggsIndex,
                                        const int& recoBjetFromHiggsIndex, const int& recoAntiBjetFromHiggsIndex,
                                        const std::vector<int>& genExtraJetIndicesIso08,
                                        const std::vector<int>& genVisJetIndices):
genBjetFromTopIndex_(genBjetFromTopIndex),
genAntiBjetFromTopIndex_(genAntiBjetFromTopIndex),
recoBjetFromTopIndex_(recoBjetFromTopIndex),
recoAntiBjetFromTopIndex_(recoAntiBjetFromTopIndex),
genBjetFromHiggsIndex_(genBjetFromHiggsIndex),
genAntiBjetFromHiggsIndex_(genAntiBjetFromHiggsIndex),
recoBjetFromHiggsIndex_(recoBjetFromHiggsIndex),
recoAntiBjetFromHiggsIndex_(recoAntiBjetFromHiggsIndex),
genExtraJetIndicesIso08_(genExtraJetIndicesIso08),
genVisJetIndices_(genVisJetIndices)
{}



bool ttbar::GenObjectIndices::uniqueGenTopMatching()const
{
    return genBjetFromTopIndex_>=0 && genAntiBjetFromTopIndex_>=0 && genBjetFromTopIndex_!=genAntiBjetFromTopIndex_;
}



bool ttbar::GenObjectIndices::uniqueRecoTopMatching()const
{
    if(!this->uniqueGenTopMatching()) return false;
    
    return recoBjetFromTopIndex_>=0 && recoAntiBjetFromTopIndex_>=0 && recoBjetFromTopIndex_!=recoAntiBjetFromTopIndex_;
}



bool ttbar::GenObjectIndices::uniqueGenHiggsMatching()const
{
    return genBjetFromHiggsIndex_>=0 && genAntiBjetFromHiggsIndex_>=0 && genBjetFromHiggsIndex_!=genAntiBjetFromHiggsIndex_;
}



bool ttbar::GenObjectIndices::uniqueRecoHiggsMatching()const
{
    if(!this->uniqueGenHiggsMatching()) return false;
    
    return recoBjetFromHiggsIndex_>=0 && recoAntiBjetFromHiggsIndex_>=0 && recoBjetFromHiggsIndex_!=recoAntiBjetFromHiggsIndex_;
}



bool ttbar::GenObjectIndices::uniqueGenMatching()const
{
    if(!this->uniqueGenTopMatching() || !this->uniqueGenHiggsMatching()) return false;
    
    return !(genBjetFromTopIndex_==genBjetFromHiggsIndex_ || genBjetFromTopIndex_==genAntiBjetFromHiggsIndex_ ||
             genAntiBjetFromTopIndex_==genBjetFromHiggsIndex_ || genAntiBjetFromTopIndex_==genAntiBjetFromHiggsIndex_);
}



bool ttbar::GenObjectIndices::uniqueRecoMatching()const
{
    if(!this->uniqueRecoTopMatching() || !this->uniqueRecoHiggsMatching() || !this->uniqueGenMatching()) return false;
    
    return !(recoBjetFromTopIndex_==recoBjetFromHiggsIndex_ || recoBjetFromTopIndex_==recoAntiBjetFromHiggsIndex_ ||
             recoAntiBjetFromTopIndex_==recoBjetFromHiggsIndex_ || recoAntiBjetFromTopIndex_==recoAntiBjetFromHiggsIndex_);
}



bool ttbar::GenObjectIndices::isCorrectPairFromTop(const int& bIndex, const int& antiBIndex)const
{
    return bIndex==recoBjetFromTopIndex_ && antiBIndex==recoAntiBjetFromTopIndex_;
}



bool ttbar::GenObjectIndices::isSwappedPairFromTop(const int& bIndex, const int& antiBIndex)const
{
    return bIndex==recoAntiBjetFromTopIndex_ && antiBIndex==recoBjetFromTopIndex_;
}



bool ttbar::GenObjectIndices::isPairFromTop(const int& bIndex, const int& antiBIndex)const
{
    return this->isCorrectPairFromTop(bIndex, antiBIndex) || this->isSwappedPairFromTop(bIndex, antiBIndex);
}



bool ttbar::GenObjectIndices::isCorrectPairFromHiggs(const int& bIndex, const int& antiBIndex)const
{
    return bIndex==recoBjetFromHiggsIndex_ && antiBIndex==recoAntiBjetFromHiggsIndex_;
}



bool ttbar::GenObjectIndices::isSwappedPairFromHiggs(const int& bIndex, const int& antiBIndex)const
{
    return bIndex==recoAntiBjetFromHiggsIndex_ && antiBIndex==recoBjetFromHiggsIndex_;
}



bool ttbar::GenObjectIndices::isPairFromHiggs(const int& bIndex, const int& antiBIndex)const
{
    return this->isCorrectPairFromHiggs(bIndex, antiBIndex) || this->isSwappedPairFromHiggs(bIndex, antiBIndex);
}



ttbar::RecoObjectIndices::RecoObjectIndices(const std::vector<int>& allLeptonIndices,
                                  const std::vector<int>& leptonIndices, const std::vector<int>& antiLeptonIndices,
                                  const int& leptonIndex, const int& antiLeptonIndex,
                                  const int& leadingLeptonIndex, const int& nLeadingLeptonIndex,
                                  const int& leptonXIndex, const int& leptonYIndex,
                                  const std::vector<int>& extraJetIndicesIso08,
                                  const std::vector<int>& jetIndices,
                                  const std::vector<int>& bjetIndices):
allLeptonIndices_(allLeptonIndices),
leptonIndices_(leptonIndices),
antiLeptonIndices_(antiLeptonIndices),
leptonIndex_(leptonIndex),
antiLeptonIndex_(antiLeptonIndex),
leadingLeptonIndex_(leadingLeptonIndex),
nLeadingLeptonIndex_(nLeadingLeptonIndex),
leptonXIndex_(leptonXIndex),
leptonYIndex_(leptonYIndex),
extraJetIndicesIso08_(extraJetIndicesIso08),
jetIndices_(jetIndices),
bjetIndices_(bjetIndices)
{}





