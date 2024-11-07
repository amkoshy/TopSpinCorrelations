#ifndef TreeHandlerTTBarPlain_h
#define TreeHandlerTTBarPlain_h

#include <vector>

class TString;
class TTree;

#include "TreeHandlerBase.h"

class VariablesBase;
class EventMetadata;
class RecoObjects;
class CommonGenObjects;
class TopGenObjects;
class KinematicReconstructionSolutions;
namespace ttbar{
    class RecoLevelWeights;
    class GenLevelWeights;
    class GenObjectIndices;
    class RecoObjectIndices;
}

// OZ 2017.03.15
// Created this by copy-paste VariablesTTBar, the idea is to hold only plain particle 4-momenta
// all comments below concerning "MVA" come from VariablesTTBar and seem to be wrong




class TopAnalysisPlain;

/// Class for handling the trees of input variables for MVA,
/// trying to identify the jets coming from (anti)b's from (anti)tops
class TreeHandlerTTBarPlain : public TreeHandlerBase{
    
public:
    
    /// Constructor for selection steps
    TreeHandlerTTBarPlain(const char* inputDir,
                          const std::vector<TString>& selectionStepsNoCategories);
    
    /// Destructor
    ~TreeHandlerTTBarPlain(){}
    
    void SetAnalysisPlainPtr(TopAnalysisPlain* ptr) { ptrTopAnalysisPlain_ = ptr; }
    
    
    
    
private:
    
    /// Fill all variables for given selection step
    virtual void fillVariables(const EventMetadata& eventMetadata, 
                               const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                               const TopGenObjects& topGenObjects,
                               const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                               const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                               const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                               const double& weight, const TString& step,
                               std::vector<VariablesBase*>& variables)const;

    /// Book branches for TTree variables (dummy method, override in inherited TreeHandler)
    virtual void bookBranches(TTree* tree, VariablesBase* const variables_)const;
    
    /// Fill branches for TTree variables
    virtual void fillBranches(TTree* tree, const std::vector<VariablesBase*>& v_variables);
    
    /// Import all branches from TTree
    //virtual void importBranches(TTree* tree, std::vector<VariablesBase*>& variables)const;
    
    // OZ 30.07.2017 pointer to access everything in the analysis class
    TopAnalysisPlain* ptrTopAnalysisPlain_;
};





#endif







