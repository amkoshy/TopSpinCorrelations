#include <iostream>
#include <cstdlib>

#include <TTree.h>
#include <TString.h>
#include <Rtypes.h>

#include "TreeHandlerTTBarPlain.h"
#include "VariablesBase.h"
#include "VariablesTTBarPlain.h"
#include "analysisStructs.h"
#include "../../common/include/analysisObjectStructs.h"
#include "../../common/include/KinematicReconstructionSolution.h"





TreeHandlerTTBarPlain::TreeHandlerTTBarPlain(const char* inputDir,
                                             const std::vector<TString>& selectionStepsNoCategories):
TreeHandlerBase("ttBar_", inputDir, selectionStepsNoCategories, new VariablesTTBarPlain())
{
    std::cout<<"--- Beginning setting up MVA tree handler for top jets assignment\n";
    std::cout<<"=== Finishing setting up MVA tree handler for top jets assignment\n\n";
}



void TreeHandlerTTBarPlain::fillVariables(const EventMetadata& eventMetadata,
                                     const RecoObjects& recoObjects, const CommonGenObjects& commonGenObjects,
                                          const TopGenObjects& topGenObjects,
                                          const KinematicReconstructionSolutions& kinematicReconstructionSolutions,
                                          const ttbar::RecoObjectIndices& recoObjectIndices, const ttbar::GenObjectIndices& genObjectIndices,
                                          const ttbar::GenLevelWeights& genLevelWeights, const ttbar::RecoLevelWeights& recoLevelWeights,
                                          const double& weight, const TString&,
                                          std::vector<VariablesBase*>& variables)const
{
    // Loop over all jet combinations and get MVA input variables
    const std::vector<VariablesBase*> v_variablesTTbar = 
            VariablesTTBarPlain::fillVariables(eventMetadata, recoObjects, commonGenObjects,topGenObjects,kinematicReconstructionSolutions, recoObjectIndices,  genObjectIndices,genLevelWeights,recoLevelWeights,weight, ptrTopAnalysisPlain_);
    
    // Fill the MVA variables
    variables.insert(variables.end(), v_variablesTTbar.begin(), v_variablesTTbar.end());
}



void TreeHandlerTTBarPlain::bookBranches(TTree* tree, VariablesBase* const variables_)const
{
    VariablesTTBarPlain* const variablesTTBar = dynamic_cast<VariablesTTBarPlain*>(variables_);
    
    //this->createBranch(tree, variablesTTBar->isTopGen_);
    this->createBranch(tree, variablesTTBar->entry_);
    //this->createBranch(tree, variablesTTBar->isKinReco_);
    
    // determine step (0 gen flag = 2, 8 rec flag = 1)
    const char& lastTreeChar = tree->GetName()[TString(tree->GetName()).Length() - 1];
    printf("TreeHandlerTTBarPlain booking branches for step %c\n", lastTreeChar);
    int flagLevel = 0;
    if(lastTreeChar == '0')
      flagLevel = 1;
    else if(lastTreeChar == '8')
      flagLevel = 2;
    else
      throw std::logic_error(TString::Format("Error in TreeHandlerTTBarPlain::bookBranches(): unknown step = %c\n", flagLevel));
    
    // >>>>>>>> reco >>>>>>>>>>>>>
    if(flagLevel == 2)
    {
      this->createBranch(tree, variablesTTBar->eventWeight_);

      this->createBranch(tree, variablesTTBar->jet_multiplicity_);
      this->createBranch(tree, variablesTTBar->jets_x_);
      this->createBranch(tree, variablesTTBar->jets_y_);
      this->createBranch(tree, variablesTTBar->jets_z_);
      this->createBranch(tree, variablesTTBar->jets_e_);

      this->createBranch(tree, variablesTTBar->bjet_multiplicity_);
      this->createBranch(tree, variablesTTBar->bjets_x_);
      this->createBranch(tree, variablesTTBar->bjets_y_);
      this->createBranch(tree, variablesTTBar->bjets_z_);
      this->createBranch(tree, variablesTTBar->bjets_e_);

      this->createBranch(tree, variablesTTBar->primvtx_multiplicity_);

      //this->createBranch(tree, variablesTTBar->met_x_);

      // top
      this->createBranch(tree, variablesTTBar->top1_p_x_);
      this->createBranch(tree, variablesTTBar->top1_p_y_);
      this->createBranch(tree, variablesTTBar->top1_p_z_);
      this->createBranch(tree, variablesTTBar->top1_p_e_);

      // antitop
      this->createBranch(tree, variablesTTBar->top2_p_x_);
      this->createBranch(tree, variablesTTBar->top2_p_y_);
      this->createBranch(tree, variablesTTBar->top2_p_z_);
      this->createBranch(tree, variablesTTBar->top2_p_e_);

      // lep1
      this->createBranch(tree, variablesTTBar->lep1_p_x_);
      this->createBranch(tree, variablesTTBar->lep1_p_y_);
      this->createBranch(tree, variablesTTBar->lep1_p_z_);
      this->createBranch(tree, variablesTTBar->lep1_p_e_);

      // lep2
      this->createBranch(tree, variablesTTBar->lep2_p_x_);
      this->createBranch(tree, variablesTTBar->lep2_p_y_);
      this->createBranch(tree, variablesTTBar->lep2_p_z_);
      this->createBranch(tree, variablesTTBar->lep2_p_e_);

      this->createBranch(tree, variablesTTBar->lepPdgId1_);
      this->createBranch(tree, variablesTTBar->lepPdgId2_);

      // MET
      this->createBranch(tree, variablesTTBar->met_pt_);
      //this->createBranch(tree, variablesTTBar->met_eta_);
      this->createBranch(tree, variablesTTBar->met_phi_);
    }

    // >>>>>>>> gen >>>>>>>>>>>>>
    if(flagLevel == 1)
    {
      this->createBranch(tree, variablesTTBar->trueLevelWeight_);

      this->createBranch(tree, variablesTTBar->gen_jet_multiplicity_);
      this->createBranch(tree, variablesTTBar->gen_jets_x_);
      this->createBranch(tree, variablesTTBar->gen_jets_y_);
      this->createBranch(tree, variablesTTBar->gen_jets_z_);
      this->createBranch(tree, variablesTTBar->gen_jets_e_);

      this->createBranch(tree, variablesTTBar->ttbar_decay_mode_);
      this->createBranch(tree, variablesTTBar->pdfWeights_);

      // top
      this->createBranch(tree, variablesTTBar->gen_top1_p_x_);
      this->createBranch(tree, variablesTTBar->gen_top1_p_y_);
      this->createBranch(tree, variablesTTBar->gen_top1_p_z_);
      this->createBranch(tree, variablesTTBar->gen_top1_p_e_);

      // antitop
      this->createBranch(tree, variablesTTBar->gen_top2_p_x_);
      this->createBranch(tree, variablesTTBar->gen_top2_p_y_);
      this->createBranch(tree, variablesTTBar->gen_top2_p_z_);
      this->createBranch(tree, variablesTTBar->gen_top2_p_e_);

      // lep1
      this->createBranch(tree, variablesTTBar->gen_lep1_p_x_);
      this->createBranch(tree, variablesTTBar->gen_lep1_p_y_);
      this->createBranch(tree, variablesTTBar->gen_lep1_p_z_);
      this->createBranch(tree, variablesTTBar->gen_lep1_p_e_);

      // lep2
      this->createBranch(tree, variablesTTBar->gen_lep2_p_x_);
      this->createBranch(tree, variablesTTBar->gen_lep2_p_y_);
      this->createBranch(tree, variablesTTBar->gen_lep2_p_z_);
      this->createBranch(tree, variablesTTBar->gen_lep2_p_e_);
    }
}



void TreeHandlerTTBarPlain::fillBranches(TTree* tree, const std::vector<VariablesBase*>& v_variables)
{
    
    for(const VariablesBase* variablesTmp : v_variables){
        const VariablesTTBarPlain* variablesTTBarTmp = dynamic_cast<const VariablesTTBarPlain*>(variablesTmp);
        if(!variablesTTBarTmp){
            std::cerr<<"ERROR in TreeHandlerTTBar::fillBranches()! variables are of wrong type, cannot typecast\n"
                     <<"...break\n"<<std::endl;
            exit(395);
        }
        
        *(dynamic_cast<VariablesTTBarPlain*>(variables_)) = *variablesTTBarTmp;
        
        tree->Fill();
    }
}
