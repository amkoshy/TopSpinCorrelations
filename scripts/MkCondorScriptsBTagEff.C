#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;
void MkCondorScriptsBTagEff(void){
  gSystem->mkdir("condor_jobs");
  gSystem->mkdir("condor_jobs/logs");
  ofstream BigJob("RunBTagEff.sh");
  
  string chan[3]={"ee","emu","mumu"}; 
  string syst1[29]={"Nominal", 
		    "JER_UP", "JER_DOWN", 
		    "JES_UP", "JES_DOWN",
		    "PU_UP", "PU_DOWN",
		    "TRIG_UP", "TRIG_DOWN",
		    "TRIG_ETA_UP", "TRIG_ETA_DOWN",
		    "LEPT_UP", "LEPT_DOWN",
		    "UNCLUSTERED_UP", "UNCLUSTERED_DOWN",
		    "KIN_UP", "KIN_DOWN",
		    "BTAG_UP", "BTAG_DOWN", "BTAG_LJET_UP", "BTAG_LJET_DOWN",
		    "BTAG_PT_UP", "BTAG_PT_DOWN", "BTAG_ETA_UP", "BTAG_ETA_DOWN",
		    "BTAG_LJET_PT_UP", "BTAG_LJET_PT_DOWN", "BTAG_LJET_ETA_UP", "BTAG_LJET_ETA_DOWN"};
  
  string syst2[52]={"JESAbsoluteStat_UP", "JESAbsoluteStat_DOWN",
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
		    "JESPileUpPtBB_UP", "JESPileUpPtBB_DOWN"};

  vector<string>systVar_arg;
  vector<string>output_arg;
  
  for(int i=0; i<3; i++){
    for(int j=0; j<29; j++){
      systVar_arg.push_back("-f ttbarsignalplustau.root -c "+chan[i]+" -s "+syst1[j]);
      output_arg.push_back("selectionRoot/BTagEff/"+syst1[j]+"/"+chan[i]);
    }
    for(int j=0; j<52; j++){
      systVar_arg.push_back("-f ttbarsignalplustau.root -c "+chan[i]+" -s "+syst2[j]);
      output_arg.push_back("selectionRoot/BTagEff/"+syst2[j]+"/"+chan[i]);
    }
    
    systVar_arg.push_back("-f ttbarsignalplustau.root -s PDF_ALPHAS_UP -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau.root -s PDF_ALPHAS_DOWN -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau.root -s MESCALE_UP -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau.root -s MESCALE_DOWN -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau.root -s MEFACSCALE_UP -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau.root -s MEFACSCALE_DOWN -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau.root -s MERENSCALE_UP -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau.root -s MERENSCALE_DOWN -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau.root -s BFRAG_UP -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau.root -s BFRAG_DOWN -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau.root -s BFRAG_CENTRAL -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau.root -s BFRAG_PETERSON -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau.root -s BSEMILEP_UP -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau.root -s BSEMILEP_DOWN -c "+chan[i]);
    systVar_arg.push_back("-f tau_psfsrscaleup.root -c "+chan[i]);
    systVar_arg.push_back("-f tau_psfsrscaledown.root -c "+chan[i]);
    systVar_arg.push_back("-f tau_psisrscaleup.root -c "+chan[i]);
    systVar_arg.push_back("-f tau_psisrscaledown.root -c "+chan[i]);
    systVar_arg.push_back("-f tau_175_massup.root -c "+chan[i]);
    systVar_arg.push_back("-f tau_169_massdown.root -c "+chan[i]);
    systVar_arg.push_back("-f tau_uetuneup.root -c "+chan[i]);
    systVar_arg.push_back("-f tau_uetunedown.root -c "+chan[i]);
    systVar_arg.push_back("-f tau_matchup.root -c "+chan[i]);
    systVar_arg.push_back("-f tau_matchdown.root -c "+chan[i]);
    systVar_arg.push_back("-f tau_powhegv2Herwig.root -c "+chan[i]);
    systVar_arg.push_back("-f tau_fromDilepton_amcatnlofxfx.root -c "+chan[i]);
    systVar_arg.push_back("-f tau_erdon.root -c "+chan[i]);
    systVar_arg.push_back("-f tau_erdonretune.root -c "+chan[i]);
    systVar_arg.push_back("-f tau_gluonmovetune.root -c "+chan[i]);
    
    
    output_arg.push_back("selectionRoot/BTagEff/PDF_ALPHAS_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/PDF_ALPHAS_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/MESCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/MESCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/MEFACSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/MEFACSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/MERENSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/MERENSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/BFRAG_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/BFRAG_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/BFRAG_CENTRAL/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/BFRAG_PETERSON/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/BSEMILEP_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/BSEMILEP_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/PSFSRSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/PSFSRSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/PSISRSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/PSISRSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/MASS_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/MASS_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/UETUNE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/UETUNE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/MATCH_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/MATCH_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/POWHEGV2HERWIG/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/AMCATNLOFXFX/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/ERDON/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/ERDONRETUNE/"+chan[i]);
    output_arg.push_back("selectionRoot/BTagEff/GLUONMOVETUNE/"+chan[i]);
  }
  
  BigJob << "#!/bin/bash" <<endl;
  BigJob << "tar -X excludefromtar.txt -zcvf CMSSW.tgz /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/" <<endl;
  for (int i=0; i<systVar_arg.size(); i++){
    
    char ShellName[600];
    sprintf(ShellName, "condor_jobs/condor_doit_%i.sh",i);
    
    char JDLName[600];
    sprintf(JDLName, "condor_jobs/AnaCondor_%i.jdl",i);
    
    ofstream jdl(JDLName);
    ofstream shel(ShellName);
    BigJob << "condor_submit " << JDLName << endl;
    
    //Make the JDL file:
    
    jdl << "Universe                = vanilla" << endl;
    jdl << "Environment             = CONDORJOBID=$(Process)" << endl;
    jdl << "notification            = Error" << endl;
    jdl << "X509UserProxy = /tmp/x509up_u587091" << endl;
    jdl << "should_transfer_files   = YES" << endl;
    jdl << "Transfer_input_files = condor_jobs/condor_doit_"<<i<<".sh" << endl;
    jdl << "Executable = /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<ShellName << endl;
    jdl << "when_to_transfer_output = ON_EXIT" << endl;
    jdl << "Output = /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/condor_jobs/logs/batch_"<<i<<"_$(cluster)_$(process).stdout" << endl;
    jdl << "Error  = /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/condor_jobs/logs/batch_"<<i<<"_$(cluster)_$(process).stderr" << endl;
    jdl << "Log    = /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/condor_jobs/logs/batch_"<<i<<"_$(cluster)_$(process).condor" << endl;
    jdl << "Queue 1" <<endl; 
    
    //Make the Shell script: 
    shel << "#!/bin/sh" << endl;
    shel << "source /cvmfs/cms.cern.ch/cmsset_default.sh" << endl;
    shel << "cp /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/CMSSW.tgz ." <<endl;
    shel << "tar -zxf CMSSW.tgz" <<endl;
    shel << "rm CMSSW.tgz" <<endl;
    shel << "ls -lhart"<< endl;
    shel << "cd home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src" <<endl;
    shel << "export SCRAM_ARCH=slc6_amd64_gcc530" <<endl;
    shel << "scram b ProjectRename" <<endl;
    shel << "eval `scramv1 runtime -sh`" <<endl;
    shel << "cd TopAnalysis/Configuration/analysis/diLeptonic" <<endl;
    shel << "./install/bin/load_Analysis "<<systVar_arg[i].c_str()<<endl;
    shel << "mkdir -p /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<output_arg[i].c_str()<<endl;
    shel << "cp "<<output_arg[i].c_str()<<"/ttbarsignalplustau.root /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<output_arg[i].c_str()<<"/ttbarsignalplustau.root" <<endl;
    shel << "rm -rf CMSSW_8_0_26_patch2" <<endl;
    
  }
}
