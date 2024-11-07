#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;
void MkCondorScriptsBTagEff_2017(void){
  const int nChan = 3;
  const int nSyst1 = 30;
  const int nSyst2 = 52;

  gSystem->mkdir("condor_jobs");
  gSystem->mkdir("condor_jobs/logs");
  ofstream BigJob("RunBTagEff_2017.sh");

  const char* env_cmsswbase = std::getenv("CMSSW_BASE");
  
  string chan[nChan]={"ee","emu","mumu"}; 
  string syst1[nSyst1]={
                       "Nominal", 
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
		       "BTAG_LJET_PT_UP", "BTAG_LJET_PT_DOWN", "BTAG_LJET_ETA_UP", "BTAG_LJET_ETA_DOWN",    "TOP_PT"
		    };
  
  string syst2[nSyst2]={"JESAbsoluteStat_UP", "JESAbsoluteStat_DOWN",
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
  
  for(int i=0; i<nChan; i++){
    for(int j=0; j<nSyst1; j++){
      systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -c "+chan[i]+" -s "+syst1[j]);
      output_arg.push_back("selectionRoot_2017/BTagEff_2017/"+syst1[j]+"/"+chan[i]);
    }
    /*    
    for(int j=0; j<nSyst2; j++){
      systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -c "+chan[i]+" -s "+syst2[j]);
      output_arg.push_back("selectionRoot_2017/BTagEff_2017/"+syst2[j]+"/"+chan[i]);
    }
    */
        
    systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -s PDF_ALPHAS_UP -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -s PDF_ALPHAS_DOWN -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -s MESCALE_UP -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -s MESCALE_DOWN -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -s MEFACSCALE_UP -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -s MEFACSCALE_DOWN -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -s MERENSCALE_UP -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -s MERENSCALE_DOWN -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -s BFRAG_UP -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -s BFRAG_DOWN -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -s BFRAG_CENTRAL -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -s BFRAG_PETERSON -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -s BSEMILEP_UP -c "+chan[i]);
    systVar_arg.push_back("-f ttbarsignalplustau_fromDilepton.root -s BSEMILEP_DOWN -c "+chan[i]);
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
    
    
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/PDF_ALPHAS_UP/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/PDF_ALPHAS_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/MESCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/MESCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/MEFACSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/MEFACSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/MERENSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/MERENSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/BFRAG_UP/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/BFRAG_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/BFRAG_CENTRAL/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/BFRAG_PETERSON/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/BSEMILEP_UP/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/BSEMILEP_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/PSFSRSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/PSFSRSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/PSISRSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/PSISRSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/MASS_UP/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/MASS_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/UETUNE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/UETUNE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/MATCH_UP/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/MATCH_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/POWHEGV2HERWIG/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/AMCATNLOFXFX/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/ERDON/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/ERDONRETUNE/"+chan[i]);
    output_arg.push_back("selectionRoot_2017/BTagEff_2017/GLUONMOVETUNE/"+chan[i]);
    
    }
  
  BigJob << "#!/bin/bash" <<endl;
  BigJob << "voms-proxy-init -voms cms -valid 999:00:00" <<endl;
  BigJob << "tar -X excludefromtar.txt -zcvf CMSSW.tgz ../../../../../../CMSSW_10_2_22" <<endl;
  BigJob << "export X509_USER_CERT=`voms-proxy-info -path`" <<endl;
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
    jdl << "Requirements            = (Distro == \"RHEL7\")" << endl;
    jdl << "request_cpus            = 1.0" << endl;
    jdl << "request_disk            = 9000000KB" << endl;
    jdl << "request_memory          = 2000MB" << endl;
    //if proxy is not specified it will use environment variable $X509_USER_CERT                                                                                                                        
    jdl << "X509UserProxy = /tmp/x509up_u607473" << endl;                                                                                                                                             
    jdl << "should_transfer_files   = YES" << endl;
    jdl << "Transfer_input_files = condor_jobs/condor_doit_"<<i<<".sh" << endl;
    jdl << "Executable = "<<env_cmsswbase<<"/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<ShellName << endl;
    jdl << "when_to_transfer_output = ON_EXIT" << endl;
    jdl << "Output = "<<env_cmsswbase<<"/src/TopAnalysis/Configuration/analysis/diLeptonic/condor_jobs/logs/batch_"<<i<<"_$(cluster)_$(process).stdout" << endl;
    jdl << "Error  = "<<env_cmsswbase<<"/src/TopAnalysis/Configuration/analysis/diLeptonic/condor_jobs/logs/batch_"<<i<<"_$(cluster)_$(process).stderr" << endl;
    jdl << "Log    = "<<env_cmsswbase<<"/src/TopAnalysis/Configuration/analysis/diLeptonic/condor_jobs/logs/batch_"<<i<<"_$(cluster)_$(process).condor" << endl;
    jdl << "Queue 1" <<endl;
    
    //Make the Shell script:                                                                                                                                                                            
    shel << "#!/bin/sh" << endl;
    shel << "pwd" << endl;
    shel << "source /cvmfs/cms.cern.ch/cmsset_default.sh" << endl;
    shel << "cp "<<env_cmsswbase<<"/src/TopAnalysis/Configuration/analysis/diLeptonic/CMSSW.tgz ." <<endl;
    shel << "tar -zxf CMSSW.tgz" <<endl;
    shel << "rm CMSSW.tgz" <<endl;
    shel << "ls -lhart"<< endl;
    shel << "cd CMSSW_10_2_22/src" <<endl;
    shel << "export SCRAM_ARCH=slc7_amd64_gcc700" <<endl;
    shel << "scram b ProjectRename" <<endl;
    shel << "eval `scramv1 runtime -sh`" <<endl;
    shel << "cd TopAnalysis/Configuration/analysis/diLeptonic" <<endl;
    shel << "./install/bin/load_Analysis "<<systVar_arg[i].c_str()<<endl;
    shel << "mkdir -p "<<env_cmsswbase<<"/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<output_arg[i].c_str()<<endl;
    shel << "cp "<<output_arg[i].c_str()<<"/ttbarsignalplustau_fromDilepton.root "<<env_cmsswbase<<"/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<output_arg[i].c_str()<<"/ttbarsignalplustau_fromDilepton.root" <<endl;
    shel << "rm -rf CMSSW_10_2_22" <<endl;    
  }
}
