#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;
void MkCondorScriptsTheoryVariations(void){
  gSystem->mkdir("condor_jobs");
  gSystem->mkdir("condor_jobs/logs");
  ofstream BigJob("RunTheoryVariations.sh");
  
  //string mode="";
  string mode="-m spinCorr";
  string chan[3]={"ee","emu","mumu"}; 

  vector<string>sampleVar_arg;
  vector<string>output_arg;
  vector<string>tree_arg;
  
  for(int i=0; i<3; i++){

    /*sampleVar_arg.push_back("-f 175_massup "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f 169_massdown "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f psfsrscaleup "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f psfsrscaledown "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f psisrscaleup "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f psisrscaledown "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f uetuneup "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f uetunedown "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f matchup "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f matchdown "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f powhegv2Herwig "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f amcatnlofxfx "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f erdon.root "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f erdonretune "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f gluonmovetune "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau_175_massup "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau_169_massdown "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau_psfsrscaleup "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau_psfsrscaledown "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau_psisrscaleup "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau_psisrscaledown "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau_uetuneup "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau_uetunedown "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau_matchup "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau_matchdown "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau_powhegv2Herwig "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau_fromDilepton_amcatnlofxfx "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau_erdon.root "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau_erdonretune "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau_gluonmovetune "+mode+" -c "+chan[i]+" --bgviatau");*/
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s PDF_ALPHAS_UP "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s PDF_ALPHAS_DOWN "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s MESCALE_UP "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s MESCALE_DOWN "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s MEFACSCALE_UP "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s MEFACSCALE_DOWN "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s MERENSCALE_UP "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s MERENSCALE_DOWN "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s BFRAG_UP "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s BFRAG_DOWN "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s BFRAG_CENTRAL "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s BFRAG_PETERSON "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s BSEMILEP_UP "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s BSEMILEP_DOWN "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s PDF_ALPHAS_UP "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s PDF_ALPHAS_DOWN "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s MESCALE_UP "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s MESCALE_DOWN "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s MEFACSCALE_UP "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s MEFACSCALE_DOWN "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s MERENSCALE_UP "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s MERENSCALE_DOWN "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s BFRAG_UP "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s BFRAG_DOWN "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s BFRAG_CENTRAL "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s BFRAG_PETERSON "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s BSEMILEP_UP "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarsignalplustau.root -s BSEMILEP_DOWN "+mode+" -c "+chan[i]+" --bgviatau");
    sampleVar_arg.push_back("-f ttbarbg.root -s PDF_ALPHAS_UP "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarbg.root -s PDF_ALPHAS_DOWN "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarbg.root -s MESCALE_UP "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarbg.root -s MESCALE_DOWN "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarbg.root -s MEFACSCALE_UP "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarbg.root -s MEFACSCALE_DOWN "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarbg.root -s MERENSCALE_UP "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarbg.root -s MERENSCALE_DOWN "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarbg.root -s BFRAG_UP "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarbg.root -s BFRAG_DOWN "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarbg.root -s BFRAG_CENTRAL "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarbg.root -s BFRAG_PETERSON "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarbg.root -s BSEMILEP_UP "+mode+" -c "+chan[i]);
    sampleVar_arg.push_back("-f ttbarbg.root -s BSEMILEP_DOWN "+mode+" -c "+chan[i]);        

    /*output_arg.push_back("selectionRoot/MASS_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/MASS_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/PSFSRSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/PSFSRSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/PSISRSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/PSISRSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/UETUNE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/UETUNE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/MATCH_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/MATCH_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/POWHEGV2HERWIG/"+chan[i]);
    output_arg.push_back("selectionRoot/AMCATNLOFXFX/"+chan[i]);
    output_arg.push_back("selectionRoot/ERDON/"+chan[i]);
    output_arg.push_back("selectionRoot/ERDONRETUNE/"+chan[i]);
    output_arg.push_back("selectionRoot/GLUONMOVETUNE/"+chan[i]);
    output_arg.push_back("selectionRoot/MASS_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/MASS_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/PSFSRSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/PSFSRSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/PSISRSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/PSISRSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/UETUNE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/UETUNE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/MATCH_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/MATCH_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/POWHEGV2HERWIG/"+chan[i]);
    output_arg.push_back("selectionRoot/AMCATNLOFXFX/"+chan[i]);
    output_arg.push_back("selectionRoot/ERDON/"+chan[i]);
    output_arg.push_back("selectionRoot/ERDONRETUNE/"+chan[i]);
    output_arg.push_back("selectionRoot/GLUONMOVETUNE/"+chan[i]);*/
    output_arg.push_back("selectionRoot/PDF_ALPHAS_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/PDF_ALPHAS_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/MESCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/MESCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/MEFACSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/MEFACSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/MERENSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/MERENSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BFRAG_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BFRAG_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BFRAG_CENTRAL/"+chan[i]);
    output_arg.push_back("selectionRoot/BFRAG_PETERSON/"+chan[i]);
    output_arg.push_back("selectionRoot/BSEMILEP_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BSEMILEP_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/PDF_ALPHAS_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/PDF_ALPHAS_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/MESCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/MESCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/MEFACSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/MEFACSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/MERENSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/MERENSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BFRAG_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BFRAG_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BFRAG_CENTRAL/"+chan[i]);
    output_arg.push_back("selectionRoot/BFRAG_PETERSON/"+chan[i]);
    output_arg.push_back("selectionRoot/BSEMILEP_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BSEMILEP_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/PDF_ALPHAS_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/PDF_ALPHAS_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/MESCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/MESCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/MEFACSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/MEFACSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/MERENSCALE_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/MERENSCALE_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BFRAG_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BFRAG_DOWN/"+chan[i]);
    output_arg.push_back("selectionRoot/BFRAG_CENTRAL/"+chan[i]);
    output_arg.push_back("selectionRoot/BFRAG_PETERSON/"+chan[i]);
    output_arg.push_back("selectionRoot/BSEMILEP_UP/"+chan[i]);
    output_arg.push_back("selectionRoot/BSEMILEP_DOWN/"+chan[i]);        

    /*tree_arg.push_back("spinCorrInput/MASS_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MASS_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/PSFSRSCALE_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/PSFSRSCALE_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/PSISRSCALE_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/PSISRSCALE_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/UETUNE_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/UETUNE_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MATCH_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MATCH_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/POWHEGV2HERWIG/"+chan[i]);
    tree_arg.push_back("spinCorrInput/AMCATNLOFXFX/"+chan[i]);
    tree_arg.push_back("spinCorrInput/ERDON/"+chan[i]);
    tree_arg.push_back("spinCorrInput/ERDONRETUNE/"+chan[i]);
    tree_arg.push_back("spinCorrInput/GLUONMOVETUNE/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MASS_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MASS_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/PSFSRSCALE_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/PSFSRSCALE_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/PSISRSCALE_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/PSISRSCALE_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/UETUNE_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/UETUNE_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MATCH_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MATCH_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/POWHEGV2HERWIG/"+chan[i]);
    tree_arg.push_back("spinCorrInput/AMCATNLOFXFX/"+chan[i]);
    tree_arg.push_back("spinCorrInput/ERDON/"+chan[i]);
    tree_arg.push_back("spinCorrInput/ERDONRETUNE/"+chan[i]);
    tree_arg.push_back("spinCorrInput/GLUONMOVETUNE/"+chan[i]);*/
    tree_arg.push_back("spinCorrInput/PDF_ALPHAS_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/PDF_ALPHAS_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MESCALE_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MESCALE_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MEFACSCALE_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MEFACSCALE_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MERENSCALE_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MERENSCALE_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BFRAG_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BFRAG_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BFRAG_CENTRAL/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BFRAG_PETERSON/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BSEMILEP_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BSEMILEP_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/PDF_ALPHAS_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/PDF_ALPHAS_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MESCALE_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MESCALE_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MEFACSCALE_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MEFACSCALE_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MERENSCALE_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MERENSCALE_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BFRAG_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BFRAG_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BFRAG_CENTRAL/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BFRAG_PETERSON/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BSEMILEP_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BSEMILEP_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/PDF_ALPHAS_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/PDF_ALPHAS_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MESCALE_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MESCALE_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MEFACSCALE_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MEFACSCALE_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MERENSCALE_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/MERENSCALE_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BFRAG_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BFRAG_DOWN/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BFRAG_CENTRAL/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BFRAG_PETERSON/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BSEMILEP_UP/"+chan[i]);
    tree_arg.push_back("spinCorrInput/BSEMILEP_DOWN/"+chan[i]);        
  }
  
  BigJob << "#!/bin/bash" <<endl;
  BigJob << "voms-proxy-init -voms cms" <<endl;
  //BigJob << "tar -X excludefromtar.txt -zcvf CMSSW.tgz /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/" <<endl;
  for (int i=0; i<sampleVar_arg.size(); i++){

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
    shel << "./install/bin/load_Analysis "<<sampleVar_arg[i].c_str()<<endl;
    shel << "mkdir -p /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<output_arg[i].c_str()<<endl;
    shel << "cp "<<output_arg[i].c_str()<<"/*.root /home/akhatiwa/depot_cms/TTBarSpinCor_fullInstallation_DESYRepo/CMSSW_8_0_26_patch2/src/TopAnalysis/Configuration/analysis/diLeptonic/"<<output_arg[i].c_str()<<"/." <<endl;
    shel << "gfal-mkdir gsiftp://cms-gridftp.rcac.purdue.edu/store/group/local/cmstop/ajeeta/"<<tree_arg[i].c_str()<<endl;
    //create a list of files created inside the spinCorrInput folder
    shel << "ls "<<tree_arg[i].c_str()<<"/*.root > filescreated.txt" <<endl;
    shel << "file=filescreated.txt" <<endl;
    shel << "mypwd=$PWD" <<endl;
    shel << "while IFS=\'\' read -r line || [[ -n \"$line\" ]]; do" <<endl;
    shel << "    gfal-copy file:////$mypwd/$line gsiftp://cms-gridftp.rcac.purdue.edu/store/group/local/cmstop/ajeeta/$line " <<endl;
    shel << "    echo \"gfal-copy file:////$mypwd/$line gsiftp://cms-gridftp.rcac.purdue.edu/store/group/local/cmstop/ajeeta/$line\" " <<endl;
    shel << "done < \"$file\" " <<endl;
    shel << "rm -rf CMSSW_8_0_26_patch2" <<endl;
    
  }
}
