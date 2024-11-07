import os
import ROOT
import uproot
import numpy as np

def get_lumi_from_samplename(samplename) :
    # 2016
    if   "2016ULpreVFP"  in samplename : lumi = 19668.0
    elif "2016ULpostVFP" in samplename : lumi = 16778.0
    elif "2016UL" in samplename : lumi = 36310.0
    
    # 2017
    elif "2017UL" in samplename : lumi = 41480.0
    
    # 2018
    elif "2018UL" in samplename : lumi = 59830.0
    
    # full Run2
    elif "fullRun2UL" in samplename : lumi = 137650.0
    
    return lumi


def get_dySF_from_samplename(samplename,channel) :
    dySF = 1.0 

    if  "2018" in samplename:
        if    channel == "ee"       : dySF = 1.12717
        elif  channel == "emu"      : dySF = 1.10941
        elif  channel == "mumu"     : dySF = 1.09194
        elif  channel == "combined" : dySF = 1.1027
    
    elif "2017" in samplename :
        if    channel == "ee"       : dySF = 1.17378
        elif  channel == "emu"      : dySF = 1.16718
        elif  channel == "mumu"     : dySF = 1.16061
        elif  channel == "combined" : dySF = 1.16462

    elif "2016postVFP" in samplename :
        if    channel == "ee"       : dySF = 1.19939
        elif  channel == "emu"      : dySF = 1.1901
        elif  channel == "mumu"     : dySF = 1.18087
        elif  channel == "combined" : dySF = 1.18599

    elif "2016preVFP" in samplename :
        if    channel == "ee"       : dySF = 1.1809
        elif  channel == "emu"      : dySF = 1.14192
        elif  channel == "mumu"     : dySF = 1.10422
        elif  channel == "combined" : dySF = 1.12537

    elif "2016" in samplename and not "VFP" in samplename: 
        if    channel == "ee"       : dySF = 1.18945
        elif  channel == "emu"      : dySF = 1.16376
        elif  channel == "mumu"     : dySF = 1.13863
        elif  channel == "combined" : dySF = 1.15266
            
    elif "fullRun2UL" in samplename :
        if    channel == "ee"       : dySF = 1.15707
        elif  channel == "emu"      : dySF = 1.14002
        elif  channel == "mumu"     : dySF = 1.12321
        elif  channel == "combined" : dySF = 1.13327
    
    return dySF 

def get_cross_section_from_samplename(samplename, channel) :
    topxsec = 830.91
    
    # w/o tau signal files
    if   "ee_ttbarsignalplustau_fromDilepton"   in samplename : xsection = topxsec * 0.10706 * 0.964261576
    elif "emu_ttbarsignalplustau_fromDilepton"  in samplename : xsection = topxsec * 0.10706 * 0.957058875
    elif "mumu_ttbarsignalplustau_fromDilepton" in samplename : xsection = topxsec * 0.10706 * 0.949909976
    
    # w tau signal files
    elif "ee_ttbarsignalviatau_fromDilepton"   in samplename : xsection = topxsec * 0.10706 * 1.029827957
    elif "emu_ttbarsignalviatau_fromDilepton"  in samplename : xsection = topxsec * 0.10706 * 1.026209047
    elif "mumu_ttbarsignalviatau_fromDilepton" in samplename : xsection = topxsec * 0.10706 * 1.022670477
    
    # backgrounds
    
    # ttbar backgrounds
    elif "bg_fromDilepton" in samplename : xsection = topxsec * 0.10706
    elif "fromLjets"       in samplename : xsection = topxsec * 0.44113
    elif "fromHadronic"    in samplename : xsection = topxsec * 0.45441
    
    # elif ("ttbar"  in samplename) and not("ttbarW" in samplename) and not("ttbarW" in samplename) : xsection = topxsec 
    
    # Single top
    elif ("single" in samplename)    and ("tw" in samplename) : xsection = 35.85 * 0.54559
    elif ("single" in samplename)    and ("_s" in samplename) : xsection = 10.32
    elif ("singletop" in samplename) and ("_t" in samplename) : xsection = 136.02

    # Single antitop
    elif ("singleantitop" in samplename) and ("_t" in samplename) : xsection = 80.95
    
    # VV
    elif "ww" in samplename : xsection = 118.7
    elif "wz" in samplename : xsection = 47.13
    elif "zz" in samplename : xsection = 16.523
    
    # Drell-Yan
    elif "1050" in samplename               : xsection = 22635.1 * get_dySF_from_samplename(samplename, channel)
    elif "50inf_amcatnlofxfx" in samplename : xsection = 0.5 * 3.*2075.14 * get_dySF_from_samplename(samplename, channel)
    elif "50inf_madgraphmlm"  in samplename : xsection = 0.5 * 3.*2075.14 * get_dySF_from_samplename(samplename, channel)
    
    # Smaller backgrounds
    elif "wtolnu"          in samplename : xsection = 61526.7
    elif "ttbarWjetstolnu" in samplename : xsection = 0.2043
    elif "ttbarWjetstoqq"  in samplename : xsection = 0.4062
    elif "ttbarZtollnunu"  in samplename : xsection = 0.2529
    elif "ttbarZtoqq"      in samplename : xsection = 0.5297
    else : xsection = topxsec
    
    return xsection

def get_lumi_weight(root_fileptr, lumi, xsection) :
    nevents = root_fileptr.Get('hNrOfEvts').GetBinContent(1)
    lumi_wt = (lumi * xsection) / nevents
    return lumi_wt


# Open txt file ready root fileptrs
with open("TUnfoldFileList_Nominal_combined.txt") as ipfile :
  allfiles = ipfile.readlines()

for file in allfiles :
  file = allfiles.strip('\n')

# // ***************
# // FillHistoArrays 
# // ***************

# // Get the histograms stored in a file and store the pointer in a TObjArray for further use
# // also pass a string here so that that string can be appended in front of the name

def FillHistoArrays(channel) :

    # TString bgfracplotsfile = "BGfracplots/"+CHAN+"/";
    # gSystem->mkdir(bgfracplotsfile.Data(), true);

    fDataArray   = []
    fRecArray    = []
    fGenArray    = []
    fVisGenArray = []
    fMigMatArray = []
    fBgArray     = []
    fBgSamplesArray = []
    fMigMatSysArray = []
    fTrueDataArray  = []
  
    for i in range(len(fVariableNames)) :

    # ****
    # DATA
    # ****
    
    # No requirement in number of background files because should have option to not use any background at all
    
    channels = []
    channels.append("ee");
    channels.append("emu");
    channels.append("mumu");
    channels.append("combined");

    # // Amandeep : kMaxNumDataFiles is set to 100 in the .h file
    # TH1D* htempdata[kMaxNumDataFiles+1];
    # TH1D* hdatabgsub_channels[4];

    data_yields          = [0,0,0,0]
    bkgsubdata_yields    = [0,0,0,0]
    hdatabgsub_channels  = []
    htempdata            = []

    if fDataFiles :
        for j in range(len(channels)) :
            hdatabgsub_channels[j]= fDataFiles.At(0).Get("hreco_" + fVariableNames[i]).Clone(Form("%sDataBgSub_%s",fVariableNames[i],channels[j]))
            hdatabgsub_channels[j].Reset()

        htempdata[0]= (TH1D*)((TH1D*)((TFile*)fDataFiles->At(0))->Get(Form("hreco_%s",fVariableNames[i].Data())))->Clone(Form("%sData",fVariableNames[i].Data()));
        htempdata[0].Reset()

        for(int j=0; j < ndatafiles; j++)
        
        htempdata[j+1]= (TH1D*) ((TFile*)fDataFiles->At(j))->Get(Form("hreco_%s",fVariableNames[i].Data()));
        std::cout << fDataFiles->At(j)->GetName() << " Integral: " << htempdata[j+1]->Integral() << std::endl;
        if( TString(fDataFiles->At(j)->GetName()).Contains("histosTUnfold_ee_") )   {hdatabgsub_channels[0]->Add(htempdata[j+1]); data_yields[0]+=htempdata[j+1]->Integral();}
        if( TString(fDataFiles->At(j)->GetName()).Contains("histosTUnfold_emu_") )  {hdatabgsub_channels[1]->Add(htempdata[j+1]); data_yields[1]+=htempdata[j+1]->Integral();}
        if( TString(fDataFiles->At(j)->GetName()).Contains("histosTUnfold_mumu_") ) {hdatabgsub_channels[2]->Add(htempdata[j+1]); data_yields[2]+=htempdata[j+1]->Integral();}
        hdatabgsub_channels[3]->Add(htempdata[j+1]);
        data_yields[3]+=htempdata[j+1]->Integral();

        htempdata[0]->Add(htempdata[j+1]);

        std::cout << "data yields ee emu mumu combined: " << data_yields[0] << " "<< data_yields[1] << " " << data_yields[2]<<" "<<data_yields[3]<<std::endl;
        fDataArray->AddLast(htempdata[0])



    # // ***************
    # // MC signal files
    # // ***************

    # // Reco
    # TH1D* htempreco[nMCfiles+1];
    # // Gen
    # TH1D* htempgen[nMCfiles+1];
    # // VisGen
    # TH1D* htempvisgen[nMCfiles+1];
    # // MigMat  
    # TH2D* htemprecoVsgen[nMCfiles+1];
    
    TH1D* hrecoplustau_channels[4];
    double sig_yields[4] = {0,0,0,0};
    
    if(fMCFiles) 
    {

      htempreco[0]  = (TH1D*)((TH1D*)((TFile*)fMCFiles->At(0))->Get(Form("hreco_%s",fVariableNames[i].Data())))->Clone(Form("%sReco",fVariableNames[i].Data()));
      htempreco[0]->Reset();
      // std::cout<<fMCFiles->At(0)->GetName()<<" Integral: "<<htempreco[0]->Integral()<<" xsec: "<<fXSectionMC.at(0)<<std::endl;

      for (int j = 0; j < 4; ++j)
      {
        hrecoplustau_channels[j]= (TH1D*)((TH1D*)((TFile*)fMCFiles->At(0))->Get(Form("hreco_%s",fVariableNames[i].Data())))->Clone(Form("%sReco_%s",fVariableNames[i].Data(),channels.at(j).Data()));
        hrecoplustau_channels[j]->Reset();
      }

      htempgen[0]    = (TH1D*)((TH1D*)((TFile*)fMCFiles->At(0))->Get(Form("hgen_%s",fVariableNames[i].Data())))->Clone(Form("%sGen",fVariableNames[i].Data()));
      htempgen[0]->Reset();

      htempvisgen[0] = (TH1D*)((TH1D*)((TFile*)fMCFiles->At(0))->Get(Form("hvisgen_%s",fVariableNames[i].Data())))->Clone(Form("%sVisGen",fVariableNames[i].Data()));
      htempvisgen[0]->Reset();

      for (int j=0; j < nMCfiles; j++)
      {
        TString MCFileName = fMCFiles->At(j)->GetName();

        // Amandeep : Check here maybe ?
        htempreco[j+1] = (TH1D*)((TFile*)fMCFiles->At(j))->Get(Form("hreco_%s",fVariableNames[i].Data()));
        htempreco[j+1]->Scale(getLumiWeight( (TFile*) fMCFiles->At(j), fXSectionMC.at(j), fLumiMC.at(j)));
        std::cout << MCFileName << " Integral: " << htempreco[j+1]->Integral() << " xsec: " << fXSectionMC.at(j) << std::endl;
        
        if( MCFileName.Contains("histosTUnfold_ee_") )   {hrecoplustau_channels[0]->Add(htempreco[j+1]); sig_yields[0]+=htempreco[j+1]->Integral();}
        if( MCFileName.Contains("histosTUnfold_emu_") )  {hrecoplustau_channels[1]->Add(htempreco[j+1]); sig_yields[1]+=htempreco[j+1]->Integral();}
        if( MCFileName.Contains("histosTUnfold_mumu_") ) {hrecoplustau_channels[2]->Add(htempreco[j+1]); sig_yields[2]+=htempreco[j+1]->Integral();}
        
        hrecoplustau_channels[3]->Add(htempreco[j+1]);
        sig_yields[3]+=htempreco[j+1]->Integral();

        htempreco[0]->Add(htempreco[j+1]);

        htempgen[j+1]    = (TH1D*)((TFile*)fMCFiles->At(j))->Get(Form("hgen_%s",fVariableNames[i].Data()));
        htempgen[j+1]->Scale(getLumiWeight( (TFile*) fMCFiles->At(j), fXSectionMC.at(j), fLumiMC.at(j)));
        htempgen[0]->Add(htempgen[j+1]);

        htempvisgen[j+1] = (TH1D*)((TFile*)fMCFiles->At(j))->Get(Form("hvisgen_%s",fVariableNames[i].Data()));
        htempvisgen[j+1]->Scale(getLumiWeight( (TFile*) fMCFiles->At(j), fXSectionMC.at(j), fLumiMC.at(j)));
        htempvisgen[0]->Add(htempvisgen[j+1]);
      }

      std::cout << "sig yields ee emu mumu combined: " << sig_yields[0] << " " << sig_yields[1] << " " << sig_yields[2] << " " << sig_yields[3] << std::endl;

      fRecArray->AddLast(htempreco[0]);
      fGenArray->AddLast(htempgen[0]);
      fVisGenArray->AddLast(htempvisgen[0]);
    }

    // ***********
    // BACKGROUNDS
    // ***********

    // Amandeep : Maybe this is even right ?
    // kMaxNumBgFiles is defined in the .h to be 100 

    TH1D* htempbg[kMaxNumBgFiles + 1];
    TH1D* htaubg_channels[4];
    TH1D* hMCtaubg_channels[4];
    TH1D* htWfrac_channels[4];
    TH1D* htaufrac_channels[4];

    double bkg_yields[4] = {0,0,0,0};
    double tau_yields[4] = {0,0,0,0};
    
    // Amandeep : commenting since we don't have taus as bgs
    // double tau_SF[4]={0,0,0,0};
    // End

    double channel_SF[4] = {1,1,1,1};
    double taufrac_SF[4] = {1,1,1,1};

    if(fBgFiles) {
      htempbg[0]= (TH1D*)((TH1D*)((TFile*)fBgFiles->At(0))->Get(Form("hreco_%s",fVariableNames[i].Data())))->Clone(Form("%sRecoBg",fVariableNames[i].Data()));
      htempbg[0]->Reset();
      //std::cout<<fBgFiles->At(0)->GetName()<<" Integral: "<<htempbg[0]->Integral()<<" xsec: "<<fXSectionBg.at(0)<<std::endl;

      for (int j = 0; j < 4; ++j)
      {
        htWfrac_channels[j]= (TH1D*)((TH1D*)((TFile*)fBgFiles->At(0))->Get(Form("hreco_%s",fVariableNames[i].Data())))->Clone(Form("%stWfrac_%s",fVariableNames[i].Data(),channels.at(j).Data()));
        htWfrac_channels[j]->Reset();
      }
      
      std::vector<TString> taubgname;
      std::vector<int> itau;
      //std::vector<TH1D*> tauhist;

      for (int j=0; j < nbgfiles; j++)
      {
        TString BgFileName = fBgFiles->At(j)->GetName();

        htempbg[j+1]       = (TH1D*)((TH1D*)((TFile*) fBgFiles->At(j))->Get(Form("hreco_%s",fVariableNames[i].Data())))->Clone(Form("%sRecoBg_%s",fVariableNames[i].Data(),BgFileName.Data()));
        double lumiweight  = getLumiWeight(((TFile*)fBgFiles->At(j)), fXSectionBg[j],  fLumiBg[j]);
        
        htempbg[j+1]->Scale(lumiweight);

        double bkgintegral, bkgerror, bkgnentries, bkgneffentries;
        bkgintegral    = htempbg[j+1]->IntegralAndError(1, htempbg[j+1]->GetNbinsX(), bkgerror);
        bkgnentries    = htempbg[j+1]->GetEntries();
        bkgneffentries = htempbg[j+1]->GetEffectiveEntries();

        std::cout<< BgFileName << " Integral: " << bkgintegral << " +/- " << bkgerror << " xsec: " << fXSectionBg.at(j) << " lumiweight: " << lumiweight << " bkgnentries: " << bkgnentries << " bkgneffentries: " << bkgneffentries << std::endl;
        
        if( BgFileName.Contains("histosTUnfold_ee_") )   {hdatabgsub_channels[0]->Add(htempbg[j+1],-1.); bkg_yields[0] += bkgintegral;}
        if( BgFileName.Contains("histosTUnfold_emu_") )  {hdatabgsub_channels[1]->Add(htempbg[j+1],-1.); bkg_yields[1] += bkgintegral;}
        if( BgFileName.Contains("histosTUnfold_mumu_") ) {hdatabgsub_channels[2]->Add(htempbg[j+1],-1.); bkg_yields[2] += bkgintegral;}
        
        hdatabgsub_channels[3]->Add(htempbg[j+1],-1.);
        bkg_yields[3] += bkgintegral;
        htempbg[0]->Add(htempbg[j+1]);

        if( BgFileName.Contains("single") ) {
          if( BgFileName.Contains("histosTUnfold_ee_") )   {htWfrac_channels[0]->Add(htempbg[j+1]);}
          if( BgFileName.Contains("histosTUnfold_emu_") )  {htWfrac_channels[1]->Add(htempbg[j+1]);}
          if( BgFileName.Contains("histosTUnfold_mumu_") ) {htWfrac_channels[2]->Add(htempbg[j+1]);}
          htWfrac_channels[3]->Add(htempbg[j+1]);
        }

      }

      for (int j = 0; j < 4; ++j)
      {
        for (int i_bin = 1; i_bin <= hdatabgsub_channels[j]->GetNbinsX(); ++i_bin)
        {
          if(hdatabgsub_channels[j]->GetBinContent( i_bin )<0) {
            std::cout<<"belowzerobin: "<<i_bin<<" "<<hdatabgsub_channels[j]->GetBinContent( i_bin )<<std::endl;
            hdatabgsub_channels[j]->SetBinContent( i_bin, 0 );
          }
        }
      }


        std::cout<<"bkg yields ee emu mumu combined: "<<bkg_yields[0]<<" "<<bkg_yields[1]<<" "<<bkg_yields[2]<<" "<<bkg_yields[3]<<std::endl;
        std::cout<<"tau yields ee emu mumu combined: "<<tau_yields[0]<<" "<<tau_yields[1]<<" "<<tau_yields[2]<<" "<<tau_yields[3]<<std::endl;
        // std::cout<<"tau SF ee emu mumu combined: "<<tau_SF[0]<<" "<<tau_SF[1]<<" "<<tau_SF[2]<<" "<<tau_SF[3]<<std::endl;


        for (unsigned int i_file = 0; i_file < itau.size(); ++i_file)
        {
          // Amandeep : Revisit this logic
          int i_chan = 3;
        
          if( taubgname.at(i_file).Contains("histosTUnfold_ee_") )   i_chan = 0;
          if( taubgname.at(i_file).Contains("histosTUnfold_emu_") )  i_chan = 1;
          if( taubgname.at(i_file).Contains("histosTUnfold_mumu_") ) i_chan = 2;

          // std::cout<<"Scaling "<< taubgname.at(i_file) <<" by "<<tau_SF[i_chan]<<" :"<<std::endl;
          // htempbg[ itau.at(i_file) ]->Scale( tau_SF[i_chan] );

          std::cout<< taubgname.at(i_file) << " Integral: " << htempbg[ itau.at(i_file) ]->Integral() << " itau: " << itau.at(i_file) << std::endl;
          bkg_yields[i_chan] += htempbg[ itau.at(i_file) ]->Integral();

          if(i_chan!=3) bkg_yields[3]+=htempbg[ itau.at(i_file) ]->Integral();
          htempbg[0]->Add(htempbg[ itau.at(i_file) ]);

          // htempbg[ itau.at(i_file) ]->Print("all");
          std::cout << "Integral: "<< htempbg[ itau.at(i_file)]->Integral() << std::endl;
        }


        for (int i_chan = 0; i_chan < 4; ++i_chan)
        {
          // tau_yields[i_chan]*=tau_SF[i_chan];
          bkgsubdata_yields[i_chan] = data_yields[i_chan] - bkg_yields[i_chan];
        }

        std::cout<<"final tau yields ee emu mumu combined: "<<tau_yields[0]<<" "<<tau_yields[1]<<" "<<tau_yields[2]<<" "<<tau_yields[3]<<std::endl;
        std::cout<<"final bkg yields ee emu mumu combined: "<<bkg_yields[0]<<" "<<bkg_yields[1]<<" "<<bkg_yields[2]<<" "<<bkg_yields[3]<<std::endl;
        std::cout<<"bkgsubdata yields ee emu mumu combined: "<<bkgsubdata_yields[0]<<" "<<bkgsubdata_yields[1]<<" "<<bkgsubdata_yields[2]<<" "<<bkgsubdata_yields[3]<<std::endl;
      }

      // Amandeep : nMCfiles was hardcoded here | Replace with channel.compare('combined') ?
      if( CHAN.CompareTo("combined") == 0 ) 
      { 
        CalculateChannelFractions(bkgsubdata_yields, sig_yields, channel_SF);
      }

      fBgArray->AddLast(htempbg[0] );

      for(int j=0;j<nbgfiles;j++)
      {
        fBgSamplesArray->AddLast(htempbg[j+1] );
      }

    bool doCorrectChannelFractions = kTRUE;

    if(fMCFiles) 
    {
      htemprecoVsgen[0] = (TH2D*)((TH2D*)((TFile*)fMCFiles->At(0))->Get(Form("hrecoVsgen_%s",fVariableNames[i].Data())))->Clone(Form("%sRespMat",fVariableNames[i].Data()));
      htemprecoVsgen[0]->Reset();

      for(int j=0; j < nMCfiles; j++)
      {
        TString MCFileName  = fMCFiles->At(j)->GetName();

        htemprecoVsgen[j+1] = (TH2D*)((TFile*)fMCFiles->At(j))->Get(Form("hrecoVsgen_%s",fVariableNames[i].Data()));
        htemprecoVsgen[j+1]->Scale(getLumiWeight( (TFile*) fMCFiles->At(j), fXSectionMC.at(j), fLumiMC.at(j)));

        // Amandeep : nMCfiles was hardcoded here | Replace with channel.compare('combined')
        if( doCorrectChannelFractions && (!fSelfConsistency) && CHAN.CompareTo("combined") == 0 ) {
          CorrectChannelFractions(htemprecoVsgen[j+1], MCFileName, channel_SF);
        }
        htemprecoVsgen[0]->Add(htemprecoVsgen[j+1]);
      }
      fMigMatArray->AddLast(htemprecoVsgen[0]);
    }

    // ***********
    // SYSTEMATICS
    // ***********
    
    // Same for systematics matrix

    if(fsetSyst)
    {
      // Migration Matrix 
      // Amandeep : implementing this a vector breaks 
      // Original :
      TH2D *htemp4[kMaxNumSysFiles];
      // Modified :
      // TH2D *htemp4[nsysfiles];

      // Over the systematics 
      for(int j=0; j < nsysfiles; j++){
        std::cout << "Filling histograms for systematics" << std::endl;

        // Amandeep : nMCfiles was hardcoded here | Replace with channel.compare('combined')
        if( CHAN.CompareTo("combined") == 0 ) 
        {
          // re-evaluate the SFs for each systematic variation
          double syst_yields[4] = {0,0,0,0};

          for(int k=0; k < nMCfiles; k++)
          {
            TString SysFileName = fSysFiles->At(j*nMCfiles+k)->GetName();
            TH1D* hrecosysttemp = (TH1D*) ((TFile*)fSysFiles->At(j*nMCfiles+k))->Get(Form("hreco_%s",fVariableNames[i].Data()));
            
            hrecosysttemp->Scale(getLumiWeight( (TFile*) fSysFiles->At(j*nMCfiles+k), fXSectionSys.at(j*nMCfiles+k), fLumiSys.at(j*nMCfiles+k)));

            if( SysFileName.Contains("histosTUnfold_ee_") )   syst_yields[0]+=hrecosysttemp->Integral();
            if( SysFileName.Contains("histosTUnfold_emu_") )  syst_yields[1]+=hrecosysttemp->Integral();
            if( SysFileName.Contains("histosTUnfold_mumu_") ) syst_yields[2]+=hrecosysttemp->Integral();
            syst_yields[3]+=hrecosysttemp->Integral();
          }

          if(fSelfConsistency) CalculateChannelFractions(sig_yields, syst_yields, channel_SF);
          else CalculateChannelFractions(bkgsubdata_yields, syst_yields, channel_SF);

        }

        htemp4[j] = (TH2D*)((TH2D*) ((TFile*)fSysFiles->At(j*nMCfiles))->Get(Form("hrecoVsgen_%s",fVariableNames[i].Data())))->Clone(Form("%sRespMatSys_%s",fVariableNames[i].Data(), fVectorOfValidSystematics[j].c_str()));
        htemp4[j]->Reset();

        // Over the channels/ files ??
        for(int k=0; k < nMCfiles; k++) {
          TString SysFileName = fSysFiles->At(j*nMCfiles+k)->GetName();          
          
          // Amandeep : Adding this here to debug
          std::cout << i  << " , " << j << " , " << k << std::endl;
          std::cout << "Sysfilename :: " << SysFileName << ", " << nMCfiles  << " , " << nsysfiles << " , " << j*nMCfiles+k << std::endl;
          // End         
          
          TH2D* htemp4temp    = (TH2D*) ((TFile*)fSysFiles->At(j*nMCfiles+k))->Get(Form("hrecoVsgen_%s",fVariableNames[i].Data()));
          htemp4temp->Scale(getLumiWeight( (TFile*) fSysFiles->At(j*nMCfiles+k), fXSectionSys.at(j*nMCfiles+k), fLumiSys.at(j*nMCfiles+k)));
          
          // Amandeep : nMCfiles was hardcoded here | Replace with channel.compare('combined')
          if( doCorrectChannelFractions && CHAN.CompareTo("combined") == 0 ) {
            // Not sure about the math here, need to follow the logic fully
            CorrectChannelFractions(htemp4temp, SysFileName, channel_SF);
          }
          // This adds all channel files of a particular syst to the htemp
          htemp4[j]->Add(htemp4temp);
        }
        fMigMatSysArray->AddLast(htemp4[j]); //consists of more entries than fMigMatArray (nsysfiles*nvariables)
      }  
    }

    if(fClosureTest){
     TH1D* htemptruedata[kMaxNumDataFiles];
     htemptruedata[0]= (TH1D*)((TH1D*)((TFile*)fDataFiles->At(0))->Get(Form("hgen_%s",fVariableNames[i].Data())))->Clone(Form("%sTrueData",fVariableNames[i].Data()));
     for(int j = 1; j < ndatafiles; j++)
     {
       htemptruedata[j]= (TH1D*) ((TFile*)fDataFiles->At(j))->Get(Form("hgen_%s",fVariableNames[i].Data()));
       htemptruedata[0]->Add(htemptruedata[j]);
     }
     fTrueDataArray->AddLast(htemptruedata[0]);
   }
    
  }

  return kTRUE;
  
}


void MatrixUnfControl::CalculateChannelFractions(double bkgsubdata_yields[], double sig_yields[], double channel_SF[]){

  for (int i_chan = 0; i_chan < 4; ++i_chan)
  {
    if( sig_yields[i_chan] > 0 ) channel_SF[i_chan] = ( bkgsubdata_yields[i_chan]/(bkgsubdata_yields[0]+bkgsubdata_yields[1]+bkgsubdata_yields[2]) ) / ( sig_yields[i_chan]/(sig_yields[0]+sig_yields[1]+sig_yields[2]) );
    else channel_SF[i_chan] = 0;
  }

  std::cout<<"channel SF ee emu mumu combined: "<<channel_SF[0]<<" "<<channel_SF[1]<<" "<<channel_SF[2]<<" "<<channel_SF[3]<<std::endl;

}

// Amandeep : There is something going on here
void MatrixUnfControl::CorrectChannelFractions(TH2D*& htemprecoVsgen, TString MCFileName, double channel_SF[]){

  std::cout<<"Correcting channel yield fractions using the following SFs (ee, emu, mumu, combined): "<<channel_SF[0]<<" "<<channel_SF[1]<<" "<<channel_SF[2]<<" "<<channel_SF[3]<<std::endl;

  //check integral and error
  double integralN, integralNerror, integralD, integralDerror;

  integralD = htemprecoVsgen->IntegralAndError(1,htemprecoVsgen->GetNbinsX(),0,htemprecoVsgen->GetNbinsY()+1,integralDerror);
  integralN = htemprecoVsgen->IntegralAndError(1,htemprecoVsgen->GetNbinsX(),1,htemprecoVsgen->GetNbinsY(),integralNerror);
  std::cout<<MCFileName<<" IntegralDB: "<<integralD<<" +/- "<<integralDerror<<" IntegralNB: "<<integralN<<" +/- "<<integralNerror<<std::endl;

  double tempSF=1;
  if( MCFileName.Contains("histosTUnfold_ee_") ) tempSF = channel_SF[0];
  if( MCFileName.Contains("histosTUnfold_emu_") ) tempSF = channel_SF[1];
  if( MCFileName.Contains("histosTUnfold_mumu_") ) tempSF = channel_SF[2];

  double delta_uf[htemprecoVsgen->GetNbinsX()] = {0};
  double deltaerr2_uf[htemprecoVsgen->GetNbinsX()] = {0};


  for(Int_t ri =1; ri <= htemprecoVsgen->GetNbinsX(); ri++)  {
    for(Int_t rj =1; rj <= htemprecoVsgen->GetNbinsY(); rj++) { 
      delta_uf[ri-1] += (1.-tempSF)*htemprecoVsgen->GetBinContent(ri,rj);
      deltaerr2_uf[ri-1] += (1.-tempSF*tempSF)*htemprecoVsgen->GetBinError(ri,rj)*htemprecoVsgen->GetBinError(ri,rj);
      htemprecoVsgen->SetBinContent(ri,rj,tempSF*htemprecoVsgen->GetBinContent(ri,rj));
      htemprecoVsgen->SetBinError(ri,rj,tempSF*htemprecoVsgen->GetBinError(ri,rj));
    }
  }

  for(Int_t ri =1; ri <= htemprecoVsgen->GetNbinsX(); ri++)  { 
      htemprecoVsgen->SetBinContent(ri,0,delta_uf[ri-1]+htemprecoVsgen->GetBinContent(ri,0));
      htemprecoVsgen->SetBinError(ri,0,sqrt(deltaerr2_uf[ri-1]+htemprecoVsgen->GetBinError(ri,0)*htemprecoVsgen->GetBinError(ri,0)) );
  }

  integralD = htemprecoVsgen->IntegralAndError(1,htemprecoVsgen->GetNbinsX(),0,htemprecoVsgen->GetNbinsY()+1,integralDerror);
  integralN = htemprecoVsgen->IntegralAndError(1,htemprecoVsgen->GetNbinsX(),1,htemprecoVsgen->GetNbinsY(),integralNerror);
  std::cout<<MCFileName<<" IntegralDA: "<<integralD<<" +/- "<<integralDerror<<" IntegralNA: "<<integralN<<" +/- "<<integralNerror<<std::endl;

}
