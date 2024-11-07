#include "ttmd_ZMadGraph.h"
#include "PlotterConfigurationHelper.h"
#include "UnfoldingIOHandler.h"
//#include "ttmd_settings.h"
#include "UnfoldingXSecHandler.h"
#include "ttmd_measur.h"
//#include <TKey.h>
#include <sstream>
#include <ttmd_fileReader.h>

int ZMadGraph::GetNMt(const double mt) const
{
  for(size_t m = 0; m < VMt.size(); m++)
    if(VMt[m] == mt)
      return m;
  throw std::runtime_error(TString::Format("Error: no mt = %f", mt).Data());
}

void ZMadGraph::ResizeStoredPredictions()
{
  const int npdf = 500;
  const int nmu = 9;
  const int nmt = VMt.size();

  //vvvmPredNom.resize(7);
  //vvvmPredCov.resize(7);
  vvvmPredMeas.resize(npdf);
  for(int pdf = 0; pdf < npdf; pdf++)
  {
    //vvvmPredNom[pdf].resize(nmu);
    //vvvmPredCov[pdf].resize(nmu);
    vvvmPredMeas[pdf].resize(nmu);
    for(int mu = 0; mu < nmu; mu++)
    {
      //vvvmPredNom[pdf][mu].resize(nmt);
      //vvvmPredCov[pdf][mu].resize(nmt);
      vvvmPredMeas[pdf][mu].resize(nmt);
    }
  }
}

ZMadGraph::ZMadGraph(PlotterConfigurationHelper *configHelper, const bool flagAG, const bool flagUseNomPDFCovMat)
{
  configHelper_= configHelper;
  if(!flagAG && VMt.size() > 3)
    throw std::logic_error("flagAG = false and VMt.size() > 3 not supported");
  ResizeStoredPredictions();
  FlagAG = flagAG;
  FlagUseNomPDFCovMat = flagUseNomPDFCovMat;
  //std::fill(v5PredAGnbins,v5PredAGnbins + sizeof(v5PredAGnbins), 0);
  std::fill(&v5PredAGnbins[0][0][0][0][0], &v5PredAGnbins[0][0][0][0][0] + 3 * 3 * 7 * 500 * 500, 0);
  if(FlagAG)
  {
    /*vvvvvHistosAG.resize(NOrder);
    for(int ord = 0; ord < NOrder; ord++)
    {
      vvvvvHistosAG[ord].resize(NOrder);
      for(size_t nmt = 0; nmt < VMt.size(); nmt++)
      {
        vvvvvHistosAG[ord][nmt].resize(VMu.size());
        for(size_t nmu = 0; nmu < VMu.size(); nmu++)
        {
          vvvvvHistosAG[ord][nmt][nmu].resize(NPDF);
          for(int npdf = 0; npdf < NPDF; npdf++)
          {
            vvvvvHistosAG[ord][nmt][nmu][npdf].resize(NMaxHistoID);
          }
        }
      }
    }*/
    vvDirAG.resize(NOrder);
    for(int ord = 0; ord < NOrder; ord++)
      vvDirAG[ord].resize(VMt.size());
    for(size_t m = 0; m < VMt.size(); m++)
    {
      if(configHelper_->gFlagMG == 4)
      {
        // Jun18
        TString strBins15 = "b15.15";
        TString strBins10 = "b10.10";
        TString strBins0j = strBins15;
        TString strBins1j = (VMt[m] == 172.5 || VMt[m] == 177.5) ? strBins10 : strBins15;
        // 15.07.18 put symlinks b10 -> b15 for two tt1j directories with b10
        strBins1j = strBins15;
        TString strBins2j = strBins10;
        //vvDirAG[1][m] = gMGTabsDir + TString::Format("/Jun18/aMCfast-tt1j-p0.001-%s-%s/conv", strBins1j.Data(), VMtName[m].Data());
        //vvDirAG[0][m] = gMGTabsDir + TString::Format("/Jun18/aMCfast-tt0j-p0.0002-%s-%s/conv", strBins0j.Data(), VMtName[m].Data());
        // "final" Sep18
        vvDirAG[1][m] = configHelper_->gMGTabsDir + TString::Format("/Sep18/tt1j-iappl2-p0.0003-%s-%s/conv", strBins1j.Data(), VMtName[m].Data());
        vvDirAG[0][m] = configHelper_->gMGTabsDir + TString::Format("/Sep18/tt0j-iappl2-p0.00005-%s-%s/conv", strBins0j.Data(), VMtName[m].Data());
        // still Jun18
        vvDirAG[2][m] = configHelper_->gMGTabsDir + TString::Format("/Jun18/aMCfast-tt2j-p0.003-%s-%s/conv", strBins2j.Data(), VMtName[m].Data());

        /*
        //vvDirAG[0][m] = gMGTabsDir + TString::Format("/aMCfast-140418-iappl2-tt0j-p0.0002-%s/conv", VMtName[m].Data());
        //vvDirAG[0][m] = gMGTabsDir + TString::Format("/aMCfast-190418-iappl2-tt0j-p0.0002-b20.20-%s/conv", VMtName[0].Data());
        if(m < 3)
        {
          // 172.5, 167.5, 177.5
          vvDirAG[0][m] = gMGTabsDir + TString::Format("/aMCfast-200418-iappl2-tt0j-p0.0002-b15.15-%s/conv", VMtName[m].Data());
          vvDirAG[1][m] = gMGTabsDir + TString::Format("/aMCfast-140418-iappl2-tt1j-p0.001-%s/conv", VMtName[m].Data());
          vvDirAG[2][m] = gMGTabsDir + TString::Format("/aMCfast-190418-iappl2-tt2j-p0.003-%s/conv", VMtName[m].Data());
        }
        else
        {
          // 170.0, 175.0
          //vvDirAG[0][m] = gMGTabsDir + TString::Format("/aMCfast-200418-iappl2-tt0j-p0.0002-b15.15-%s/conv", VMtName[m].Data());
          //vvDirAG[1][m] = gMGTabsDir + TString::Format("/aMCfast-250418-iappl2-tt1j-p0.001-%s/conv", VMtName[m].Data());
          //vvDirAG[2][m] = gMGTabsDir + TString::Format("/aMCfast-250418-iappl2-tt2j-p0.003-%s/conv", VMtName[m].Data());
          vvDirAG[0][m] = vvDirAG[0][0];
          vvDirAG[1][m] = vvDirAG[1][0];
          vvDirAG[2][m] = vvDirAG[2][0];
        }
        //for(int ord = 0; ord < NOrder; ord++)
        //  printf("ord = %d m = %d dir = %s\n", ord, m, vvDirAG[ord][m].Data());
        //*/
      }
      else if (configHelper_->gFlagMG == 5)
      {
        vvDirAG[0][m] = configHelper_->gMGTabsDir + TString::Format("/aMCfast-230518-iappl2-tt0j-p0.001-b20.20-%s/conv", VMtName[m].Data());
        vvDirAG[1][m] = configHelper_->gMGTabsDir + TString::Format("/aMCfast-230518-iappl2-tt1j-p0.005-b20.20-%s/conv", VMtName[m].Data());
        //vvDirAG[2][m] = gMGTabsDir + TString::Format("/aMCfast-190418-iappl2-tt2j-p0.003-%s/conv", VMtName[m].Data());
      }
      else
        throw std::logic_error(TString::Format("Error: not supported gFlagMG = %d", configHelper_->gFlagMG).Data());
    }
  }
  else
  {
    vvFile = utils::ttmd::CreateVectorPtr2D<TFile>(VMt.size(), NOrder);
    vvvHistos = utils::ttmd::CreateVectorPtr3D<TH1D>(NOrder, VMt.size(), NMaxHistoID * NWeightOffset);
    if(configHelper_->gFlagMG == 4)
    {
      for(size_t m = 0; m < VMt.size(); m++)
      {
          vvFile[0][m] = TFile::Open(configHelper_->gMGTabsDir + TString::Format("/MADatNLO-070418-iappl0-tt0j-p0.0002-%s.root", VMtName[m].Data()), "read");
          vvFile[1][m] = TFile::Open(configHelper_->gMGTabsDir + TString::Format("/MADatNLO-080418-iappl0-tt1j-p0.001-ptj15-%s.root", VMtName[m].Data()), "read");
          vvFile[2][m] = TFile::Open(configHelper_->gMGTabsDir + TString::Format("/MADatNLO-090418-iappl0-tt2j-p0.003-ptj15-%s.root", VMtName[0].Data()), "read");
      }
    }
    else
      throw std::logic_error(TString::Format("Error: not supported gFlagMG = %d for FlagAG = %d", configHelper_->gFlagMG, FlagAG).Data());
  }
}

ZMadGraph::~ZMadGraph()
{
  printf("timeReadPred: %.2fs  ", timeReadPred);
  printf("timeReadFiles: %.2fs\n", timeReadFiles);
  printf("timeProcess: %.2fs\n", timeProcess);
  printf("timeAddVar: %.2fs\n", timeAddVar);
  //printf("OZ: in ZMadGraph::~ZMadGraph()\n");
  //delete[] vHistos;
  for(auto& vFile : vvFile)
    for(auto& file : vFile)
      if(file)
        file->Close();
  /*for(auto& v1 : vvvvvHistosAG)
    for(auto& v2 : v1)
      for(auto& v3 : v2)
        for(auto& v4 : v3)
          for(auto& v5 : v4)
            delete v5;*/
}

int ZMadGraph::GetOffsetHistoPtjThreshold(PlotterConfigurationHelper *configHelper__, const double pt)
{
  static std::vector<double> vec;
  if(vec.size() == 0)
  {
    if(configHelper__->gFlagMG == 4)
      vec = { 30.0,40.0,50.0,75.0,100.0,150.0 };
    else if(configHelper__->gFlagMG == 5)
      vec = { 30.0 };
    else
      throw std::logic_error(TString::Format("Error: not supported gFlagMG = %d", configHelper__->gFlagMG).Data());
  }
  for(size_t i = 0; i < vec.size(); i++)
    if(vec[i] == pt)
      return i + 1;
  throw std::logic_error(TString::Format("Error: not found pt jet threshold = %f", pt).Data());
}

int ZMadGraph::GetOffsetHistoPtttThreshold(PlotterConfigurationHelper *configHelper__, const double pt)
{
  static std::vector<double> vec;
  if(vec.size() == 0)
  {
    if(configHelper__->gFlagMG == 4)
      vec = { 30.0,40.0,50.0,75.0,100.0,150.0 };
    else if(configHelper__->gFlagMG == 5)
      vec = { 50.0 };
    else
      throw std::logic_error(TString::Format("Error: not supported gFlagMG = %d", configHelper__->gFlagMG).Data());
  }
  if(pt < 0.0)
    return vec.size() - 1;
  for(size_t i = 0; i < vec.size(); i++)
    if(vec[i] == pt)
      return i;
  throw std::logic_error(TString::Format("Error: not found pt(ttbar) threshold = %f", pt).Data());
}

const TH1D* ZMadGraph::ReadPred(const ZPredHisto& predHisto)
{
  if(FlagAG)
  {
    clock_t tt = clock();
    // access histogram
    const int nmt = GetNMt(predHisto.PredSetup.Mt);
    const int nPDF = predHisto.PredSetup.PDF + predHisto.PredSetup.PDFmem - ZPredSetup::CT14nlo;
    TH1D*& h = vvvvvHistosAG[predHisto.Order][nmt][predHisto.PredSetup.NMu][nPDF][predHisto.HistoID];
    if(h)
      return h;

    // read histograms
    assert(predHisto.Order < vvDirAG.size());
    assert(nmt < vvDirAG[predHisto.Order].size());
    //assert(predHisto.PredSetup.NMu < vvDirAG[predHisto.Order][nmt].size());
    assert(vvDirAG[predHisto.Order][nmt] != "");
    // 16.07.18 support alternative scales
    TString _fileName = TString::Format("%s/h%d-mr%.1f-mf%.1f.txt", vvDirAG[predHisto.Order][nmt].Data(), predHisto.HistoID, predHisto.PredSetup.Mur, predHisto.PredSetup.Muf);
    if(predHisto.PredSetup.Mur == 0.0)
    {
      _fileName = TString::Format("%s/h%d-mr%.1f-mf%.1f.txt", vvDirAG[predHisto.Order][nmt].Data(), predHisto.HistoID, 1.0, 1.0);
      _fileName.ReplaceAll("/conv/", TString::Format("-mu%.0f/conv/", predHisto.PredSetup.Muf));
    }
    const TString fileName = _fileName;
    //const int nLines = ReadNumberOfLines(fileName);
    const int nColumns = utils::ttmd::ReadNumberOfColumns(fileName);
    const int nbins = nColumns - 1;

    //static TH1D* hStatic = NULL;
    //if(!hStatic)
    //  hStatic = new TH1D("", "", nbins, 0.0, nbins);
    //for(int i = 0; i < nLines; i++)
      //vvvvvHistosAG[predHisto.Order][nmt][predHisto.PredSetup.NMu][i][predHisto.HistoID] = (TH1D*) hStatic->Clone();
    //  vvvvvHistosAG[predHisto.Order][nmt][predHisto.PredSetup.NMu][i][predHisto.HistoID] = new TH1D("", "", nbins, 0.0, nbins);

    static TH1D* h0 = NULL;
    if(!h0)
      h0 = new TH1D("", "", nbins, 0.0, nbins);
    //h = new TH1D("", "", nbins, 0.0, nbins);
    h = h0;
    h->Reset();
    //const int lineStart = 1;
    const int lineStart = nPDF + 1;

    /*for(int pdf = 0; pdf < nLines; pdf++)
      for (int col = 0; col < nbins; col++)
        vvvvvHistosAG[predHisto.Order][nmt][predHisto.PredSetup.NMu][pdf][predHisto.HistoID]->SetBinContent(col + 1, 1.0);*/

    //int pdf = 0;
    // open file
    utils::ttmd::FileExists(fileName, 1);

    /*std::ifstream file(fileName.Data());
    std::string line;

    // skip until lineStart lines
    for(int l = 1; l < lineStart; l++)
      getline(file, line);

    //while (1)
    {
      // read new line
      getline(file, line);
      //if (true == file.eof()) break;
      //if (line.at(0) == '#' ) continue; //ignore comments
      line.erase(line.find_last_not_of(" \n\r\t")+1); // trim trailing whitespaces
      line.erase(line.find_last_not_of(" ")+1); // trim starting whitespaces
      //line = line.c_str() + line.find_first_of(' ') + 1;
      std::stringstream sline(line);

      double val = 0.0;
      for (int col = 0; col < nbins; col++)
      {
        sline >> val;
        h->SetBinContent(col + 1, val);
      }

      // increment bin
      //pdf++;
    }

    //if(pdf != nLines)
    //  throw std::runtime_error("Error: wrong number of read lines");

    file.close();*/

    clock_t t = clock();
    TextFile* fastFile = NULL;
    fastFile = new TextFile(fileName.Data());
    /*auto it = mFileAG.find(fileName);
    if(it == mFileAG.end())
    {
      fastFile = new TextFile(fileName.Data());
      mFileAG[fileName] = fastFile;
    }
    else
      fastFile = it->second;*/

    //std::string str = fastFile->ReadNthLine(lineStart + 1);
    std::string str = fastFile->ReadNthLine(lineStart - 1);
    //fastFile->close();
    delete fastFile;
    // replace multiple spaces by one space
    //std::string::iterator new_end = std::unique(str.begin(), str.end(), BothAreSpaces);
    //str.erase(new_end, str.end());
    // remove leading space if present
    //if(str[0] == ' ')
    //  str = str.c_str() + 1;
    std::vector<std::string> vStr = split(str, ' ');
    for (int col = 0; col < nbins; col++)
    {
      double val = atof(vStr[col].c_str());
      h->SetBinContent(col + 1, val);
    }
    t = clock() - t;
    timeReadFiles += ((double)t)/CLOCKS_PER_SEC;

    tt = clock() - tt;
    timeReadPred += ((double)tt)/CLOCKS_PER_SEC;

    return h;
  }
  else
  {
    // check order
    assert(predHisto.Order < vvFile.size());

    // get and check mt
    const int nmt = GetNMt(predHisto.PredSetup.Mt);
    assert(nmt < vvFile[predHisto.Order].size());

    // get histo id
    int histoID = GetFullHistoID(predHisto);
    //assert(fullID < vHistos.size());
    TH1D* h = vvvHistos[predHisto.Order][nmt][histoID];
    if(h)
      return h;

    // check that file exists
    assert(vvFile[predHisto.Order][nmt]);

    TString nameHistoFull = TString::Format("h%d", histoID);
    vvFile[predHisto.Order][nmt]->GetObject(nameHistoFull, h);

    assert(h);
    h->SetDirectory(0);
    // set another histogram name
    h->SetName(TString::Format("h%dh%d%s", predHisto.Order, nmt, nameHistoFull.Data()));
    vvvHistos[predHisto.Order][nmt][histoID] = h;

    return h;
  }
}

float*ZMadGraph::ReadPredFast(const ZPredHisto& predHisto, int* nbinsPtr)
{
  clock_t tt = clock();
  // access prediction
  const int nmt = GetNMt(predHisto.PredSetup.Mt);
  const int nPDF = predHisto.PredSetup.PDF + predHisto.PredSetup.PDFmem - ZPredSetup::CT14nlo;
  char& nbins = v5PredAGnbins[predHisto.Order][nmt][predHisto.PredSetup.NMu][nPDF][predHisto.HistoID];
  float* res = v5PredAG[predHisto.Order][nmt][predHisto.PredSetup.NMu][nPDF][predHisto.HistoID];
  // ApplGrid or histogram mode
  if(FlagAG)
  {
    if(nbins != 0)
    {
      // everything is ready
    }
    else
    {
      // read and store prediction
      assert(predHisto.Order < vvDirAG.size());
      assert(nmt < vvDirAG[predHisto.Order].size());
      //assert(predHisto.PredSetup.NMu < vvDirAG[predHisto.Order][nmt].size());
      assert(vvDirAG[predHisto.Order][nmt] != "");

      // 16.07.18 support alternative scales
      //const TString fileName = TString::Format("%s/h%d-mr%.1f-mf%.1f.txt", vvDirAG[predHisto.Order][nmt].Data(), predHisto.HistoID, predHisto.PredSetup.Mur, predHisto.PredSetup.Muf);
      TString _fileName = TString::Format("%s/h%d-mr%.1f-mf%.1f.txt", vvDirAG[predHisto.Order][nmt].Data(), predHisto.HistoID, predHisto.PredSetup.Mur, predHisto.PredSetup.Muf);
      if(predHisto.PredSetup.Mur == 0.0)
      {
        _fileName = TString::Format("%s/h%d-mr%.1f-mf%.1f.txt", vvDirAG[predHisto.Order][nmt].Data(), predHisto.HistoID, 1.0, 1.0);
        _fileName.ReplaceAll("/conv/", TString::Format("-mu%.0f/conv/", predHisto.PredSetup.Muf));
      }
      const TString fileName = _fileName;
      //printf("fileName: %s\n", fileName.Data());

      //const int nLines = ReadNumberOfLines(fileName);
      const int nColumns = utils::ttmd::ReadNumberOfColumns(fileName);
      nbins = nColumns - 1;
      assert(nbins > 0);
      const int lineStart = nPDF + 1;

      // open file
      utils::ttmd::FileExists(fileName, 1);

      /*std::ifstream file(fileName.Data());
      std::string line;

      // skip until lineStart lines
      for(int l = 1; l < lineStart; l++)
        getline(file, line);

      //while (1)
      {
        // read new line
        getline(file, line);
        //if (true == file.eof()) break;
        //if (line.at(0) == '#' ) continue; //ignore comments
        line.erase(line.find_last_not_of(" \n\r\t")+1); // trim trailing whitespaces
        line.erase(line.find_last_not_of(" ")+1); // trim starting whitespaces
        //line = line.c_str() + line.find_first_of(' ') + 1;
        std::stringstream sline(line);

        double val = 0.0;
        for (int col = 0; col < nbins; col++)
        {
          sline >> val;
          h->SetBinContent(col + 1, val);
        }

        // increment bin
        //pdf++;
      }

      //if(pdf != nLines)
      //  throw std::runtime_error("Error: wrong number of read lines");

      file.close();*/

      clock_t t = clock();
      TextFile* fastFile = NULL;
      fastFile = new TextFile(fileName.Data(), lineStart + 1);
      /*auto it = mFileAG.find(fileName);
      if(it == mFileAG.end())
      {
        fastFile = new TextFile(fileName.Data());
        mFileAG[fileName] = fastFile;
        printf("mFileAG size: %lu\n", mFileAG.size());
      }
      else
        fastFile = it->second;*/

      std::string str = fastFile->ReadNthLine(lineStart - 1);
      delete fastFile;
      // replace multiple spaces by one space
      //std::string::iterator new_end = std::unique(str.begin(), str.end(), BothAreSpaces);
      //str.erase(new_end, str.end());
      // remove leading space if present
      //if(str[0] == ' ')
      //  str = str.c_str() + 1;
      std::vector<std::string> vStr = split(str, ' ');
      for (int col = 0; col < nbins; col++)
      {
        double val = atof(vStr[col].c_str());
        res[col] = val;
      }
      t = clock() - t;
      timeReadFiles += ((double)t)/CLOCKS_PER_SEC;
    }
  }
  else
  {
    // check order
    assert(predHisto.Order < vvFile.size());

    // get and check mt
    const int nmt = GetNMt(predHisto.PredSetup.Mt);
    assert(nmt < vvFile[predHisto.Order].size());

    // get MG histo id
    int histoID = GetFullHistoID(predHisto);

    // check that file exists
    assert(vvFile[predHisto.Order][nmt]);

    TString nameHistoFull = TString::Format("h%d", histoID);
    TH1D* h = NULL;
    vvFile[predHisto.Order][nmt]->GetObject(nameHistoFull, h);
    assert(h);
    //h->SetDirectory(0);
    // set another histogram name
    //h->SetName(TString::Format("h%dh%d%s", predHisto.Order, nmt, nameHistoFull.Data()));
    //vvvHistos[predHisto.Order][nmt][histoID] = h;
    nbins = h->GetNbinsX();
    for(int b = 0; b < nbins; b++)
      res[b] = h->GetBinContent(b + 1);
  }

  if(nbinsPtr)
    *nbinsPtr = nbins;

  tt = clock() - tt;
  timeReadPred += ((double)tt)/CLOCKS_PER_SEC;

  return res;
}

std::vector<ZPredSetup> ZMadGraph::GetPredGroup(const TString& name)
{
  std::vector<ZPredSetup> vPred;
  if(name == "CT14PDFUnc")
  {
    vPred.resize(28 * 2);
    //int offset = 0;
    for(int i = 1; i <= 28 * 2; i++)
    {
      vPred[i - 1] = ZPredSetup(ZPredSetup::CT14nlo, i, 1.0, 1.0, MtNom);
      //offset = i;
    }
    assert(ZPredSetup::CT14nlo + offset + 1 == ZPredSetup::CT14nlo_as_0111);
  }
  if(name == "MMHT2014PDFUnc")
  {
    vPred.resize(25 * 2);
    //int offset = 0;
    for(int i = 1; i <= 25 * 2; i++)
    {
      vPred[i - 1] = ZPredSetup(ZPredSetup::MMHT2014nlo68cl, i, 1.0, 1.0, MtNom);
      //offset = i;
    }
    assert(ZPredSetup::MMHT2014nlo68cl + offset + 1 == ZPredSetup::MMHT2014nlo_asmzlargerange);
  }
  if(name == "HERAPDF20PDFUnc")
  {
    vPred.resize(14 * 2);
    //int offset = 0;
    for(int i = 1; i <= 14 * 2; i++)
    {
      vPred[i - 1] = ZPredSetup(ZPredSetup::HERAPDF20_NLO_EIG, i, 1.0, 1.0, MtNom);
      //offset = i;
    }
    assert(ZPredSetup::HERAPDF20_NLO_EIG + offset + 1 == ZPredSetup::HERAPDF20_NLO_VAR);
  }
  if(name == "NNPDF31PDFUnc")
  {
    vPred.resize(100);
    //int offset = 0;
    for(int i = 1; i <= 100; i++)
    {
      vPred[i - 1] = ZPredSetup(ZPredSetup::NNPDF31_nlo_as_0118, i, 1.0, 1.0, MtNom);
      //offset = i;
    }
    assert(ZPredSetup::NNPDF31_nlo_as_0118 + offset + 1 == ZPredSetup::NNPDF31_nlo_as_0120);
  }
  if(name == "ABMP16PDFUnc" || name == "ABMP16as0.118PDFUnc")
  {
    vPred.resize(29);
    //int offset = 0;
    for(int i = 1; i <= 29; i++)
    {
      vPred[i - 1] = ZPredSetup(ZPredSetup::ABMP16als118_5_nlo, i, 1.0, 1.0, MtNom);
      //offset = i;
    }
    assert(ZPredSetup::ABMP16als118_5_nlo + offset + 1 == ZPredSetup::ABMP16als122_5_nlo);
  }
  if(name == "ABMP16as0.114PDFUnc")
  {
    vPred.resize(29);
    //int offset = 0;
    for(int i = 1; i <= 29; i++)
    {
      vPred[i - 1] = ZPredSetup(ZPredSetup::ABMP16als114_5_nlo, i, 1.0, 1.0, MtNom);
      //offset = i;
    }
    assert(ZPredSetup::ABMP16als114_5_nlo + offset + 1 == ZPredSetup::ABMP16als118_5_nlo);
  }
  if(name == "ABMP16as0.122PDFUnc")
  {
    vPred.resize(29);
    //int offset = 0;
    for(int i = 1; i <= 29; i++)
    {
      vPred[i - 1] = ZPredSetup(ZPredSetup::ABMP16als122_5_nlo, i, 1.0, 1.0, MtNom);
      //offset = i;
    }
    assert(ZPredSetup::ABMP16als122_5_nlo + offset + 1 == ZPredSetup::JR14NLO08VF);
  }
  if(name == "JR14PDFUnc")
  {
    vPred.resize(38);
    //int offset = 0;
    for(int i = 1; i <= 38; i++)
    {
      vPred[i - 1] = ZPredSetup(ZPredSetup::JR14NLO08VF, i, 1.0, 1.0, MtNom);
      //offset = i;
    }
    assert(ZPredSetup::JR14NLO08VF + offset + 1 == ZPredSetup::CJ15nlo);
  }
  if(name == "CJ15PDFUnc")
  {
    vPred.resize(48);
    //int offset = 0;
    for(int i = 1; i <= 48; i++)
    {
      vPred[i - 1] = ZPredSetup(ZPredSetup::CJ15nlo, i, 1.0, 1.0, MtNom);
      //offset = i;
    }
    //assert(ZPredSetup::CJ15nlo + offset + 1 == ZPredSetup::CJ15nlo);
  }
  else if(name == "mu")
  {
    vPred.resize(7);
    vPred[0] = ZPredSetup(ZPredSetup::CT14nlo, 0, 2.0, 1.0, MtNom);
    //vPred[0].Title = "#mu_{r}"
    vPred[1] = ZPredSetup(ZPredSetup::CT14nlo, 0, 0.5, 1.0, MtNom);
    vPred[2] = ZPredSetup(ZPredSetup::CT14nlo, 0, 1.0, 2.0, MtNom);
    vPred[3] = ZPredSetup(ZPredSetup::CT14nlo, 0, 1.0, 0.5, MtNom);
    vPred[4] = ZPredSetup(ZPredSetup::CT14nlo, 0, 2.0, 2.0, MtNom);
    vPred[5] = ZPredSetup(ZPredSetup::CT14nlo, 0, 0.5, 0.5, MtNom);
    vPred[6] = ZPredSetup(ZPredSetup::CT14nlo, 0, 0.0, 12.0, MtNom);
  }
  else if(name == "as")
  {
    vPred.resize(2);
    vPred[0] = ZPredSetup(ZPredSetup::CT14nlo_as_0113, 0, 1.0, 1.0, MtNom);
    vPred[1] = ZPredSetup(ZPredSetup::CT14nlo_as_0123, 0, 1.0, 1.0, MtNom);
  }
  else if(name == "mt")
  {
    vPred.resize(2);
    vPred[0] = ZPredSetup(ZPredSetup::CT14nlo, 0, 1.0, 1.0, 167.5);
    vPred[1] = ZPredSetup(ZPredSetup::CT14nlo, 0, 1.0, 1.0, 177.5);
  }
  return vPred;
}

ZPredSetup ZMadGraph::GetPred(const TString& name)
{
  ZPredSetup pred;
  if(name == "CT14")
    pred = ZPredSetup(ZPredSetup::CT14nlo, 0, 1.0, 1.0, MtNom);
  else if(name == "MMHT2014")
    pred = ZPredSetup(ZPredSetup::MMHT2014nlo68cl, 0, 1.0, 1.0, MtNom);
  else if(name == "HERAPDF20")
    pred = ZPredSetup(ZPredSetup::HERAPDF20_NLO_EIG, 0, 1.0, 1.0, MtNom);
  else if(name == "NNPDF31")
    pred = ZPredSetup(ZPredSetup::NNPDF31_nlo_as_0118, 0, 1.0, 1.0, MtNom);
  else if(name == "ABMP16" || name == "ABMP16as0.118")
    pred = ZPredSetup(ZPredSetup::ABMP16als118_5_nlo, 0, 1.0, 1.0, MtNom);
  else if(name == "ABMP16as0.114")
    pred = ZPredSetup(ZPredSetup::ABMP16als114_5_nlo, 0, 1.0, 1.0, MtNom);
  else if(name == "ABMP16as0.122")
    pred = ZPredSetup(ZPredSetup::ABMP16als122_5_nlo, 0, 1.0, 1.0, MtNom);
  else if(name == "JR14")
    pred = ZPredSetup(ZPredSetup::JR14NLO08VF, 0, 1.0, 1.0, MtNom);
  else if(name == "CJ15")
    pred = ZPredSetup(ZPredSetup::CJ15nlo, 0, 1.0, 1.0, MtNom);
  else
    throw std::runtime_error(TString::Format("Error: should not be here"));
  return pred;
}

ZPredSetup::ZPredSetupPDF ZMadGraph::GetPDFForCovMatrix(ZPredSetup::ZPredSetupPDF pdfIn)
{
  if(pdfIn == ZPredSetup::CT14nlo ||
     pdfIn == ZPredSetup::CT14nlo_as_0111 || pdfIn == ZPredSetup::CT14nlo_as_0113 || pdfIn == ZPredSetup::CT14nlo_as_0116 ||
     pdfIn == ZPredSetup::CT14nlo_as_0118 || pdfIn == ZPredSetup::CT14nlo_as_0120 || pdfIn == ZPredSetup::CT14nlo_as_0123)
    return ZPredSetup::CT14nlo;
  if(pdfIn == ZPredSetup::MMHT2014nlo68cl || pdfIn == ZPredSetup::MMHT2014nlo_asmzlargerange)
    return ZPredSetup::MMHT2014nlo68cl;
  if(pdfIn == ZPredSetup::NNPDF31_nlo_as_0116 || pdfIn == ZPredSetup::NNPDF31_nlo_as_0118 || pdfIn == ZPredSetup::NNPDF31_nlo_as_0120)
    return ZPredSetup::NNPDF31_nlo_as_0118;
  if(pdfIn == ZPredSetup::HERAPDF20_NLO_EIG ||
     pdfIn == ZPredSetup::HERAPDF20_NLO_ALPHAS_110 || pdfIn == ZPredSetup::HERAPDF20_NLO_ALPHAS_113 || pdfIn == ZPredSetup::HERAPDF20_NLO_ALPHAS_116 ||
     pdfIn == ZPredSetup::HERAPDF20_NLO_ALPHAS_118 || pdfIn == ZPredSetup::HERAPDF20_NLO_ALPHAS_122 || pdfIn == ZPredSetup::HERAPDF20_NLO_ALPHAS_126)
    return ZPredSetup::HERAPDF20_NLO_EIG;
  if(pdfIn == ZPredSetup::ABMP16als114_5_nlo || pdfIn == ZPredSetup::ABMP16als118_5_nlo || pdfIn == ZPredSetup::ABMP16als122_5_nlo)
    return pdfIn;
  if(pdfIn == ZPredSetup::JR14NLO08VF || pdfIn == ZPredSetup::CJ15nlo)
    return pdfIn;
  throw std::runtime_error(TString::Format("Error: cannot find PDF for cov. matrix for input PDF = %d", pdfIn));
}

std::vector<TH1D*> ZMadGraph::GetPredictionWithUnc(const UnfoldingIOHandler* unfoldingIOHandler, const ZPredSetup& predNom, const std::vector<ZPredSetup>& vPredUnc, const TString& uncType, ZMeasur* th, const bool flagNorm, const bool flagDivideBW, const double rescale)
{
  // uncType: pdfeig, pdfrep, mu, as
  assert(uncType != "pdfrep" || rescale == 1.0);
  assert(!uncType.BeginsWith("pdfeig") || (vPredUnc.size() % 2) == 0);
  assert(uncType != "as" || vPredUnc.size() == 2);
  assert(uncType != "mt" || vPredUnc.size() == 2);
  // size 6 for mu unc if skipping mu12
  assert(uncType != "mu" || (vPredUnc.size() == 7 || vPredUnc.size() == 6));

  std::vector<TH1D*> result(vPredUnc.size());
  //static double timeAddVar = 0.0;
  UnfoldingXSecHandler* unfoldingXSecHandler = new UnfoldingXSecHandler;
  if(predNom.Mt > 0.0)
  {
    TH1D* h = unfoldingXSecHandler->GetPredictionHisto(unfoldingIOHandler, this, predNom, flagNorm, flagDivideBW);
    th->SetNom(h);
  }

  if(uncType == "pdfrep")
    th->ResizeRep(vPredUnc.size());

  for(size_t i = 0; i < vPredUnc.size(); i++)
  {
    TH1D* h = unfoldingXSecHandler->GetPredictionHisto(unfoldingIOHandler, this, vPredUnc[i], flagNorm, flagDivideBW);
    //TH1D* h = unfoldingXSecHandler->GetPredictionHisto(unfoldingIOHandler, this, predNom, flagNorm, flagDivideBW);
    // 20.04.18 (and similar below) for faster dealing with var/sys names in maps inside ZMeasur
    TString var;
    if(uncType.BeginsWith("pdfeig") || uncType == "pdfsymeig")
      var = TString::Format("%zu", i);
    else
      var = TString::Format("%s%zu", uncType.Data(), i);
    TString sys;
    int nvar = -1;
    if(uncType == "pdfrep")
    {
      th->AddRep(i, h);
    }
    else
    {
      if(uncType.BeginsWith("pdfeig") || uncType == "as" || uncType == "mt")
      {
        if(uncType.BeginsWith("pdfeig"))
          sys = TString::Format("%zu", i / 2 + 1);
        else
          sys = TString::Format("%s%zu", uncType.Data(), i / 2 + 1);
        nvar = 2;
        if(uncType == "pdfeig90")
          utils::ttmd::ScaleRelativeHistoVar(h, 1.0 / 1.64, th->GetNom());
        else if(uncType == "pdfeig" || uncType == "as" || uncType == "mt")
          ;
        else
          throw std::logic_error(TString::Format("Error in ReadXfitterTheory(): unknown uncType = %s", uncType.Data()));
      }
      else if(uncType == "pdfsymeig")
      {
        //sys = TString::Format("%s%zu", uncType.Data(), i);
        sys = TString::Format("%zu", i);
        nvar = 1;
      }
      else if(uncType == "mu")
      {
        sys = TString::Format("mu");
        nvar = 7;
      }
      else
        throw std::logic_error(TString::Format("Error in ReadXfitterTheory(): unknown uncType = %s", uncType.Data()));

      if(rescale != 1.0)
        utils::ttmd::ScaleRelativeHistoVar(h, rescale, th->GetNom());
      clock_t t = clock();
      th->AddVar(var, h, nvar, sys);
      t = clock() - t;
      timeAddVar += ((double)t)/CLOCKS_PER_SEC;
      result[i] = h;
    }
  }
  clock_t t = clock();
  th->Process();
  t = clock() - t;
  timeProcess += ((double)t)/CLOCKS_PER_SEC;

  return result;
}

std::vector<ZMeasur*> ZMadGraph::PredMeasVec(const UnfoldingIOHandler* unfoldingIOHandler, const std::pair<int, int> pdf, const std::pair<double, double> mu, const std::vector<double>& vMt)
{
  std::vector<ZMeasur*> vMeas(vMt.size());
  for(size_t i = 0; i < vMt.size(); i++)
    vMeas[i] = PredMeas(unfoldingIOHandler, pdf, mu, vMt[i]);
  return vMeas;
}

std::vector<TH2D*> ZMadGraph::PredCovVec(const UnfoldingIOHandler* unfoldingIOHandler, const std::pair<int, int> pdf, const std::pair<double, double> mu, const std::vector<double>& vMt)
{
  std::vector<TH2D*> vCov(vMt.size());
  for(size_t i = 0; i < vMt.size(); i++)
    vCov[i] = PredMeas(unfoldingIOHandler, pdf, mu, vMt[i])->GetCovFullMatrix();
  return vCov;
}



int ZMadGraph::GetFullHistoID(const ZPredHisto& predHisto) const
{
  int nWeight = -1;
  int nScales = predHisto.PredSetup.GetNScales(predHisto.PredSetup.Mur, predHisto.PredSetup.Muf);
  if(nScales == 1)
    nWeight = predHisto.PredSetup.PDF + predHisto.PredSetup.PDFmem;
  else
  {
    assert(predHisto.PredSetup.PDF == ZPredSetup::CT14nlo && predHisto.PredSetup.PDFmem == 0);
    nWeight = nScales;
  }
  // int nmt = ZMadGraph::GetNMt(predHisto.PredSetup.Mt);
  /*if(nmt == 2)
  {
    printf("nmt\n");
  }*/

  //
  assert(nWeight > 0 && nWeight < NWeightOffset);
  assert(predHisto.HistoID < NMaxHistoID);
  assert(predHisto.Order < NOrder);
  assert(ZMadGraph::GetNMt(predHisto.PredSetup.Mt) < VMt.size());
  int histoID = predHisto.HistoID * NWeightOffset + nWeight;
  return histoID;
  //int fullID = nmt * NOrder * NMaxHistoID * NWeightOffset + predHisto.Order * NMaxHistoID * NWeightOffset + histoID;
  //int fullID = (predHisto.Order + nmt * NOrder) * (NWeightOffset * NMaxHistoID) + histoID;
  //int fullID = (predHisto.Order * ZMadGraph::VMt.size() + nmt) * (ZMadGraph::NWeightOffset * ZMadGraph::NMaxHistoID) + histoID;
  //return fullID;
}

int ZPredSetup::GetNScales(const double mr, const double mf)
{
  if(mr == 1.0 && mf == 1.0) return 1;
  if(mr == 2.0 && mf == 1.0) return 2;
  if(mr == 0.5 && mf == 1.0) return 3;
  if(mr == 1.0 && mf == 2.0) return 4;
  if(mr == 2.0 && mf == 2.0) return 5;
  if(mr == 0.5 && mf == 2.0) return 6;
  if(mr == 1.0 && mf == 0.5) return 7;
  if(mr == 2.0 && mf == 0.5) return 8;
  if(mr == 0.5 && mf == 0.5) return 9;
  throw std::runtime_error(TString::Format("Error: not supported scales mr = %f mf = %f", mr, mf));
}
