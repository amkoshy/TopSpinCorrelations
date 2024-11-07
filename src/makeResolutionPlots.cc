#include "TChain.h"
#include "TProfile.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include <TCanvas.h>
#include "utils.h"
#include <vector>
#include <cassert>
#include <TLegend.h>

/*class Var
{
  public:
    TString name;

};*/

void DrawGraphAsTH2(const TGraph* g)
{
  const int nb = 10;
  TH2F* h = new TH2F("", "", nb, g->GetXaxis()->GetXmin(), g->GetXaxis()->GetXmax(), nb, g->GetYaxis()->GetXmin(), g->GetYaxis()->GetXmax());
  h->GetXaxis()->SetTitle(g->GetXaxis()->GetTitle());
  h->GetYaxis()->SetTitle(g->GetYaxis()->GetTitle());
  //double max = h->GetMaximum();
  //printf("maximum: %f\n", max);
  for(int n = 0; n < g->GetN(); n++)
    h->Fill(g->GetX()[n], g->GetY()[n], 1.0);
  //h->SetMaximum(10.0);
  h->Draw("same colz");
}

void MoveStatBox()
{
  gPad->Update();
  TPaveStats *s = (TPaveStats*) gPad->GetPrimitive("stats");
  s->SetX1NDC(0.12);
  s->SetY1NDC(0.55);
  s->SetX2NDC(0.5);
  s->SetY2NDC(0.88);
}

int main(int argc, char** argv)
{
  bool flagHighStat = 1;
  //bool flagHighStat = 0;

  std::vector<TString> vFileName;

  vFileName.push_back("./mergedPlainTree_2018/Nominal/emu/emu_ttbarsignalplustau_fromDilepton.root");
  vFileName.push_back("./mergedPlainTree_2018/Nominal/ee/ee_ttbarsignalplustau_fromDilepton.root");
  vFileName.push_back("./mergedPlainTree_2018/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton.root");

  vFileName.push_back("./mergedPlainTree_2017/Nominal/emu/emu_ttbarsignalplustau_fromDilepton_PSweights.root");
  vFileName.push_back("./mergedPlainTree_2017/Nominal/ee/ee_ttbarsignalplustau_fromDilepton_PSweights.root");
  vFileName.push_back("./mergedPlainTree_2017/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton_PSweights.root");

  vFileName.push_back("./mergedPlainTree_2016/Nominal/emu/emu_ttbarsignalplustau.root");
  vFileName.push_back("./mergedPlainTree_2016/Nominal/ee/ee_ttbarsignalplustau.root");
  vFileName.push_back("./mergedPlainTree_2016/Nominal/mumu/mumu_ttbarsignalplustau.root");

  /*
  vFileName.push_back("./plainTree_2018/Nominal/emu/emu_ttbarsignalplustau_fromDilepton_0.root");
  vFileName.push_back("./plainTree_2018/Nominal/ee/ee_ttbarsignalplustau_fromDilepton_0.root");
  vFileName.push_back("./plainTree_2018/Nominal/mumu/mumu_ttbarsignalplustau_fromDilepton_0.root");*/

  TChain* chGen = new TChain("plainTree_gen_step0");
  TChain* chRec = new TChain("plainTree_rec_step8");
  TChain* chRecKR9 = new TChain("plainTree_recSimple_step7L");
  //TChain* chRecSimple = new TChain("plainTree_rec_SimpleKinRec_step8");
  //TChain* chRecSimpleMlbCut = new TChain("plainTree_rec_SimpleKinRecMlbCut_step8");
  //TChain* chRecSimpleMlbCut150 = new TChain("plainTree_rec_SimpleKinRecMlbCut150_step8");
  //TChain* chRecSimpleMW = new TChain("plainTree_rec_SimpleKinRecMW_step8");
  //TChain* chRecSimpleMWE = new TChain("plainTree_rec_SimpleKinRecMWE_step8");
  //std::vector<TChain*> vCh = { chGen, chRec, chRecSimple, chRecSimpleMW };
  //std::vector<TChain*> vCh = { chGen, chRec, chRecSimple, chRecSimpleMlbCut, chRecSimpleMlbCut150, chRecSimpleMW, chRecSimpleMWE };
  // std::vector<TChain*> vCh = { chGen, chRec, chRecSimpleMlbCut };
  std::vector<TChain*> vCh = { chGen, chRec, chRecKR9};

  const TString prefixPlot = "kinrecoPlots/200423/plot-";
  std::vector<TString> vVarName;
  //std::vector<TString> vVarNameGen;
  std::vector<TH2D*> vVarHRange;

  std::vector<std::vector<float> > vVarVal;
  std::vector<std::vector<int> > vVarValInt;
  //std::vector<std::vector<TGraph*> > vVarGraph2d(2);
  //std::vector<std::vector<TProfile*> > vVarProfile(2);
  //std::vector<std::vector<TH1D*> > vVarResidual(2);
  std::vector<std::vector<TGraph*> > vVarGraph2d(6);
  std::vector<std::vector<TProfile*> > vVarProfile(6);
  std::vector<std::vector<TH1D*> > vVarResidual(6);
  std::vector<std::vector<TH2F*> > vVarTH2F(6);
  std::vector<std::pair<double,double>> extremeVals;

  //std::vector<TString> vVarName;

  // kinreco suffices
  /*std::vector<TString> vKRName   = { "full", "vis", "visE", "visllbb",  "ll", "ll+MET", "ll+METtot", "ll+METtot(cor)", "llbar", "llbarS",   "MW",   "MWE",
                                     "visE[mlb]", "ll+MET[mlb]", "ll+METtot[mlb]" };
  std::vector<int> vKRChain   = {         1,     2,      2,         2,     2,        2,           2,                2,       2,        2,      4,       5,
                                               3,             3,                3 };
  std::vector<TString> vKRVarExt = {     "", "KR1",  "KR2",     "KR3", "KR4",    "KR5",       "KR6",            "KR9",   "KR7",    "KR8", "KRMW", "KRMWE",
                                           "KR2",         "KR5",            "KR6" };
  std::vector<int> vKRActive  = {         1,     0,      0,         0,     0,        1,           1,                0,       0,        0,      0,       0,
                                               0,             1,                1 };*/

  std::vector<TString> vKRName ;
  std::vector<size_t> vKRChain;
  std::vector<TString> vKRVarExt;
  std::vector<int> vKRActive;


  // full
  vKRName.push_back("full KR");
  vKRChain.push_back(1);
  vKRVarExt.push_back("");
  vKRActive.push_back(1);

  // loose
  vKRName.push_back("loose KR");
  vKRChain.push_back(2);
  vKRVarExt.push_back("KR9");
  vKRActive.push_back(1);

  //for(int mlb = 2; mlb <= 4; mlb++)
  /*for(int mlb = 2; mlb <= 3; mlb++)
  {
    //TString mlbStr = (mlb == 2) ? "[no M(lb) cut]" : "";
    //if(mlb == 4)
    //  mlbStr = "[M(lb)>150GeV]";
    //TString mlbStr = "";

    // vis
    //vKRName.push_back("visible" + mlbStr);
    //vKRChain.push_back(mlb);
    //vKRVarExt.push_back("KR2");
    //vKRActive.push_back(1);

    // ll
    //vKRName.push_back("ll" + mlbStr);
    //vKRChain.push_back(mlb);
    //vKRVarExt.push_back("KR5");
    //vKRActive.push_back(0);

    // loose
    //vKRName.push_back("loose" + mlbStr);
    //vKRChain.push_back(mlb);
    //vKRVarExt.push_back("KR6");
    //vKRActive.push_back(1);

    // loose2
    //vKRName.push_back("loose KR" + mlbStr);
    vKRChain.push_back(mlb);
    vKRVarExt.push_back("KR9");
    vKRActive.push_back(1);
  }*/

  // ytt
  vVarName.push_back("ytt");
  vVarHRange.push_back(new TH2D("", "", 1, -3.0, 3.0, 1, -3.0, 3.0));
  vVarHRange.back()->GetXaxis()->SetTitle("y(t#bar{t})^{gen}");
  vVarHRange.back()->GetYaxis()->SetTitle("y(t#bar{t})^{rec}");
  utils::ttmd::SetHistoAxisFonts(vVarHRange.back());

  vVarVal.push_back(std::vector<float>(vKRName.size() + 1));
  vVarValInt.push_back(std::vector<int>(vKRName.size() + 1));
  std::vector<double> vBins_ytt;
  extremeVals.push_back(std::pair<double,double> (-3.0,3.0));
  if(flagHighStat)
    vBins_ytt.push_back(-3.0);
  vBins_ytt.push_back(-2.5);
  vBins_ytt.push_back(-2.2);
  vBins_ytt.push_back(-2.0);
  for(int i = 1; i <= 20; i++)
    vBins_ytt.push_back(-2.0 + i * 0.2);
  assert(vBins_ytt.back() == 2.0);
  vBins_ytt.push_back(2.2);
  vBins_ytt.push_back(2.5);
  if(flagHighStat)
    vBins_ytt.push_back(3.0);
  for(size_t s = 0; s < vKRName.size(); s++)
  {
    TString s_var_name = "ytt-" + std::to_string(s);
    vVarGraph2d[0].push_back(new TGraph);
    //vVarProfile[0].push_back(new TProfile("", "", 25, -2.5, 2.5, "s"));
    vVarProfile[0].push_back(new TProfile(s_var_name, "", vBins_ytt.size() - 1, &vBins_ytt[0], "s"));

    vVarResidual[0].push_back(new TH1D("", "", 40, -2.0, 2.0));
    vVarResidual[0].back()->GetXaxis()->SetTitle("y(t#bar{t})^{gen} - y(t#bar{t})^{rec}");

    vVarTH2F[0].push_back(new TH2F("","", 40, -3.0, 3.0, 40, -3.0, 3.0));
    vVarTH2F[0].back()->GetXaxis()->SetTitle("y(t#bar{t})^{gen}");
    vVarTH2F[0].back()->GetYaxis()->SetTitle("y(t#bar{t})^{rec}");
    utils::ttmd::SetHistoAxisFonts(vVarTH2F[0].back());
  }

  // mtt
  vVarName.push_back("mtt");
  vVarHRange.push_back(new TH2D("", "", 1, 0.0, 3000.0, 1, 0.0, 3000.0));
  vVarHRange.back()->GetXaxis()->SetTitle("M(t#bar{t})^{gen} [GeV]");
  vVarHRange.back()->GetYaxis()->SetTitle("M(t#bar{t})^{rec} [GeV]");
  utils::ttmd::SetHistoAxisFonts(vVarHRange.back());

  vVarVal.push_back(std::vector<float>(vKRName.size() + 1));
  vVarValInt.push_back(std::vector<int>(vKRName.size() + 1));
  extremeVals.push_back(std::pair<double,double> (0.,3000.));
  std::vector<double> vBins_mtt;
  vBins_mtt.push_back(300.0);
  vBins_mtt.push_back(350.0);
  vBins_mtt.push_back(400.0);
  for(int i = 1; i <= 10; i++)
    vBins_mtt.push_back(400.0 + i * 30.0);
  assert(vBins_mtt.back() == 700.0);
  for(int i = 1; i <= 6; i++)
    vBins_mtt.push_back(700.0 + i * 50.0);
  assert(vBins_mtt.back() == 1000.0);
  if(flagHighStat)
  {
    vBins_mtt.push_back(1100.0);
    vBins_mtt.push_back(1300.0);
  }
  vBins_mtt.push_back(1500.0);
  if(flagHighStat)
  {
    vBins_mtt.push_back(1700.0);
    vBins_mtt.push_back(2000.0);
    vBins_mtt.push_back(2500.0);
  }
  vBins_mtt.push_back(3000.0);
  for(size_t s = 0; s < vKRName.size(); s++)
  {
    TString s_var_name = "mtt-" + std::to_string(s);
    vVarGraph2d[1].push_back(new TGraph);
    vVarProfile[1].push_back(new TProfile(s_var_name, "", vBins_mtt.size() - 1, &vBins_mtt[0], "s"));
    vVarResidual[1].push_back(new TH1D("", "", 40, -1000.0, 1000.0));
    vVarResidual[1].back()->GetXaxis()->SetTitle("M(t#bar{t})^{gen} - M(t#bar{t})^{rec}");

    vVarTH2F[1].push_back(new TH2F("","", 50, 0.0, 2000.0, 50, 0.0, 2000.0));
    vVarTH2F[1].back()->GetXaxis()->SetTitle("M(t#bar{t})^{gen} [GeV]");
    vVarTH2F[1].back()->GetYaxis()->SetTitle("M(t#bar{t})^{rec} [GeV]");
    utils::ttmd::SetHistoAxisFonts(vVarTH2F[1].back());

  }

  // nj
  vVarName.push_back("njdefIso04Pt30");
  //vVarNameGen.push_back("njdefIso04Pt30");
  vVarHRange.push_back(new TH2D("", "", 1, 0.0, 6.0, 1, 0.0, 5.0));
  vVarHRange.back()->GetXaxis()->SetTitle("N_{jet}^{gen}");
  vVarHRange.back()->GetYaxis()->SetTitle("N_{jet}^{rec}");
  extremeVals.push_back(std::pair<double,double> (0.0,6.0));
  utils::ttmd::SetHistoAxisFonts(vVarHRange.back());
  vVarVal.push_back(std::vector<float>(vKRName.size() + 1));
  vVarValInt.push_back(std::vector<int>(vKRName.size() + 1));
  std::vector<double> vBins_nj;
  vBins_nj.push_back(0.0);
  vBins_nj.push_back(1.0);
  vBins_nj.push_back(2.0);
  vBins_nj.push_back(3.0);
  vBins_nj.push_back(4.0);
  vBins_nj.push_back(5.0);
  //vBins_nj.push_back(6.0);
  //vBins_nj.push_back(7.0);
  //vBins_nj.push_back(8.0);
  //vBins_nj.push_back(10.0);
  for(size_t s = 0; s < vKRName.size(); s++)
  {
    TString s_var_name = "nj-" + std::to_string(s);
    vVarGraph2d[2].push_back(new TGraph);
    vVarProfile[2].push_back(new TProfile(s_var_name, "", vBins_nj.size() - 1, &vBins_nj[0], "s"));
    vVarResidual[2].push_back(new TH1D("", "", 10, -5.0, 5.0));
    vVarResidual[2].back()->GetXaxis()->SetTitle("N_{jet}^{gen} - N_{jet}^{rec}");

    vVarTH2F[2].push_back(new TH2F("","", 6, 0.0, 6.0, 6, 0.0, 6.0));
    vVarTH2F[2].back()->GetXaxis()->SetTitle("N_{jet}^{gen}");
    vVarTH2F[2].back()->GetYaxis()->SetTitle("N_{jet}^{rec}");
    utils::ttmd::SetHistoAxisFonts(vVarTH2F[2].back());
  }

  // ptt
  vVarName.push_back("ptt");
  vVarHRange.push_back(new TH2D("", "", 1, 0.0, 700.0, 1, 0.0, 700.0));
  vVarHRange.back()->GetXaxis()->SetTitle("p_{T}(t)^{gen} [GeV]");
  vVarHRange.back()->GetYaxis()->SetTitle("p_{T}(t)^{rec} [GeV]");
  utils::ttmd::SetHistoAxisFonts(vVarHRange.back());
  extremeVals.push_back(std::pair<double,double> (0.,700.));
  vVarVal.push_back(std::vector<float>(vKRName.size() + 1));
  vVarValInt.push_back(std::vector<int>(vKRName.size() + 1));
  std::vector<double> vBins_ptt;
  for(int i = 0; i <= 15; i++)
    vBins_ptt.push_back(i * 20.0);
  assert(vBins_ptt.back() == 300.0);
  for(int i = 1; i <= 10; i++)
    vBins_ptt.push_back(300.0 + i * 40.0);
  assert(vBins_ptt.back() == 700.0);
  for(size_t s = 0; s < vKRName.size(); s++)
  {
    TString s_var_name = "ptt-" + std::to_string(s);
    vVarGraph2d[3].push_back(new TGraph);
    vVarProfile[3].push_back(new TProfile(s_var_name, "", vBins_ptt.size() - 1, &vBins_ptt[0], "s"));
    vVarResidual[3].push_back(new TH1D("", "", 40, -400.0, 400.0));
    vVarResidual[3].back()->GetXaxis()->SetTitle("p_{T}(t)^{gen} - p_{T}(t)^{rec}");

    vVarTH2F[3].push_back(new TH2F("","", 40, 0.0, 700.0, 40, 0.0, 700.0));
    vVarTH2F[3].back()->GetXaxis()->SetTitle("p_{T}(t)^{gen} [GeV]");
    vVarTH2F[3].back()->GetYaxis()->SetTitle("p_{T}(t)^{rec} [GeV]");
    utils::ttmd::SetHistoAxisFonts(vVarTH2F[3].back());
  }

  // pttt
  vVarName.push_back("pttt");
  vVarHRange.push_back(new TH2D("", "", 1, 0.0, 500.0, 1, 0.0, 500.0));
  vVarHRange.back()->GetXaxis()->SetTitle("p_{T}(t#bar{t})^{gen} [GeV]");
  vVarHRange.back()->GetYaxis()->SetTitle("p_{T}(t#bar{t})^{rec} [GeV]");
  utils::ttmd::SetHistoAxisFonts(vVarHRange.back());
  extremeVals.push_back(std::pair<double,double> (0.0,500.0));
  vVarVal.push_back(std::vector<float>(vKRName.size() + 1));
  vVarValInt.push_back(std::vector<int>(vKRName.size() + 1));
  std::vector<double> vBins_pttt;
  for(int i = 0; i <= 15; i++)
    vBins_pttt.push_back(i * 20.0);
  assert(vBins_pttt.back() == 300.0);
  for(int i = 1; i <= 5; i++)
    vBins_pttt.push_back(300.0 + i * 40.0);
  assert(vBins_pttt.back() == 500.0);
  for(size_t s = 0; s < vKRName.size(); s++)
  {
    TString s_var_name = "pttt-" + std::to_string(s);
    vVarGraph2d[4].push_back(new TGraph);
    vVarProfile[4].push_back(new TProfile(s_var_name, "", vBins_pttt.size() - 1, &vBins_pttt[0], "s"));
    vVarResidual[4].push_back(new TH1D("", "", 40, -400.0, 400.0));
    vVarResidual[4].back()->GetXaxis()->SetTitle("p_{T}(t#bar{t})^{gen} - p_{T}(t#bar{t})^{rec}");

    vVarTH2F[4].push_back(new TH2F("","", 50, 0.0, 500.0, 50, 0.0, 500.0));
    vVarTH2F[4].back()->GetXaxis()->SetTitle("p_{T}(t#bar{t})^{gen} [GeV]");
    vVarTH2F[4].back()->GetYaxis()->SetTitle("p_{T}(t#bar{t})^{rec} [GeV]");
    utils::ttmd::SetHistoAxisFonts(vVarTH2F[4].back());
  }


  // yt
  vVarName.push_back("yt");
  vVarHRange.push_back(new TH2D("", "", 1, -3.0, 3.0, 1, -3.0, 3.0));
  vVarHRange.back()->GetXaxis()->SetTitle("y(t)^{gen} [GeV]");
  vVarHRange.back()->GetYaxis()->SetTitle("y(t)^{rec} [GeV]");
  utils::ttmd::SetHistoAxisFonts(vVarHRange.back());
  extremeVals.push_back(std::pair<double,double> (-3.0,3.0));
  vVarVal.push_back(std::vector<float>(vKRName.size() + 1));
  vVarValInt.push_back(std::vector<int>(vKRName.size() + 1));
  std::vector<double> vBins_yt;
  //vBins_yt.push_back(-3.0);
  vBins_yt.push_back(-2.5);
  vBins_yt.push_back(-2.2);
  vBins_yt.push_back(-2.0);
  for(int i = 1; i <= 20; i++)
    vBins_yt.push_back(-2.0 + i * 0.2);
  assert(vBins_yt.back() == 2.0);
  vBins_yt.push_back(2.2);
  vBins_yt.push_back(2.5);
  //vBins_yt.push_back(3.0);
  for(size_t s = 0; s < vKRName.size(); s++)
  {
    TString s_var_name = "yt-" + std::to_string(s);
    vVarGraph2d[5].push_back(new TGraph);
    vVarProfile[5].push_back(new TProfile(s_var_name, "", vBins_yt.size() - 1, &vBins_yt[0], "s"));
    vVarResidual[5].push_back(new TH1D("", "", 25, -2.5, 2.5));
    vVarResidual[5].back()->GetXaxis()->SetTitle("y(t)^{gen} - y(t)^{rec}");

    vVarTH2F[5].push_back(new TH2F("","", 40, -3.0, 3.0, 40, -3.0, 3.0));
    vVarTH2F[5].back()->GetXaxis()->SetTitle("y(t)^{gen} [GeV]");
    vVarTH2F[5].back()->GetYaxis()->SetTitle("y(t)^{rec} [GeV]");
    utils::ttmd::SetHistoAxisFonts(vVarTH2F[5].back());
  }

  // gen level
  std::vector<TString> vVarGenName;
  std::vector<TH2D*> vVarGenHRange;
  std::vector<std::vector<float> > vVarGenVal;
  std::vector<std::vector<TGraph*> > vVarGenGraph2d;
  std::vector<TProfile*> vVarGenProfile;
  // ttbare
  vVarGenName.push_back("ttbare");
  vVarGenHRange.push_back(new TH2D("", "", 1, 0.0, 4000.0, 1, 0.0, 4000.0));
  vVarGenVal.push_back(std::vector<float>(2));
  vVarGenGraph2d.push_back(std::vector<TGraph*>{ new TGraph });
  vVarGenProfile.push_back(new TProfile("", "", 80, 0.0, 4000.0, "s"));
  // ttbarz
  vVarGenName.push_back("ttbarz");
  vVarGenHRange.push_back(new TH2D("", "", 1, -4000.0, 4000.0, 1, -4000.0, 4000.0));
  vVarGenVal.push_back(std::vector<float>(2));
  vVarGenGraph2d.push_back(std::vector<TGraph*>{ new TGraph });
  vVarGenProfile.push_back(new TProfile("", "", 80, -4000.0, 4000.0, "s"));

  std::vector<int> eventCounter(vCh.size(), -1);
  for(size_t ch = 0; ch < vCh.size(); ch++)
  {
    std::cout << "Loading ch " << ch << std::endl;
    for(size_t f = 0; f < vFileName.size(); f++)
      vCh[ch]->Add(vFileName[f]);

    if(ch == 1 || ch == 2)
    {
      std::vector<double> vMat = { 0.25, 0.5, 1.0 };
      // calculate correct and wrong jet assignments
      printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
      for(size_t m = 0; m < vMat.size(); m++)
      {
        printf("ch = %lu  matched %f:\n", ch, vMat[m]);
        TString strMat = TString::Format("%f", vMat[m]);
        TString strMatchedPt = "matchedPt";
        TString strMatchedR = "matchedR";
        if(ch == 1)
        {
          strMatchedPt += "FullKR";
          strMatchedR += "FullKR";
        }
        printf("  -> %.2f", ((double)vCh[ch]->Draw("", strMatchedPt + " < " + strMat + " && " + strMatchedR + " < " + strMat)/(double)vCh[ch]->GetEntries()));
        if(ch == 1)
          printf("[%.2f]", ((double)vCh[ch]->Draw("", strMatchedPt + "2" + " < " + strMat + " && " + strMatchedR + "2" + " < " + strMat)/(double)vCh[ch]->GetEntries()));
        printf("\n");
      }
      printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
    }

    vCh[ch]->SetBranchStatus("*", 0);
    vCh[ch]->SetBranchStatus("eventCounter", 1);
    vCh[ch]->SetBranchAddress("eventCounter", &eventCounter[ch]);
    for(size_t v = 0; v < vVarName.size(); v++)
    {
      if(ch == 0)
      {
        // gen level
        vCh[ch]->SetBranchStatus(vVarName[v], 1);
        if(vVarName[v].BeginsWith("nj"))
          vCh[ch]->SetBranchAddress(vVarName[v], &vVarValInt[v][0]);
        else
          vCh[ch]->SetBranchAddress(vVarName[v], &vVarVal[v][0]);
      }
      else
      {
        for(size_t vv = 0; vv < vKRChain.size(); vv++)
        {
          if(vKRChain[vv] != ch)
            continue;
          if(vVarName[v].BeginsWith("nj"))
          {
            vCh[ch]->SetBranchStatus(vVarName[v], 1);
            vCh[ch]->SetBranchAddress(vVarName[v], &vVarValInt[v][1 + vv]);
          }

          //else if (vVarName[v] == "mtt" || vVarName[v] == "ytt" || vVarName[v] == "pttt" || vVarName[v] == "minmlb")
          else
          {
            vCh[ch]->SetBranchStatus(vVarName[v] + vKRVarExt[vv], 1);
            vCh[ch]->SetBranchAddress(vVarName[v] + vKRVarExt[vv], &vVarVal[v][1 + vv]);
          }
          //else std::cout << "WARNING: skipping var " << vVarName[v] << " in branch " << ch << std::endl;
        }
      }
    }
  }

  std::cout << "files loaded!!" << std::endl;

  // mll
  chGen->SetBranchStatus("mll", 1);
  float mll;
  chGen->SetBranchAddress("mll", &mll);
  TProfile* prof_mtt_mll = new TProfile("", "", 60, 0.0, 3000.0, "s");

  long maxEvents = 0;
  if(argc >= 2)
    maxEvents = atoi(argv[1]);

  std::cout << "Filling Plots ..." << std::endl;

  std::vector<int> vEntry(vCh.size() - 1, 0);
  std::vector<bool> vFlagRecUsed(vCh.size() - 1, true);
  for(int e0 = 0; e0 < (maxEvents ? maxEvents : chGen->GetEntries()); e0++)
  {
    //printf(" >>>>> e0 = %d\n", e0);
    vCh[0]->GetEntry(e0);
    //
    float mtt = vVarVal[1][0];
    //printf("%f %f\n", mtt, mll);
    prof_mtt_mll->Fill(mtt, mll);
    //
    /*for(size_t v = 0; v < vVarName.size(); v++)
    {
      vVarGenGraph2d[v][0]->SetPoint(vVarGenGraph2d[v][0]->GetN(), vVarGenVal[v][0], vVarGenVal[v][1]);
      vVarGenProfile[v]->Fill(vVarGenVal[v][0], vVarGenVal[v][0] - vVarGenVal[v][1]);
    }*/

    for(size_t ch = 1; ch < vCh.size(); ch++)
    {
      if(vFlagRecUsed[ch - 1])
      {
        vCh[ch]->GetEntry(vEntry[ch - 1]);
        vEntry[ch - 1]++;
        vFlagRecUsed[ch - 1] = false;
      }
      if(eventCounter[0] == eventCounter[ch])
      {
        vFlagRecUsed[ch - 1] = true;
        for(size_t vv = 0; vv < vKRChain.size(); vv++)
        {
          if(vKRChain[vv] != ch)
            continue;
          for(size_t v = 0; v < vVarName.size(); v++)
          {
            //std::cout << vVarName[v] << std::endl;
            double x = (vVarName[v].BeginsWith("nj")) ? vVarValInt[v][0] : vVarVal[v][0];
            double y = (vVarName[v].BeginsWith("nj")) ? vVarValInt[v][1 + vv] : vVarVal[v][1 + vv];
            if(vVarName[v].BeginsWith("nj"))
            {
              // printf("ch = %lu x,y = %d,%d\n", ch, (int)x, (int)y);
            }
            //printf("ch = %lu  vv = %lu  v = %lu  x,y = %f %f\n", ch, vv, v, x, y);
            vVarGraph2d[v][vv]->SetPoint(vVarGraph2d[v][vv]->GetN(), x, y);
            vVarResidual[v][vv]->Fill(y - x);
            vVarTH2F[v][vv]->Fill(x, y);

            double min_val = extremeVals.at(v).first;
            double max_val = extremeVals.at(v).second;

            if (!TMath::IsNaN(x) && !TMath::IsNaN(y)){
                if (x < min_val || y < min_val || x > max_val || y > max_val) continue;
                vVarProfile[v][vv]->Fill(x, y);
            }
          }
        }
      }
    }
  }

  std::cout << "Creating VarGenGraphs ..." << std::endl;

  // gen ttbar e,z
  for(size_t v = 0; v < vVarGenName.size(); v++)
  {
    for(size_t g = 0; g < vVarGenGraph2d[v].size(); g++)
    {
      printf("graph2d var = %s  graph = %lu  npoints = %d\n", vVarGenName[v].Data(), g, vVarGenGraph2d[v][g]->GetN());
      TCanvas* c = new TCanvas("", "", 600, 600);
      c->cd();
      vVarGenHRange[v]->Draw();
      vVarGenGraph2d[v][g]->Draw("p");
      TPaveText* pt = new TPaveText(0.15, 0.85, 0.25, 0.88);
      pt->AddText(vKRName[g]);
      pt->Draw();
      vVarGenHRange[v]->Draw("axis same");
      //MoveStatBox();
      c->SaveAs(TString::Format("%s-gen-graph-%lu.pdf", (prefixPlot + vVarGenName[v]).Data(), g));
      //c->SaveAs(TString::Format("%s-gen-graph-%lu.png", (prefixPlot + vVarGenName[v]).Data(), g));
    }
    // profile
    printf("profile var = %s  profile =  npoints = %.0f\n", vVarGenName[v].Data(), vVarGenProfile[v]->GetEntries());
    TCanvas* c = new TCanvas("", "", 600, 600);
    c->cd();
    //vVarHRange[v]->Draw();
    vVarGenProfile[v]->Fit("pol1", "w");
    vVarGenProfile[v]->Draw();
    //vVarHRange[v]->Draw("axis same");
    //MoveStatBox();
    c->SaveAs(TString::Format("%s-gen-profile.pdf", (prefixPlot + vVarGenName[v]).Data()));
    //c->SaveAs(TString::Format("%s-gen-profile.png", (prefixPlot + vVarGenName[v]).Data()));
  }

  std::cout << std::endl << "Plotting VarGraphs ..." << std::endl;

  gStyle->SetOptFit(1);
  utils::ttmd::SetOZStyle();
  gStyle->SetOptStat("mruo");
  for(size_t v = 0; v < vVarName.size(); v++)
  {
    for(size_t g = 0; g < vVarGraph2d[v].size(); g++)
    {
      if(g != 2 && vVarName[v].Contains("nocor")) continue;
      std::cout << std::endl << std::endl;
      printf("graph2d var = %s  graph = %lu  npoints = %d\n", vVarName[v].Data(), g, vVarGraph2d[v][g]->GetN());
      TCanvas* c = new TCanvas("", "", 600, 600);
      c->cd();

      //std::cout << vVarGraph2d[v][g]->GetMean(1) << "    " << vVarGraph2d[v][g]->GetMean(2) << std::endl;
      //vVarHRange[v]->Draw();
      //vVarGraph2d[v][g]->Draw();
      //vVarGraph2d[v][g]->Draw("colz");
      //DrawGraphAsTH2(vVarGraph2d[v][g]);
      //vVarHRange[v]->Draw("axis same");
      //MoveStatBox();
      //c->SaveAs(TString::Format("%s-graph-%lu.pdf", (prefixPlot + vVarName[v]).Data(), g));
      gStyle->SetOptStat(0);
      vVarTH2F[v][g]->Draw("COLZ");
      TAxis* xaxis = vVarTH2F[v][g]->GetXaxis();
      utils::ttmd::SetAxisFonts(xaxis,63,15);
      TAxis* yaxis = vVarTH2F[v][g]->GetYaxis();
      utils::ttmd::SetAxisFonts(yaxis,63,15);
      TAxis* zaxis = vVarTH2F[v][g]->GetZaxis();
      utils::ttmd::SetAxisFonts(zaxis,63,15);

      utils::ttmd::SaveCanvas(c, TString::Format("%s-graph-%lu", (prefixPlot + vVarName[v]).Data(), g));
      //c->SaveAs(TString::Format("%s-graph-%lu.png", (prefixPlot + vVarName[v]).Data(), g));
      // profile
      gStyle->SetOptStat("mruo");
      printf("profile var = %s  profile = %lu  npoints = %.0f\n", vVarName[v].Data(), g, vVarProfile[v][g]->GetEntries());
      c = new TCanvas("", "", 600, 600);
      c->cd();
      //vVarHRange[v]->Draw();
      vVarProfile[v][g]->Fit("pol1", "w");
      vVarProfile[v][g]->Draw();
      //vVarHRange[v]->Draw("axis same");
      //MoveStatBox();
      c->SaveAs(TString::Format("%s-profile-%lu.pdf", (prefixPlot + vVarName[v]).Data(), g));
      //c->SaveAs(TString::Format("%s-profile-%lu.png", (prefixPlot + vVarName[v]).Data(), g));

      // residual
      c = new TCanvas("", "", 600, 600);
      c->cd();
      vVarResidual[v][g]->Fit("gaus");
      vVarResidual[v][g]->SetMarkerColor(1);
      vVarResidual[v][g]->SetLineColor(1);
      //vVarResidual[v][g]->Draw("e0");
      vVarResidual[v][g]->Draw();
      utils::ttmd::SaveCanvas(c, TString::Format("%s-residual-%lu", (prefixPlot + vVarName[v]).Data(), g));
    }
  }

  // mtt vs mll
  TCanvas* c = new TCanvas("", "", 600, 600);
  c->cd();
  prof_mtt_mll->Fit("pol1");
  prof_mtt_mll->Draw();
  utils::ttmd::SaveCanvas(c, prefixPlot + "mttmll");

  // one plot with mean and RMS
  std::vector<std::vector<TH1D*> > vResMean(vVarGraph2d.size());
  std::vector<std::vector<TH1D*> > vResRMS(vVarGraph2d.size());
  TH2D* hr[6];
  hr[0] = new TH2D("hr0", "", 1, -3.0, 3.0, 1, -1.0, 1.0);
  hr[0]->GetXaxis()->SetTitleOffset(1.2);
  hr[0]->GetYaxis()->SetTitleOffset(1.7);
  hr[0]->GetXaxis()->SetTitle("y(t#bar{t})^{gen}");
  hr[0]->GetYaxis()->SetTitle("mean & RMS [y(t#bar{t})^{gen} - y(t#bar{t})^{rec}]");
  utils::ttmd::SetHistoAxisFonts(hr[0]);
  hr[1] = new TH2D("hr1", "", 1, 300.0, 3000.0, 1, -100.0, 1000.0);
  hr[1]->GetXaxis()->SetTitle("M(t#bar{t})^{gen} [GeV]");
  hr[1]->GetYaxis()->SetTitle("mean & RMS [M(t#bar{t})^{gen} - M(t#bar{t})^{rec}] [GeV]");
  hr[1]->GetXaxis()->SetTitleOffset(1.2);
  hr[1]->GetYaxis()->SetTitleOffset(1.7);
  utils::ttmd::SetHistoAxisFonts(hr[1]);
  hr[2] = new TH2D("hr2", "", 1, 0.0, 10.0, 1, -1.0, 5.0);
  hr[2]->GetXaxis()->SetTitle("N_{jet}^{gen} [GeV]");
  hr[2]->GetYaxis()->SetTitle("mean & RMS [N_{jet}^{gen} - N_{jet}^{rec}]");
  hr[2]->GetXaxis()->SetTitleOffset(1.2);
  hr[2]->GetYaxis()->SetTitleOffset(1.7);
  utils::ttmd::SetHistoAxisFonts(hr[2]);
  hr[3] = new TH2D("hr3", "", 1, 0.0, 700.0, 1, -100.0, 350.0);
  hr[3]->GetXaxis()->SetTitle("p_{T}(t)^{gen} [GeV]");
  hr[3]->GetYaxis()->SetTitle("mean & RMS [p_{T}(t)^{gen} - p_{T}(t)^{rec}]");
  hr[3]->GetXaxis()->SetTitleOffset(1.2);
  hr[3]->GetYaxis()->SetTitleOffset(1.7);
  utils::ttmd::SetHistoAxisFonts(hr[3]);
  hr[4] = new TH2D("hr4", "", 1, 0.0, 500.0, 1, -100.0, 350.0);
  hr[4]->GetXaxis()->SetTitle("p_{T}(t#bar{t})^{gen} [GeV]");
  hr[4]->GetYaxis()->SetTitle("mean & RMS [p_{T}(t#bar{t})^{gen} - p_{T}(t#bar{t})^{rec}]");
  hr[4]->GetXaxis()->SetTitleOffset(1.2);
  hr[4]->GetYaxis()->SetTitleOffset(1.7);
  utils::ttmd::SetHistoAxisFonts(hr[4]);
  hr[5] = new TH2D("hr5", "", 1, -3.0, 3.0, 1, -1.0, 1.5);
  hr[5]->GetXaxis()->SetTitle("y(t)^{gen}");
  hr[5]->GetYaxis()->SetTitle("mean & RMS [y(t)^{gen} - y(t)^{rec}]");
  hr[5]->GetXaxis()->SetTitleOffset(1.2);
  hr[5]->GetYaxis()->SetTitleOffset(1.7);
  utils::ttmd::SetHistoAxisFonts(hr[5]);
  double markerSize = 1.0;
  utils::ttmd::SetOZStyle();
  for(size_t s = 0; s < vVarProfile.size(); s++)
  {
    TCanvas* c = new TCanvas("", "", 600, 600);
    c->SetMargin(0.13, 0.05, 0.10, 0.03);
    TLegend* leg = NULL;
    if(s == 0)
      leg = new TLegend(0.4, 0.12, 0.94, 0.37);
    else if(s == 1)
      leg = new TLegend(0.15, 0.68, 0.70, 0.93);
    else
      leg = new TLegend(0.15, 0.68, 0.70, 0.93);
    leg->SetNColumns(3);
    //leg->SetFillColor(0);
    c->cd();
    if(s == 1)
    {
      gPad->SetLogx();
      hr[s]->GetXaxis()->SetMoreLogLabels();
    }
    hr[s]->Draw();
    vResMean[s].resize(vVarProfile[s].size());
    vResRMS[s].resize(vVarProfile[s].size());
    int lastColor = 1;
    for(size_t g = 0; g < vVarProfile[s].size(); g++)
    {
      if(vKRActive[g] == 0)
        continue;
      //if(g == 3 || g == 4 || g == 6) continue;
      //if(g == 4 || g == 6) continue;
      TString vResMean_name = "ResMean-" + std::to_string(s) + "-" + std::to_string(g);
      TString vResRMS_name = "ResRMS-" + std::to_string(s) + "-" + std::to_string(g);

      vResMean[s][g] = vVarProfile[s][g]->ProjectionX(vResMean_name);
      vResMean[s][g]->SetLineColor(utils::ttmd::GetDistinctColor(lastColor));
      vResMean[s][g]->SetMarkerColor(utils::ttmd::GetDistinctColor(lastColor));
      vResMean[s][g]->SetMarkerStyle(24);
      vResMean[s][g]->SetMarkerSize(markerSize);
      vResRMS[s][g] = vVarProfile[s][g]->ProjectionX(vResRMS_name);
      vResRMS[s][g]->SetLineColor(utils::ttmd::GetDistinctColor(lastColor));
      vResRMS[s][g]->SetMarkerColor(utils::ttmd::GetDistinctColor(lastColor));
      lastColor++;
      vResRMS[s][g]->SetMarkerStyle(20);
      double meanAv = 0.0;
      double meanAvErr = 0.0;
      int goodBins = 0;
      vResRMS[s][g]->SetMarkerSize(markerSize);

      vVarProfile[s][g]->Print("all");
      for(int b = 0; b < vResMean[s][g]->GetNbinsX(); b++)
      {
        if(!TMath::IsNaN(vResMean[s][g]->GetBinContent(b + 1)))
        {
          vResMean[s][g]->SetBinContent(b + 1, vVarProfile[s][g]->GetBinCenter(b + 1) - vVarProfile[s][g]->GetBinContent(b + 1));
          vResRMS[s][g]->SetBinContent(b + 1, vVarProfile[s][g]->GetBinError(b + 1));
          meanAv += TMath::Abs(vResMean[s][g]->GetBinContent(b + 1));
          meanAvErr += TMath::Abs(vResRMS[s][g]->GetBinContent(b + 1));
          goodBins++;
        }
        else
        {
          //vResMean[s][g]->SetBinContent(b + 1, vVarProfile[s][g]->GetBinCenter(b + 1) - vVarProfile[s][g]->GetBinContent(b + 1));
          //vResRMS[s][g]->SetBinContent(b + 1, vVarProfile[s][g]->GetBinError(b + 1));
          vResMean[s][g]->SetBinContent(b + 1, -10000.0);
          vResRMS[s][g]->SetBinContent(b + 1, -10000.0);
        }
      }

      if(g == 0 || s < 3 || s == 4)
      {
          std::cout << "meanAv: " << meanAv << std::endl;
          std::cout << "NbinsX(mean): " << vResMean[s][g]->GetNbinsX() << std::endl;
          std::cout << "NbinsX(RMS): " << vResRMS[s][g]->GetNbinsX() << std::endl;
          std::cout << "good bins: " << goodBins << std::endl;

          //vResMean[s][g]->Print("all");
          //vResRMS[s][g]->Print("all");

        vResMean[s][g]->Draw("p hist same");
        vResRMS[s][g]->Draw("p hist same");
        //vResMean[s][g]->Draw("cp hist same");
        //vResRMS[s][g]->Draw("cp hist same");
        leg->AddEntry((TObject*)(NULL), vKRName[g], "");
        //leg->AddEntry(vResMean[s][g], "mean", "p");
        if(s == 0)
        {
          leg->AddEntry(vResMean[s][g], TString::Format("mean [%.2f]", meanAv / goodBins), "p");
          leg->AddEntry(vResRMS[s][g], TString::Format("RMS [%.2f]", meanAvErr / goodBins), "p");
        }
        else if(s == 1)
        {
          leg->AddEntry(vResMean[s][g], TString::Format("mean [%.0f]", meanAv / goodBins), "p");
          leg->AddEntry(vResRMS[s][g], TString::Format("RMS [%.0f]", meanAvErr / goodBins), "p");
        }
        else
        {
          leg->AddEntry(vResMean[s][g], TString::Format("mean [%.0f]", meanAv / goodBins), "p");
          leg->AddEntry(vResRMS[s][g], TString::Format("RMS [%.0f]", meanAvErr / goodBins), "p");
        }
      }
    }
    leg->Draw();
    hr[s]->Draw("axis same");
    TString fileName = prefixPlot + "summary-" + vVarName[s];
    utils::ttmd::SaveCanvas(c, fileName);
  }

  return 0;
}
