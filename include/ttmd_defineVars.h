#ifndef defineVars_h
#define defineVars_h

#include "ttmd_vars.h"
#include "TreePlain.h"
#include "utils.h"
#include "PlotterConfigurationHelper.h"
//#include "ttmd_settings.h"


class ZVarYttSigned: public ZVar
{
  public:
  ZVarYttSigned(const int mode = 0, const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "visytts";
      else
	Expression = "ytts";
      TitleTxt = "y(ttbar)";
      if(mode)
      {
	Expression += TString::Format("KR%d", mode);
        //Expression += gvKRSuffix[mode];
        //TitleTxt += TString::Format("[KR%d]", mode);
        //Title += TString::Format("_{KR%d}", mode);
      }
    }


    /*virtual double Get(const TreePlain* tree) const
    {
      double var = tree->Vars.ytt;
      return var;
    }*/
};

/*class ZVarYtt2jSigned: public ZVar
{
  public:
    ZVarYtt2jSigned(const double ptMin = 30.0)
    {
      if(ptMin == 30.0)
        Expression = "ytt2js";
      else
        Expression = TString::Format("ytt2j%.0fs", ptMin);
      TitleTxt = TString::Format("y(ttbar2jets%.0f)", ptMin);
      Title = TString::Format("|y(t#bar{t}2j%.0f)|", ptMin);
      Units = "";
    }
};

class ZVarYtt2j: public ZVar
{
  public:
    ZVarYtt2j(const double ptMin = 30.0)
    {
      if(ptMin == 30.0)
        Expression = "ytt2js";
      else
        Expression = TString::Format("ytt2j%.0fs", ptMin);
      TitleTxt = TString::Format("|y(ttbar2jets%.0f)|", ptMin);
      Title = TString::Format("|y(t#bar{t}2j%.0f)|", ptMin);
      Units = "";
      NDigits = 2;
    }

    virtual float Get(const TreePlain* tree)
    {
      return GetAbs(tree);
    }

    virtual ZVar* clone()
    {
      return (ZVar*) (new ZVarYtt2j(*this));
    }
};

class ZVarYttajSigned: public ZVar
{
  public:
    ZVarYttajSigned()
    {
      Expression = "yttajs";
      TitleTxt = "y(ttbarjets%.0f)";
    }
};

class ZVarYttaj: public ZVar
{
  public:
    ZVarYttaj()
    {
      Expression = "yttajs";
      TitleTxt = "|y(ttbarjets%.0f)|";
      var->Title = "|y(t#bar{t}j)|";
      var->Units = "";
    }

    virtual float Get(const TreePlain* tree)
    {
      return GetAbs(tree);
    }

    virtual ZVar* clone()
    {
      return (ZVar*) (new ZVarYtt2j(*this));
    }
};*/

class ZVarYtt: public ZVar
{
  public:
    ZVarYtt(const int mode = 0)
    {
      Expression = "ytt";
      TitleTxt = "y(ttbar)";
      //TitleTxt = "|y(ttbar)|";
      //Title = "y(t#bar{t})";
      Title = "|y(t#bar{t})|";
      Units = "";
      NDigits = 2;
      if(mode)
      {
        Expression += TString::Format("KR%d", mode);
        //Expression += gvKRSuffix[mode];
        if(mode != 9)
        {
          TitleTxt += TString::Format("[KR%d]", mode);
          Title += TString::Format("_{KR%d}", mode);
        }
      }
    }

    virtual float Get(const TreePlain* tree)
    {
      return GetAbs(tree);
    }

    virtual ZVar* clone()
    {
      return (ZVar*) (new ZVarYtt(*this));
    }
};

class ZVarMtt: public ZVar
{
  public:
    ZVarMtt(const int mode = 0, const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "vismtt";
      else
	Expression = "mtt";
      TitleTxt = "M(ttbar)";
      Title = "M(t#bar{t})";
      Units = "GeV";
      NDigits = 0;
      if(mode)
      {
        //Expression += gvKRSuffix[mode];
        Expression += TString::Format("KR%d", mode);
        if(mode != 9)
        {
          TitleTxt += TString::Format("[KR%d]", mode);
          Title += TString::Format("_{KR%d}", mode);
        }
      }
    }
};

class ZVarMll: public ZVar
{
  public:
    ZVarMll()
    {
      Expression = "vismll";
      TitleTxt = "mll";
      Units = "GeV";
      Title = TString("M(l#bar{l})");
    }
};

class ZVarMbb: public ZVar
{
  public:
    ZVarMbb()
    {
      Expression = "vismbb";
      TitleTxt = "mbb";
      Units = "GeV";
      Title = TString("M(b#bar{b})");
    }
};

class ZVarMllbb: public ZVar
{
  public:
    ZVarMllbb()
    {
      Expression = "vismllbb";
      TitleTxt = "mllbb";
      Units = "GeV";
      Title = TString("M(l#bar{l}b#bar{b})");
    }
};

class ZVarMllbbmet: public ZVar
{
  public:
    ZVarMllbbmet()
    {
      Expression = "vismllbbmet";
      TitleTxt = "mllbbmet";
      Units = "GeV";
      Title = TString("M(l#bar{l}b#bar{b}MET)");
    }
};

/*class ZVarMtt2j: public ZVar
{
  public:
    ZVarMtt2j(const double ptMin = 30.0)
    {
      if(ptMin == 30.0)
        Expression = TString::Format("mttj");
      else
        Expression = TString::Format("mttj%.0f", ptMin);
      TitleTxt = TString::Format("M(ttbar2jets%.0f)", ptMin);
      Title = TString::Format("M(t#bar{t}2j)_{p_{T}(j) > %.0f GeV}", ptMin);
      Units = "GeV";
      NDigits = 0;
    }
};*/

class ZVarRhoj: public ZVar
{
  public:
    ZVarRhoj(const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const TString& leadStr = "2")
    {
      Expression = utils::ttmd::GetNameExtraJetVar("rhoj", def, iso, pt, leadStr);
      TitleTxt = TString::Format("#rho(ttbar%sjets%.0f)", leadStr.Data(), pt);
      Title = TString::Format("#rho(t#bar{t}%sj)", leadStr.Data());
      TitleExtension = TString::Format("p_{T}(j) > %.0f GeV", pt);
      Units = "";
      NDigits = 2;
    }
};

class ZVarYttjSigned: public ZVar
{
  public:
    ZVarYttjSigned(const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const TString& leadStr = "2")
    {
      Expression = utils::ttmd::GetNameExtraJetVar("yttj", def, iso, pt, leadStr);
      TitleTxt = TString::Format("y(ttbar%sjets%.0f)", leadStr.Data(), pt);
      Title = TString::Format("y(t#bar{t}%sj)", leadStr.Data());
      TitleExtension = TString::Format("p_{T}(j) > %.0f GeV", pt);
      Units = "";
      NDigits = 2;
    }
};

class ZVarYttj: public ZVar
{
  public:
    ZVarYttj(const TString& def = "def", const double iso = 0.4, const double pt = 30.0, const TString& leadStr = "2")
    {
      Expression = utils::ttmd::GetNameExtraJetVar("yttj", def, iso, pt, leadStr);
      TitleTxt = TString::Format("|y(ttbar%sjets%.0f)|", leadStr.Data(), pt);
      Title = TString::Format("|y(t#bar{t}%sj)|", leadStr.Data());
      TitleExtension = TString::Format("p_{T}(j) > %.0f GeV", pt);
      Units = "";
      NDigits = 2;
    }

    virtual float Get(const TreePlain* tree)
    {
      return GetAbs(tree);
    }

    virtual ZVar* clone()
    {
      return (ZVar*) (new ZVarYttj(*this));
    }
};

class ZVarMttaj: public ZVar
{
  public:
    ZVarMttaj()
    {
      Expression = TString::Format("mttaj");
      TitleTxt = "M(ttbaraj)";
      Title = "M(t#bar{t}aj)";
      Units = "GeV";
    }
};

class ZVarPttt: public ZVar
{
  public:
    ZVarPttt(const int mode = 0, const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "vispttt";
      else
	Expression = "pttt";
      Title = "p_{T}(t#bar{t})";
      TitleTxt = "pT(ttbar)";
      Units = "GeV";
      NDigits = 0;
      NDiv = 305;
      if(mode)
      {
        Expression += TString::Format("KR%d", mode);
        //Expression += gvKRSuffix[mode];
        if(mode != 9)
        {
          TitleTxt += TString::Format("[KR%d]", mode);
          Title += TString::Format("_{KR%d}", mode);
        }
      }
    }
};

class ZVarPtll: public ZVar
{
  public:
    ZVarPtll()
    {
      Expression = "visptll";
      Title = "p_{T}(l#bar{l})";
      TitleTxt = "pT(llbar)";
      Units = "GeV";
    }
};

class ZVarPtbb: public ZVar
{
  public:
    ZVarPtbb()
    {
      Expression = "visptbb";
      Title = "p_{T}(b#bar{b})";
      TitleTxt = "pT(bbbar)";
      Units = "GeV";
    }
};

class ZVarDetatt: public ZVar
{
  public:
    ZVarDetatt()
    {
      Expression = "detatt";
      TitleTxt = "Delta_eta(ttbar)";
      NDiv = 305;
    }
};

class ZVarDytt: public ZVar
{
  public:
    ZVarDytt(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "visdytt";
      else
	Expression = "dytt";
      TitleTxt = "Delta_y(ttbar)";
    }
};

class ZVarDetall: public ZVar
{
  public:
    ZVarDetall()
    {
      Expression = "visdetall";
      TitleTxt = "Delta_#eta(llbar)";
    }
};

class ZVarDphitt: public ZVar
{
  public:
    ZVarDphitt(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "visdphitt";
      else
	Expression = "dphitt";
      TitleTxt = "Delta_phi(ttbar)";
    }
};

class ZVarDphill: public ZVar
{
  public:
    ZVarDphill()
    {
      Expression = "visdphill";
      TitleTxt = "Delta_phi(llbar)";
    }
};

/*class ZVarYtLSigned: public ZVar
{
    virtual double Get(const ZTree* tree) const
    {
      double pt1 = (tree->LVt1()).Pt();
      double pt2 = (tree->LVt2()).Pt();
      if(pt1 > pt2)
        return (tree->LVt1()).Rapidity();
      else
        return (tree->LVt2()).Rapidity();
    }

  public:
    ZVarYtLSigned()
    {
      Expression = "ytls";
    }
};

class ZVarYtSSigned: public ZVar
{
    virtual double Get(const ZTree* tree) const
    {
      double pt1 = (tree->LVt1()).Pt();
      double pt2 = (tree->LVt2()).Pt();
      if(pt1 < pt2)
        return (tree->LVt1()).Rapidity();
      else
        return (tree->LVt2()).Rapidity();
    }

  public:
    ZVarYtSSigned()
    {
      Expression = "ytss";
    }
};*/

class ZVarYtSigned: public ZVar
{
  public:
    ZVarYtSigned(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "visyts";
      else
	Expression = "yts";
      TitleTxt = "y(t)";
      //TitleTxt = "|y(t)|";
      //Title = "|y(t)|";
      Title = "y(t)";
    }
};

class ZVarYatSigned: public ZVar
{
  public:
    ZVarYatSigned(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "visyats";
      else
	Expression = "yats";
      TitleTxt = "y(tbar)";
      Title = "y(#bar{t})";
    }
};

class ZVarYt: public ZVar
{
  public:
    ZVarYt()
    {
      Expression = "yt";
      //TitleTxt = "|y(t)|";
      TitleTxt = "y(t)";
      Title = "|y(t)|";
    }

    virtual float Get(const TreePlain* tree)
    {
      return GetAbs(tree);
    }

    virtual ZVar* clone()
    {
      return (ZVar*) (new ZVarYt(*this));
    }
};

class ZVarYtLead: public ZVar
{
  public:
    ZVarYtLead(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "visytLead";
      else
	Expression = "ytLead";
      TitleTxt = "y(t) leading";
      Title = "y(t) leading";
    }
};

class ZVarYtNLead: public ZVar
{
  public:
    ZVarYtNLead(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "visytNLead";
      else
	Expression = "ytNLead";
      TitleTxt = "y(t) trailing";
      Title = "y(t) trailing";
    }
};

class ZVarEtal: public ZVar
{
  public:
    ZVarEtal()
    {
      Expression = "visetal";
      TitleTxt = "eta(l)";
      Title = "#eta(l)";
    }
};

class ZVarEtaal: public ZVar
{
  public:
    ZVarEtaal()
    {
      Expression = "visetaal";
      TitleTxt = "eta(lbar)";
      Title = "#eta(#bar{l})";
    }
};

class ZVarEtalLead: public ZVar
{
  public:
    ZVarEtalLead()
    {
      Expression = "visetalLead";
      TitleTxt = "eta(l) leading";
      Title = "#eta(l) leading";
    }
};

class ZVarEtalNLead: public ZVar
{
  public:
    ZVarEtalNLead()
    {
      Expression = "visetalNLead";
      TitleTxt = "eta(l) trailing";
      Title = "#eta(l) trailing";
    }
};

class ZVarEtabLead: public ZVar
{
  public:
    ZVarEtabLead()
    {
      Expression = "visetabLead";
      TitleTxt = "eta(b) leading";
      Title = "#eta(b) leading";
    }
};

class ZVarEtabNLead: public ZVar
{
  public:
    ZVarEtabNLead()
    {
      Expression = "visetabNLead";
      TitleTxt = "eta(b) trailing";
      Title = "#eta(b) trailing";
    }
};

class ZVarPtt: public ZVar
{
  public:
    ZVarPtt(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "visptt";
      else
	Expression = "ptt";
      TitleTxt = "pT(t)";
    }
};

class ZVarPtat: public ZVar
{
  public:
    ZVarPtat(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "visptat";
      else
	Expression = "ptat";
      TitleTxt = "pT(tbar)";
    }
};

class ZVarPttLead: public ZVar
{
  public:
    ZVarPttLead(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "vispttLead";
      else
	Expression = "pttLead";
      TitleTxt = "pT(t) leading";
    }
};

class ZVarPttNLead: public ZVar
{
  public:
    ZVarPttNLead(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "vispttNLead";
      else
	Expression = "pttNLead";
      TitleTxt = "pT(t) trailing";
    }
};

class ZVarPtl: public ZVar
{
  public:
    ZVarPtl()
    {
      Expression = "visptl";
      TitleTxt = "pT(l)";
    }
};

class ZVarPtal: public ZVar
{
  public:
    ZVarPtal()
    {
      Expression = "visptal";
      TitleTxt = "pT(lbar)";
    }
};

class ZVarPtlLead: public ZVar
{
  public:
    ZVarPtlLead()
    {
      Expression = "visptlLead";
      TitleTxt = "pT(l) leading";
    }
};

class ZVarPtlNLead: public ZVar
{
  public:
    ZVarPtlNLead()
    {
      Expression = "visptlNLead";
      TitleTxt = "pT(l) trailing";
    }
};

class ZVarPtbLead: public ZVar
{
  public:
    ZVarPtbLead()
    {
      Expression = "visptbLead";
      TitleTxt = "pT(b) leading";
    }
};

class ZVarPtbNLead: public ZVar
{
  public:
    ZVarPtbNLead()
    {
      Expression = "visptbNLead";
      TitleTxt = "pT(b) trailing";
    }
};

class ZVarPttTTRestFrame: public ZVar
{
  public:
    ZVarPttTTRestFrame(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "vispttTTRestFrame";
      else
	Expression = "pttTTRestFrame";
      TitleTxt = "pT(t) (ttbar RF)";
    }
};

class ZVarPtatTTRestFrame: public ZVar
{
  public:
    ZVarPtatTTRestFrame(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "visptatTTRestFrame";
      else
	Expression = "ptatTTRestFrame";
      TitleTxt = "pT(tbar) (ttbar RF)";
    }
};

/*class ZVarPttL: public ZVar
{
    virtual double Get(const ZTree* tree) const
    {
      double pt1 = (tree->LVt1()).Pt();
      double pt2 = (tree->LVt2()).Pt();
      if(pt1 > pt2)
        return (tree->LVt1()).Pt();
      else
        return (tree->LVt2()).Pt();
    }

  public:
    ZVarPttL()
    {
      Expression = "pttl";
    }
};

class ZVarPttS: public ZVar
{
    virtual double Get(const ZTree* tree) const
    {
      double pt1 = (tree->LVt1()).Pt();
      double pt2 = (tree->LVt2()).Pt();
      if(pt1 < pt2)
        return (tree->LVt1()).Pt();
      else
        return (tree->LVt2()).Pt();
    }

  public:
    ZVarPttS()
    {
      Expression = "ptts";
    }
};*/

class ZVarNj: public ZVar
{
  public:
    ZVarNj(const TString& def = "def", const double iso = 0.4, const double pt = 30.0)
    {
      //Title = TString("Jets");
      //Title = TString("N_{j}");
      Title = TString("N_{jet}");
      //TitleExtension = TString::Format("p_{T}(j) > %.0f GeV", pt);
      TitleExtension = "";
      TitleTxt = TString::Format("N(jets%.0f)", pt);
      Expression = utils::ttmd::GetNameExtraJetVar("nj", def, iso, pt);
      shortExpression = "nj";
      Units = "";
      NDigits = 1;
    }

    virtual float Get(const TreePlain* tree)
    {
      float var = (float) *(boost::get<std::shared_ptr<int>>(zTrees[tree]));
      return var;
    }

    virtual ZVar* clone()
    {
      return (ZVar*) (new ZVarNj(*this));
    }
};

class ZVarNj_inside: public ZVar
{
  public:
    ZVarNj_inside(const double iso = 0.4, const double pt = 30.0)
    {
      //Title = TString("Jets");
      //Title = TString("N_{j}");
      Title = TString("N_{jet}^{inside}");
      //TitleExtension = TString::Format("p_{T}(j) > %.0f GeV", pt);
      TitleExtension = "";
      TitleTxt = TString::Format("N(jets_inside%.0f)", pt);
      Expression = utils::ttmd::GetNameExtraJetVar("njinetattbar", "", iso, pt);
      shortExpression = "njinetattbar";
      Units = "";
      NDigits = 1;
    }

    virtual float Get(const TreePlain* tree)
    {
      float var = (float) *(boost::get<std::shared_ptr<int>>(zTrees[tree]));
      return var;
    }

    virtual ZVar* clone()
    {
      return (ZVar*) (new ZVarNj_inside(*this));
    }
};

// for control plots

// b-tagged jet multiplicity
class ZVarNbj: public ZVar
{
  public:
    ZVarNbj()
    {
      Title = TString("Jets b-tagged");
      TitleTxt = "nbj";
      Expression = "nbj";
      Units = "";
      NDigits = 0;
    }

    virtual float Get(const TreePlain* tree)
    {
      float var = (float) *(boost::get<std::shared_ptr<int>>(zTrees[tree]));
      return var;
    }

    virtual ZVar* clone()
    {
      return (ZVar*) (new ZVarNbj(*this));
    }
};

// all jet multiplicity (but with jet-lepton cleaning cut)
class ZVarNaj: public ZVar
{
  public:
    ZVarNaj()
    {
      Title = TString("All jets");
      TitleTxt = "naj";
      Expression = "naj";
      Units = "";
      NDigits = 0;
    }

    virtual float Get(const TreePlain* tree)
    {
      float var = (float) *(boost::get<std::shared_ptr<int>>(zTrees[tree]));
      return var;
    }

    virtual ZVar* clone()
    {
      return (ZVar*) (new ZVarNaj(*this));
    }
};

class ZVarNvtx: public ZVar
{
  public:
    ZVarNvtx()
    {
      Expression = "nvtx";
      TitleTxt = "nvtx";
      Units = "";
      NDigits = 0;
      Title = TString("N_{vtx}");
      //NDiv = 305;
    }

    virtual float Get(const TreePlain* tree)
    {
      float var = (float) *(boost::get<std::shared_ptr<int>>(zTrees[tree]));
      return var;
    }

    virtual ZVar* clone()
    {
      return (ZVar*) (new ZVarNvtx(*this));
    }
};

class ZVarMet: public ZVar
{
  public:
    ZVarMet()
    {
      Expression = "met";
      TitleTxt = "met";
      Units = "GeV";
      NDigits = 0;
      Title = TString("E_{T}^{miss}");
    }
};

class ZVarPtDaug: public ZVar
{
  public:
    ZVarPtDaug(const TString daug)
    {
      assert(daug == "j1" || daug == "j2"
             || daug == "ej1" || daug == "ej2");
      Expression = "pt" + daug;
      TitleTxt = "pt" + daug;
      Units = "GeV";
      NDigits = 0;
      Title = TString("p_{T}(" + daug + ")");
    }
};

class ZVarEtaDaug: public ZVar
{
  public:
    ZVarEtaDaug(const TString daug)
    {
      assert(daug == "j1" || daug == "j2"
             || daug == "ej1" || daug == "ej2");
      Expression = "eta" + daug;
      TitleTxt = "eta" + daug;
      Units = "";
      NDigits = 1;
      Title = TString("#eta(" + daug + ")");
    }
};

class ZVarRatioPttMtt: public ZVar
{
  public:
    ZVarRatioPttMtt(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "visratioPttMtt";
      else
	Expression = "ratioPttMtt";
      TitleTxt = "ratioPttMtt";
      Units = "";
      Title = "p_{T}(t)/M(t#bar{t})";
    }
};

class ZVarRatioPtttMtt: public ZVar
{
  public:
    ZVarRatioPtttMtt(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "visratioPtttMtt";
      else
	Expression = "ratioPtttMtt";
      TitleTxt = "ratioPtttMtt";
      Units = "";
      Title = "p_{T}(t#bar{t})/M(t#bar{t})";
    }
};

class ZVarRatioPtbLeadt: public ZVar
{
  public:
    ZVarRatioPtbLeadt()
    {
      Expression = "visratioPtbLeadt";
      TitleTxt = "ratioPtttMtt";
      Units = "";
      Title = "p_{T}(b leading)/P_{T}(t)";
    }
};

class ZVarRatioPtbNLeadt: public ZVar
{
  public:
    ZVarRatioPtbNLeadt()
    {
      Expression = "visratioPtbNLeadt";
      TitleTxt = "ratioPtbNLeadt";
      Units = "";
      Title = "p_{T}(b trailing)/P_{T}(t)";
    }
};

class ZVarRatioPtlt: public ZVar
{
  public:
    ZVarRatioPtlt()
    {
      Expression = "visratioPtlt";
      TitleTxt = "visratioPtlt";
      Units = "";
      Title = "p_{T}(l)/P_{T}(t)";
    }
};

class ZVarRecoPartonMomFraction: public ZVar
{
  public:
    ZVarRecoPartonMomFraction(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "visrecoPartonMomFraction";
      else
	Expression = "recoPartonMomFraction";
      TitleTxt = "x1";
      Units = "";
      Title = "log(x_{1})";
    }
};

class ZVarRecoAntipartonMomFraction: public ZVar
{
  public:
    ZVarRecoAntipartonMomFraction(const bool DoParticle = false)
    {
      if(DoParticle)
	Expression = "visrecoAntipartonMomFraction";
      else
	Expression = "recoAntipartonMomFraction";
      TitleTxt = "x2";
      Units = "";
      Title = "log(x_{2})";
    }
};

#endif
