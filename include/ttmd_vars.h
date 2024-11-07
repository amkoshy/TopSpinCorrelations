#ifndef ZVars_h
#define ZVars_h

#include <vector>
#include "TreePlain.h"
#include <TH1D.h>
#include <exception>
#include <map>
#include <cassert>

class ZVar
{
  protected:
    std::map<const TreePlain*, boost::variant<std::shared_ptr<int>, std::shared_ptr<float>, std::vector<std::shared_ptr<float>> >> zTrees;
    TString shortExpression;
  public:
    std::vector<std::vector<double> > Bins;
    int NActiveBinnings = -1;
    TString Expression;
    TString ExpressionCut;
    TString Title;
    TString TitleExtension;
    TString TitleTxt;
    TString Units;
    int NDigits = 1;
    int NDiv = -1;

    /*TString& ExpressionShort()
    {
      return shortExpression;
    }*/
    TString ExpressionShort() const
    {
      if(shortExpression == "")
        return Expression;
      else
        return shortExpression;
    }

    /*TString TitleSplitted(const int nParts) const
    {
      // nParts >= 0 returns part from nParts-th to (nParts+1)-th ';'
      // nParts = -1 returns all title but ';' are separated by '; '
      // nParts = -2 returns all after 1st ';'
      printf("TitleSplitted nParts = %d\n", nParts);
      printf("Title: %s\n", Title.Data());

      if(nParts == -1)
      {
        int nPart = 0;
        TString str = this->TitleSplitted(nPart);
        while(true)
        {
          nPart++;
          TString strSuffix = this->TitleSplitted(nPart);
          printf("strSuffix: %s\n", strSuffix.Data());
          if(strSuffix != "")
            str += "; " + strSuffix;
          else
            break;
        }
        printf("str: %s\n", str.Data());
        return str;
      }

      if(nParts == -2)
      {
        TString strAll = this->TitleSplitted(-1);
        TString strBeginning = this->TitleSplitted(0);
        if(strAll.Length() == strBeginning.Length())
          return "";
        TString str = strAll.Data() + strBeginning.Length() + 1;
        printf("str: %s\n", str.Data());
        return str;
      }

      int pos = 0;
      //printf("TitleSplitted  nParts = %d  this->Title.Data() = %s\n", nParts, this->Title.Data());
      for(int p = 0; p < nParts + 1; p++)
      {
        int posCurrent = TString(this->Title.Data() + pos).First(';');
        //printf("posCurrent = %d\n", posCurrent);
        if(posCurrent == -1)
        {
          if(nParts == 0)
            return "";
          else
            return TString(this->Title.Data() + pos);
        }
        if(p == nParts)
          return TString(this->Title.Data() + pos, posCurrent);
        pos += posCurrent + 1;
        //printf("pos = %d\n", pos);
      }
      return "---";
    }*/

    static TString TitleWithExtensionCombined(std::vector<ZVar*> vVars)
    {
      TString xTitle;
      TString xTitleSuf;
      for(size_t i = 0; i < vVars.size(); i++)
      {
        if(i > 0)
          xTitle += ",";
        xTitle += vVars[i]->Title;
        TString suf = vVars[i]->TitleExtension;
        assert(suf == "" || i == 0 || suf == xTitleSuf);
        if(i == 0)
          xTitleSuf = suf;
      }
      if(xTitleSuf != "")
        xTitle += "; " + xTitleSuf;
      return xTitle;
    }

    static TString TitleFull(std::vector<ZVar*> vVars)
    {
      TString xsecTitle = "[" + ZVar::TitleWithExtensionCombined(vVars) + "]";
      if(vVars.size() == 3 && vVars[0]->Expression.Contains("nj"))
      {
        //xsecTitle += TString::Format(" (%lu N_{j} bins)", vVars[0]->BinsC().size() - 1);
        int nb = vVars[0]->Bins[0].size() - 1;
        if(nb == 2)
          xsecTitle.ReplaceAll("N_{jet}", "N_{jet}^{0,1+}");
        else if(nb == 3)
          xsecTitle.ReplaceAll("N_{jet}", "N_{jet}^{0,1,2+}");
        else if(nb == 4)
          xsecTitle.ReplaceAll("N_{jet}", "N_{jet}^{0,1,2,3+}");
      }
      xsecTitle.ReplaceAll("|y(t#bar{t})|", "y(t#bar{t})");
      xsecTitle.ReplaceAll("|y(t)|", "y(t)");
      return xsecTitle;
    }

    int NBins(const int n) const { return Bins[n].size() - 1; }
    // TODO fix inconsisteny with BinsF(), BinsC()
    //int NBinsC() const { return Bins[1].size() - 1; }
    //int NBinsF() const { return Bins[0].size() - 1; }
    int NBinnings() const { return (NActiveBinnings < 0) ? Bins.size() : NActiveBinnings; }
    const std::vector<double>& BinsF() const { return Bins[1]; }
    const std::vector<double>& BinsC() const { return Bins[0]; }
    
    TString GetXTitle() const
    {
      TString str = this->Title;
      if(Units != "")
        str += TString::Format(" [%s]", Units.Data());
      if(this->TitleExtension != "")
        str += "; " + this->TitleExtension;
      return str;
    }
    
    TString GetYTitleXSec() const
    {
      TString str = TString::Format("d#sigma/d%s", this->Title.Data());
      if(Units != "")
        str += TString::Format(" [pb %s^{-1}]", Units.Data());
      else
        str += TString::Format(" [pb]");
      if(this->TitleExtension != "")
        str += "; " + this->TitleExtension;
      return str;
    }

    TString GetYTitleXSecNorm() const
    {
      TString str = TString::Format("1/#sigma d#sigma/d%s", this->Title.Data());
      if(Units != "")
        str += TString::Format(" [%s^{-1}]", Units.Data());
      if(this->TitleExtension != "")
        str += "; " + this->TitleExtension;
      return str;
    }
    
    TH1D* CreateTH1D(const int binning, const TString name = "", const TString title = "") const
    {
      TH1D* h = new TH1D(name, title, Bins[binning].size() - 1, &Bins[binning][0]);
      h->GetXaxis()->SetTitle(GetXTitle());
      //h->GetYaxis()->SetTitle("Events");
      return h;
    }
    
    TString GetIneq(const int binning, const int b) const
    {
      TString str = TString::Format("%.*f < %s < %.*f", NDigits, Bins[binning][b], Title.Data(), NDigits, Bins[binning][b + 1]);
      if(Units != "")
        str += " " + Units;
      //printf("str: %s\n", str.Data());
      return str;
    }

    void RegisterTree(TreePlain* tree)
    {
      if(zTrees.find(tree) != zTrees.end())
        throw std::logic_error(TString::Format("Error in ZVar::RegisterTree(): another call for tree = %s", "tree.Data()"));
      tree->InitVar(Expression, ExpressionCut, zTrees);
    }

    void DeregisterTree(TreePlain* tree)
    {
      if(zTrees.find(tree) == zTrees.end())
        throw std::logic_error(TString::Format("Error in ZVar::DeregisterTree(): no registered tree = %s", "tree.Data()"));
      if (zTrees[tree].type() == typeid(std::vector<std::shared_ptr<float>>))
	boost::get<std::vector<std::shared_ptr<float>>>(zTrees[tree]).clear();
      zTrees.erase(tree);
    }

    virtual float Get(const TreePlain* tree)
    {

      if (zTrees[tree].type() == typeid(std::vector<std::shared_ptr<float>>)) {
	float var1 = *(boost::get<std::vector<std::shared_ptr<float>>>(zTrees[tree])[0]);
	float var2 = *(boost::get<std::vector<std::shared_ptr<float>>>(zTrees[tree])[1]);

	float var = 0.;
	if(var2 != 0.)
	  var = var1/var2;
	else
	  var = 0;

        //printf("Exp = %s var1 = %f\n", Expression.Data(), var1);
        //printf("Exp = %s var2 = %f\n", Expression.Data(), var2);
	//printf("Exp = %s var = %f\n", Expression.Data(), var);

	return var;
      }
      else if (zTrees[tree].type() == typeid(std::shared_ptr<float>)) {

	float var = *(boost::get<std::shared_ptr<float>>(zTrees[tree]));

	if (Expression == "recoPartonMomFraction" || Expression == "visrecoPartonMomFraction" ||
	    Expression == "recoAntipartonMomFraction" || Expression == "visrecoAntipartonMomFraction")
	  {
	    var = TMath::Log10(var);
	  }


	//printf("Exp = %s var = %f\n", Expression.Data(), var);

	return var;

      }
      else if (zTrees[tree].type() == typeid(std::shared_ptr<int>)) {

	int var = *(boost::get<std::shared_ptr<int>>(zTrees[tree]));

	//printf("Exp = %s var = %d\n", Expression.Data(), var);

	return var;
      }
      else {
	throw std::logic_error(TString::Format("Error in ZVar::Get(): invalid typeid"));
      }

    }

    //virtual ZVar* clone() const = 0;
    // TODO is there better solution, without replicating code into each derived class which needs to be copied?
    // (currently needed to do so only in derived classes which redefine Get())
    virtual ZVar* clone()
    {
      return new ZVar(*this);
    }

  protected:
    float GetAbs(const TreePlain* tree)
    {
      float var = ZVar::Get(tree);
      var = TMath::Abs(var);
      return var;
    }
};

#endif
