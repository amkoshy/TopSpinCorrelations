#ifndef TTBAR_measur_h
#define TTBAR_measur_h

#include <map>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>

#include "PlotterConfigurationHelper.h"

class ZMeasur
{
  public:
    ZMeasur(PlotterConfigurationHelper *configHelper);
    ~ZMeasur();

    bool DoEigenvectors = false;
    bool DoEnvelope = true;
    bool FlagShapeOnly = false;

    void SetNom(const TH1D* h);
    void SetCovStatMatrix(const TH2D* h);
    void AddVar(const TString& var, const TH1D* h, int nVarInSys = 0, TString sys = "");

    // TODO integrate this with eigenvector approach
    void AddRep(const int nrep, const TH1D* h);
    void ResizeRep(const int nrep) { zHRep.resize(nrep); }
    //std::vector<TH1D*> GetVHRep() { return zHRep; }
    //TH2D* GetCovRepMatrix() const;

    void Process();

    TH1D* GetNom() const;
    //const TH1D* GetVar(const TString& varName);

    std::pair<TH1D*, TH1D*> GetSystUD() const;
    std::pair<TH1D*, TH1D*> GetTotalUD() const;
    TH2D* GetCovSysMatrix() const { return zHCovSys; }
    TH2D* GetCovStatMatrix() const { return zHCovStat; }
    TH2D* GetCovFullMatrix() const;
    void SetCovFullMatrix(TH2D* h);
    TH2D* GetCovFullMatrixEnvelope() const;
    double Chi2(const TH1D* h, int skip = 0) const;
    double Chi2(const TH1D* h, const std::vector<int>& vSkip) const;
    const std::map<TString, TH1D*>& GetEigenvectors() const { return zHEigens; }

private:
    TH1D* zHNom;
    TH2D* zHCovStat;
    std::map<TString, TH1D*> zHEigens;
    std::map<TString, std::pair<TH1D*, TH1D*> > zHSys;

    std::vector<TH1D*> zHRep;
    TH2D* zHCovRep = NULL;

    TH2D* zHCovSys;
    TH2D* zHCovSysEnvelope;
    TH2D* zHCovFull;
    TH2D* zHCovFullEnvelope;
    std::pair<TH1D*, TH1D*> zHSystUD;
    std::pair<TH1D*, TH1D*> zHTotalUD;

    void AddVarDiff(const TString& var, TH1D* h, int nVarInSys = 0, TString sys = "");
    void AddVarToSys(const TString& sys, const TH1D* h);

    static int nCallAddVar;

    PlotterConfigurationHelper *configHelper_;

};

#endif // TTBAR_measurement_h
