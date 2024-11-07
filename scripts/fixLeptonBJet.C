#include <TFile.h>
#include <TH1.h>
void fixLeptonBJet(char* filename)
{
    TH1D *lb;
    TFile f(filename, "UPDATE");
    f.GetObject("HypLeptonBAjetMass", lb);
    if (lb) {
        cout << "HypLeptonBAjetMass found, saving as HypLeptonBjetMass\n";
        lb->SetName("HypLeptonBjetMass");
        lb->Write();
    }
    f.Close();
}