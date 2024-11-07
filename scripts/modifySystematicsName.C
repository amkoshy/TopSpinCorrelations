#include <TString.h>
#include <TFile.h>
#include <TObjString.h>
#include <iostream>

using namespace std;

void modifySystematicsName(const char *filename, const char *newSystematics) {
    TFile *f = TFile::Open(filename, "update");
    if (!f) {
        cout << "cant open " << filename << endl;
        return;
    }
    f->cd("writeNTuple");
    TObjString *sys = new TObjString(newSystematics);
    sys->Write("systematicsName");
    f->Write();
    f->Close();
    delete f;
}
