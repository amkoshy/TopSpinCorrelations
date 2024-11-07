#include <string>
#include <vector>

#include "../../common/include/utils.h"



using namespace std;

void printGenXsec(string quantity, string variation) {

    auto fileReader = RootFileReader::getInstance();
    string filename = "selectionRoot/closure/emu/emu_ttbarsignalplustau_fakerun_" + variation + ".root";
    
    TH1 *unbinned; vector<double> bins;
    if (quantity == "pt") {
        unbinned = fileReader->GetClone<TH1D>(filename.c_str(), "VisGenToppT");
        unbinned->Add(fileReader->Get<TH1D>(filename.c_str(), "VisGenAntiToppT"));    
        bins =  {0,80,130,200,300,400};
    } else if (quantity == "y") {
        unbinned = fileReader->GetClone<TH1D>(filename.c_str(), "VisGenTopRapidity");
        unbinned->Add(fileReader->Get<TH1D>(filename.c_str(), "VisGenAntiTopRapidity"));
        bins = {-5, -2.5, -1.3, -0.8, -0.4, 0.0, 0.4, 0.8, 1.3, 2.5, 5};
    } else {
        cerr << "Wrong quantity\n";
        std::exit(1);
    }
    unbinned->Scale(0.5);
    
    TH1 *binned = unbinned->Rebin(bins.size()-1, "bh", &bins[0]);
    for (Int_t i=1; i<(int)bins.size(); ++i) {
        double binWidth = binned->GetBinLowEdge(i+1)-binned->GetBinLowEdge(i);
        binned->SetBinContent(i, binned->GetBinContent(i)/binWidth);
    }
    
    double genscale = 1/binned->Integral("width");
    binned->Scale(genscale);
    int start_bin = 1, stop_bin = binned->GetNbinsX(); 
    if (quantity == "y") { ++start_bin; --stop_bin; }       //remove outer bins in y
    for (int i = start_bin; i <= stop_bin; ++i)
        cout << i << " " << binned->GetBinContent(i) << "\n";
    delete binned; delete unbinned;
}
    
    
int main() {
    cout << "Pt variations\n";
    for (string i : { "0.002", "0.004", "0.006", "0.008", "-0.002", "-0.004", "-0.006", "-0.008" })
        printGenXsec("pt", "pttop" + i);

    cout << "\n\n\nRap variations\n";
    for (string i : { "0.080", "0.160", "0.240", "0.320", "-0.080", "-0.160", "-0.240", "-0.320" })
        printGenXsec("y", "ytop" + i);
        
}
