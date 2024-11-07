#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <utility>
#include <algorithm>
#include <cassert>
#include <string>

#include <TROOT.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TGaxis.h>
#include <TPad.h>
#include <TExec.h>
#include <TLine.h>
#include <TMath.h>
#include <sys/stat.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TMatrixD.h>
#include <TSystem.h>
#include <Rtypes.h>
#include <TString.h>
#include <TFile.h>

#include "utils.h"
#include "../../common/include/utils.h"
#include "../../common/include/RootFileReader.h"
#include "ttmd_fileReader.h"


void utils::addAndDelete_or_Assign(TH1*& addToThis, TH1* addThis)
{
    if (!addToThis) addToThis = addThis;
    else {
        addToThis->Add(addThis);
        delete addThis;
    }
}

void utils::drawRatio(TH1* histNumerator, TH1* histDenominator, const TH1* uncband,
               const Double_t& ratioMin, const Double_t& ratioMax, TString ratioType, TH1D* hist, TGraphAsymmErrors* graph)
{
    // check that histos have the same binning
    if(histNumerator->GetNbinsX()!=histDenominator->GetNbinsX()){
        std::cout << "error when calling drawRatio - histos have different number of bins" << std::endl;
        std::cout << "building ratio plot of " << histNumerator->GetName();
        std::cout << " and " << histDenominator->GetName() << std::endl;
        return;
    }

    // create ratio of uncertainty band
    TH1 *band = nullptr;
    if (uncband) {
        band = (TH1*)uncband->Clone("band");
        band->Divide(band);
    }


    // create ratio
    TH1* ratio = (TH1*)histNumerator->Clone();
    ratio->Divide(histDenominator);
    ratio->SetLineColor(1);
    ratio->SetLineStyle(histNumerator->GetLineStyle());
    ratio->SetMarkerColor(histNumerator->GetMarkerColor());
    ratio->SetMarkerStyle(histNumerator->GetMarkerStyle());
    ratio->SetMarkerSize(histNumerator->GetMarkerSize());
    // calculate error for ratio
        for(int bin=1; bin<=histNumerator->GetNbinsX(); bin++){
            ratio->SetBinError(bin, histNumerator->GetBinError(bin)/histDenominator->GetBinContent(bin) );
        }


    // get some values from old pad
    Int_t    logx = gPad->GetLogx();
    Double_t left = 0.20;
    Double_t right = 0.05;
    Double_t bottom = gPad->GetBottomMargin();

    // y:x size ratio for canvas
    double canvAsym = 4./3.;
    // ratio size of pad with plot and pad with ratio
    double ratioSize = 0.32;


    // create graph ratio
    TGraphAsymmErrors* ratioGR = 0;
    if(graph)
    {
        ratioGR = (TGraphAsymmErrors*)graph->Clone();
        double *y = ratioGR->GetY();
        double *yUP = ratioGR->GetEYhigh();
        double *yDOWN = ratioGR->GetEYlow();
        int n = ratioGR->GetN();
        if(n!=histDenominator->GetNbinsX())printf("Error in utils::drawRatio: number of graph points are not equal to number of bins in histo.");
        for (int i=0;i<n;i++) {
            double content = histDenominator->GetBinContent(i+1);
            y[i] = y[i]/content;
            yUP[i] = yUP[i]/content;
            yDOWN[i] = yDOWN[i]/content;
        }

        ratioSize = 0.42;
    }


    // change old pad
    gPad->SetBottomMargin(ratioSize);
    gPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(0);
    gPad->SetFillColor(10);
//     // create new pad for ratio plot
    TPad *rPad;
    rPad = new TPad("rPad","",0,0,1,ratioSize+0.001);
    rPad->SetFillStyle(0);
    rPad->SetFillColor(0);
    rPad->SetBorderSize(0);
    rPad->SetBorderMode(0);
    rPad->Draw();
    rPad->cd();
    rPad->SetLogy(0);//log
    rPad->SetLogx(logx);
    rPad->SetTicky(1);
    // configure ratio plot
    double scaleFactor = 1./(canvAsym*ratioSize);
    ratio->SetStats(kFALSE);
    ratio->SetTitle("");
    ratio->SetName("ratio");
    ratio->SetMaximum(ratioMax);
    ratio->SetMinimum(ratioMin);
    ratio->SetLineWidth(1);
    // configure axis of ratio plot
    ratio->SetTitleSize(histNumerator->GetXaxis()->GetTitleSize()*scaleFactor*1.9,"X");
    ratio->SetTitleOffset(histNumerator->GetXaxis()->GetTitleOffset()*1,"X");
    ratio->SetLabelSize(histNumerator->GetXaxis()->GetLabelSize()*scaleFactor*1.4,"X");
    //ratio->SetTitle(histNumerator->GetXaxis()->GetTitle(),"X");
    //ratio->GetXaxis()->SetNdivisions(histNumerator->GetNdivisions());
    ratio->GetYaxis()->CenterTitle();
    if(ratioType == "cp"){
        ratio->GetYaxis()->SetTitle("#frac{N_{Data}}{N_{MC}}");
    }
    else if(ratioType == "xsec"){
        ratio->GetYaxis()->SetTitle("#frac{Theory}{Data}");
    }

    ratio->GetYaxis()->SetTitleSize(histNumerator->GetYaxis()->GetTitleSize()*scaleFactor);
    ratio->SetTitleOffset(histNumerator->GetYaxis()->GetTitleOffset()/scaleFactor*0.8,"Y");
    ratio->GetYaxis()->SetLabelSize(histNumerator->GetYaxis()->GetLabelSize()*scaleFactor);
    ratio->SetLabelOffset(histNumerator->GetYaxis()->GetLabelOffset()*3.3,"Y");
    ratio->GetYaxis()->SetTickLength(0.03);
    ratio->GetYaxis()->SetNdivisions(405);
    ratio->GetXaxis()->SetRange(histNumerator->GetXaxis()->GetFirst(), histNumerator->GetXaxis()->GetLast());
    // delete axis of initial plot
    histNumerator->GetXaxis()->SetLabelSize(0);
    histNumerator->GetXaxis()->SetTitleSize(0);
    histDenominator->GetXaxis()->SetLabelSize(0);
    histDenominator->GetXaxis()->SetTitleSize(0);
    // draw ratio plot
    if(hist){
        hist->SetMaximum(ratioMax);
        hist->SetMinimum(ratioMin);
        hist->SetTitle(0);
        hist->GetXaxis()->SetLabelSize(0.08);
        hist->GetXaxis()->SetTitleSize(0.08);
        hist->GetYaxis()->SetLabelSize(histNumerator->GetYaxis()->GetLabelSize()*scaleFactor);
        hist->SetLabelOffset(histNumerator->GetYaxis()->GetLabelOffset()*3.3,"Y");
        hist->GetYaxis()->SetTickLength(0.03);

        hist->GetYaxis()->SetNdivisions(405);
        hist->Draw();
        for(int i=1;i<hist->GetXaxis()->GetNbins();++i)
        {
            double xLine = hist->GetBinLowEdge(i+1);
            TLine ln(xLine,ratioMin,xLine,ratioMax);
            ln.DrawClone("same");
        }
        ratio->DrawClone("p e1 X0 same");
        if(graph)ratioGR->Draw("P,SAME");
    }
    else
    {
        ratio->DrawClone("p e1 X0");
        if(graph)ratioGR->Draw("P,SAME");
    }
    rPad->SetTopMargin(0.0);
    rPad->SetBottomMargin(bottom*scaleFactor);
    rPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->RedrawAxis();
    //draw grid
    rPad->SetGrid(0,1);

    // draw a horizontal lines on a given histogram
    // a) at 1
//     Double_t xmin = ratio->GetXaxis()->GetXmin();
//     Double_t xmax = ratio->GetXaxis()->GetXmax();
//     TString height = ""; height += 1;
//     TF1 *f = new TF1("f", height, xmin, xmax);
//     f->SetLineStyle(1);
//     f->SetLineWidth(1);
//     f->SetLineColor(kBlack);
//     f->Draw("L same");
//     // b) at upper end of ratio pad
//     TString height2 = ""; height2 += ratioMax;
//     TF1 *f2 = new TF1("f2", height2, xmin, xmax);
//     f2->SetLineStyle(1);
//     f2->SetLineWidth(1);
//     f2->SetLineColor(kBlack);
//     f2->Draw("L same");
}

// OZ 15.01.2017 for PlotterForTUnfold, copied from Jenya Korol's code
void utils::drawRatio2d_Eban(TH1* histNumerator, TH1* histDenominator, const Double_t& ratioMin, const Double_t& ratioMax, TString ratioType ,TGraphAsymmErrors* graph1 ,TGraphAsymmErrors* graph2 )
    {
    // OZ to eliminate unused var warning and not change the function prototype
    graph1 = graph1;

    // check that histos have the same binning
    if(histNumerator->GetNbinsX()!=histDenominator->GetNbinsX()){
        std::cout << "error when calling drawRatio - histos have different number of bins" << std::endl;
        std::cout << "building ratio plot of " << histNumerator->GetName();
        std::cout << " and " << histDenominator->GetName() << std::endl;
        return;
    }

    (histDenominator->GetXaxis())->SetTitle((histNumerator->GetXaxis())->GetTitle());

    // create ratio
    TH1* ratio = (TH1*)histNumerator->Clone();
    ratio->Divide(histDenominator);
    ratio->SetLineColor(1);
    ratio->SetLineStyle(histDenominator->GetLineStyle());
    ratio->SetMarkerColor(histDenominator->GetMarkerColor());
    ratio->SetMarkerStyle(histDenominator->GetMarkerStyle());
    ratio->SetMarkerSize(histDenominator->GetMarkerSize());
    // calculate error for ratio
        for(int bin=1; bin<=histNumerator->GetNbinsX(); bin++){
            ratio->SetBinError(bin, histNumerator->GetBinError(bin)/histDenominator->GetBinContent(bin) );
            //ratio->SetBinError(bin, histDenominator->GetBinError(bin)*histNumerator->GetBinContent(bin)/pow(histDenominator->GetBinContent(bin),2) );
        }

    // get some values from old pad
    Int_t    logx = gPad->GetLogx();
    Double_t left = 0.20;
    Double_t right = 0.05;
    Double_t bottom = gPad->GetBottomMargin();

    // y:x size ratio for canvas
    double canvAsym = 4./3.;
    // ratio size of pad with plot and pad with ratio
    double ratioSize = 0.32;



    // change old pad
    gPad->SetBottomMargin(ratioSize);
    gPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(0);
    gPad->SetFillColor(10);


    // create new pad for ratio plot
    TPad *rPad;
    rPad = new TPad("rPad","",0,0,1,ratioSize+0.001);
    rPad->SetFillStyle(0);
    rPad->SetFillColor(0);
    rPad->SetBorderSize(0);
    rPad->SetBorderMode(0);
    rPad->Draw();
    rPad->cd();
    rPad->SetLogy(0);//log
    rPad->SetLogx(logx);
    rPad->SetTicky(1);
    // configure ratio plot
    double scaleFactor = 1./(canvAsym*ratioSize);
    double yOffset = histDenominator->GetYaxis()->GetTitleOffset()*0.6;
    yOffset = histDenominator->GetYaxis()->GetTitleOffset()*0.3;
    double yTitelSize = histDenominator->GetYaxis()->GetTitleSize()*scaleFactor;
    double yLabelSize = histDenominator->GetYaxis()->GetLabelSize()*scaleFactor;
    double xTitelSize = 1.2*histDenominator->GetXaxis()->GetTitleSize()*scaleFactor;
    double xLabelSize = 1.2*histDenominator->GetXaxis()->GetLabelSize()*scaleFactor;

    // delete axis of initial plot
    histNumerator->GetXaxis()->SetLabelSize(0);
    histNumerator->GetXaxis()->SetTitleSize(0);
    histDenominator->GetXaxis()->SetLabelSize(0);
    histDenominator->GetXaxis()->SetTitleSize(0);


    ratio->SetStats(kFALSE);
    ratio->SetTitle("");
    ratio->SetName("ratio");
    ratio->SetMaximum(ratioMax);
    ratio->SetMinimum(ratioMin);
    ratio->SetLineWidth(1);
    ratio->GetYaxis()->CenterTitle();
    if(ratioType == "cp"){
        ratio->GetYaxis()->SetTitle("#frac{Data}{MC}");
    }
    else if(ratioType == "xsec"){
        ratio->GetYaxis()->SetTitle("#frac{Theory}{Data}");
    }
    ratio->SetTitleOffset(yOffset*1.6,"Y");
    ratio->SetTitleSize(yTitelSize*1.6,"Y");
    ratio->SetLabelSize(yLabelSize,"Y");
    ratio->GetYaxis()->SetTickLength(0.03);
    ratio->GetYaxis()->SetNdivisions(405);
    ratio->SetTitleSize(xTitelSize*1.2,"X");
    ratio->SetLabelSize(xLabelSize,"X");
    ratio->GetXaxis()->SetTitle(histDenominator->GetXaxis()->GetTitle());

//     ratio->SetLabelOffset(histNumerator->GetYaxis()->GetLabelOffset()*3.3,"Y");

//     ratio->GetXaxis()->SetRange(histNumerator->GetXaxis()->GetFirst(), histNumerator->GetXaxis()->GetLast());
    //ratio->GetXaxis()->SetNdivisions(histNumerator->GetNdivisions());

    // draw ratio plot

    ratio->DrawClone("p e1 X0");
    graph2->Draw("2");
    rPad->SetTopMargin(0.0);
    rPad->SetBottomMargin(bottom*scaleFactor);
    rPad->SetRightMargin(right);
    rPad->SetLeftMargin(left);
    gPad->RedrawAxis();
    //draw grid
    rPad->SetGrid(0,1);
}

// OZ 15.01.2017 for PlotterForTUnfold, copied from Jenya Korol's code
void utils::drawRatio2d(TH1* histNumerator, TH1* histDenominator,
               const Double_t& ratioMin, const Double_t& ratioMax,
               TString ratioType, TH1D* hist, TGraphAsymmErrors* graph)
{
    // check that histos have the same binning
    if(histNumerator->GetNbinsX()!=histDenominator->GetNbinsX()){
        std::cout << "error when calling drawRatio - histos have different number of bins" << std::endl;
        std::cout << "building ratio plot of " << histNumerator->GetName();
        std::cout << " and " << histDenominator->GetName() << std::endl;
        return;
    }

    (histDenominator->GetXaxis())->SetTitle((histNumerator->GetXaxis())->GetTitle());

    // create ratio
    TH1* ratio = (TH1*)histNumerator->Clone();
    ratio->Divide(histDenominator);
    ratio->SetLineColor(1);
    ratio->SetLineStyle(histDenominator->GetLineStyle());
    ratio->SetMarkerColor(histDenominator->GetMarkerColor());
    ratio->SetMarkerStyle(histDenominator->GetMarkerStyle());
    ratio->SetMarkerSize(histDenominator->GetMarkerSize());
    // calculate error for ratio
        for(int bin=1; bin<=histNumerator->GetNbinsX(); bin++){
            ratio->SetBinError(bin, histNumerator->GetBinError(bin)/histDenominator->GetBinContent(bin) );
            //ratio->SetBinError(bin, histDenominator->GetBinError(bin)*histNumerator->GetBinContent(bin)/pow(histDenominator->GetBinContent(bin),2) );
        }

    // get some values from old pad
    Int_t    logx = gPad->GetLogx();
    Double_t left = 0.20;
    Double_t right = 0.05;
    Double_t bottom = gPad->GetBottomMargin();

    // y:x size ratio for canvas
    double canvAsym = 4./3.;
    // ratio size of pad with plot and pad with ratio
    double ratioSize = 0.32;


    // create graph ratio
    TGraphAsymmErrors* ratioGR = 0;
    if(graph)
    {
        ratioGR = (TGraphAsymmErrors*)graph->Clone();
        double *y = ratioGR->GetY();
        double *yUP = ratioGR->GetEYhigh();
        double *yDOWN = ratioGR->GetEYlow();
        int n = ratioGR->GetN();
        if(n!=histDenominator->GetNbinsX())printf("Error in utils::drawRatio: number of graph points are not equal to number of bins in histo.");
        for (int i=0;i<n;i++) {
            double content = histDenominator->GetBinContent(i+1);
            y[i] = y[i]/content;
            yUP[i] = yUP[i]/content;
            yDOWN[i] = yDOWN[i]/content;
        }

        ratioSize = 0.42;
    }


    // change old pad
    gPad->SetBottomMargin(ratioSize);
    gPad->SetRightMargin(right);
    gPad->SetLeftMargin(left);
    gPad->SetBorderMode(0);
    gPad->SetBorderSize(0);
    gPad->SetFillColor(10);


    // create new pad for ratio plot
    TPad *rPad;
    rPad = new TPad("rPad","",0,0,1,ratioSize+0.001);
    rPad->SetFillStyle(0);
    rPad->SetFillColor(0);
    rPad->SetBorderSize(0);
    rPad->SetBorderMode(0);
    rPad->Draw();
    rPad->cd();
    rPad->SetLogy(0);//log
    rPad->SetLogx(logx);
    rPad->SetTicky(1);
    // configure ratio plot
    double scaleFactor = 1./(canvAsym*ratioSize);
    double yOffset = histDenominator->GetYaxis()->GetTitleOffset()*0.6;
    if(!hist)yOffset = histDenominator->GetYaxis()->GetTitleOffset()*0.3;
    double yTitelSize = histDenominator->GetYaxis()->GetTitleSize()*scaleFactor;
    double yLabelSize = histDenominator->GetYaxis()->GetLabelSize()*scaleFactor;
    double xTitelSize = 1.2*histDenominator->GetXaxis()->GetTitleSize()*scaleFactor;
    double xLabelSize = 1.2*histDenominator->GetXaxis()->GetLabelSize()*scaleFactor;

    // delete axis of initial plot
    histNumerator->GetXaxis()->SetLabelSize(0);
    histNumerator->GetXaxis()->SetTitleSize(0);
    histDenominator->GetXaxis()->SetLabelSize(0);
    histDenominator->GetXaxis()->SetTitleSize(0);


    ratio->SetStats(kFALSE);
    ratio->SetTitle("");
    ratio->SetName("ratio");
    ratio->SetMaximum(ratioMax);
    ratio->SetMinimum(ratioMin);
    ratio->SetLineWidth(1);
    ratio->GetYaxis()->CenterTitle();
    if(ratioType == "cp"){
        ratio->GetYaxis()->SetTitle("#frac{Data}{MC}");
    }
    else if(ratioType == "xsec"){
        ratio->GetYaxis()->SetTitle("#frac{Theory}{Data}");
    }
    ratio->SetTitleOffset(yOffset*1.6,"Y");
    ratio->SetTitleSize(yTitelSize*1.6,"Y");
    ratio->SetLabelSize(yLabelSize,"Y");
    ratio->GetYaxis()->SetTickLength(0.03);
    ratio->GetYaxis()->SetNdivisions(405);
    ratio->SetTitleSize(xTitelSize*1.2,"X");
    ratio->SetLabelSize(xLabelSize,"X");
    ratio->GetXaxis()->SetTitle(histDenominator->GetXaxis()->GetTitle());

//     ratio->SetLabelOffset(histNumerator->GetYaxis()->GetLabelOffset()*3.3,"Y");

//     ratio->GetXaxis()->SetRange(histNumerator->GetXaxis()->GetFirst(), histNumerator->GetXaxis()->GetLast());
    //ratio->GetXaxis()->SetNdivisions(histNumerator->GetNdivisions());

    // draw ratio plot
    if(hist){
        hist->SetMaximum(ratioMax);
        hist->SetMinimum(ratioMin);
        hist->SetTitle("");
        hist->GetYaxis()->CenterTitle();
        hist->GetYaxis()->SetTitle(ratio->GetYaxis()->GetTitle());
        hist->SetTitleOffset(yOffset,"Y");
        hist->SetTitleSize(yTitelSize,"Y");
        hist->SetLabelSize(yLabelSize,"Y");
        hist->GetYaxis()->SetTickLength(0.03);
        hist->GetYaxis()->SetNdivisions(405);
        hist->SetTitleSize(xTitelSize,"X");
        hist->SetLabelSize(xLabelSize,"X");

        //hist->GetXaxis()->SetLabelSize(0.08);
        //hist->GetXaxis()->SetTitleSize(0.08);

        hist->Draw();
        for(int i=1;i<hist->GetXaxis()->GetNbins();++i)
        {
            double xLine = hist->GetBinLowEdge(i+1);
            TLine ln(xLine,ratioMin,xLine,ratioMax);
            ln.DrawClone("same");
        }
        ratio->DrawClone("p e1 X0 same");
        if(graph)ratioGR->Draw("P,SAME");
    }
    else
    {
        ratio->DrawClone("p e1 X0");
        if(graph)ratioGR->Draw("P,SAME");
    }
    rPad->SetTopMargin(0.0);
    rPad->SetBottomMargin(bottom*scaleFactor);
    rPad->SetRightMargin(right);
    rPad->SetLeftMargin(left);
    gPad->RedrawAxis();
    //draw grid
    rPad->SetGrid(0,1);
}

void utils::rebin2d(TH2D*& histOut,TH2D* histIn,TString /*name*/,TString xAxisName_,TString yAxisName_,const int rebinX_,const int rebinY_, const int nbinX_, const double* x_binsArr_, const int nbinY_, const double* y_binsArr_)
{
        if(rebinX_||rebinY_){
            histOut = (TH2D* )histIn->Rebin2D(rebinX_,rebinY_,histIn->GetName());
        }
        else
        {
            //histOut = new TH2D(name,"",nbinX_,x_binsArr_,nbinY_,y_binsArr_);
            histOut->GetXaxis()->SetTitle(xAxisName_);
            histOut->GetYaxis()->SetTitle(yAxisName_);
            double binWx = (histIn->GetXaxis())->GetBinWidth(1);
            double binWy = (histIn->GetYaxis())->GetBinWidth(1);
            for(int ix=0;ix<nbinX_;ix++){
                for(int iy=0;iy<nbinY_;iy++){
                    int binx1 = (histIn->GetXaxis())->FindBin(x_binsArr_[ix]+0.1*binWx);
                    int binx2 = (histIn->GetXaxis())->FindBin(x_binsArr_[ix+1]-0.1*binWx);
                    int biny1 = (histIn->GetYaxis())->FindBin(y_binsArr_[iy]+0.1*binWy);
                    int biny2 = (histIn->GetYaxis())->FindBin(y_binsArr_[iy+1]-0.1*binWy);
                    double content=histIn->Integral(binx1,binx2,biny1,biny2);
                    histOut->SetBinContent(ix+1,iy+1,content);
                }
            }
        }
}



TString utils::numToString(double val)
{
    std::stringstream ss;
    ss.str("");
    ss  << val;
    return (TString)ss.str();
}


TString utils::makeBinTitle(TString axisName,double x1,double x2)
{
    std::stringstream ss;
    ss.str("");
    ss  << axisName << " [" << x1 << " : " << x2 << "] ";
    return (TString)ss.str();

}

TString utils::makeTitleBins(TString plotNameUnits,std::vector<double>& v_bin,int underflow, int overflow)
{
    std::stringstream ss;
    ss.str("");
    ss  << plotNameUnits << " bins: [";
    if(underflow == 1) ss << "underflow  ";
    for(int i=0; i<(int)v_bin.size(); ++i){
        double bin = v_bin.at(i);
        if(i==0)ss << bin;
        else ss << ",  " << bin;
        //if(i!=0)ss << "  ";
        //ss << "[" << v_bin.at(i) << "  " << v_bin.at(i+1) << "]";

    }
    if(overflow == 1) ss << "  overflow";
    ss  << "]";
    return (TString)ss.str();

}



void utils::readLineToVector(const TString& file, const TString& keyWord,std::vector<double>& outVector)
{
    std::ifstream tempStream(file.Data(), std::ifstream::in);
    if (!tempStream.good()) {
        std::cerr<<"Error in utils::readLineToVector! Cannot find file with name: "<< file <<"\n...break\n"<<std::endl;
        exit(12);
    }
    while(tempStream.good()){
        std::string line;
        getline(tempStream, line);
        line.erase(0, line.find_first_not_of(" \t"));
        if (line.size() == 0 || line[0] == '#') continue;
        std::vector<TString> vWord;
        std::string word;
        for (std::stringstream ss(line); ss >> word; ){
            vWord.push_back(word);
        }
        if(vWord.at(0) == keyWord.Data())
        {
            vWord.erase(vWord.begin());
            for(auto word: vWord)outVector.push_back(word.Atof());
        }
    }
}



void utils::cat(const TString& file)
{
    std::ifstream tempStream(file.Data(), std::ifstream::in);
    if (!tempStream.good()) {
        std::cerr<<"Error in utils::readLineToVector! Cannot find file with name: "<< file <<"\n...break\n"<<std::endl;
        exit(12);
    }
    while(tempStream.good()){
        std::string line;
        getline(tempStream, line);
        std::cout << line << std::endl;
    }
}



// std::vector<double> utils::addUpDownVect(const std::vector<double>& a,const std::vector<double>& b,
//                                    const double scale, const int precision)
// {
//     std::vector<double> result;
//     for(size_t i=0;i<a.size();i++){
//         double c = 0;
//         if(a.at(i)*b.at(i) > 0)
//         {}
//         result.push_back(scale*(a.at(i)+b.at(i)));
//
//     }
//     if(precision>0)utils::setPrecision(result,precision);
//     return result;
// }



void utils::fillSysUpDown(const std::vector<double>& vDiffUp,const std::vector<double>& vDiffDown,
                                     std::vector<double>& vSystErrPos, std::vector<double>& vSystErrNeg)
{
                                        if(vSystErrPos.size()<1){
                                            vSystErrPos.assign((int)vDiffUp.size(),0);
                                            vSystErrNeg.assign((int)vDiffDown.size(),0);
                                        }
                                        else{
                                            for(size_t i=0;i<vDiffUp.size();i++){
                                                double up = vDiffUp.at(i);
                                                double down = vDiffDown.at(i);
                                                if(up*down>0){
                                                    if(up*up>down*down){
                                                        vSystErrPos.at(i) = sqrt(pow(vSystErrPos.at(i),2) + (up>=0)*pow(up,2));
                                                        vSystErrNeg.at(i) = sqrt(pow(vSystErrNeg.at(i),2) + (up<0)*pow(up,2));
                                                    }
                                                    else{
                                                        vSystErrPos.at(i) = sqrt(pow(vSystErrPos.at(i),2) + (down>=0)*pow(down,2));
                                                        vSystErrNeg.at(i) = sqrt(pow(vSystErrNeg.at(i),2) + (down<0)*pow(down,2));
                                                    }
                                                }
                                                else{
                                                    vSystErrPos.at(i) = sqrt(pow(vSystErrPos.at(i),2) + (up>=0)*pow(up,2));
                                                    vSystErrNeg.at(i) = sqrt(pow(vSystErrNeg.at(i),2) + (up<0)*pow(up,2));
                                                    vSystErrPos.at(i) = sqrt(pow(vSystErrPos.at(i),2) + (down>=0)*pow(down,2));
                                                    vSystErrNeg.at(i) = sqrt(pow(vSystErrNeg.at(i),2) + (down<0)*pow(down,2));
                                                }
                                            }
                                        }
}

TString utils::GetNameFromUpDownVariation(TString name_with_up)
{
    TString result = name_with_up;
    if (name_with_up.Contains("_UP")) result.ReplaceAll("_UP", "");
    else if (name_with_up.Contains("_DOWN")) result.ReplaceAll("_DOWN", "");

    return result;
}



void utils::addSqrVect(std::vector<double>& a, const std::vector<double>& b)
{
    if(a.size()<1)a.assign((int)b.size(),0);
    for(size_t i=0;i<a.size();i++){
        a.at(i) = sqrt(a.at(i)*a.at(i)+b.at(i)*b.at(i));

    }
}



std::vector<double> utils::addVect(const std::vector<double>& a,const std::vector<double>& b,
                                   const double scale, const int precision)
{
    std::vector<double> result;
    for(size_t i=0;i<a.size();i++)result.push_back(scale*(a.at(i)+b.at(i)));
    if(precision>0)utils::setPrecision(result,precision);
    return result;
}



std::vector<double> utils::diffVect(const std::vector<double>& a,const std::vector<double>& b,
                                    const double scale, const int precision)
{
    std::vector<double> result;
    for(size_t i=0;i<a.size();i++)result.push_back(scale*(a.at(i)-b.at(i)));
    if(precision>0)utils::setPrecision(result,precision);
    return result;
}



std::vector<double> utils::divideVect(const std::vector<double>& a,const std::vector<double>& b,
                                      const double scale, const int precision)
{
    std::vector<double> result;
    for(size_t i=0;i<a.size();i++)result.push_back(scale*(a.at(i)/b.at(i)));
    if(precision>0)utils::setPrecision(result,precision);
    return result;
}



std::vector<double> utils::multiplyVect(const std::vector<double>& a,const std::vector<double>& b,
                                      const double scale, const int precision)
{
    std::vector<double> result;
    for(size_t i=0;i<a.size();i++)result.push_back(scale*(a.at(i)*b.at(i)));
    if(precision>0)utils::setPrecision(result,precision);
    return result;
}


std::vector<double> utils::relativeDiff(const std::vector<double>& a,const std::vector<double>& b,
                                        const double scale, const int precision)
{
    std::vector<double> result;
    for(size_t i=0;i<a.size();i++)result.push_back(scale*((a.at(i)-b.at(i))/b.at(i)));
    if(precision>0)utils::setPrecision(result,precision);
    return result;
}



void utils::setPrecision(std::vector<double>& a,int precision)
{
    for(size_t i=0;i<a.size();i++){
       a.at(i) = round(a.at(i) * pow10(precision))/pow10(precision);
    }
}


TCanvas* utils::setCanvas()
{
    TCanvas* canvas = new TCanvas("","",600,600);
    canvas->Clear();
    canvas->SetName("");
    canvas->SetTitle("");
    return canvas;
}


TLegend* utils::setLegend(const double x1,const double y1,const double x2,const double y2)
{
    TLegend* legend = new TLegend(x1,y1,x2,y2);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetX1NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength() - 0.25);
    legend->SetY1NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength() - 0.05 - legend->GetBorderSize()*0.04);
    legend->SetX2NDC(1.0 - gStyle->GetPadRightMargin() - gStyle->GetTickLength());
    legend->SetY2NDC(1.0 - gStyle->GetPadTopMargin()  - gStyle->GetTickLength());
    legend->Clear();
    return legend;
}


const std::string utils::DATA_PATH_DILEPTONIC()
{
std::string result(common::CMSSW_BASE());
result.append("/src/TopAnalysis/Configuration/analysis/diLeptonic/data");
return result;
}


std::vector<std::pair<TString, TString> > utils::nameStepPairs(const char* filename,
                                                             const char* objectNameFragment,
                                                             const std::vector<TString>& selectedSteps)
{
    std::vector<std::pair<TString, TString> > result;
    std::vector<TString> v_step;
    std::stringstream ss_step;
    std::vector<TString> v_unselectedStep;
    std::stringstream ss_unselected;

    // Search for all steps
    RootFileReader* fileReader(RootFileReader::getInstance());
    const std::vector<TString>& v_treeName = fileReader->findObjects(filename, objectNameFragment, false);
    for(std::vector<TString>::const_iterator i_objectName = v_treeName.begin(); i_objectName != v_treeName.end(); ++i_objectName){
        const TString& step = utils::extractSelectionStepAndJetCategory(*i_objectName);

        // Reject steps in case only selected steps should be chosen
        if(selectedSteps.size() && std::find(selectedSteps.begin(), selectedSteps.end(), step) == selectedSteps.end()){
            if(std::find(v_unselectedStep.begin(), v_unselectedStep.end(), step) == v_unselectedStep.end()){
                ss_unselected<<step<<", ";
                v_unselectedStep.push_back(step);
            }
            continue;
        }

        // Add selection steps
        if(std::find(v_step.begin(), v_step.end(), step) == v_step.end()){
            ss_step<<step<<", ";
            v_step.push_back(step);
        }
        result.push_back(std::make_pair(*i_objectName, step));
    }

    std::cout<<"Found chosen selection steps:\n"<<ss_step.str()<<std::endl;
    if(selectedSteps.size() && ss_unselected.str()!="") std::cout<<"Found rejected selection steps:\n"<<ss_unselected.str()<<std::endl;
    return result;
}



TString utils::stepName(const TString& stepShort, const int& category)
{
    std::ostringstream result;
    result<<"_step"<<stepShort;
    if(category>=0){
        result<<"_cate"<<category;
    }
    return result.str().c_str();
}



TString utils::stepName(const TString& stepShort, const std::vector<int>& v_category)
{
    std::ostringstream result;
    result<<"_step"<<stepShort;
    std::vector<int> v_clone = v_category;
    std::sort(v_clone.begin(), v_clone.end());
    for(const auto& category : v_clone){
        result<<"_cate"<<category;
    }
    return result.str().c_str();
}



TString utils::extractSelectionStep(const TString& name)
{
    TString result(name);
    // Remove e.g. things like .root ??
    result = helper::searchFragmentByToken(result, "step", ".");
    // Further separation by "_" to find string containing step
    result = helper::searchFragmentByToken(result, "step", "_");
    //std::cout<<"The extracted selection step is (step/histogram name): "<<result<<" / "<<name<<std::endl;
    return result;
}



TString utils::extractJetCategory(const TString& name)
{
    TString result(name);
    // Remove e.g. things like .root ??
    result = helper::searchFragmentByToken(result, "cate", ".", true);
    // Further separation by "_" to find string containing step
    result = helper::searchFragmentByToken(result, "cate", "_", true);
    //std::cout<<"The extracted category bin is (bin/histogram name): "<<result<<" / "<<name<<std::endl;
    return result;
}



TString utils::extractSelectionStepAndJetCategory(const TString& name)
{
    const TString step(extractSelectionStep(name));
    const TString cate(extractJetCategory(name));
    const TString result = step+cate;
    return result;
}



TString utils::helper::searchFragmentByToken(const TString& name, const TString& searchPattern, const TString& token, const bool allowMultiple)
{
    TString result;
    TObjArray* a_nameFragment = TString(name).Tokenize(token);
    bool alreadyFound(false);
    for(Int_t iFragment= 0; iFragment < a_nameFragment->GetEntriesFast(); ++iFragment){
        const TString& fragment = a_nameFragment->At(iFragment)->GetName();
        if(fragment.Contains(searchPattern)){
            if(!allowMultiple && alreadyFound){
                std::cerr<<"ERROR in searchFragmentByToken()! Ambiguous string, name contains search pattern \""<<searchPattern
                         <<"\" in more than one fragment: "<<name
                         <<"\n...cannot extract fragment, so break\n";
                exit(33);
            }
            alreadyFound = true;
            result.Append("_").Append(fragment);
        }
    }
    return result;
}



TString utils::assignFolder(const char* baseDir, const TString& channel, const TString& systematic)
{
    TString path("");

    // Create all subdirectories contained in baseDir
    TObjArray* a_subDir = TString(baseDir).Tokenize("/");
    for(Int_t iSubDir = 0; iSubDir < a_subDir->GetEntriesFast(); ++iSubDir){
        const TString& subDir = a_subDir->At(iSubDir)->GetName();
        path.Append(subDir);
        path.Append("/");
        gSystem->MakeDirectory(path);
    }

    // Create subdirectories for systematic and channel
    path.Append(systematic);
    path.Append("/");
    gSystem->MakeDirectory(path);
    path.Append(channel);
    path.Append("/");
    gSystem->MakeDirectory(path);

    return path;
}



TString utils::accessFolder(const char* baseDir, const TString& channel,
                            const TString& systematic, const bool allowNonexisting)
{
    // Build directory path
    TString path(baseDir);
    path.Append("/");
    path.Append(systematic);
    path.Append("/");
    path.Append(channel);
    path.Append("/");

    // Check if directory really exists
    if(!gSystem->OpenDirectory(path)){
        if(allowNonexisting){
            // It is allowed to request a folder which does not exist, so return empty string silently
            return "";
        }
        else{
            std::cerr<<"ERROR! Request to access directory is not possible, because it does not exist. Directory name: "<<path
                     <<"\n...break\n"<<std::endl;
            exit(237);
        }
    }

    return path;
}



void utils::setPlotNames(const std::string& nameListFile, std::vector<std::vector<TString>>& varNames)
{
       std::ifstream nameListStream(nameListFile.data(), std::ifstream::in);
       if (!nameListStream.good()) {
        std::cerr<<"Error in utils::::setVarNames! Cannot find NameList with name: "<<nameListFile<<"\n...break\n"<<std::endl;
        exit(12);
       }

    // Loop over all variables in NameList
    while(nameListStream.good()){
        // Read NameList-File
        std::string line;
        getline(nameListStream, line);
        //remove leading whitespace
        line.erase(0, line.find_first_not_of(" \t"));

        //skip empty lines and/or comments
        if (line.size() == 0 || line[0] == '#') continue;

        // Loop over words in NameList-File line
        std::vector<TString> vect;
        std::string word;
        for (std::stringstream ss(line); ss >> word; ){
            vect.push_back(word);
        }
	varNames.push_back(vect);
    }
}

void utils::selectIndicesFromFlags(std::vector<int>& v_indices, const std::vector<bool>& v_flags){

  std::vector<int> result;

  for(const int index : v_indices){

    if(index >= int(v_flags.size()))
    {
      std::ostringstream ss;
      ss << "utils::selectIndicesFromFlags -- index (" << index << ") out-of-range in vector of flags";

      throw std::runtime_error(ss.str());
    }

    if(v_flags.at(index)){

      result.emplace_back(index);
    }
  }

  v_indices.clear();
  v_indices = result;
}


std::vector<bool> utils::flagsForPassingCutInPtRange(const std::vector<int>& v_ints, const int cutValue, const std::vector<LV>& v_p4, const float min_pt, const float max_pt){

  std::vector<bool> result;

  if((min_pt >= 0.) and (max_pt >= 0.) and (max_pt <= min_pt)){

    std::ostringstream ss;
    ss << "utils::flagsForPassingCutInPtRange -- logic error: inconsistent values for min-pT (" << min_pt << ") and max-pT (" << max_pt << ")";

    throw std::runtime_error(ss.str());
  }

  if(v_ints.size() != v_p4.size()){

    std::ostringstream ss;
    ss << "utils::flagsForPassingCutInPtRange -- logic error: vectors of ints and 4-momenta differ in size (" << v_ints.size() << " != " << v_p4.size() << ")";

    throw std::runtime_error(ss.str());
  }

  for(uint idx=0; idx<v_p4.size(); ++idx){

    const auto& p4 = v_p4.at(idx);

    bool ignore_cut(false);
    if((min_pt >= 0.) or (max_pt >= 0.)){

      if(((min_pt >= 0.) and (p4.pt() < min_pt)) or ((max_pt >= 0.) and (p4.pt() > max_pt))){

        ignore_cut = true;
      }
    }

    const bool passes(ignore_cut or (v_ints.at(idx) >= cutValue));
    result.emplace_back(passes);
  }

  return result;
}

std::vector<bool> utils::flagsForPassingJetVetoMap(const utils::JetVetoMap* vetoMap, const std::vector<LV>& v_p4, const float& cutValue){

  std::vector<bool> result;

  for(uint idx=0; idx<v_p4.size(); ++idx){
    const auto& p4 = v_p4.at(idx);
    result.emplace_back(vetoMap->veto(p4, cutValue));
  }

  return result;
}


// JetVetoMap --------------------------
utils::JetVetoMap::JetVetoMap(const std::string& name, const bool verbose): name_(name), verbose_(verbose){}

utils::JetVetoMap::~JetVetoMap(){}

void utils::JetVetoMap::initMap(utils::JetVetoMap::Map& map, const std::string& tfile_path, const std::string& th2_key, const std::string& th2_format){
  map = utils::JetVetoMap::Map(tfile_path, th2_key, th2_format, verbose_);
}

void utils::JetVetoMap::initMap(const std::string& tfile_path, const std::string& th2_key, const std::string& th2_format){

  this->initMap(map_, tfile_path, th2_key, th2_format);
}

utils::JetVetoMap::Map::Map(const std::string& tfile_path, const std::string& th2_key, const std::string& th2_format, const bool verbose){

  TFile* tfile = TFile::Open(tfile_path.c_str());
  if((tfile == nullptr) or tfile->IsZombie())
    {
      throw std::runtime_error("JetVetoMap::Map::Map -- failed to open TFile under path "+tfile_path);
    }
  else if(verbose)
    {
      std::cout << "JetVetoMap::Map::Map -- opened input TFile under path = " << tfile_path << std::endl;
    }

  th2_ = dynamic_cast<TH2*>(tfile->Get(th2_key.c_str()));
  if(th2_ == nullptr)
    {
      throw std::runtime_error("JetVetoMap::Map::Map -- target TH2 object (key="+th2_key+") not found in input file "+tfile_path);
    }
  else if(verbose)
    {
      std::cout << "JetVetoMap::Map::Map -- accessed TH2* object, key = " << th2_key << std::endl;
    }

  th2_->SetDirectory(0);

  tfile->Close();

  if(verbose) std::cout << "JetVetoMap::Map::Map -- accessed TH2* object, th2_format = " << th2_format << std::endl;
  if(verbose) std::cout << "JetVetoMap::Map::Map -- accessed TH2* object, th2_format_ = " << th2_format_ << std::endl;

  th2_format_ = th2_format;
  if(verbose) std::cout << "JetVetoMap::Map::Map -- accessed TH2* object, th2_format_ = " << th2_format_ << std::endl;
  if(   (th2_format_ != "phi_vs_abseta")
	&& (th2_format_ != "phi_vs_eta"   )
	&& (th2_format_ != "abseta_vs_phi")
	&& (th2_format_ !=    "eta_vs_phi")
	)
    {
      throw std::runtime_error("JetVetoMap::Map::Map -- unrecognized key for TH2 format: "+th2_format_);
    }
}

void utils::JetVetoMap::Map::getValues(float& binContent, const LV& p4) const {

  float valX(0.), valY(0.);
  if     (th2_format_ == "phi_vs_abseta"){ valX = p4.Phi()       ; valY = fabs(p4.Eta()); }
  else if(th2_format_ == "phi_vs_eta"   ){ valX = p4.Phi()       ; valY =      p4.Eta() ; }
  else if(th2_format_ == "abseta_vs_phi"){ valX = fabs(p4.Eta()); valY = p4.Phi()       ; }
  else if(th2_format_ ==    "eta_vs_phi"){ valX =      p4.Eta() ; valY = p4.Phi()       ; }
  else
    {
      throw std::runtime_error("JetVetoMap::Map::getValues -- unrecognized key for TH2 format: "+th2_format_);
    }

  if(th2_ == nullptr){

    throw std::runtime_error("JetVetoMap::Map::getValues -- null pointer to TH2 object");
  }

  int binX(th2_->GetXaxis()->FindBin(valX));
  if(binX > th2_->GetXaxis()->GetNbins()){
    binX = th2_->GetXaxis()->GetNbins();
  }
  else if(binX == 0){
    binX = 1;
  }

  int binY(th2_->GetYaxis()->FindBin(valY));
  if(binY > th2_->GetYaxis()->GetNbins()){
    binY = th2_->GetYaxis()->GetNbins();
  }
  else if(binY == 0){
    binY = 1;
  }

  binContent = th2_->GetBinContent(binX, binY);
}

float utils::JetVetoMap::veto(const LV& p4, const float& cutValue) const {

  float val(0.);
  map_.getValues(val, p4);
  bool passesCut = val < cutValue ? true : false;

  if(verbose_){
    std::cout << std::endl;
    std::cout << name_ << "::veto -- jet eta              = " << p4.Eta() << std::endl;
    std::cout << name_ << "::veto -- jet phi              = " << p4.Phi() << std::endl;
    std::cout << name_ << "::veto -- jet passes veto     = " << passesCut << std::endl;
    std::cout << name_ << "::veto -- value TH2 bin content   = " << val << std::endl;
  }

  return passesCut;
}
// ------------------------------------------------------


void utils::ttmd::HCovToFArrays(const TH1* h, const TH2* hCov, float* d, float* dCovDiag, float* dCovUndiag)
{
  int n = h->GetNbinsX();
  int c = 0;
  for(int b = 0; b < n; b++)
  {
    d[b] = h->GetBinContent(b + 1);
    dCovDiag[b] = hCov->GetBinContent(b + 1, b + 1);
    for(int bb = 0; bb < b; bb++)
      dCovUndiag[c++] = hCov->GetBinContent(b + 1, bb + 1);
  }
}



void utils::ttmd::FArraysToHCov(const int n, const float* d, const float* dCovDiag, const float* dCovUndiag, TH1* h, TH2* hCov)
{
  int c = 0;
  for(int b = 0; b < n; b++)
  {
    h->SetBinContent(b + 1, d[b]);
    h->SetBinError(b + 1, TMath::Sqrt(dCovDiag[b]));
    hCov->SetBinContent(b + 1, b + 1, dCovDiag[b]);
    //hCov->SetBinError(b + 1, b + 1, 0.0);
    if(dCovUndiag)
    {
      for(int bb = 0; bb < b; bb++)
      {
        hCov->SetBinContent(b + 1, bb + 1, dCovUndiag[c]);
        hCov->SetBinContent(bb + 1, b + 1, dCovUndiag[c]);
        c++;
      }
    }
  }
}



std::pair<double, double> utils::ttmd::GetMinMaxFromTH1(const TH1* h)
{
  double min = h->GetBinContent(1);
  double max = h->GetBinContent(1);
  for(int b = 1; b <= h->GetNbinsX(); b++)
  {
    double content = h->GetBinContent(b);
    if(min > content)
      min = content;
    if(max < content)
      max = content;
  }
  return std::pair<double, double>(min, max);
}



std::pair<double, double> utils::ttmd::GetMinMaxFromTH1DVector(const std::vector<TH1D*>& vH)
{
  std::pair<double, double> result;
  for(unsigned int h = 0; h < vH.size(); h++)
  {
    std::pair<double, double> current = GetMinMaxFromTH1(vH[h]);
    if(result.first > current.first)
      result.first = current.first;
    if(result.second < current.second)
      result.second = current.second;
  }
  return result;
}



std::pair<double, double> utils::ttmd::GetMinMaxFromTH1DVector2D(const std::vector<std::vector<TH1D*> >& vH)
{
  std::pair<double, double> result;
  for(unsigned int h = 0; h < vH.size(); h++)
  {
    std::pair<double, double> current = GetMinMaxFromTH1DVector(vH[h]);
    result.first = std::min(result.first, current.first);
    result.second = std::max(result.second, current.second);
  }
  return result;
}



TH2D* utils::ttmd::MakeProbMatrix(const TH2D* hResp)
{
  TH2D* hProb = (TH2D*) hResp->Clone();
  for(int g = 0; g <= hResp->GetNbinsY() + 1; g++)
  {
    double tot = hResp->Integral(0, hResp->GetNbinsX() + 1, g, g);
    if(tot == 0.0)
      tot = 1.0;
    for(int r = 0; r <= hResp->GetNbinsX() + 1; r++)
      hProb->SetBinContent(r, g, hResp->GetBinContent(r, g) / tot);
  }
  return hProb;
}



TH2D* utils::ttmd::MakeCorrMatrix(const TH2D* hCov)
{
  TH2D* hCorr = (TH2D*) hCov->Clone();
  for(int i = 0; i < hCov->GetNbinsX(); i++)
    for(int j = 0; j < hCov->GetNbinsY(); j++)
    {
      double norm = TMath::Sqrt(hCov->GetBinContent(i + 1, i + 1) * hCov->GetBinContent(j + 1, j + 1));
      hCorr->SetBinContent(i + 1, j + 1, hCov->GetBinContent(i + 1, j + 1) / norm);
    }
  return hCorr;
}



TH2D* utils::ttmd::MakeCovMatrixFromCorr(const TH2D* hCor, const TH1D* hVec)
{
  TH2D* hCov = (TH2D*) hCor->Clone();
  for(int i = 0; i < hCov->GetNbinsX(); i++)
    for(int j = 0; j < hCov->GetNbinsY(); j++)
    {
      double cov = hVec->GetBinError(i + 1) * hVec->GetBinError(j + 1) * hCor->GetBinContent(i + 1, j + 1);
      //printf("cov = %f   %f   %f\n", cov, hCor->GetBinContent(i + 1, j + 1), hVec->GetBinError(i + 1));
      hCov->SetBinContent(i + 1, j + 1, cov);
    }
  return hCov;
}



void utils::ttmd::DivideSquareWH(const TCanvas* fCanvas, const int n, int& w, int& h)
{
  if (fCanvas->GetWindowWidth() > fCanvas->GetWindowHeight()) {
    w = TMath::Ceil(TMath::Sqrt(n));
    h = TMath::Floor(TMath::Sqrt(n));
    if (w*h < n) w++;
  } else {
    h = TMath::Ceil(TMath::Sqrt(n));
    w = TMath::Floor(TMath::Sqrt(n));
    if (w*h < n) h++;
  }
}



void utils::ttmd::SetErrorsToZero(TH1* h)
{
  for(int b = 1; b <= h->GetNbinsX(); b++)
    h->SetBinError(b, 0.0);
}



double utils::ttmd::Chi2(const TH1* hData, const TH2* covmat, const TH1* hGen, int skip/* = -1*/)
{
  //hGen->Print("all");
  //covmat->Print("all");
  //hData->Print("all");
  //throw;
  int n = hData->GetNbinsX();
  if(skip)
    n -= 1;
  TMatrixD res(1, n);
  for(int i = 0; i < n; i++)
    res(0, i) = hData->GetBinContent(i + 1) - ((hGen) ? hGen->GetBinContent(i + 1) : 0.0);
  //res.Print("all");
  //hData->Print("all");
  //hGen->Print("all");
  TMatrixD resT = res;
  resT.T();
  TMatrixD cov(n, n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
      cov(i, j) = covmat->GetBinContent(i + 1, j + 1);
  cov.Invert();
  double chi2 = (res * cov * resT)(0, 0);
  //hData->Print("all");
  //hGen->Print("all");
  //printf("CHI2 %.6f\n", chi2);
  return chi2;
}



double utils::ttmd::Chi2(const TH1D* hData, const TH2D* covmat, const TH1* hGen, const std::vector<int>& vSkip)
{
  if(vSkip.size() == 1 && vSkip[0] == -1)
    return Chi2(hData, covmat, hGen, 1);
  assert(hData->GetNbinsX() == hGen->GetNbinsX());
  assert(hData->GetNbinsX() == covmat->GetNbinsX());
  assert(hData->GetNbinsX() == covmat->GetNbinsY());
  assert(hData->GetNbinsX() > vSkip.size());
  assert(hGen->GetNbinsX() > vSkip.size());
  assert(covmat->GetNbinsX() > vSkip.size());
  TH1D* hDataNew = new TH1D("", "", hData->GetNbinsX() - vSkip.size(), 0.0, 1.0);
  TH1D* hGenNew = new TH1D("", "", hData->GetNbinsX() - vSkip.size(), 0.0, 1.0);
  TH2D* covmatNew = new TH2D("", "", covmat->GetNbinsX() - vSkip.size(), 0.0, 1.0, covmat->GetNbinsY() - vSkip.size(), 0.0, 1.0);
  int bxNew = 0;
  for(int bx = 1; bx <= hData->GetNbinsX(); bx++)
  {
    if(std::find(vSkip.begin(), vSkip.end(), bx) != vSkip.end())
      continue;
    bxNew++;
    hDataNew->SetBinContent(bxNew, hData->GetBinContent(bx));
    hGenNew->SetBinContent(bxNew, hGen->GetBinContent(bx));
    int byNew = 0;
    for(int by = 1; by <= hData->GetNbinsX(); by++)
    {
      if(std::find(vSkip.begin(), vSkip.end(), by) != vSkip.end())
        continue;
      byNew++;
      covmatNew->SetBinContent(bxNew, byNew, covmat->GetBinContent(bx, by));
    }
    assert(byNew == hData->GetNbinsX() - vSkip.size());
  }
  assert(bxNew == hData->GetNbinsX() - vSkip.size());
  double chi2 = Chi2(hDataNew, covmatNew, hGenNew);
  delete hDataNew;
  delete hGenNew;
  delete covmatNew;
  return chi2;
}



TH1D* utils::ttmd::Convolute(const TH1* hGen, const TH2* hResp)
{
  TH1D* hConv = new TH1D("", "", hResp->GetNbinsX(), 0.5, hResp->GetNbinsX() + 0.5);
  for(int bx = 0; bx < hResp->GetNbinsX(); bx++)
  {
    double sum = 0;
    for(int by = 0; by < hResp->GetNbinsY(); by++)
      sum += hResp->GetBinContent(by + 1, bx + 1) * hGen->GetBinContent(by + 1);
    hConv->SetBinContent(bx + 1, sum);
  }
  return hConv;
}



TH2D* utils::ttmd::MakeCovTH2(const TH1* hData)
{
  int n = hData->GetNbinsX();
  TH2D* covmat = new TH2D("", "", n, hData->GetXaxis()->GetXbins()->GetArray(), n, hData->GetXaxis()->GetXbins()->GetArray());
  for(int i = 0; i < n; i++)
    covmat->SetBinContent(i + 1, i + 1, TMath::Power(hData->GetBinError(i + 1), 2.0));
  return covmat;
}



TGraphAsymmErrors* utils::ttmd::DrawAsGraph(const std::pair<TH1D*, TH1D*> pairUD, const TH1D* h, const bool flagUnc/* = true*/, const double offset/* = 0.5*/, const TString& option/* = "ZP0"*/)
{
  TGraphAsymmErrors* g = new TGraphAsymmErrors;
  for(int b = 0; b < h->GetNbinsX(); b++)
  {
    double x = h->GetBinLowEdge(b + 1) + offset * h->GetBinWidth(b + 1);
    double y = h->GetBinContent(b + 1);
    double uncLowY = pairUD.second ? -1 * pairUD.second->GetBinContent(b + 1) : h->GetBinError(b + 1);
    double uncHighY = pairUD.first ? pairUD.first->GetBinContent(b + 1) : h->GetBinError(b + 1);
    g->SetPoint(b, x, y);
    g->SetPointError(b, 0.0, 0.0, uncLowY, uncHighY);
  }
  g->SetLineColor(h->GetLineColor());
  g->SetMarkerColor(h->GetMarkerColor());
  g->SetMarkerStyle(h->GetMarkerStyle());
  g->SetLineStyle(h->GetLineStyle());
  g->SetMarkerSize(h->GetMarkerSize());
  TString drawOption = option;
  if(!flagUnc)
    drawOption += "X";
  g->Draw(drawOption);
  return g;
}



TGraphAsymmErrors* utils::ttmd::DrawAsGraph(const TH1D* h, const bool flagUnc/* = true*/, const double offset/* = 0.5*/, const TString& option/* = "ZP0"*/)
{
  const std::pair<TH1D*, TH1D*> pairUD = std::make_pair<TH1D*, TH1D*>(NULL, NULL);
  return DrawAsGraph(pairUD, h, flagUnc, offset, option);
}



void utils::ttmd::AddVarInQuadrature(double& u, double&d, const double& var)
{
  if(var > 0.0)
    u = TMath::Sqrt(TMath::Power(u, 2.0) + TMath::Power(var, 2.0));
  else
    d = -1 * TMath::Sqrt(TMath::Power(d, 2.0) + TMath::Power(var, 2.0));
}



void utils::ttmd::AddVarInQuadrature(double& u, double&d, const std::vector<double>& var)
{
  double max = *std::max_element(var.begin(), var.end());
  if(max > 0.0)
    u = TMath::Sqrt(TMath::Power(u, 2.0) + TMath::Power(max, 2.0));
  double min = *std::min_element(var.begin(), var.end());
  if(min < 0.0)
    d = -1 * TMath::Sqrt(TMath::Power(d, 2.0) + TMath::Power(min, 2.0));
}



void utils::ttmd::AddVarInQuadrature(std::pair<TH1D*, TH1D*> hud, const TH1D* hv)
{
  double* u = hud.first->GetArray();
  double* d = hud.second->GetArray();
  const double* v = hv->GetArray();
  for(int b = 1; b <= hv->GetNbinsX(); b++)
    utils::ttmd::AddVarInQuadrature(u[b], d[b], v[b]);
}



void utils::ttmd::AddVarInQuadrature(std::pair<TH1D*, TH1D*> hud, std::pair<TH1D*, TH1D*> hv)
{
  assert(*std::min_element(hud.first->GetArray(), hud.first->GetArray() + hud.first->GetNbinsX()) >= 0.0);
  assert(*std::min_element(hv.first->GetArray(), hv.first->GetArray() + hv.first->GetNbinsX()) >= 0.0);
  assert(*std::max_element(hud.second->GetArray(), hud.second->GetArray() + hud.second->GetNbinsX()) <= 0.0);
  assert(*std::max_element(hv.second->GetArray(), hv.second->GetArray() + hv.second->GetNbinsX()) <= 0.0);

  utils::ttmd::AddVarInQuadrature(hud, hv.first);
  utils::ttmd::AddVarInQuadrature(hud, hv.second);
}



void utils::ttmd::AddVarToEnvelope(double& u, double&d, const double& var)
{
  if(var > u)
    u = var;
  else if(var < d)
    d = var;
}



void utils::ttmd::AddVarToEnvelope(double& u, double&d, const std::vector<double>& var)
{
  double max = *std::max_element(var.begin(), var.end());
  if(max > u)
    u = max;
  double min = *std::min_element(var.begin(), var.end());
  if(min < d)
    d = min;
}



void utils::ttmd::AddVarToEnvelope(std::pair<TH1D*, TH1D*> hud, const TH1D* hv)
{
  double* u = hud.first->GetArray();
  double* d = hud.second->GetArray();
  const double* v = hv->GetArray();
  for(int b = 1; b <= hv->GetNbinsX(); b++)
    utils::ttmd::AddVarToEnvelope(u[b], d[b], v[b]);
}



void utils::ttmd::AddVarToUncTH1D(std::vector<std::vector<TH1D*> >& vSysUnc, std::vector<std::vector<TH1D*> > vVarUnc)
{
  if(vSysUnc.size() != 2)
    throw std::runtime_error(TString::Format("Error in AddVarUncToSysUnc(): vSysUnc.size() != 2").Data());
  if(vSysUnc[0].size() != vVarUnc[0].size())
    throw std::runtime_error(TString::Format("Error in AddVarUncToSysUnc(): vVarUnc.size() = %ld != %ld = vSysUnc[0].size()", vVarUnc[0].size(), vSysUnc[0].size()).Data());
  for(unsigned int a = 0; a < vSysUnc[0].size(); a++)
  {
    for(int v = 0; v < 2; v++)
      if(!vSysUnc[v][a])
      {
        //printf("creating vSysUnc[%d][%d]\n", v, a);
        vSysUnc[v][a] = (TH1D*)vVarUnc[v][a]->Clone();
        vSysUnc[v][a]->Reset();
      }
    for(int b = 0; b < vSysUnc[0][a]->GetNbinsX(); b++)
    {
      double u = vSysUnc[0][a]->GetBinContent(b + 1);
      double d = vSysUnc[1][a]->GetBinContent(b + 1);
      std::vector<double> var(5);
      for(unsigned int v = 0; v < vVarUnc.size(); v++)
        var.push_back(vVarUnc[v][a]->GetBinContent(b + 1));
      utils::ttmd::AddVarInQuadrature(u, d, var);
      vSysUnc[0][a]->SetBinContent(b + 1, u);
      vSysUnc[1][a]->SetBinContent(b + 1, d);
    }
  }
}



// Jacobian_ij = (delta_ij * I - a_i) / I^2, I = sum_i a_i
void utils::ttmd::Normalise(TH1* h, TH2* covmat)
{
  int n = h->GetNbinsX();
  if(n != covmat->GetNbinsX() || n != covmat->GetNbinsY())
  {
    printf("Error in Normalise() incosistent input h->GetNbinsX() = %d covmat->GetNbinsX() = %d covmat->GetNbinsY() = %d\n", h->GetNbinsX(), covmat->GetNbinsX(), covmat->GetNbinsY());
    throw;
  }
  double integral = h->Integral();
  double integral2 = integral * integral;

  TMatrixD matrixG(n, n);
  TMatrixD matrixCov(n, n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
    {
      matrixCov(i, j) = covmat->GetBinContent(i + 1, j + 1);
      if(i == j)
        matrixG(i, i) = (integral - h->GetBinContent(i + 1)) / integral2;
      else
        matrixG(i, j) = -1 * h->GetBinContent(i + 1) / integral2;
    }
  TMatrixD matrixGT = matrixG;
  matrixGT.T();
  TMatrixD res = matrixG * matrixCov * matrixGT;

  for(int i = 0; i < n; i++)
  {
    h->SetBinContent(i + 1, h->GetBinContent(i + 1) / integral);
    h->SetBinError(i + 1, TMath::Sqrt(res(i, i)));
    for(int j = 0; j < n; j++)
      covmat->SetBinContent(i + 1, j + 1, res(i, j));
  }
}



bool utils::ttmd::IsEqual(const double val1, const double val2, const double eps/* = 1e-6*/)
{
  return (TMath::Abs(val1 - val2) < eps || TMath::Abs((val1 - val2) / val1) < eps);
}


//std::vector<TString> utils::ttmd::gTCanvasExt = { "pdf", "png" };
std::vector<TString> utils::ttmd::gTCanvasExt = {"pdf", "eps"};
void utils::ttmd::SaveCanvas(const TVirtualPad* c, const TString& baseName)
{
  CheckFile(baseName);
  bool flagWasEmpty = false;
  if(utils::ttmd::gTCanvasExt.size() == 0)
  {
    utils::ttmd::gTCanvasExt.push_back("pdf");
    flagWasEmpty = true;
  }
  for(unsigned int i = 0; i < utils::ttmd::gTCanvasExt.size(); i++)
    c->SaveAs(baseName + "." + utils::ttmd::gTCanvasExt[i]);
  if(flagWasEmpty)
    utils::ttmd::gTCanvasExt.clear();
}



void utils::ttmd::ScaleAxisFonts(TAxis* axis, const double scale)
{
  axis->SetTitleSize(scale * axis->GetTitleSize());
  axis->SetLabelSize(scale * axis->GetLabelSize());
}



void utils::ttmd::ScaleHistoFonts(TH1* h, const double scale)
{
  TAxis* xaxis = h->GetXaxis();
  utils::ttmd::ScaleAxisFonts(xaxis, scale);
  TAxis* yaxis = h->GetYaxis();
  utils::ttmd::ScaleAxisFonts(yaxis, scale);
}



void utils::ttmd::ScaleHistoFonts(TGraph* h, const double scale)
{
  TAxis* xaxis = h->GetXaxis();
  utils::ttmd::ScaleAxisFonts(xaxis, scale);
  TAxis* yaxis = h->GetYaxis();
  utils::ttmd::ScaleAxisFonts(yaxis, scale);
}



void utils::ttmd::SetAxisFonts(TAxis* axis, const int font/* = 63*/, const int size/* = 20 */)
{
  axis->SetTitleFont(font);
  axis->SetTitleSize(size);
  axis->SetLabelFont(font);
  axis->SetLabelSize(size);
}



void utils::ttmd::SetHistoAxisFonts(TH1* h, const int font/* = 63*/)
{
  TAxis* xaxis = h->GetXaxis();
  SetAxisFonts(xaxis, font);
  TAxis* yaxis = h->GetYaxis();
  SetAxisFonts(yaxis, font);
}



void utils::ttmd::SetHistoAxisFonts(TGraph* h, const int font/* = 63*/)
{
  TAxis* xaxis = h->GetXaxis();
  SetAxisFonts(xaxis, font);
  TAxis* yaxis = h->GetYaxis();
  SetAxisFonts(yaxis, font);
}



void utils::ttmd::AdjustDividedCanvasHistograms(TH2D* hr, const int p, const int npads, const int w/*, const int h*/, const double scaleOffset/* = 1.0*/)
{
  if((p % w) == 0 && (npads - p) <= w)
  {
    ScaleHistoFonts(hr, 0.71);
    hr->GetXaxis()->SetLabelOffset(hr->GetXaxis()->GetLabelOffset() + 0.04 * scaleOffset);
    hr->GetXaxis()->SetTitleOffset(hr->GetYaxis()->GetTitleOffset() + 0.33 * scaleOffset);
    hr->GetYaxis()->SetTitleOffset(hr->GetYaxis()->GetTitleOffset() + 0.42 * scaleOffset);
  }
  if((p % w) != 0)
    utils::ttmd::ScaleAxisFonts(hr->GetYaxis(), 0.0);
  if((npads - p) > w)
    utils::ttmd::ScaleAxisFonts(hr->GetXaxis(), 0.0);
}



void utils::ttmd::SetOZStyle()
{
  //return;
  gStyle->SetOptStat(000000000);
  gStyle->SetTitle(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  //gStyle->SetPadGridX(1);
  //gStyle->SetPadGridY(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetEndErrorSize(5);
  gStyle->SetNdivisions(310);
  TGaxis::SetMaxDigits(4);
  SetNicePalette();
  // 31.10.17
  //printf("GetHatchesSpacing: %f\n", gStyle->GetHatchesSpacing());
  //gStyle->SetHatchesSpacing(0.60);
  //printf("GetHatchesLineWidth: %f\n", gStyle->GetHatchesLineWidth());
  //gStyle->SetHatchesLineWidth(1.5);
  // 26.10.18 extra line styles
  //gStyle->SetLineStyleString(11, "25 10");
  //gStyle->SetLineStyleString(12, "25 10 6 10");
  //gStyle->SetLineStyleString(13, "25 10 6 10 6 10");
  //gStyle->SetLineStyleString(14, "25 10 6 10 6 10 6 10");
  //gStyle->SetLineStyleString(15, "6 6");
  //gStyle->SetLineStyleString(16, "15 10");
  gStyle->SetLineStyleString(11, "17 7");
  gStyle->SetLineStyleString(12, "17 7 4 7");
  gStyle->SetLineStyleString(13, "17 7 4 7 4 7");
  gStyle->SetLineStyleString(14, "17 7 4 7 4 7 4 7");
  gStyle->SetLineStyleString(15, "4 4");
  gStyle->SetLineStyleString(16, "11 7");
}



void utils::ttmd::SetBWPalette(int ncol/* = 999*/)
{
  const Int_t NCont = ncol;

  const UInt_t Number = 2;
  Double_t Red[Number]   = { 0.00, 1.00};
  Double_t Green[Number] = { 0.00, 1.00};
  Double_t Blue[Number]  = { 0.00, 1.00};
  Double_t Stops[Number] = { 0.00, 1.00};    TColor::CreateGradientColorTable(Number, Stops, Red, Green, Blue, NCont);
  gStyle->SetNumberContours(NCont);
}



void utils::ttmd::SetNicePalette(int ncol/* = 999*/)
{
  const Int_t NRGBs = 5;
  const Int_t NCont = ncol;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}



void utils::ttmd::SetNicePaletteEBW()
{
  const Int_t NRGBs = 6;
  const Int_t NCont = 999;

  Double_t stops[NRGBs] = { 0.0,  1.0000001/NCont, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 1.00, 0.00,            0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 1.00, 0.00,            0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 1.00, 0.51,            1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  //    printf("NCont: %d\n", NCont);
  gStyle->SetNumberContours(NCont);
}



/*int GetNiceColor(int i, int n)
{
  SetNicePalette(n);
  return gStyle->GetColorPalette(i);
}*/



int utils::ttmd::GetDistinctColor(int i)
{
#define NCOLSXYZ 13
  static int cols[NCOLSXYZ] = {kBlack, kRed, kBlue, kMagenta, kGreen + 2, kYellow + 1, kAzure + 4, kSpring + 4, kOrange + 2, kRed - 7, kBlue - 9, kRed + 3, kViolet - 7};
  static std::vector<TString> pdf = {"CT14", "MMHT2014", "NNPDF31", "HERAPDF20", "ABMP16", "JR14", "CJ15"};
  if(i >= NCOLSXYZ)
  {
    printf("Error in GetDistinctColor(): color for i = %d not available, returning color for %d\n", i, NCOLSXYZ - 1);
    i = NCOLSXYZ - 1;
  }
  return cols[i];
}



int utils::ttmd::GetDistinctColor(const TString& pdf)
{
#define NPDFSXYZ 7
  static std::vector<TString> pdfs = {"CT14", "MMHT2014", "NNPDF31", "HERAPDF20", "ABMP16", "JR14", "CJ15"};
  for(size_t i = 0; i < NPDFSXYZ; i++)
    if(pdfs[i] == pdf)
      return GetDistinctColor(i);
  printf("Error in GetDistinctColor(): color for pdf = %s not available, returning color 1\n", pdf.Data());
  return GetDistinctColor(1);
}



int utils::ttmd::CheckFile(const TString& fileName, const bool flagCreateDir/* = true*/)
{
  //printf("CheckFile fileName = %s\n", fileName.Data());
  if(fileName == "")
    throw std::runtime_error("Error: checking empty fileName");

  if(FileExists(fileName))
    return 0;
  Ssiz_t last = fileName.Last('/');
  if(last == kNPOS) // current directory
    return 1;
  TString dir = TString(fileName.Data(), last);
  struct stat sb;
  if (stat(dir.Data(), &sb) == 0 && S_ISDIR(sb.st_mode))
    return 1;
  else if(flagCreateDir)
  {
    printf("CheckFile(): creating dir %s\n", dir.Data());
    const int dir_err = system(TString::Format("mkdir -p %s", dir.Data()));
    if (-1 == dir_err)
      throw(TString::Format("Error in CheckFile(): cannot create directory %s\n", dir.Data()));
    return 2;
  }
  throw std::runtime_error("utils::ttmd::CheckFile: should not be here!");
  return -999;
}



bool utils::ttmd::FileExists(const TString& name, const bool flagThrow/* = false*/)
{
  struct stat buffer;
  bool ret = (stat (name.Data(), &buffer) == 0);
  if(!ret && flagThrow)
    throw std::runtime_error(TString::Format("Error: cannot access file %s", name.Data()).Data());
  return ret;
}

double utils::ttmd::DVectorAppendNWidth(std::vector<double>& vec, const double min, const int n, const double width)
{
  for(int i = 0; i < n; i++)
    vec.push_back(min + i * width);
  return vec.back() + width;
}



TGraphAsymmErrors* utils::ttmd::DrawBand(const std::pair<TH1D*, TH1D*>& hPairUD, const TH1D* h)
{
  TGraphAsymmErrors* g = new TGraphAsymmErrors;
  TH1D* hd = (TH1D*) h->Clone();
  hd->SetLineWidth(0);
  hd->SetFillStyle(0);
  TH1D* hu = (TH1D*) h->Clone();
  hu->SetLineWidth(0);
  hu->SetFillStyle(0);
  for(int b = 0; b < h->GetNbinsX(); b++)
  {
    double xMin = h->GetBinLowEdge(b + 1);
    double xMax = h->GetBinLowEdge(b + 2);
    double x = xMin + (xMax - xMin) / 2;
    double y = h->GetBinContent(b + 1);
    //double unc = h->GetBinError(b + 1);
    g->SetPoint(b, x, y);
    g->SetPointError(b, x - xMin, xMax - x, -1 * hPairUD.second->GetBinContent(b + 1), hPairUD.first->GetBinContent(b + 1));
    //g->SetPointError(b, x - xMin, xMax - x, 1000, 1000);
    hd->SetBinContent(b + 1, y + hPairUD.second->GetBinContent(b + 1));
    hu->SetBinContent(b + 1, y + hPairUD.first->GetBinContent(b + 1));
  }
  g->SetLineColor(1);
  g->SetLineStyle(1);
  //g->SetFillColor(6);
  //g->Print("all");
  //g->SetFillStyle(3001);
  g->SetFillStyle(3004);
  //g->SetFillStyle(3006);
  //g->SetFillStyle(3354);
  //g->SetFillStyle(1);
  g->Draw("2");
  hu->Draw("same hist");
  hd->Draw("same hist");
  return g;
}



void utils::ttmd::ScaleRelativeHistoVar(TH1* hvar, const double scale, const TH1* hNom, const bool flagKeepArea/* = false*/)
{
  if(flagKeepArea)
    hvar->Scale(hNom->Integral() / hvar->Integral());
  assert(hvar);
  assert(hNom);
  hvar->Add(hNom, -1.0);
  hvar->Scale(scale);
  hvar->Add(hNom, 1.0);
}



double utils::ttmd::GetMeanDeviationFromRef(const TH1* h, const TH1* hRef)
{
  double dev2 = 0.0;
  for(int b = 0; b < h->GetNbinsX(); b++)
  {
    if(hRef)
      //dev2 += TMath::Power((h->GetBinContent(b + 1) / hRef->GetBinContent(b + 1) - 1.0), 2.0);
      dev2 += TMath::Abs((h->GetBinContent(b + 1) / hRef->GetBinContent(b + 1)) - 1.0);
    else
      //dev2 += TMath::Power(h->GetBinContent(b + 1), 2.0);
      dev2 += TMath::Abs(h->GetBinContent(b + 1));
  }
  //double dev = TMath::Sqrt(dev2 / h->GetNbinsX());
  double dev = dev2 / h->GetNbinsX();
  return dev;
}



double utils::ttmd::GetMeanDeviationFromRef(const TH1* h1, const TH1* h2, const TH1* hRef)
{
  double dev1 = utils::ttmd::GetMeanDeviationFromRef(h1, hRef);
  double dev2 = utils::ttmd::GetMeanDeviationFromRef(h2, hRef);
  //double dev = TMath::Sqrt((dev1 * dev1 + dev2 * dev2) / 2.0);
  double dev = (dev1 + dev2) / 2.0;
  return dev;
}



double utils::ttmd::GetBinContentForValueNoUF(TH2* h, const double x, const double y, int* counterPtr)
{
  const int bin = h->FindBin(x, y);
  int binx = -1;
  int biny = -1;
  int binz = -1;
  h->GetBinXYZ(bin, binx, biny, binz);

  bool flagOut = false;
  if(binx < 1)
  {
    flagOut = true;
    binx = 1;
  }
  else if(binx > h->GetNbinsX())
  {
    flagOut = true;
    binx = h->GetNbinsX();
  }
  if(biny < 1)
  {
    flagOut = true;
    biny = 1;
  }
  else if(biny > h->GetNbinsY())
  {
    flagOut = true;
    biny = h->GetNbinsY();
  }
  if(counterPtr && flagOut)
  {
    (*counterPtr)++;
  }

  const double val = h->GetBinContent(binx, biny);
  return val;
}



void utils::ttmd::ErFgets()
{
  throw std::runtime_error("fgets returned NULL");
}



int utils::ttmd::ReadNumberOfLines(const TString& fileName)
{
  // open file
  FileExists(fileName, 1);
  std::ifstream file(fileName.Data());
  std::string line;
  int nLine = 0;
  while (1)
  {
    // read new line
    getline(file, line);
    if (true == file.eof()) break;
    nLine++;
  }
  return nLine;
}



// this code will work only for separator == single space, and no leading/trailing spaces
int utils::ttmd::ReadNumberOfColumns(const TString& fileName, const int lineStart)
{
  // open file
  //if(fileName == ".././mgtabs/aMCfast-200418-iappl2-tt0j-p0.0002-b15.15-mt1725/conv/h105-mr1.0-mf1.0.txt")
  //  printf("dupa\n");
  FileExists(fileName, 1);
  std::ifstream file(fileName.Data());
  std::string line;

  // skip lineStart lines
  for(int l = 1; l < lineStart; l++)
    getline(file, line);
  if(!file.is_open() || !file.good())
    throw std::runtime_error(TString::Format("Error: can not open file %s", fileName.Data()));

  // read new line
  getline(file, line);
  // TODO why sometimes line is empty?
  // tmp: repeat until it is not empty
  /*while(line == "")
  {
    file = std::ifstream(fileName.Data());
    getline(file, line);
  }*/
  //if(line == "")
  //  printf("dupa\n");

  // count columns in line
  std::vector<std::string> vStr = split(line, ' ');
  int nColumns = vStr.size();
  // somehow the code below does not work sometimes
  /*int nColumns = 0;
  line.erase(line.find_last_not_of(" \n\r\t")+1); // trim trailing whitespaces
  std::stringstream sline(line);
  while (sline.good())
  {
    std::string dummy;
    sline >> dummy;
    nColumns++;
  }*/
  //if(nColumns < 1)
  //  printf("dupa\n");
  file.close();
  return nColumns;
}



void utils::ttmd::ReadVectorKFactor(std::vector<double>& v, const TString& fileName, const int column, const int lineStart, const bool flagMakeNorm)
{
  size_t b = 0;
  // open file
  FileExists(fileName, 1);
  std::ifstream file(fileName.Data());
  std::string line;

  // skip lineStart lines
  for(int l = 1; l < lineStart; l++)
    getline(file, line);

  while (1)
  {
    // read new line
    getline(file, line);
    if (true == file.eof()) break;
    if (line.at(0) == '#' ) continue; //ignore comments
    line.erase(line.find_last_not_of(" \n\r\t")+1); // trim trailing whitespaces
    line.erase(line.find_last_not_of(" ")+1); // trim starting whitespaces
    //line = line.c_str() + line.find_first_of(' ') + 1;
    std::stringstream sline(line);

    // count columns in line
    int nColumns = 0;
    while (sline.good())
    {
      std::string dummy;
      sline >> dummy;
      nColumns++;
    }

    // check that the number of columns is not smaller than the requested column
    if(column > nColumns)
      throw std::runtime_error("Error in ReadHistoKFactor(): no column = " + std::to_string(column) + " in file " + fileName);

    // read kfactor value
    sline.clear();
    sline.seekg(0);
    sline.str(line);
    double val = 0.0;
    for (int col = 0; col < column; col++)
      sline >> val;

    // store kfactor value
    v[b] = val;

    // increment bin
    b++;
  }

  //printf("b = %d  v.size() = %d\n", b, v.size());
  if(b == (v.size() - 1) && flagMakeNorm)
  {
    double sum = 0.0;
    // does not work?
    //std::for_each(v.rbegin(), v.rend() - 1, [&](int n) { sum += n; });
    for(size_t i = 0; i < b; i++)
      sum += v[i];
    //printf("sum = %f\n", sum);
    v[b] = 1.0 - sum;
    b++;
  }

  if(b != v.size())
    throw std::runtime_error("Error in ReadVectorKFactor(): wrong number of read bins");

  if(flagMakeNorm)
  {
    double sum = 0.0;
    for(size_t i = 0; i < v.size(); i++)
      sum += v[i];
    //printf("sum = %f\n", sum);
    for(size_t i = 0; i < v.size(); i++)
      v[i] /= sum;
  }

  file.close();
}



void utils::ttmd::ReadHistoKFactor(TH1D* h, const TString& fileName, const int column, const int lineStart, const bool flagMakeNorm)
{
  std::vector<double> v(h->GetNbinsX());
  ReadVectorKFactor(v, fileName, column, lineStart, flagMakeNorm);
  for(size_t i = 0; i < v.size(); i++)
    h->SetBinContent(i + 1, v[i]);
}



double utils::ttmd::Rhoj(const double mttj, const int njet, const double ptjMin)
{
  const double mt0 = 160.0;
  double val = (2 * mt0 * TMath::Sqrt(1 + njet * ptjMin / mt0)) / mttj;
  return val;
}



TString utils::ttmd::GetNameExtraJetVar(const TString& varName, const TString& def, const double iso, const double pt, const TString& leadStr)
{
  TString str = TString::Format("%s%sIso%.1fPt%.0f", varName.Data(), def.Data(), iso, pt);
  if(leadStr != "")
    str = TString::Format("%sLead%s", str.Data(), leadStr.Data());

  // remove dot in isolation radius, e.g. 0.4 -> 04
  str.ReplaceAll(".", "");

  return str;
}



TString utils::ttmd::Execute(const TString& cmd)
{
  //printf("Execute: %s\n", cmd.Data());
  char buffer[1024];
  TString result = "";
  FILE* pipe = popen(cmd.Data(), "r");
  if (!pipe) throw std::runtime_error("popen() failed!");
  try {
    while (!feof(pipe)) {
      if (fgets(buffer, 1024, pipe) != NULL)
        result += buffer;
    }
  } catch (...) {
    pclose(pipe);
    throw;
  }
  pclose(pipe);
  //printf("result: %s\n", result.Data());
  return result;
}



TH1D* utils::ttmd::CloneTH1D(const TH1D* h, const double content)
{
  TH1D* hClone = (TH1D*) h->Clone();
  for(int b = 0; b < hClone->GetNbinsX(); b++)
    hClone->SetBinContent(b + 1, content);
  return hClone;
}



int utils::ttmd::GetDistinctMarkerStyle(int i)
{
#define NSTYLESSXYZ 7
  static int styles[NSTYLESSXYZ] = {24, 25, 26, 32, 27, 28, 30};
  if(i >= NSTYLESSXYZ)
  {
    printf("Error in GetDistinctMarkerStyle(): color for i = %d not available, returning style for %d\n", i, NSTYLESSXYZ - 1);
    i = NSTYLESSXYZ - 1;
  }
  return styles[i];
}



int utils::ttmd::GetDistinctLineStyle(int i)
{
#define NSTYLESSXYZ 7
  static int styles[NSTYLESSXYZ] = { 1, 11, 12, 13, 14, 15, 16 };
  if(i >= NSTYLESSXYZ)
  {
    printf("Error in GetDistinctLineStyle(): style for i = %d not available, returning style for %d\n", i, NSTYLESSXYZ - 1);
    i = NSTYLESSXYZ - 1;
  }
  return styles[i];
}



std::vector<double> utils::ttmd::MakeSubbins(const std::vector<double>& v, const int nsubbins)
{
  std::vector<double> vNewBins;
  for(size_t b = 0; b < v.size() - 1; b++)
  {
    assert(nsubbins > 1);
    double min = v[b];
    double max = v[b + 1];
    double width = (max - min) / nsubbins;
    for(int bb = 0; bb < nsubbins; bb++)
      vNewBins.push_back(min + width * bb);
  }
  vNewBins.push_back(v.back());
  return vNewBins;
}



const TString utils::ttmd::TexAdoptVarTitle(const TString& str)
{
  TString res = str;
  res.ReplaceAll("t#bar{t}", "{\\rm t}#bar{\\rm t}");
  res.ReplaceAll("t,#bar{t}", "{\\rm t},#bar{\\rm t}");
  res.ReplaceAll("(t)", "({\\rm t})");
  res.ReplaceAll("N_{j}", "N_{\\rm jet}");
  res.ReplaceAll("N_{jet}", "N_{\\rm jet}");
  res.ReplaceAll("#", "\\");
  res.ReplaceAll("p_{T}", "p_{\\rm T}");
  return res;
}



const TString utils::ttmd::GetDigitTimesFormat(const double val, const int nDig)
{
  assert(val > 0.0);
  int n = 0;
  double valCopyPlus = val;
  while(valCopyPlus < 1.0)
  {
    valCopyPlus = valCopyPlus * 10.0;
    n--;
  }
  double valCopyMinus = val;
  while(valCopyMinus > 10.0)
  {
    valCopyMinus = valCopyMinus / 10.0;
    n++;
  }
  TString str = TString::Format("%.*f \\times 10^{%d}", nDig, (n > 0) ? valCopyMinus : valCopyPlus, n);
  return str;
}



const TString utils::ttmd::TexAdoptBinLabel(const TString& str)
{
  TString res = str;
  res.ReplaceAll("$-0.5$--$0.5$", "$0$");
  res.ReplaceAll("$0.5$--$8.5$", "$\\ge 1$");
  res.ReplaceAll("$0.5$--$1.5$", "$1$");
  res.ReplaceAll("$1.5$--$8.5$", "$\\ge 2$");
  res.ReplaceAll("$1.5$--$2.5$", "$2$");
  res.ReplaceAll("$2.5$--$8.5$", "$\\ge 3$");
  return res;
}



const TString utils::ttmd::TexAdoptVarYAxisTitle(const TString& str)
{
  TString title = TexAdoptVarTitle(str);
  title.ReplaceAll("d", "{\\rm d}");
  title.ReplaceAll("ra{\\rm d}", "{\\rm rad}");
  title = "\\frac{1}{\\sigma(\\ttbar)} \\frac{" + title;
  title.ReplaceAll("/", "}{");
  title.ReplaceAll("[", "} [");
  title.ReplaceAll("pb", "");
  title.ReplaceAll("GeV", "{\\rm GeV}");
  title.ReplaceAll("rad", "{\\rm rad}");
  // 20.11.18 avoid empty brackets
  title.ReplaceAll("[]", "");
  return title;
}



const TString utils::ttmd::TexXsecTitle(const TString& str, TString& strExtraBins)
{
  TString strXSec = str;
  strExtraBins = "";
  strXSec.ReplaceAll("-xsec8", "");
  if(strXSec == "njmttyttdefIso04Pt30-b2-mtt3")
  {
    //strExtraBins = " (with two \\nj bins)";
    strXSec = "njmttytttwo";
  }
  else if(strXSec == "njmttyttdefIso04Pt30-b3-mtt3")
  {
    //strExtraBins = " (with three \\nj bins)";
    strXSec = "njmttyttthree";
  }
  strXSec = "\\" + strXSec;
  return strXSec;
}



const TString utils::ttmd::TexXsecLabel(const TString& str, const int nj)
{
  TString strXSecLabel = str;
  strXSecLabel.ReplaceAll("\\", "");
  if(strXSecLabel.BeginsWith("nj"))
  {
    if(nj == 2)
      strXSecLabel.ReplaceAll("nj", "nj2");
    else if(nj == 3)
      strXSecLabel.ReplaceAll("nj", "nj3");
    else if(nj == 4)
      strXSecLabel.ReplaceAll("nj", "nj4");
  }
  return strXSecLabel;
}



// Jacobian_ij = (delta_ij - delta_i(j+-12)) * a[i+-12] / s^2, s = a[i] + a[i+-12]
void utils::ttmd::DoRjHisto(TH1* h, TH2* covmat)
{
  assert(h->GetNbinsX() == 24);
  const int offset = 12;
  int n = h->GetNbinsX();
  if(covmat && (n != covmat->GetNbinsX() || n != covmat->GetNbinsY()))
  {
    printf("Error in TransformHisto() incosistent input h->GetNbinsX() = %d covmat->GetNbinsX() = %d covmat->GetNbinsY() = %d\n", h->GetNbinsX(), covmat->GetNbinsX(), covmat->GetNbinsY());
    throw;
  }
  std::vector<double> vSum(n);
  std::vector<double> vSum2(n);
  std::vector<double> vPairIndex(n);
  for(int i = 0; i < n; i++)
  {
    if(i < offset)
    {
      vSum[i] = h->GetBinContent(i + 1) + h->GetBinContent(i + 1 + offset);
      vPairIndex[i] = i + offset;
    }
    else
    {
      vPairIndex[i] = i;
      vSum[i] = h->GetBinContent(i + 1) + h->GetBinContent(i + 1 - offset);
    }
    vSum2[i] = vSum[i] * vSum[i];
  }

  TMatrixD matrixG(n, n);
  TMatrixD matrixCov(n, n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j < n; j++)
    {
      if(covmat)
        matrixCov(i, j) = covmat->GetBinContent(i + 1, j + 1);
      else
      {
        if(i == j)
          matrixCov(i, j) = h->GetBinError(i + 1) * h->GetBinError(i + 1);
        else
          matrixCov(i, j) = 0.0;
      }
      if(i == j)
        matrixG(i, j) = h->GetBinContent(vPairIndex[i] + 1) / vSum2[i];
      else if(abs(i - j) == offset)
        matrixG(i, j) = -1 * h->GetBinContent(vPairIndex[i] + 1) / vSum2[i];
      else
        matrixG(i, j) = 0;
    }
  TMatrixD matrixGT = matrixG;
  matrixGT.T();
  TMatrixD res = matrixG * matrixCov * matrixGT;

  for(int i = 0; i < n; i++)
  {
    h->SetBinContent(i + 1, h->GetBinContent(i + 1) / vSum[i]);
    h->SetBinError(i + 1, TMath::Sqrt(res(i, i)));
    if(covmat)
      for(int j = 0; j < n; j++)
        covmat->SetBinContent(i + 1, j + 1, res(i, j));
  }
}

void styleUtils::setResultLegendStyle(TLegend* leg, const bool /*result*/)
{
    double x1 = 0.7, y1 = 0.5;
    double height = 0.2+0.155, width = 0.2;

    leg->SetX1NDC(x1);
    leg->SetY1NDC(y1);
    leg->SetX2NDC(x1 + width);
    leg->SetY2NDC(y1 + height);

    leg->SetTextFont(42);
    leg->SetTextAlign(12);
    leg->SetTextSize(0.04);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
}



void styleUtils::setHHStyle(TStyle& HHStyle)
{
    //const int fontstyle=42;
    HHStyle.SetPalette(1);

    // ==============
    //  Canvas
    // ==============

//     HHStyle.SetCanvasBorderMode(0);
//     HHStyle.SetCanvasColor(kWhite);
//     HHStyle.SetCanvasDefH(600); //Height of canvas
//     HHStyle.SetCanvasDefW(600); //Width of canvas
//     HHStyle.SetCanvasDefX(0);   //Position on screen
//     HHStyle.SetCanvasDefY(0);

    // ==============
    //  Pad
    // ==============

//     HHStyle.SetPadBorderMode(0);
//     // HHStyle.SetPadBorderSize(Width_t size = 1);
//     HHStyle.SetPadColor(kWhite);
//     HHStyle.SetPadGridX(false);
//     HHStyle.SetPadGridY(false);
//     HHStyle.SetGridColor(0);
//     HHStyle.SetGridStyle(3);
//     HHStyle.SetGridWidth(1);

    // ==============
    //  Frame
    // ==============

//     HHStyle.SetFrameBorderMode(0);
//     HHStyle.SetFrameBorderSize(1);
//     HHStyle.SetFrameFillColor(0);
//     HHStyle.SetFrameFillStyle(0);
//     HHStyle.SetFrameLineColor(1);
//     HHStyle.SetFrameLineStyle(1);
//     HHStyle.SetFrameLineWidth(1);

    // ==============
    //  Histo
    // ==============

    HHStyle.SetErrorX(0.0);
    //HHStyle.SetEndErrorSize(0);

    // HHStyle.SetHistFillColor(1);
    // HHStyle.SetHistFillStyle(0);
    HHStyle.SetHistLineColor(1);
    HHStyle.SetHistLineStyle(0);
    HHStyle.SetHistLineWidth(1);
    // HHStyle.SetLegoInnerR(Float_t rad = 0.5);
    // HHStyle.SetNumberContours(Int_t number = 20);

    // HHStyle.SetErrorMarker(20);
    HHStyle.SetMarkerStyle(20);

    // ==============
    //  Fit/function
    // ==============

    HHStyle.SetOptFit(1);
    HHStyle.SetFitFormat("5.4g");
    HHStyle.SetFuncColor(2);
    HHStyle.SetFuncStyle(1);
    HHStyle.SetFuncWidth(1);

    // ==============
    //  Date
    // ==============

//     HHStyle.SetOptDate(0);
//     // HHStyle.SetDateX(Float_t x = 0.01);
//     // HHStyle.SetDateY(Float_t y = 0.01);
//
//     // =====================
//     //  Statistics Box
//     // =====================
//
//     HHStyle.SetOptFile(0);
//     HHStyle.SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
//     HHStyle.SetStatColor(kWhite);
//     HHStyle.SetStatFont(fontstyle);
//     HHStyle.SetStatFontSize(0.025);
//     HHStyle.SetStatTextColor(1);
//     HHStyle.SetStatFormat("6.4g");
//     HHStyle.SetStatBorderSize(1);
//     HHStyle.SetStatH(0.1);
//     HHStyle.SetStatW(0.15);
//     // HHStyle.SetStatStyle(Style_t style = 1001);
//     // HHStyle.SetStatX(Float_t x = 0);
//     // HHStyle.SetStatY(Float_t y = 0);
//
//     // ==============
//     //  Margins
//     // ==============
//
//     HHStyle.SetPadTopMargin(0.1);
//     HHStyle.SetPadBottomMargin(0.15);
//     HHStyle.SetPadLeftMargin(0.20);
//     HHStyle.SetPadRightMargin(0.05);
//
//     // ==============
//     //  Global Title
//     // ==============
//
//     HHStyle.SetOptTitle(0);
//     HHStyle.SetTitleFont(fontstyle);
//     HHStyle.SetTitleColor(1);
//     HHStyle.SetTitleTextColor(1);
//     HHStyle.SetTitleFillColor(10);
//     HHStyle.SetTitleFontSize(0.05);
//     // HHStyle.SetTitleH(0); // Set the height of the title box
//     // HHStyle.SetTitleW(0); // Set the width of the title box
//     // HHStyle.SetTitleX(0); // Set the position of the title box
//     // HHStyle.SetTitleY(0.985); // Set the position of the title box
//     // HHStyle.SetTitleStyle(Style_t style = 1001);
//     // HHStyle.SetTitleBorderSize(2);
//
//     // ==============
//     //  Axis titles
//     // ==============
//
//     HHStyle.SetTitleColor(1, "XYZ");
//     HHStyle.SetTitleFont(fontstyle, "XYZ");
//     HHStyle.SetTitleSize(0.04, "XYZ");
//     // HHStyle.SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
//     // HHStyle.SetTitleYSize(Float_t size = 0.02);
//     HHStyle.SetTitleXOffset(1.25);
//     HHStyle.SetTitleYOffset(1.6);
//     // HHStyle.SetTitleOffset(1.1, "Y"); // Another way to set the Offset
//
//     // ==============
//     //  Axis Label
//     // ==============
//
//     //HHStyle.SetLabelColor(1, "XYZ");
//     HHStyle.SetLabelFont(fontstyle, "XYZ");
//     HHStyle.SetLabelOffset(0.007, "XYZ");
//     HHStyle.SetLabelSize(0.04, "XYZ");
//
//     // ==============
//     //  Axis
//     // ==============
//
//     HHStyle.SetAxisColor(1, "XYZ");
//     HHStyle.SetStripDecimals(kTRUE);
//     HHStyle.SetTickLength(0.03, "XYZ");
//     HHStyle.SetNdivisions(510, "XYZ");
//     HHStyle.SetPadTickX(1);  // To get tick marks on the opposite side of the frame
//     HHStyle.SetPadTickY(1);
//
//     // Change for log plots:
//     HHStyle.SetOptLogx(0);
//     HHStyle.SetOptLogy(0);
//     HHStyle.SetOptLogz(0);
//
//     // ==============
//     //  Text
//     // ==============
//
//     HHStyle.SetTextAlign(11);
//     HHStyle.SetTextAngle(0);
//     HHStyle.SetTextColor(1);
//     HHStyle.SetTextFont(fontstyle);
//     HHStyle.SetTextSize(0.05);
//
//     // =====================
//     //  Postscript options:
//     // =====================
//
//     HHStyle.SetPaperSize(20.,20.);
//     // HHStyle.SetLineScalePS(Float_t scale = 3);
//     // HHStyle.SetLineStyleString(Int_t i, const char* text);
//     // HHStyle.SetHeaderPS(const char* header);
//     // HHStyle.SetTitlePS(const char* pstitle);
//
//     // HHStyle.SetBarOffset(Float_t baroff = 0.5);
//     // HHStyle.SetBarWidth(Float_t barwidth = 0.5);
//     // HHStyle.SetPaintTextFormat(const char* format = "g");
//     // HHStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
//     // HHStyle.SetTimeOffset(Double_t toffset);
//     // HHStyle.SetHistMinimumZero(kTRUE);
}


uint utils::helper::GetNumberOfNonZeroBins(TH1& inputHisto){
    int nBins = 0;
    for (uint i=0; i< (uint) inputHisto.GetNbinsX(); i++){
        if(inputHisto.GetBinContent(i)>0.){
            nBins = nBins + 1;
        }
    }
    return nBins;
}

std::vector<double> utils::GetVectorFromTStringData(TString vector_tstring){
    std::vector<double> output_vector = {};

    if (vector_tstring != ""){
        TObjArray* input_vector  = vector_tstring.Tokenize(",");
        int data_size = int(input_vector->GetEntries());
        for (int i=0; i<data_size; i++){
            output_vector.push_back(std::stod((((TObjString*)input_vector->At(i))->GetString()).Data()));
        }
    }
    else{
        std::cerr << "ERROR in utils::GetVectorFromTStringData -->> input text is empty!!!" << std::endl;
        exit(1);
    }
    std::cout << "Getting warnings" << std::endl;

    return output_vector;
}
