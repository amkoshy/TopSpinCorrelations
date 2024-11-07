///////////////////////////////////////////////
////
////    Execute this macro via:   root -l  macros/compareClosure.C++
////
///////////////////////////////////////////////

#include <TString.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TH1.h>

#include <vector>
#include "iostream"



void plot(TString var, TString channel, std::vector<TString> slopes)
{

    if(var != "pt" && var != "y" && 
       channel != "emu" && channel != "ee" && channel != "mumu" &&
       slopes.size()==0) {
           std::cout<<"Not valid parameters!"<<std::endl;
           exit(89);
    }
    
    double yMax = 1.0, yMin = 0.;
    bool logY = 0;
    TString xaxis = "", title = "", legendTitle = "";
    if(var =="pt"){
        yMax = 1.e-2;
        yMin = 1.e-5;
        logY = 1;
        xaxis = "p_{T}^{t}";
        title = "Top quark p_{T}";
        legendTitle = "max(0.1, (1+|p_{T}^{t}-1|*s) *((1+|p_{T}^{#bar{t}}-1|*s))";
    } else if (var == "y"){
        yMax = 0.6;
        yMin = 1.e-3;
        logY = 0;
        xaxis = "y^{t}";
        title = "Top quark y";
        legendTitle = "max(0.1, (1+|y^{t}-1|*s) *((1+|y^{#bar{t}}-1|*s))";
    }
    TString yaxis = TString("#frac{1}{#sigma} #frac{d #sigma}{d"+xaxis+"}");

    std::vector<TGraphErrors*> gr;
    TLegend *leg = new TLegend(0.65, 0.6, 0.9, 0.9);
    leg->SetHeader(legendTitle);
    leg->SetTextSize(0.03);

    for(unsigned int i=0; i<slopes.size(); i++){
        TString file = TString("Closure/closure_result_"+var+"_"+channel+"_"+slopes.at(i)+".txt");
        TGraphErrors *tmp = new TGraphErrors(file , "%lg %lg %lg %lg");
        if(tmp->GetN()==0) {
            std::cout<<"File: "<<file<<"  not plottable"<<std::endl;
            continue;
        }

        tmp->SetMarkerStyle(20+i);
        tmp->SetMarkerColor(1+i);
        tmp->GetXaxis()->SetTitle(xaxis);
        tmp->GetYaxis()->SetTitle(yaxis);
        tmp->SetTitle(title);

        leg->AddEntry(tmp, "s="+slopes.at(i), "p");
        gr.push_back(tmp);
    }
    TCanvas *c = new TCanvas();
    for(unsigned int i=0; i<gr.size();i++)
    {
        if(i==0){
            gr.at(i)->GetHistogram()->SetMaximum(yMax);
            gr.at(i)->GetHistogram()->SetMinimum(yMin);
            gr.at(i)->Draw("ap");
        } else {
            gr.at(i)->Draw("p,same");
        }
    }
    leg->Draw("same");
    c->SetLogy(logY);
    c->Print("Closure/closureComparison_"+var+"_"+channel+".eps");
    c->Print("Closure/closureComparison_"+var+"_"+channel+".C");

    delete c;
    delete leg;
}


void compareClosure()
{
    const TString yslope_Array[] = {"0.320", "0.240", "0.160", "0.080", "-0.080", "-0.160", "-0.240", "-0.320"};
    const TString ptslope_Array[] = {"0.008", "0.006", "0.004", "0.002", "-0.002", "-0.004", "-0.006", "-0.008"};

    std::vector<TString> ptslope (ptslope_Array, ptslope_Array + sizeof(ptslope_Array)/sizeof(ptslope_Array[0]));
    std::vector<TString> yslope (yslope_Array, yslope_Array + sizeof(yslope_Array)/sizeof(yslope_Array[0]));

    TString channel[] = {"emu", "ee", "mumu"};

    for (unsigned int i=0; i<3; i++){
        plot(TString("pt"), channel[i], ptslope);
        plot(TString("y"), channel[i], yslope );
    }

}