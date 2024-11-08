
#include <iostream>
#include <fstream>

#include "TString.h"

#include "Output.h"



Output::Output(const TString& outputType):
nLines(0),
type(outputType)
{
    
    if(outputType.EqualTo("sample"))
    {
        headerLine_.push_back("entries");
        headerLine_.push_back("weight");
        headerLine_.push_back("file");
        headerLine_.push_back("name");
    }
    else if(outputType.EqualTo("xsec"))
    {
        
    }
    else if(outputType.EqualTo("syst"))
    {
        
    }
    else if(outputType.EqualTo("tau"))
    {
        
    }
    else
    {
        std::cout << "[Output]:  ERROR: \"outputType\" does not match to any existing types.";
    }
    
}



void Output::addLine(TString value)
{
    for(const auto& name : headerLine_) {
        if(name!=headerLine_.at(0))value= " ";
        this->add(name,value);
    }
    
}



void Output::add(const TString& name, TString value)
{
    if(value=="")value="...";
    mainBodyS_[name].push_back(value);
    if(name==headerLine_.at(0))nLines++;//for samples
}



void Output::add(const TString& name, std::vector<double> vect)
{
    auto tempPair = std::make_pair(name, vect);
    mainBodyD_.push_back(tempPair);
    nLines++;
}



void Output::save(const TString& saveFile)
{
    
    std::ofstream file_out(saveFile);
    
    if(type.EqualTo("sample")){
        for(const auto& s : headerLine_) file_out << s << "\t\t";
        file_out << std::endl;
        for(int i=0;i<nLines;++i)
        {
            for(const auto& s : headerLine_)
            {
                file_out << mainBodyS_[s].at(i) << "\t\t";
            }
            file_out << std::endl;
        }
    }
    
    if(type.EqualTo("xsec") || type.EqualTo("syst") || type.EqualTo("tau")){
        for (mapDouble::iterator i = mainBodyD_.begin(); i != mainBodyD_.end(); ++i) {
            file_out << i->first << " ";
            for(const auto& s :  i->second )
            {
                file_out << s << " ";
            }
            file_out << std::endl;
        }
    }
    
    file_out.close();
    
    std::cout << "Info in <Output::save>  " << saveFile << " " << "has been created !!!" << std::endl;
    
}

// void Output::addPar(const TString& parName)
// {
//     
// }



