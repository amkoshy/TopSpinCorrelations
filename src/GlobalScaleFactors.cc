#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstdlib>

#include <TString.h>

#include "GlobalScaleFactors.h"
#include "LuminosityScaleFactors.h"
#include "DyScaleFactors.h"
#include "Samples.h"
#include "../../common/include/sampleHelpers.h"
#include "../../common/include/RootFileReader.h"











GlobalScaleFactors::GlobalScaleFactors(const std::vector<Channel::Channel>& v_channel,
                                       const std::vector<Systematic::Systematic>& v_systematic,
                                       const double& luminosityInInversePb,
                                       const double& luminosityUncertaintyRelative,
                                       const bool dyCorrection):
luminosityInInversePb_(luminosityInInversePb),
luminosityUncertaintyRelative_(luminosityUncertaintyRelative),
luminosityScaleFactors_(0),
dyScaleFactors_(0),
rootFileReader_(RootFileReader::getInstance())
{
    std::cout<<"--- Beginning to set up global scale factors\n";
    
    // Need to check which samples are needed for calculating the requested scale factors
    // Drell-Yan scale factors require Samples for channels "ee" "emu" "mumu"
    std::vector<Channel::Channel> v_channelForCorrections = v_channel;
    if(dyCorrection){
        if(std::find(v_channel.begin(), v_channel.end(), Channel::ee) == v_channel.end()) v_channelForCorrections.push_back(Channel::ee);
        if(std::find(v_channel.begin(), v_channel.end(), Channel::emu) == v_channel.end()) v_channelForCorrections.push_back(Channel::emu);
        if(std::find(v_channel.begin(), v_channel.end(), Channel::mumu) == v_channel.end()) v_channelForCorrections.push_back(Channel::mumu);
    }
    Samples scalingSamples("FileLists", v_channelForCorrections, v_systematic, 0);
    
    // Produce luminosity scale factors
    luminosityScaleFactors_ = new LuminosityScaleFactors(scalingSamples, luminosityInInversePb_, luminosityUncertaintyRelative_, rootFileReader_);
    scalingSamples.setGlobalWeights(this);
    
    // Produce Drell-Yan scale factors
    if(dyCorrection){
        dyScaleFactors_ = new DyScaleFactors(scalingSamples, rootFileReader_);
        scalingSamples.setGlobalWeights(this);
    }
    
    std::cout<<"=== Finishing to set up global scale factors\n\n";
}



std::pair<SystematicChannelFactors, bool> GlobalScaleFactors::scaleFactors(const Samples& samples, const TString& step,
                                                                           const bool dyCorrection)const
{
    bool anyCorrectionApplied(false);
    SystematicChannelFactors systematicChannelFactors;
    const SystematicChannelFactors& luminosityScaleFactors = luminosityScaleFactors_->scaleFactorMap();
    
    for(const auto& systematicChannelSamples : samples.getSystematicChannelSamples()){
        const Systematic::Systematic& systematic(systematicChannelSamples.first);
        for(const auto& channelSample : systematicChannelSamples.second){
            const Channel::Channel& channel(channelSample.first);
            const std::vector<Sample>& samples(channelSample.second);
            for(size_t iSample = 0; iSample < samples.size(); ++iSample){
                const Sample& sample(samples.at(iSample));
                double weight(-999.);
                if(sample.sampleType() != Sample::data){
                    weight = luminosityScaleFactors.at(systematic).at(channel).at(iSample);
                    
                    const int dyStatus = (dyCorrection && dyScaleFactors_) ? dyScaleFactors_->applyScaleFactor(weight, step, sample, systematic) : 0;
                    if(dyStatus == 1) anyCorrectionApplied = true;
                }
                else weight = 1; //data weight!!!
                systematicChannelFactors[systematic][channel].push_back(weight);
            }
        }
    }
    
    return std::make_pair(systematicChannelFactors, anyCorrectionApplied);
}



std::pair<SystematicChannelFactors, bool> GlobalScaleFactors::scaleFactors(const Samples& samples, const TString& step)const
{
    return this->scaleFactors(samples, step, dyScaleFactors_);
}



double GlobalScaleFactors::luminosityInInversePb()const
{
    return luminosityInInversePb_;
}







