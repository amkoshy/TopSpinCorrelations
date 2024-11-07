import ROOT
from array import array

outDir = "TriggerSFs/";

filesAndHistos = {
    "2016": {
        "ee": ["trigger/94X_2016/tth_dileptonic_2DscaleFactors_MC_PU-lep_SFs_withSysts_2016BCDEFGH_24-06-20_updated_axis_limits.root", "h_DoubleEl_OR__X__allMET_el0_pt_vs_el1_pt_withSysts", "2016 dielectron trigger scale factors", "p_{T}^{e_{leading}} [GeV]", "p_{T}^{e_{trailing}} [GeV]"],
        "emu":  ["trigger/94X_2016/tth_dileptonic_2DscaleFactors_MC_PU-lep_SFs_withSysts_2016BCDEFGH_24-06-20_updated_axis_limits.root", "h_EMu_OR__X__allMET_mu0_pt_vs_el0_pt_withSysts", "2016 electron-muon trigger scale factors", "p_{T}^{#mu} [GeV]", "p_{T}^{e} [GeV]"],
        "mumu": ["trigger/94X_2016/tth_dileptonic_2DscaleFactors_MC_PU-lep_SFs_withSysts_2016BCDEFGH_24-06-20_updated_axis_limits.root", "h_DoubleMu_OR__X__allMET_mu0_pt_vs_mu1_pt_withSysts", "2016 dimuon trigger scale factors", "p_{T}^{#mu_{leading}} [GeV]", "p_{T}^{#mu_{trailing}} [GeV]"],
    },
    "2017": {
        "ee": ["trigger/94X_2017/tth_dileptonic_2DscaleFactors_MC_PU-lep_SFs_withSysts_2017BCDEF_24-06-20_updated_axis_limits.root", "h_DoubleEl_OR__X__allMET_el0_pt_vs_el1_pt_withSysts", "2017 dielectron trigger scale factors", "p_{T}^{e_{leading}} [GeV]", "p_{T}^{e_{trailing}} [GeV]"],
        "emu":  ["trigger/94X_2017/tth_dileptonic_2DscaleFactors_MC_PU-lep_SFs_withSysts_2017BCDEF_24-06-20_updated_axis_limits.root", "h_EMu_OR__X__allMET_mu0_pt_vs_el0_pt_withSysts", "2017 electron-muon trigger scale factors", "p_{T}^{#mu} [GeV]", "p_{T}^{e} [GeV]"],
        "mumu": ["trigger/94X_2017/tth_dileptonic_2DscaleFactors_MC_PU-lep_SFs_withSysts_2017BCDEF_24-06-20_updated_axis_limits.root", "h_DoubleMu_OR__X__allMET_mu0_pt_vs_mu1_pt_withSysts", "2017 dimuon trigger scale factors", "p_{T}^{#mu_{leading}} [GeV]", "p_{T}^{#mu_{trailing}} [GeV]"],
    },
    "2018": {
        "ee": ["trigger/102X_2018/tth_dileptonic_2DscaleFactors_MC_PU-lep_SFs_withSysts_2018ABCD_24-06-20_updated_axis_limits.root", "h_DoubleEl_OR__X__allMET_el0_pt_vs_el1_pt_withSysts", "2018 dielectron trigger scale factors", "p_{T}^{e_{leading}} [GeV]", "p_{T}^{e_{trailing}} [GeV]"],
        "emu":  ["trigger/102X_2018/tth_dileptonic_2DscaleFactors_MC_PU-lep_SFs_withSysts_2018ABCD_24-06-20_updated_axis_limits.root", "h_EMu_OR__X__allMET_mu0_pt_vs_el0_pt_withSysts", "2018 electron-muon trigger scale factors", "p_{T}^{#mu} [GeV]", "p_{T}^{e} [GeV]"],
        "mumu": ["trigger/102X_2018/tth_dileptonic_2DscaleFactors_MC_PU-lep_SFs_withSysts_2018ABCD_24-06-20_updated_axis_limits.root", "h_DoubleMu_OR__X__allMET_mu0_pt_vs_mu1_pt_withSysts", "2018 dimuon trigger scale factors", "p_{T}^{#mu_{leading}} [GeV]", "p_{T}^{#mu_{trailing}} [GeV]"],
    },
}





def set_plot_style():
    # NRGBs = 5
    # NCont = 255
    #
    # stops = [ 0.00, 0.34, 0.61, 0.84, 1.00 ]
    # red   = [ 0.00, 0.00, 0.87, 1.00, 0.51 ]
    # green = [ 0.00, 0.81, 1.00, 0.20, 0.00 ]
    # blue  = [ 0.51, 1.00, 0.12, 0.00, 0.00 ]

    NRGBs = 9
    NCont = 255
    stops = [0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000]

    # #Bird
    # red   = [ 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764]
    # green = [ 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832]
    # blue  = [ 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539]
    #
    #Light Temperature
    red   = [  31./255.,  71./255., 123./255., 160./255., 210./255., 222./255., 214./255., 199./255., 183./255.]
    green = [  40./255., 117./255., 171./255., 211./255., 231./255., 220./255., 190./255., 132./255.,  65./255.]
    blue  = [ 234./255., 214./255., 228./255., 222./255., 210./255., 160./255., 105./255.,  60./255.,  34./255.]
    #
    # #Rainbow
    # red   = [  0./255.,   5./255.,  15./255.,  35./255., 102./255., 196./255., 208./255., 199./255., 110./255.]
    # green = [  0./255.,  48./255., 124./255., 192./255., 206./255., 226./255.,  97./255.,  16./255.,   0./255.]
    # blue  = [ 99./255., 142./255., 198./255., 201./255.,  90./255.,  22./255.,  13./255.,   8./255.,   2./255.]
    #
    # #Pastel
    # red   = [ 180./255., 190./255., 209./255., 223./255., 204./255., 228./255., 205./255., 152./255.,  91./255.]
    # green = [  93./255., 125./255., 147./255., 172./255., 181./255., 224./255., 233./255., 198./255., 158./255.]
    # blue  = [ 236./255., 218./255., 160./255., 133./255., 114./255., 132./255., 162./255., 220./255., 218./255.]
    #
    # #Cool
    # red   = [  33./255.,  31./255.,  42./255.,  68./255.,  86./255., 111./255., 141./255., 172./255., 227./255.]
    # green = [ 255./255., 175./255., 145./255., 106./255.,  88./255.,  55./255.,  15./255.,   0./255.,   0./255.]
    # blue  = [ 255./255., 205./255., 202./255., 203./255., 208./255., 205./255., 203./255., 206./255., 231./255.]

    s = array('d', stops)
    r = array('d', red)
    g = array('d', green)
    b = array('d', blue)
    ROOT.TColor.CreateGradientColorTable(NRGBs, s, r, g, b, NCont)
    ROOT.gStyle.SetNumberContours(NCont)

def setStyle(myStyle):
    myStyle.SetCanvasBorderMode(0)
    myStyle.SetCanvasColor(ROOT.kWhite)
    myStyle.SetCanvasDefH(600) #Height of canvas
    myStyle.SetCanvasDefW(600) #Width of canvas
    myStyle.SetCanvasDefX(0);  #Position on screen
    myStyle.SetCanvasDefY(0)
    myStyle.SetPadBorderMode(0)
    myStyle.SetPadColor(ROOT.kWhite)
    myStyle.SetPadGridX(False)
    myStyle.SetPadGridY(False)
    myStyle.SetGridColor(0)
    myStyle.SetGridStyle(3)
    myStyle.SetGridWidth(1)
    myStyle.SetFrameBorderMode(0)
    myStyle.SetFrameBorderSize(1)
    myStyle.SetFrameFillColor(0)
    myStyle.SetFrameFillStyle(0)
    myStyle.SetFrameLineColor(1)
    myStyle.SetFrameLineStyle(1)
    myStyle.SetFrameLineWidth(1)
    myStyle.SetPadTopMargin(0.1)
    myStyle.SetPadBottomMargin(0.15)
    myStyle.SetPadLeftMargin(0.18)
    myStyle.SetPadRightMargin(0.15)
    myStyle.SetPaperSize(20.,20.)
    myStyle.SetOptStat(0)
    # #myStyle.SetTitleX(0.25) # to set the x start position in NDC
    myStyle.SetTitleY(0.96) # to set the y position in NDC [0,1]
    # #myStyle.SetTitleW(0.5) #to set the title box width
    # #myStyle.SetTitleH(0.05)


def getEff(channel, year):
    fileO,histO,title,xtitle,ytitle = filesAndHistos[year][channel]
    # file = ROOT.TFile("../common/data/triggerSummary_"+channel+"_FullRun2016_ReReco_TightMu_TightEle_FineBins_24July17.root")
    file = ROOT.TFile("../common/data/"+fileO,'OPEN')

    setStyle(ROOT.gStyle)
    set_plot_style()

    beff = file.Get(histO)

    text = ROOT.TPaveText(0.1,0.80,0.3, 0.9, "NDC")
    text.SetFillColor(0)
    text.SetFillStyle(0)
    text.SetBorderSize(0)
    text.SetTextSize(0.05)
    text.SetTextFont(42)


    channelLatex = ""
    if(channel == "mumu"):
        channelLatex = "#mu#mu"
    elif (channel == "emu"):
        channelLatex = "e#mu"
    elif (channel == "ee"):
        channelLatex = "ee"
    text.AddText(channelLatex);

    ROOT.gSystem.mkdir(outDir)
    ROOT.gSystem.mkdir(outDir+"/"+year+"/")
    c = ROOT.TCanvas("", "", 0,0, 1000, 1000)

    ROOT.gStyle.SetPaintTextFormat("4.2f")

    # xtitle = "|#eta_{lepton-1}|"
    # ytitle = "|#eta_{lepton-2}|"
    # title = ""
    # if(channel == "mumu"):
    #     title = "Dimuon trigger scale factors"
    # elif (channel == "emu"):
    #     title = "Electron-muon trigger scale factors"
    # elif (channel == "ee"):
    #     title = "Dielectron trigger scale factors"

    beff.GetXaxis().SetTitle(xtitle)
    beff.GetYaxis().SetTitle(ytitle)
    beff.GetYaxis().SetTitleOffset(1.15*beff.GetYaxis().GetTitleOffset())
    beff.GetXaxis().SetLabelFont(42)
    beff.GetXaxis().SetLabelOffset(0.008)
    beff.GetXaxis().SetLabelSize(0.035)
    beff.GetXaxis().SetTitleSize(0.045)
    beff.GetXaxis().SetNdivisions(-506, ROOT.kTRUE)
    beff.GetYaxis().SetLabelFont(42)
    beff.GetYaxis().SetLabelOffset(0.008)
    beff.GetYaxis().SetLabelSize(0.035)
    beff.GetYaxis().SetTitleSize(0.045)
    beff.GetYaxis().SetNdivisions(-506, ROOT.kTRUE)
    beff.SetTitle(title)
    beff.GetZaxis().SetRangeUser(0.95,1.05)
    beff.SetMarkerSize(1.4)
    beff.Draw("colz,text,E")
    c.SaveAs(outDir+"/"+year+"/"+channel+"_trigsf.pdf")

    file.Close()





def getTriggerSFs():
    for year in filesAndHistos:
        for channel in filesAndHistos[year]:
            getEff(channel, year)


def main():
    getTriggerSFs()


if __name__ == "__main__":
    main()
