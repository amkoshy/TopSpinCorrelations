import ROOT
import argparse
from array import array

outDir = "bTagEfficiencies/";
valid_years = ["2016preVFP", "2016postVFP", "2017", "2018"]

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
    myStyle.SetPalette(1);
    myStyle.SetCanvasBorderMode(0)
    myStyle.SetCanvasColor(ROOT.kWhite)
    myStyle.SetCanvasDefH(600) #Height of canvas
    myStyle.SetCanvasDefW(600) #Width of canvas
    myStyle.SetCanvasDefX(0)   #Position on screen
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
    myStyle.SetPadTopMargin(0.1);
    myStyle.SetPadBottomMargin(0.1)
    myStyle.SetPadLeftMargin(0.11)
    myStyle.SetPadRightMargin(0.15)
    myStyle.SetPaperSize(20.,20.)
    myStyle.SetOptStat(0)
    # myStyle.SetTitleX(0.25) // to set the x start position in NDC
    myStyle.SetTitleY(0.96) # to set the y position in NDC [0,1]
    # myStyle.SetTitleW(0.5) //to set the title box width
    # myStyle.SetTitleH(0.05)


def getEff(year, channel):
    if year=="2016preVFP":
        file = ROOT.TFile("BTagEff_"+year+"/Nominal/"+channel+"/"+"ttbarsignalplustau_fromDilepton.root","OPEN")
    if year=="2016postVFP":
        file = ROOT.TFile("BTagEff_"+year+"/Nominal/"+channel+"/"+"ttbarsignalplustau_fromDilepton.root","OPEN")
    if year=="2017":
        file = ROOT.TFile("BTagEff_"+year+"/Nominal/"+channel+"/"+"ttbarsignalplustau_fromDilepton.root","OPEN")
    if year=="2018":
        file = ROOT.TFile("BTagEff_"+year+"/Nominal/"+channel+"/"+"ttbarsignalplustau_fromDilepton.root","OPEN")

        print(file)

    setStyle(ROOT.gStyle)
    #set_plot_style()

    beff = file.Get("beff2D")
    ceff = file.Get("ceff2D")
    leff = file.Get("leff2D")

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
    text.AddText(channelLatex)

    ROOT.gSystem.mkdir(outDir)
    ROOT.gSystem.mkdir(outDir+"/"+year+"/")
    c = ROOT.TCanvas("", "", 0,0, 1000, 1000)

    ROOT.gStyle.SetPaintTextFormat("4.2f")

    xtitle = "p_{T}^{b-jet} [GeV]"
    ytitle = "|#eta_{b-jet}|"
    title = "Tagging efficiency of b-jets"
    beff.GetXaxis().SetTitle(xtitle)
    beff.GetYaxis().SetTitle(ytitle)
    beff.GetYaxis().SetTitleOffset(1.35*beff.GetYaxis().GetTitleOffset())
    beff.SetTitle(title)
    beff.GetZaxis().SetRangeUser(0,1.0)
    beff.SetMarkerSize(1.4)
    # beff.Draw("colz,text,E")
    beff.Draw("colz,text")
    text.Draw("same")
    c.SetLogx()
    c.SaveAs(outDir+"/"+year+"/"+channel+"_bJetEff.pdf")

    xtitle = "p_{T}^{c-jet} [GeV]"
    ytitle = "|#eta_{c-jet}|"
    title = "Mistag rate of c-jets"
    ceff.GetXaxis().SetTitle(xtitle)
    ceff.GetYaxis().SetTitle(ytitle)
    ceff.GetYaxis().SetTitleOffset(1.35*ceff.GetYaxis().GetTitleOffset())
    ceff.SetTitle(title)
    ceff.GetZaxis().SetRangeUser(0, 1.0)
    ceff.SetMarkerSize(1.4)
    # ceff.Draw("colz,text,E")
    ceff.Draw("colz,text")
    text.Draw("same")
    c.SaveAs(outDir+"/"+year+"/"+channel+"_cJetEff.pdf")

    xtitle = "p_{T}^{l-jet} [GeV]"
    ytitle = "|#eta_{l-jet}|"
    title = "Mistag rate of l-jets"
    leff.GetXaxis().SetTitle(xtitle)
    leff.GetYaxis().SetTitle(ytitle)
    leff.GetYaxis().SetTitleOffset(1.35*leff.GetYaxis().GetTitleOffset())
    leff.GetZaxis().SetRangeUser(0, 1.0)
    leff.SetTitle(title)
    leff.SetMarkerSize(1.4)
    # leff.Draw("colz,text,E")
    leff.Draw("colz,text")
    text.Draw("same")
    c.SaveAs(outDir+"/"+year+"/"+channel+"_lJetEff.pdf")

    file.Close()






def getBtagEfficiencies(year):
    channel = ["ee", "emu", "mumu"]
    for ch in channel:
            getEff(year,ch)

def main(year):
    getBtagEfficiencies(year)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-y', '--year', dest='year', action='store', default="", help='Year to process', required=True)
    args = parser.parse_args()

    if not args.year in valid_years:
        raise ValueError('Year is not valid: ' + str(args.year))

    main(args.year)
