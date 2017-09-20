import ROOT
import numpy as n
import ppmUtils as utils

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(False);

def SavePlot (tcanvas, name="plot", extensions=["png","pdf"]):
    utils.DrawCmsText(canvas,text="CMS Preliminary",textFont=62,textSize=0.035)
    for extension in extensions:
        tcanvas.SaveAs(extension+"s/"+name+"."+extension)

def GetHist (tfile, name, overflow = True, title=None, color=None):
    hist = tfile.Get(name)
    #hist.SetMarkerStyle(20)
    #hist.SetMarkerSize(0.9)
    hist.SetLineWidth(3)
    if overflow:
        utils.PutOverflowInLastBin(hist)
        err0 = hist.GetBinError(0)
        err1 = hist.GetBinError(1)
        hist.SetBinContent(0,hist.GetBinContent(0)+hist.GetBinContent(1))
        hist.SetBinError(0,ROOT.TMath.Sqrt(err0**2 + err1**2))
    if title != None:
        hist.SetTitle(title)
    if color != None:
        hist.SetLineColor(color)
    hist.SetTitleOffset(1.6, "y")
    return hist

def Normalize (hist):
    integral = hist.Integral()
    hist.Scale(1./integral)
    return integral

tf = ROOT.TFile("../src/merged10x.root")

# TProfiles
p_pt_EA = GetHist(tf, "p_pt_EA", overflow=False, title="#DeltaR=0.3 EA Corrections by p_{T} of Candidate;Muon p_{T} (GeV);Average Absolute EA Correction (GeV)")
p_pt_dB = GetHist(tf, "p_pt_dB", overflow=False, title="#DeltaR=0.3 #Delta#beta Corrections by p_{T} of Candidate;Muon p_{T} (GeV);Average Absolute #Delta#beta Correction (GeV)")
p_npv_EA = GetHist(tf, "p_npv_EA", overflow=False, title="#DeltaR=0.3 Corrections by N_{PV};Number of Pileup Vertices;Average Absolute Correction (GeV)",color=ROOT.kRed)
p_npv_dB = GetHist(tf, "p_npv_dB", overflow=False)
p_npv_rhos_all = GetHist(tf, "p_npv_rho_all", overflow=False, title="#rho as Function of N_{PV};Number of Pileup Vertices;Average #rho (GeV)",color=ROOT.kRed)
p_npv_rhos_ctr = GetHist(tf, "p_npv_rho_ctr", overflow=False)

h_npv_over_EA = GetHist(tf, "h_npv_over_EA", overflow=False, title="Fraction of Overcorrections by N_{PV} (Neut Iso < Correction);Number of Pileup Vertices;Overcorrected Fraction",color=ROOT.kRed)
h_npv_over_dB = GetHist(tf, "h_npv_over_dB", overflow=False)

########
# Mini #
########

h_all_rhos_ctr = GetHist(tf,"h_all_rhos_ctr", title="Central #rho Distribution;#rho (GeV);Fraction")
Normalize(h_all_rhos_ctr)

h_all_EAs = GetHist(tf, "h_all_EAs")
integral = Normalize(h_all_EAs)
h_all_EAs.SetTitle("All Mini Effective Areas ("+str(integral)+" muons);Effective Area;Fraction")

h_median_EA = GetHist(tf, "h_median_EA", overflow=False)
h_median_dB = GetHist(tf, "h_median_dB", overflow=False, title="Median Corrections as Function of N_{PV};N_{PV};Median Absolute MiniIso Correction (GeV)", color=ROOT.kRed)


h_all_EAcorrs = GetHist(tf, "h_all_EAcorrs")
h_all_dBcorrs = GetHist(tf, "h_all_dBcorrs")
integral = int(Normalize(h_all_EAcorrs))
h_all_EAcorrs.SetTitle("All MiniIso Effective Area Corrections ("+str(integral)+" muons);Absolute Effective Area Correction (GeV);Fraction")
integral = int(Normalize(h_all_dBcorrs))
h_all_dBcorrs.SetTitle("All MiniIso #Delta#beta Corrections ("+str(integral)+" muons);Absolute #Delta#beta Correction (GeV);Fraction")

h_median_EA_rel = GetHist(tf, "h_median_EA_rel", overflow=False)
h_median_dB_rel = GetHist(tf, "h_median_dB_rel", overflow=False, title="Median Corrections as Function of N_{PV};N_{PV};Median Relative MiniIso Correction", color=ROOT.kRed)

h_all_EAcorrs_rel = GetHist(tf, "h_all_EAcorrs_rel")
h_all_dBcorrs_rel = GetHist(tf, "h_all_dBcorrs_rel")
integral = int(Normalize(h_all_EAcorrs_rel))
h_all_EAcorrs_rel.SetTitle("All Relative MiniIso Effective Area Corrections ("+str(integral)+" muons);Effective Area Correction;Fraction")
h_all_EAcorrs_rel.Scale(1./integral)
integral = int(Normalize(h_all_dBcorrs_rel))
h_all_dBcorrs_rel.SetTitle("All Relative MiniIso #Delta#beta Corrections ("+str(integral)+" muons);#Delta#beta Correction;Fraction")

h_diff_0p8_EA = GetHist(tf, "h_diff_0p8_EA", title="Abs Neut MiniIso - Effective Area Correction;Difference (GeV);Fraction")
h_diff_1p3_EA = GetHist(tf, "h_diff_1p3_EA", color=ROOT.kRed)
h_diff_2p0_EA = GetHist(tf, "h_diff_2p0_EA", color=ROOT.kGreen)
h_diff_2p2_EA = GetHist(tf, "h_diff_2p2_EA", color=ROOT.kCyan)
h_diff_2p5_EA = GetHist(tf, "h_diff_2p5_EA", color=ROOT.kMagenta)
Normalize(h_diff_0p8_EA)
h_diff_0p8_EA.SetMaximum(0.5)
Normalize(h_diff_1p3_EA)
Normalize(h_diff_2p0_EA)
Normalize(h_diff_2p2_EA)
Normalize(h_diff_2p5_EA)

h_diff_0p8_dB = GetHist(tf, "h_diff_0p8_dB", title="Abs Neut MiniIso - #Delta#beta Correction;Difference (GeV);Fraction")
h_diff_1p3_dB = GetHist(tf, "h_diff_1p3_dB", color=ROOT.kRed)
h_diff_2p0_dB = GetHist(tf, "h_diff_2p0_dB", color=ROOT.kGreen)
h_diff_2p2_dB = GetHist(tf, "h_diff_2p2_dB", color=ROOT.kCyan)
h_diff_2p5_dB = GetHist(tf, "h_diff_2p5_dB", color=ROOT.kMagenta)
Normalize(h_diff_0p8_dB)
h_diff_0p8_dB.SetMaximum(0.5)
Normalize(h_diff_1p3_dB)
Normalize(h_diff_2p0_dB)
Normalize(h_diff_2p2_dB)
Normalize(h_diff_2p5_dB)


############
# Standard #
############

h_all_rhos_all = GetHist(tf, "h_all_rhos_all", title="Full-Detector #rho Distribution;#rho (GeV);Fraction")
Normalize(h_all_rhos_all)

h_all_EAs_st = GetHist(tf, "h_all_EAs_st")

h_median_EA_st = GetHist(tf, "h_median_EA_st", overflow=False)
h_median_dB_st = GetHist(tf, "h_median_dB_st", overflow=False, title="Median Corrections as Function of N_{PV};N_{PV};Median Absolute #DeltaR = 0.3 Iso Correction (GeV)", color=ROOT.kRed)

h_all_EAcorrs_st = GetHist(tf, "h_all_EAcorrs_st")
h_all_dBcorrs_st = GetHist(tf, "h_all_dBcorrs_st", )
integral = int(Normalize(h_all_dBcorrs_st))
h_all_dBcorrs_st.SetTitle("All #DeltaR = 0.3 Iso #Delta#beta Corrections ("+str(integral)+" muons);#Delta#beta Correction (GeV);Fraction")
integral = int(Normalize(h_all_EAcorrs_st))
h_all_EAcorrs_st.SetTitle("All #DeltaR = 0.3 Iso Effective Area Corrections ("+str(integral)+" muons);Effective Area Correction (GeV);Fraction")

h_median_EA_rel_st = GetHist(tf, "h_median_EA_rel_st", overflow=False)
h_median_dB_rel_st = GetHist(tf, "h_median_dB_rel_st", title="Median Corrections as Function of N_{PV};N_{PV};Median Relative #DeltaR = 0.3 Iso Correction",color=ROOT.kRed, overflow=False)

h_all_EAcorrs_rel_st = GetHist(tf, "h_all_EAcorrs_rel_st")
h_all_dBcorrs_rel_st = GetHist(tf, "h_all_dBcorrs_rel_st")
integral = int(Normalize(h_all_EAcorrs_rel_st))
h_all_EAcorrs_rel_st.SetTitle("All #DeltaR = 0.3 Iso Effective Area Corrections ("+str(integral)+" muons);Effective Area Correction;Fraction")
integral = int(Normalize(h_all_dBcorrs_rel_st))
h_all_dBcorrs_rel_st.SetTitle("All #DeltaR = 0.3 Iso #Delta#beta Corrections ("+str(integral)+" muons);#Delta#beta Correction;Fraction")

h_diff_0p8_EA_st = GetHist(tf, "h_diff_0p8_EA_st", title="Abs Neut #DeltaR=0.3 Iso - Effective Area Correction;Difference (GeV);Fraction")
h_diff_1p3_EA_st = GetHist(tf, "h_diff_1p3_EA_st", color=ROOT.kRed)
h_diff_2p0_EA_st = GetHist(tf, "h_diff_2p0_EA_st", color=ROOT.kGreen)
h_diff_2p2_EA_st = GetHist(tf, "h_diff_2p2_EA_st", color=ROOT.kCyan)
h_diff_2p5_EA_st = GetHist(tf, "h_diff_2p5_EA_st", color=ROOT.kMagenta)
Normalize(h_diff_0p8_EA_st)
h_diff_0p8_EA_st.SetMaximum(0.3)
Normalize(h_diff_1p3_EA_st)
Normalize(h_diff_2p0_EA_st)
Normalize(h_diff_2p2_EA_st)
Normalize(h_diff_2p5_EA_st)

h_diff_0p8_dB_st = GetHist(tf, "h_diff_0p8_dB_st", title="Abs Neut #DeltaR=0.3 Iso - #Delta#beta Correction;Difference (GeV);Fraction")
h_diff_1p3_dB_st = GetHist(tf, "h_diff_1p3_dB_st", color=ROOT.kRed)
h_diff_2p0_dB_st = GetHist(tf, "h_diff_2p0_dB_st", color=ROOT.kGreen)
h_diff_2p2_dB_st = GetHist(tf, "h_diff_2p2_dB_st", color=ROOT.kCyan)
h_diff_2p5_dB_st = GetHist(tf, "h_diff_2p5_dB_st", color=ROOT.kMagenta)
Normalize(h_diff_0p8_dB_st)
h_diff_0p8_dB_st.SetMaximum(0.3)
Normalize(h_diff_1p3_dB_st)
Normalize(h_diff_2p0_dB_st)
Normalize(h_diff_2p2_dB_st)
Normalize(h_diff_2p5_dB_st)


#################

##############
# Make Plots #
##############

#################

canvas = ROOT.TCanvas()
canvas.SetCanvasSize(700,700)
canvas.SetTicks(1,2)
pads = [canvas]
pads[0].SetLeftMargin(0.12)
pads[0].SetTopMargin(0.12)
pads[0].SetRightMargin(0.12)
pads[0].cd()

tl = ROOT.TLegend(0.20,0.75,0.40,0.85)

tl.AddEntry(h_median_dB,"Median #Delta#beta")
tl.AddEntry(h_median_EA,"Median Effective Area")
h_median_dB.Draw()
h_median_EA.Draw("same")
tl.Draw()
SavePlot(canvas, "median_corrs")

canvas.SetLogy(True)

h_all_EAcorrs.Draw("hist")
SavePlot(canvas, "all_EAcorrs")
h_all_dBcorrs.Draw("hist")
SavePlot(canvas, "all_dBcorrs")

#Relative 

canvas.SetLogy(False)
tl.Clear()
tl.AddEntry(h_median_dB_rel,"Median #Delta#beta")
tl.AddEntry(h_median_EA_rel,"Median Effective Area")
h_median_dB_rel.Draw()
h_median_EA_rel.Draw("same")
tl.Draw()
SavePlot(canvas, "median_corrs_rel")
canvas.SetLogy(True)

h_all_EAcorrs_rel.Draw("hist")
SavePlot(canvas, "all_EAcorrs_rel")
h_all_dBcorrs_rel.Draw("hist")
SavePlot(canvas, "all_dBcorrs_rel")

############
# Standard #
############

canvas.SetLogy(False)
tl.Clear()
tl.AddEntry(h_median_dB_st,"Median #Delta#beta")
tl.AddEntry(h_median_EA_st,"Median Effective Area")
h_median_dB_st.Draw()
h_median_EA_st.Draw("same")
tl.Draw()
SavePlot(canvas, "median_corrs_st")
canvas.SetLogy(True)

h_all_EAcorrs_st.Draw()
SavePlot(canvas, "all_EAcorrs_st")
h_all_dBcorrs_st.Draw()
SavePlot(canvas, "all_dBcorrs_st")

#Relative 

canvas.SetLogy(False)
tl.Clear()
tl.AddEntry(h_median_dB_rel_st,"Median #Delta#beta")
tl.AddEntry(h_median_EA_rel_st,"Median Effective Area")
h_median_dB_rel_st.Draw()
h_median_EA_rel_st.Draw("same")
tl.Draw()
SavePlot(canvas,"median_corrs_rel_st")
canvas.SetLogy(True)

h_all_EAcorrs_rel_st.Draw("hist")
SavePlot(canvas, "all_EAcorrs_rel_st")
h_all_dBcorrs_rel_st.Draw("hist")
SavePlot(canvas, "all_dBcorrs_rel_st")

# rhos and TProfiles

canvas.SetLogy(False)

p_pt_EA.Draw()
SavePlot(canvas,"p_pt_EA")
p_pt_dB.Draw()
SavePlot(canvas, "p_pt_dB")
tl.Clear()
tl.AddEntry(p_npv_EA,"Effective Area")
tl.AddEntry(p_npv_dB,"#Delta#beta")
p_npv_EA.SetMaximum(3)
p_npv_EA.Draw()
p_npv_dB.Draw("same")
tl.Draw()
SavePlot(canvas, "p_npv")
h_all_rhos_ctr.Draw("hist")
SavePlot(canvas, "rhos_ctr")
h_all_rhos_all.Draw("hist")
SavePlot(canvas, "rhos_all")

p_npv_rhos_all.SetMinimum(0)
tl.Clear()
tl.AddEntry(p_npv_rhos_all,"Full-Detector #rho")
tl.AddEntry(p_npv_rhos_ctr,"Central #rho")
p_npv_rhos_all.Draw()
p_npv_rhos_ctr.Draw("same")
tl.Draw()
SavePlot(canvas, "p_npv_rhos")

h_npv_over_EA.SetMaximum(0.8)
h_npv_over_EA.SetMinimum(0.4)
tl.Clear()
tl.AddEntry(h_npv_over_EA,"Effective Area")
tl.AddEntry(h_npv_over_dB,"#Delta#beta")
h_npv_over_EA.Draw()
h_npv_over_dB.Draw("same")
tl.Draw()
SavePlot(canvas, "overcorrection")

tl.Clear()
tl.AddEntry(h_diff_0p8_EA,"|#eta| < 0.8")
tl.AddEntry(h_diff_1p3_EA,"0.8 < |#eta| < 1.3")
tl.AddEntry(h_diff_2p0_EA,"1.3 < |#eta| < 2.0")
tl.AddEntry(h_diff_2p2_EA,"2.0 < |#eta| < 2.2")
tl.AddEntry(h_diff_2p5_EA,"2.2 < |#eta| < 2.5")
h_diff_0p8_EA.Draw()
h_diff_1p3_EA.Draw("same")
h_diff_2p0_EA.Draw("same")
h_diff_2p2_EA.Draw("same")
h_diff_2p5_EA.Draw("same")
tl.Draw()
SavePlot(canvas, "diff_EA")

tl.Clear()
tl.AddEntry(h_diff_0p8_dB,"|#eta| < 0.8")
tl.AddEntry(h_diff_1p3_dB,"0.8 < |#eta| < 1.3")
tl.AddEntry(h_diff_2p0_dB,"1.3 < |#eta| < 2.0")
tl.AddEntry(h_diff_2p2_dB,"2.0 < |#eta| < 2.2")
tl.AddEntry(h_diff_2p5_dB,"2.2 < |#eta| < 2.5")
h_diff_0p8_dB.Draw()
h_diff_1p3_dB.Draw("same")
h_diff_2p0_dB.Draw("same")
h_diff_2p2_dB.Draw("same")
h_diff_2p5_dB.Draw("same")
tl.Draw()
SavePlot(canvas, "diff_dB")

tl.Clear()
tl.AddEntry(h_diff_0p8_EA_st,"|#eta| < 0.8")
tl.AddEntry(h_diff_1p3_EA_st,"0.8 < |#eta| < 1.3")
tl.AddEntry(h_diff_2p0_EA_st,"1.3 < |#eta| < 2.0")
tl.AddEntry(h_diff_2p2_EA_st,"2.0 < |#eta| < 2.2")
tl.AddEntry(h_diff_2p5_EA_st,"2.2 < |#eta| < 2.5")
h_diff_0p8_EA_st.Draw()
h_diff_1p3_EA_st.Draw("same")
h_diff_2p0_EA_st.Draw("same")
h_diff_2p2_EA_st.Draw("same")
h_diff_2p5_EA_st.Draw("same")
tl.Draw()
SavePlot(canvas, "diff_EA_st")

tl.Clear()
tl.AddEntry(h_diff_0p8_dB_st,"|#eta| < 0.8")
tl.AddEntry(h_diff_1p3_dB_st,"0.8 < |#eta| < 1.3")
tl.AddEntry(h_diff_2p0_dB_st,"1.3 < |#eta| < 2.0")
tl.AddEntry(h_diff_2p2_dB_st,"2.0 < |#eta| < 2.2")
tl.AddEntry(h_diff_2p5_dB_st,"2.2 < |#eta| < 2.5")
h_diff_0p8_dB_st.Draw()
h_diff_1p3_dB_st.Draw("same")
h_diff_2p0_dB_st.Draw("same")
h_diff_2p2_dB_st.Draw("same")
h_diff_2p5_dB_st.Draw("same")
tl.Draw()
SavePlot(canvas, "diff_dB_st")
