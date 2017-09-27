import ROOT
import numpy as n
import ppmUtils as utils

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(False);

NumberOfNtuples = 269

def SavePlot (tcanvas, name="plot", extensions=["png","pdf"]):
    utils.DrawCmsText(canvas,text="CMS Preliminary",textFont=62,textSize=0.035)
    for extension in extensions:
        tcanvas.SaveAs("RecoNpv/"+extension+"s/"+name+"."+extension)

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

tf = ROOT.TFile("../src/output_RecoNpv/hists.root")

h_all_rhos_ctr = GetHist(tf,"h_all_rhos_ctr", title="#rho_{Central-Neutral};#rho (GeV);Fraction", color=ROOT.kRed)
Normalize(h_all_rhos_ctr)
h_all_rhos_all = GetHist(tf,"h_all_rhos_all", title="#rho_{All};#rho (GeV);Fraction")
Normalize(h_all_rhos_all)

p_npv_rho_ctr = GetHist(tf,"p_npv_rho_ctr",title="#rho by Reco N_{PV};Reco N_{PV};Mean #rho (GeV)",color=ROOT.kRed,overflow=False)
p_npv_rho_all = GetHist(tf,"p_npv_rho_all",title="#rho by Reco N_{PV};Reco N_{PV};Mean #rho (GeV)",overflow=False)

p_npv_eff_ctr = GetHist(tf,"p_npv_eff_ctr",title="MiniRelIso < 0.2 Eff;Reco N_{PV};Efficiency",color=ROOT.kRed,overflow=False)
p_npv_eff_all = GetHist(tf,"p_npv_eff_all",title="MiniRelIso < 0.2 Eff;Reco N_{PV};Efficiency",overflow=False)

h_all_EAs_mini = GetHist(tf, "h_all_EAs_mini",title="All Eff Areas;Eff Area;Fraction",color=ROOT.kRed)
h_all_EAs_fc = GetHist(tf, "h_all_EAs_fc",title="All Eff Ares;Eff Area;Fraction")
Normalize(h_all_EAs_mini)
Normalize(h_all_EAs_fc)

h_median_bug_mini = GetHist(tf, "h_median_bug_mini", overflow=False,color=ROOT.kRed, title="Mini;Reco N_{PV};Median (GeV)")
h_median_unbug_mini = GetHist(tf, "h_median_unbug_mini", overflow=False, title="Mini;Reco N_{PV};Median (GeV)")
h_median_bug_fc = GetHist(tf, "h_median_bug_fc",overflow=False,color=ROOT.kRed,title="Fixed-Cone;Reco N_{PV};Median (GeV)")
h_median_unbug_fc = GetHist(tf,"h_median_unbug_fc",overflow=False,title="Fixed-Cone;Reco N_{PV};Median (GeV)")
# rescale by factor equal to number of samples merged
h_median_bug_mini.Scale(1./NumberOfNtuples)
h_median_unbug_mini.Scale(1./NumberOfNtuples)
h_median_bug_fc.Scale(1./NumberOfNtuples)
h_median_unbug_fc.Scale(1./NumberOfNtuples)

h_all_bug_mini = GetHist(tf, "h_all_bug_mini",color=ROOT.kRed,title="Mini;Correction (GeV);Fraction")
h_all_nobug_mini = GetHist(tf, "h_all_nobug_mini",title="Mini;Correction (GeV);Fraction")
h_all_bug_fc = GetHist(tf, "h_all_bug_fc",color=ROOT.kRed,title="Fixed-Cone;Correction (GeV);Fraction")
h_all_nobug_fc = GetHist(tf, "h_all_nobug_fc",title="Fixed-Cone;Correction (GeV);Fraction")
Normalize(h_all_bug_mini)
Normalize(h_all_nobug_mini)
Normalize(h_all_bug_fc)
Normalize(h_all_nobug_fc)

p_pt_bug_mini=GetHist(tf,"p_pt_bug_mini",title="Mini;Muon p_{T} (GeV);Mean Correction (GeV)",color=ROOT.kRed,overflow=False)
p_pt_nobug_mini=GetHist(tf,"p_pt_nobug_mini",title="Mini;Muon p_{T} (GeV);Mean Correction (GeV)",overflow=False)
p_pt_bug_fc=GetHist(tf,"p_pt_bug_fc",title="Fixed-Cone;Muon p_{T} (GeV);Mean Correction (GeV)",color=ROOT.kRed,overflow=False)
p_pt_nobug_fc=GetHist(tf,"p_pt_nobug_fc",title="Fixed-Cone;Muon p_{T} (GeV);Mean Correction (GeV)",overflow=False)

p_npv_incl_bug_mini=GetHist(tf,"p_npv_incl_bug_mini",title="Mini;Reco N_{PV};Mean (GeV)",overflow=False,color=ROOT.kRed)
p_npv_incl_nobug_mini=GetHist(tf,"p_npv_incl_nobug_mini",title="Mini;Reco N_{PV};Mean (GeV)",overflow=False)
p_npv_incl_neut_mini=GetHist(tf,"p_npv_incl_neut_mini",title="Mini;Reco N_{PV};Mean (GeV)",color=ROOT.kGreen,overflow=False)

p_npv_0p8_bug_mini=GetHist(tf,"p_npv_0p8_bug_mini",title="Mini (|#eta|<0.8);Reco N_{PV};Mean (GeV)",overflow=False,color=ROOT.kRed)
p_npv_0p8_nobug_mini=GetHist(tf,"p_npv_0p8_nobug_mini",title="Mini (|#eta|<0.8);Reco N_{PV};Mean (GeV)",overflow=False)
p_npv_0p8_neut_mini=GetHist(tf,"p_npv_0p8_neut_mini",title="Mini (|#eta|<0.8);Reco N_{PV};Mean (GeV)",color=ROOT.kGreen,overflow=False)

p_npv_1p3_bug_mini=GetHist(tf,"p_npv_1p3_bug_mini",title="Mini (0.8<|#eta|<1.3);Reco N_{PV};Mean (GeV)",overflow=False,color=ROOT.kRed)
p_npv_1p3_nobug_mini=GetHist(tf,"p_npv_1p3_nobug_mini",title="Mini (0.8<|#eta|<1.3);Reco N_{PV};Mean (GeV)",overflow=False)
p_npv_1p3_neut_mini=GetHist(tf,"p_npv_1p3_neut_mini",title="Mini (0.8<|#eta|<1.3);Reco N_{PV};Mean (GeV)",color=ROOT.kGreen,overflow=False)

p_npv_2p0_bug_mini=GetHist(tf,"p_npv_2p0_bug_mini",title="Mini (1.3<|#eta|<2.0);Reco N_{PV};Mean (GeV)",overflow=False,color=ROOT.kRed)
p_npv_2p0_nobug_mini=GetHist(tf,"p_npv_2p0_nobug_mini",title="Mini (1.3<|#eta|<2.0);Reco N_{PV};Mean (GeV)",overflow=False)
p_npv_2p0_neut_mini=GetHist(tf,"p_npv_2p0_neut_mini",title="Mini (1.3<|#eta|<2.0);Reco N_{PV};Mean (GeV)",color=ROOT.kGreen,overflow=False)

p_npv_2p2_bug_mini=GetHist(tf,"p_npv_2p2_bug_mini",title="Mini (2.0<|#eta|<2.2);Reco N_{PV};Mean (GeV)",overflow=False,color=ROOT.kRed)
p_npv_2p2_nobug_mini=GetHist(tf,"p_npv_2p2_nobug_mini",title="Mini (2.0<|#eta|<2.2);Reco N_{PV};Mean (GeV)",overflow=False)
p_npv_2p2_neut_mini=GetHist(tf,"p_npv_2p2_neut_mini",title="Mini (2.0<|#eta|<2.2);Reco N_{PV};Mean (GeV)",color=ROOT.kGreen,overflow=False)

p_npv_2p5_bug_mini=GetHist(tf,"p_npv_2p5_bug_mini",title="Mini (2.2<|#eta<2.5);Reco N_{PV};Mean (GeV)",overflow=False,color=ROOT.kRed)
p_npv_2p5_nobug_mini=GetHist(tf,"p_npv_2p5_nobug_mini",title="Mini (2.2<|#eta<2.5);Reco N_{PV};Mean (GeV)",overflow=False)
p_npv_2p5_neut_mini=GetHist(tf,"p_npv_2p5_neut_mini",title="Mini (2.2<|#eta<2.5);Reco N_{PV};Mean (GeV)",color=ROOT.kGreen,overflow=False)

p_npv_incl_bug_fc=GetHist(tf,"p_npv_incl_bug_fc",title="Fixed-Cone;Reco N_{PV};Mean (GeV)",overflow=False,color=ROOT.kRed)
p_npv_incl_nobug_fc=GetHist(tf,"p_npv_incl_nobug_fc",title="Fixed-Cone;Reco N_{PV};Mean (GeV)",overflow=False)
p_npv_incl_neut_fc=GetHist(tf,"p_npv_incl_neut_fc",title="Fixed-Cone;Reco N_{PV};Mean (GeV)",color=ROOT.kGreen,overflow=False)

p_npv_0p8_bug_fc=GetHist(tf,"p_npv_0p8_bug_fc",title="Fixed-Cone (|#eta|<0.8);Reco N_{PV};Mean (GeV)",overflow=False,color=ROOT.kRed)
p_npv_0p8_nobug_fc=GetHist(tf,"p_npv_0p8_nobug_fc",title="Fixed-Cone (|#eta|<0.8);Reco N_{PV};Mean (GeV)",overflow=False)
p_npv_0p8_neut_fc=GetHist(tf,"p_npv_0p8_neut_fc",title="Fixed-Cone (|#eta|<0.8);Reco N_{PV};Mean (GeV)",color=ROOT.kGreen,overflow=False)

p_npv_1p3_bug_fc=GetHist(tf,"p_npv_1p3_bug_fc",title="Fixed-Cone (0.8<|#eta|<1.3);Reco N_{PV};Mean (GeV)",overflow=False,color=ROOT.kRed)
p_npv_1p3_nobug_fc=GetHist(tf,"p_npv_1p3_nobug_fc",title="Fixed-Cone (0.8<|#eta|<1.3);Reco N_{PV};Mean (GeV)",overflow=False)
p_npv_1p3_neut_fc=GetHist(tf,"p_npv_1p3_neut_fc",title="Fixed-Cone (0.8<|#eta|<1.3);Reco N_{PV};Mean (GeV)",color=ROOT.kGreen,overflow=False)

p_npv_2p0_bug_fc=GetHist(tf,"p_npv_2p0_bug_fc",title="Fixed-Cone (1.3<|#eta|<2.0);Reco N_{PV};Mean (GeV)",overflow=False,color=ROOT.kRed)
p_npv_2p0_nobug_fc=GetHist(tf,"p_npv_2p0_nobug_fc",title="Fixed-Cone (1.3<|#eta|<2.0);Reco N_{PV};Mean (GeV)",overflow=False)
p_npv_2p0_neut_fc=GetHist(tf,"p_npv_2p0_neut_fc",title="Fixed-Cone (1.3<|#eta|<2.0);Reco N_{PV};Mean (GeV)",color=ROOT.kGreen,overflow=False)

p_npv_2p2_bug_fc=GetHist(tf,"p_npv_2p2_bug_fc",title="Fixed-Cone (2.0<|#eta|<2.2);Reco N_{PV};Mean (GeV)",overflow=False,color=ROOT.kRed)
p_npv_2p2_nobug_fc=GetHist(tf,"p_npv_2p2_nobug_fc",title="Fixed-Cone (2.0<|#eta|<2.2);Reco N_{PV};Mean (GeV)",overflow=False)
p_npv_2p2_neut_fc=GetHist(tf,"p_npv_2p2_neut_fc",title="Fixed-Cone (2.0<|#eta|<2.2);Reco N_{PV};Mean (GeV)",color=ROOT.kGreen,overflow=False)

p_npv_2p5_bug_fc=GetHist(tf,"p_npv_2p5_bug_fc",title="Fixed-Cone (2.2<|#eta|<2.5);Reco N_{PV};Mean (GeV)",overflow=False,color=ROOT.kRed)
p_npv_2p5_nobug_fc=GetHist(tf,"p_npv_2p5_nobug_fc",title="Fixed-Cone (2.2<|#eta|<2.5);Reco N_{PV};Mean (GeV)",overflow=False)
p_npv_2p5_neut_fc=GetHist(tf,"p_npv_2p5_neut_fc",title="Fixed-Cone (2.2<|#eta|<2.5);Reco N_{PV};Mean (GeV)",color=ROOT.kGreen,overflow=False)

p_npv_incl_corr_bug_mini = GetHist(tf,"p_npv_incl_corr_bug_mini",title="Mini;Reco N_{PV};Mean (GeV)",color=ROOT.kRed,overflow=False)
p_npv_incl_corr_nobug_mini=GetHist(tf,"p_npv_incl_corr_nobug_mini",title="Mini;Reco N_{PV};Mean (GeV)",overflow=False)
p_npv_incl_corr_bug_fc = GetHist(tf,"p_npv_incl_corr_bug_fc",title="Fixed-Cone;Reco N_{PV};Mean (GeV)",color=ROOT.kRed,overflow=False)
p_npv_incl_corr_nobug_fc=GetHist(tf,"p_npv_incl_corr_nobug_fc",title="Fixed-Cone;Reco N_{PV};Mean (GeV)",overflow=False)


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

ROOT.gStyle.SetLegendBorderSize(0) # borderless legend
tl = ROOT.TLegend(0.15,0.70,0.60,0.85)
tl.SetMargin(0.1)

h_all_rhos_ctr.Draw("hist")
SavePlot(canvas,"rho_ctr")

h_all_rhos_all.Draw("hist")
SavePlot(canvas,"rho_all")

p_npv_rho_all.SetMinimum(0.0)
p_npv_rho_all.Draw("")
p_npv_rho_ctr.Draw("same")
tl.AddEntry(p_npv_rho_ctr,"#rho_{CN}")
tl.AddEntry(p_npv_rho_all,"#rho_{All}")
tl.Draw()
SavePlot(canvas,"rho_npv")

tl2 = ROOT.TLegend(0.15,0.15,0.60,0.30)
tl2.SetMargin(0.1)
p_npv_eff_all.SetMinimum(0.9)
p_npv_eff_all.Draw("")
p_npv_eff_ctr.Draw("same")
tl2.AddEntry(p_npv_rho_ctr,"#rho_{CN}")
tl2.AddEntry(p_npv_rho_all,"#rho_{All}")
tl2.Draw()
SavePlot(canvas,"eff_npv")

pads[0].SetLeftMargin(0.14)
p_npv_eff_all.Divide(p_npv_eff_ctr)
p_npv_eff_all.SetLineColor(ROOT.kBlack)
p_npv_eff_all.SetTitle("Efficiency Ratio;Reco N_{PV};Eff_{All} / Eff_{CN}")
p_npv_eff_all.GetYaxis().SetTitleOffset(1.9)
p_npv_eff_all.SetMinimum(0.98)
p_npv_eff_all.SetMaximum(1.02)
p_npv_eff_all.Draw("hist")
SavePlot(canvas,"eff_npv_ratio")
pads[0].SetLeftMargin(0.12)

h_all_EAs_mini.Draw("hist")
h_all_EAs_fc.Draw("same hist")
tl.Clear()
tl.AddEntry(h_all_EAs_mini,"Mini")
tl.AddEntry(h_all_EAs_fc,"Fixed-Cone")
tl.Draw()
SavePlot(canvas,"EAs")

h_median_unbug_mini.SetMinimum(0.0)
h_median_unbug_mini.Draw("")
h_median_bug_mini.Draw("same")
tl.Clear()
tl.AddEntry(h_median_bug_mini,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(h_median_unbug_mini,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"median_mini")

h_median_unbug_fc.SetMinimum(0.0)
h_median_unbug_fc.Draw("")
h_median_bug_fc.Draw("same")
tl.Clear()
tl.AddEntry(h_median_bug_fc,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(h_median_unbug_fc,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"median_fc")

h_all_bug_mini.Draw("")
h_all_nobug_mini.Draw("same")
tl.Clear()
tl.AddEntry(h_all_bug_mini,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(h_all_nobug_mini,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"all_mini")

h_all_bug_fc.Draw("")
h_all_nobug_fc.Draw("same")
tl.Clear()
tl.AddEntry(h_all_bug_fc,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(h_all_nobug_fc,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"all_fc")

p_pt_nobug_mini.SetMinimum(0.0)
p_pt_nobug_mini.Draw("")
p_pt_bug_mini.Draw("same")
tl.Clear()
tl.AddEntry(p_pt_bug_mini,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_pt_nobug_mini,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"pt_mini")

p_pt_nobug_fc.SetMinimum(0.0)
p_pt_nobug_fc.Draw("")
p_pt_bug_fc.Draw("same")
tl.Clear()
tl.AddEntry(p_pt_bug_fc,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_pt_nobug_fc,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"pt_fc")

p_npv_incl_neut_mini.SetMinimum(0.0)
p_npv_incl_neut_mini.Draw("")
p_npv_incl_nobug_mini.Draw("same")
p_npv_incl_bug_mini.Draw("same")
tl.Clear()
tl.AddEntry(p_npv_incl_neut_mini,"Sum(Neutrals)")
tl.AddEntry(p_npv_incl_bug_mini,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_npv_incl_nobug_mini,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"npv_incl_mini")

p_npv_0p8_neut_mini.SetMinimum(0.0)
p_npv_0p8_neut_mini.Draw("")
p_npv_0p8_nobug_mini.Draw("same")
p_npv_0p8_bug_mini.Draw("same")
tl.Clear()
tl.AddEntry(p_npv_0p8_neut_mini,"Sum(Neutrals)")
tl.AddEntry(p_npv_0p8_bug_mini,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_npv_0p8_nobug_mini,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"npv_0p8_mini")

p_npv_1p3_neut_mini.SetMinimum(0.0)
p_npv_1p3_neut_mini.Draw("")
p_npv_1p3_nobug_mini.Draw("same")
p_npv_1p3_bug_mini.Draw("same")
tl.Clear()
tl.AddEntry(p_npv_1p3_neut_mini,"Sum(Neutrals)")
tl.AddEntry(p_npv_1p3_bug_mini,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_npv_1p3_nobug_mini,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"npv_1p3_mini")

p_npv_2p0_neut_mini.SetMinimum(0.0)
p_npv_2p0_neut_mini.Draw("")
p_npv_2p0_nobug_mini.Draw("same")
p_npv_2p0_bug_mini.Draw("same")
tl.Clear()
tl.AddEntry(p_npv_2p0_neut_mini,"Sum(Neutrals)")
tl.AddEntry(p_npv_2p0_bug_mini,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_npv_2p0_nobug_mini,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"npv_2p0_mini")

p_npv_2p2_neut_mini.SetMinimum(0.0)
p_npv_2p2_neut_mini.Draw("")
p_npv_2p2_nobug_mini.Draw("same")
p_npv_2p2_bug_mini.Draw("same")
tl.Clear()
tl.AddEntry(p_npv_2p2_neut_mini,"Sum(Neutrals)")
tl.AddEntry(p_npv_2p2_bug_mini,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_npv_2p2_nobug_mini,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"npv_2p2_mini")

p_npv_2p5_neut_mini.SetMinimum(0.0)
p_npv_2p5_neut_mini.Draw("")
p_npv_2p5_nobug_mini.Draw("same")
p_npv_2p5_bug_mini.Draw("same")
tl.Clear()
tl.AddEntry(p_npv_2p5_neut_mini,"Sum(Neutrals)")
tl.AddEntry(p_npv_2p5_bug_mini,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_npv_2p5_nobug_mini,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"npv_2p5_mini")

p_npv_incl_neut_fc.SetMinimum(0.0)
p_npv_incl_neut_fc.Draw("")
p_npv_incl_nobug_fc.Draw("same")
p_npv_incl_bug_fc.Draw("same")
tl.Clear()
tl.AddEntry(p_npv_incl_neut_fc,"Sum(Neutrals)")
tl.AddEntry(p_npv_incl_bug_fc,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_npv_incl_nobug_fc,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"npv_incl_fc")

p_npv_0p8_neut_fc.SetMinimum(0.0)
p_npv_0p8_neut_fc.Draw("")
p_npv_0p8_nobug_fc.Draw("same")
p_npv_0p8_bug_fc.Draw("same")
tl.Clear()
tl.AddEntry(p_npv_0p8_neut_fc,"Sum(Neutrals)")
tl.AddEntry(p_npv_0p8_bug_fc,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_npv_0p8_nobug_fc,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"npv_0p8_fc")

p_npv_1p3_neut_fc.SetMinimum(0.0)
p_npv_1p3_neut_fc.Draw("")
p_npv_1p3_nobug_fc.Draw("same")
p_npv_1p3_bug_fc.Draw("same")
tl.Clear()
tl.AddEntry(p_npv_1p3_neut_fc,"Sum(Neutrals)")
tl.AddEntry(p_npv_1p3_bug_fc,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_npv_1p3_nobug_fc,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"npv_1p3_fc")

p_npv_2p0_neut_fc.SetMinimum(0.0)
p_npv_2p0_neut_fc.Draw("")
p_npv_2p0_nobug_fc.Draw("same")
p_npv_2p0_bug_fc.Draw("same")
tl.Clear()
tl.AddEntry(p_npv_2p0_neut_fc,"Sum(Neutrals)")
tl.AddEntry(p_npv_2p0_bug_fc,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_npv_2p0_nobug_fc,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"npv_2p0_fc")

p_npv_2p2_neut_fc.SetMinimum(0.0)
p_npv_2p2_neut_fc.Draw("")
p_npv_2p2_nobug_fc.Draw("same")
p_npv_2p2_bug_fc.Draw("same")
tl.Clear()
tl.AddEntry(p_npv_2p2_neut_fc,"Sum(Neutrals)")
tl.AddEntry(p_npv_2p2_bug_fc,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_npv_2p2_nobug_fc,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"npv_2p2_fc")

p_npv_2p5_neut_fc.SetMinimum(0.0)
p_npv_2p5_neut_fc.Draw("")
p_npv_2p5_nobug_fc.Draw("same")
p_npv_2p5_bug_fc.Draw("same")
tl.Clear()
tl.AddEntry(p_npv_2p5_neut_fc,"Sum(Neutrals)")
tl.AddEntry(p_npv_2p5_bug_fc,"#rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_npv_2p5_nobug_fc,"#rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"npv_2p5_fc")

p_npv_incl_corr_bug_mini.SetMinimum(0.0)
p_npv_incl_corr_bug_mini.Draw("")
p_npv_incl_corr_nobug_mini.Draw("same")
tl.Clear()
tl.AddEntry(p_npv_incl_corr_bug_mini,"Sum(Neutrals) - #rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_npv_incl_corr_nobug_mini,"Sum(Neutrals) - #rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"corr_mini")

p_npv_incl_corr_bug_fc.SetMinimum(0.0)
p_npv_incl_corr_bug_fc.Draw("")
p_npv_incl_corr_nobug_fc.Draw("same")
tl.Clear()
tl.AddEntry(p_npv_incl_corr_bug_fc,"Sum(Neutrals) - #rho_{CN} x EAConstant(#eta)")
tl.AddEntry(p_npv_incl_corr_nobug_fc,"Sum(Neutrals) - #rho_{All} x EAConstant(#eta)")
tl.Draw()
SavePlot(canvas,"corr_fc")
