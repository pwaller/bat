#! /usr/bin/env python

import numpy as np
import ROOT as R
import numpy as NP

from array import array
from math import sqrt
from rootwait import wait
from sys import argv

from scipy.optimize import bisect

from array import array
import re


from rootwait import wait

import ROOT as R

def mkarray(data): return array('d', data)

def get_mass(i):
    lo, hi, n = 400, 3000, 200
    delta = (hi-lo)/n
    return lo + delta*i

def get_quantiles(h, quantiles=None):
    """
    Return the value on the x axis of h for given quantiles
    """
    quantiles = quantiles or [0.004677, 0.15729, 0.5, 0.842701, 0.995322]
                             # erfc(2), erfc(1), 0.5, erf(1), erf(2)
    quantiles = array('d', quantiles)
    result = array('d', [0]*len(quantiles))
    h.GetQuantiles(len(quantiles), result, quantiles)
    return list(result)

def get_limit(filename):
    """
    Return the distribution of limits over pseudoexperiments in the root file 
    """
    f = R.TFile(filename)
    if f.IsZombie(): return None
    t = f.ensemble_test
    if not t: return None
    t.SetEstimate(t.GetEntries())
    
    
    t.Draw("95quantile_marginalized_1 >> hist", "", "para goff")
    h = R.gDirectory.Get("hist")
    quantiles = get_quantiles(h)
    h.SetDirectory(None)
    h.Delete()
    
    regex = re.compile("zprime_ensembles_mass([0-9]+)_run[0-9]+_gg")
    index = int(regex.search(filename).groups()[0])
    
    return get_mass(index), quantiles


def make_uncertaintyband(data):
    """
    Given quantiles, return a TGraphAsymmErrors representing the 2sig and 1sig
    uncertainty bands, and a TGraph for the median value.
    """
    masses, values = zip(*data)
    masses = mkarray(masses)
    sig2dpos, sig1dpos, med, sig1upos, sig2upos = map(mkarray, zip(*values))
    zeros = mkarray([0]*len(masses))
    
    sig2d = mkarray(+ np.array(med) - np.array(sig2dpos))
    sig2u = mkarray(- np.array(med) + np.array(sig2upos))
    sig1d = mkarray(+ np.array(med) - np.array(sig1dpos))
    sig1u = mkarray(- np.array(med) + np.array(sig1upos))
    
    g2sig = R.TGraphAsymmErrors(len(masses), masses, med, zeros, zeros, sig2d, sig2u)
    g1sig = R.TGraphAsymmErrors(len(masses), masses, med, zeros, zeros, sig1d, sig1u)
    exptd = R.TGraph(len(masses), masses, med)
    
    g2sig.SetFillColor(R.kYellow)
    g1sig.SetFillColor(R.kGreen)
    exptd.SetLineStyle(2)
    
    positions = sig2dpos, sig1dpos, med, sig1upos, sig2upos
    
    return g2sig, g1sig, exptd #, positions
    
def plot_data_limit(input_file):

    with open(input_file) as fd:
        content = [x.split() for x in fd.read().split("\n") if x]
        
    content = [x[1:] for x in content if x[0] == "limit"]

    data = [map(float, x) for x in content]
    data = [(mass*1000, x) for mass, x in data]
    data.sort()
        
    xs, ys = map(mkarray, zip(*data))
    
    c = R.TCanvas()
    c.SetLogy()
    
    g_excl_xs = g = R.TGraph(len(xs), xs, ys)
    g.SetMarkerStyle(7)
    g.GetYaxis().SetRangeUser(5e-4, 1.5)
    
    return g_excl_xs

def make_theory_functions(theory_xs, kvalues=[0.01, 0.03, 0.05, 0.1, 0.15, 0.2]):
    fs = []
    for k in kvalues:
        fnew = theory_xs.Clone()
        fs.append(fnew)
        fnew.SetParameter(5, k)
    return fs

def constrain_function(fn):
    """
    SetRange(xlo, xhi) given ylo, yhi
    """
    if not R.gPad.func(): return
    
    lo, hi = R.gPad.GetUymin(), R.gPad.GetUymax()
    if R.gPad.GetLogy():
        lo, hi = 10**lo, 10**hi
    
    lox = bisect(lambda x: fn(x) - lo, 0, 1e4)
    hix = bisect(lambda x: fn(x) - hi, 0, 1e4)
    if lox > hix:
        lox, hix = hix, lox
    #d = (hix - lox) / 10
    #lox += d
    #hix -= d
    fn.SetRange(lox, hix)

def get_values_in_plane(theory_xs, graph):
    spline = R.TSpline5("spline", graph)
    base_km = theory_xs.GetParameter(5)
    
    def compute_km_theoryintersect_at_mass(measured_xs, mass):
        measured_point = measured_xs(mass)
        theory_point = theory_xs.Eval(mass)
        if measured_point < 0: return 0
        km = base_km * sqrt(measured_point / theory_point)
        return km
    
    data = [(mass, compute_km_theoryintersect_at_mass(spline.Eval, mass))
            for mass in NP.linspace(400, 3000, 200)]
    return data

def graph_set_points(g, data):        
    g.Set(len(data))
    for i, (x, y) in enumerate(data):
        g.SetPoint(i, x, y)

def graph_to_theory_plane(theory_xs, g):
    if isinstance(g, R.TGraphAsymmErrors):
        g1, g2 = split_asym_graph_errors(g)        
        graph_set_points(g1, get_values_in_plane(theory_xs, g1))
        graph_set_points(g2, get_values_in_plane(theory_xs, g2))
        new_g = g.Clone()
        asym_graph_from_split(new_g, g1, g2)
        return new_g
    
    new_g = g.Clone()    
    values = get_values_in_plane(theory_xs, g)
    new_g.Set(len(values))
    for i, (x, y) in enumerate(values):
        new_g.SetPoint(i, x, y)
    return new_g

def asym_graph_from_split(new_g, glow, gup):
    n = glow.GetN()
    new_g.Set(n)
    for i in xrange(n):
        x, ylo, yup = glow.GetX()[i], glow.GetY()[i], gup.GetY()[i]
        mid = (ylo + yup) / 2
        new_g.SetPoint(i, x, mid)
        new_g.SetPointError(i, 0, 0, mid - ylo, yup - mid)

def split_asym_graph_errors(asym_graph):
    """
    Split a TGraphAsymmErrors into two `TGraph's
    """
    xs, ys = asym_graph.GetX(), asym_graph.GetY()
    ylos, yups = asym_graph.GetEYlow(), asym_graph.GetEYhigh()
    n = asym_graph.GetN()
    data = [(xs[i], ys[i] + yups[i], ys[i] - ylos[i]) for i in xrange(n)]
    xs, ys1, ys2 = map(mkarray, zip(*data))
    return R.TGraph(n, xs, ys1), R.TGraph(n, xs, ys2)

def get_text_pixel_width(pad, text):
    return pad.XtoAbsPixel(text.GetXsize() + pad.AbsPixeltoX(0))
    
def get_tf1_text_pad_angle(pad, func, text, start_x):
    
    from math import hypot, atan2, degrees
    text_width = get_text_pixel_width(pad, text)
    
    start_y = func.Eval(start_x)
    
    def hypotenuse_at(new_x, target_length=text_width):
        """
        Calculate hypotenuse at new_x in pixel co-ordinates
        """
        new_y = func.Eval(new_x)
        px = pad.XtoPixel(pad.XtoPad(start_x)) - pad.XtoPixel(pad.XtoPad(new_x))
        py = pad.YtoPixel(pad.YtoPad(start_y)) - pad.YtoPixel(pad.YtoPad(new_y))
        hyp = hypot(px, py)
        if new_x < start_x: hyp *= -1
        ret = hyp - target_length
        return ret
    
    final_x = bisect(hypotenuse_at, start_x, pad.PadtoX(pad.GetX2()))
    final_y = func.Eval(final_x)
    
    px = pad.XtoPixel(pad.XtoPad(final_x)) - pad.XtoPixel(pad.XtoPad(start_x))
    py = pad.YtoPixel(pad.YtoPad(final_y)) - pad.YtoPixel(pad.YtoPad(start_y))
    
    fn = text._keep_alive = func.Clone()
    fn.SetLineWidth(2)
    fn.SetLineColor(R.kWhite) #kBlue)
    fn.SetRange(start_x - (final_x - start_x)*0.1, final_x + (final_x - start_x)*0.5)
    fn.Draw("same")
    
    return degrees(atan2(-py, px))

def main():
    from sys import argv
    
    # Load limits from files specified on command line
    data = [x for x in [get_limit(filename) for filename in argv[1:]] if x]
    data.sort()
        
    g_2sig, g_1sig, g_exptd = make_uncertaintyband(data)
    
    c = R.TCanvas()
    
    # Expectation
    mg = R.TMultiGraph("excl", "")
    mg.Add(g_2sig, "e3")
    mg.Add(g_1sig, "e3")
    
    g_exptd.SetMarkerStyle(7)
    mg.Add(g_exptd, "CP")
    
    # Data
    g_excl_xs = plot_data_limit("full_log_no_ensemble")
    mg.Add(g_excl_xs, "CP")
    
    c.SetLogy()
    mg.Draw("A")
    mg.GetXaxis().SetTitle("m_{G} [TeV]")
    mg.GetYaxis().SetTitle("#sigma #times Br(G #rightarrow #gamma#gamma)")
    mg.GetYaxis().SetRangeUser(1e-3, 1.1)
    mg.Draw("A")
        
    k_factor = "1.75*"
    theory_xs = R.TF1("powerlaw", k_factor + "[0]*x^([1]+[2]*log(x))*expo(3)*([5]/0.03)^2", 100, 3000)
    theory_xs.SetLineColor(R.kRed)
    theory_xs.SetLineWidth(1)
        
    base_km = 0.1
    # Parameterization obtained from fit
    theory_xs.SetParameters(2.49412e-04, -3.19858e+01, 2.31110e+00, 1.20494e+02, -5.99267e-03, base_km)
    
    logx = 1
    if logx:
        c.SetLogx()
        mg.GetXaxis().SetRangeUser(370, 3090)
        mg.GetXaxis().SetMoreLogLabels()
        
    c.Modified(); c.Update()
    thfuncs = make_theory_functions(theory_xs)
    keep = []
    label_ypos = 3e-2
    for f in thfuncs:
        constrain_function(f)
        f.SetNpx(300)
        f.Draw("same")
        label_xpos = bisect(lambda x: f.Eval(x) - label_ypos, 0, 1e4)
        label = R.TLatex(label_xpos, label_ypos, "k/#bar{M}_{pl} = %s" % (f.GetParameter(5)))
        label.SetTextFont(label.GetTextFont()+1)
        label.SetTextSize(16)
        
        #from math import asin, pi, log, log10
        #angle = asin((log10(label_ypos) - log10(f.Eval(label_xpos + 100))) / 100)
        #angle *= 180 / pi
        
        angle = get_tf1_text_pad_angle(c, f, label, label_xpos)
        label.SetTextAngle(angle)
        label.SetTextAlign(12)
        
        label.Draw()
        keep.append(label)
        #R.TLatex(
        
    l = R.TLegend(0.572, 0.7, 0.81, 0.89)
    l.SetLineColor(R.kWhite); l.SetFillColor(R.kWhite)
    l.SetTextFont(l.GetTextFont()+1)
    l.SetTextSize(18)
    l.AddEntry(g_excl_xs, "Observed 95% C.I.", "lp")
    l.AddEntry(g_exptd, "Expected Median", "l")
    l.AddEntry(g_1sig, "Expected #pm 1#sigma", "f")
    l.AddEntry(g_2sig, "Expected #pm 2#sigma", "f")
    l.AddEntry(theory_xs, "Theoretical prediction", "l")
    l.Draw()
        
    text = R.TLatex(0.744, 0.64, "#sqrt{s} = 7 TeV")
    text.SetTextFont(text.GetTextFont()+1)
    text.SetTextSize(18)
    text.SetNDC()
    text.Draw()
    text2 = R.TLatex(0.702, 0.57, "#int L dt = 4.91 fb^{-1}")
    text2.SetTextFont(text2.GetTextFont()+1)
    text2.SetTextSize(18)
    text2.SetNDC()
    text2.Draw()
    
    keepalive = [l, text, text2]
    
    c.Update()
    import IPython; IPython.embed()
    #wait()
    return
    
    c1 = R.TCanvas()
    
    mg_plane = R.TMultiGraph("exclplane", "")
    
    
    # Expectation
    g_2sig_plane = graph_to_theory_plane(theory_xs, g_2sig)
    mg_plane.Add(g_2sig_plane, "e3")
    g_1sig_plane = graph_to_theory_plane(theory_xs, g_1sig)
    mg_plane.Add(g_1sig_plane, "e3")
    g_exptd_plane = graph_to_theory_plane(theory_xs, g_exptd)
    g_exptd_plane.SetMarkerStyle(7)
    mg_plane.Add(g_exptd_plane, "CP")
    
    # Data
    g_excl_plane = graph_to_theory_plane(theory_xs, g_excl_xs)
    mg_plane.Add(g_excl_plane, "CP")
    
    mg_plane.Draw("A")
    mg_plane.GetXaxis().SetRangeUser(350, 2500)
    mg_plane.GetYaxis().SetRangeUser(1e-3, 0.2)
    mg_plane.GetXaxis().SetTitle("m_{G} [TeV]")
    mg_plane.GetYaxis().SetTitle("k/#bar{M}_{pl}")
    mg_plane.Draw("A")
    
    l = R.TLegend(0.11, 0.7, 0.48, 0.89)
    l.SetLineColor(R.kWhite); l.SetFillColor(R.kWhite)
    l.SetTextFont(l.GetTextFont()+1)
    l.SetTextSize(18)
    l.AddEntry(g_excl_plane, "Data Observed", "lp")
    l.AddEntry(g_exptd_plane, "Expected", "l")
    l.AddEntry(g_1sig_plane, "Expected #pm 1#sigma", "f")
    l.AddEntry(g_2sig_plane, "Expected #pm 2#sigma", "f")
    l.Draw()
    
    text = R.TLatex(0.11, 0.65, "#sqrt{s} = 7 TeV")
    text.SetTextFont(text.GetTextFont()+1)
    text.SetTextSize(18)
    text.SetNDC()
    text.Draw()
    text2 = R.TLatex(0.11, 0.55, "#int L dt = 4.91 fb^{-1}")
    text2.SetTextFont(text2.GetTextFont()+1)
    text2.SetTextSize(18)
    text2.SetNDC()
    text2.Draw()
    
    mg_plane.GetYaxis().SetTitleOffset(mg_plane.GetYaxis().GetTitleOffset()*1.1)
    
    c.Modified()
    c.Update()
    c1.Modified()
    c1.Update()
    
    fs = []
            
    """
    xs, ys = map(mkarray, zip(*data))
    g_excl_plane = g = R.TGraph(len(xs), xs, ys)
    g.SetMarkerStyle(7)
    g.GetXaxis().SetRangeUser(600, 2400)
    
    """
    
    #g.Draw("ACP")
    #for f in fs: f.Draw("same")
    #g_excl_xs.Draw("CP")
    
    
    #wait()
    import IPython; IPython.embed()
    
    

if __name__ == "__main__":
    main()
    
