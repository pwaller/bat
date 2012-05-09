#! /usr/bin/env python

import numpy as np

from array import array
import re


from rootwait import wait

import ROOT as R

def get_mass(i):
    lo, hi, n = 400, 3000, 200
    delta = (hi-lo)/n
    return lo + delta*i

def get_quantiles(h, quantiles=None):
    quantiles = quantiles or [0.05, 0.32, 0.5, 0.68, 0.95]
    input = array('d', quantiles)
    result = array('d', [0]*len(quantiles))
    h.GetQuantiles(len(input), result, input)
    return list(result)

def get_limit(filename):
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
    
    
    """
    t.Draw(
        "median_marginalized_1:"
        "5quantile_marginalized_1:"
        "16quantile_marginalized_1:"
        "84quantile_marginalized_1:"
        "95quantile_marginalized_1",
        "", "para goff")
        
    def get_mean(n):
        v = t.GetVal(n)
        values = [v[i] for i in xrange(t.GetEntries())]
        return sum(values) / len(values)
    return get_mass(index), tuple(get_mean(i) for i in xrange(5))
    """

def mkarray(data): return array('d', data)

def make_uncertaintyband(data):
    
    masses, values = zip(*data)
    masses = mkarray(masses)
    sig2d, sig1d, med, sig1u, sig2u = map(mkarray, zip(*values))
    zeros = mkarray([0]*len(masses))
    
    sig2d = mkarray(+ np.array(med) - np.array(sig2d))
    sig2u = mkarray(- np.array(med) + np.array(sig2u))
    sig1d = mkarray(+ np.array(med) - np.array(sig1d))
    sig1u = mkarray(- np.array(med) + np.array(sig1u))
    
    g2sig = R.TGraphAsymmErrors(len(masses), masses, med, zeros, zeros, sig2d, sig2u)
    g1sig = R.TGraphAsymmErrors(len(masses), masses, med, zeros, zeros, sig1d, sig1u)
    exptd = R.TGraph(len(masses), masses, med)
    
    g2sig.SetFillColor(R.kYellow)
    g1sig.SetFillColor(R.kGreen)
    exptd.SetLineStyle(2)
    
    return g2sig, g1sig, exptd
    
import ROOT as R
import numpy as NP

from array import array
from math import sqrt
from rootwait import wait
from sys import argv

def mkarray(x): return array('d', x)

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
    #g.Draw("ACP")
    
    k_factor = "1.75*"
    fn = R.TF1("powerlaw", k_factor + "[0]*x^([1]+[2]*log(x))*expo(3)*([5]/0.03)^2", 100, 3000)
    fn.SetLineColor(R.kRed)
    
    base_km = 0.1
    fn.SetParameters(2.49412e-04, -3.19858e+01, 2.31110e+00, 1.20494e+02, -5.99267e-03, base_km)
    
    
    fs = []
    for i in [0.01, 0.03, 0.05, 0.1, 0.15, 0.2]:
        fnew = fn.Clone()
        fs.append(fnew)
        fnew.SetParameter(5, i)
        #fnew.Draw("same")
    
    c1 = R.TCanvas()
    
    
    spline = R.TSpline5("spline", g)
    
    def compute_km_theoryintersect_at_mass(kmfunc, mass):
        measured_point = kmfunc(mass)
        theory_point = fn.Eval(mass)
        km = base_km * sqrt(measured_point / theory_point)
        return km
    
    data = [(mass, compute_km_theoryintersect_at_mass(spline.Eval, mass))
            for mass in NP.linspace(400, 3000, 200)]
        
    xs, ys = map(mkarray, zip(*data))
    g_excl_plane = g = R.TGraph(len(xs), xs, ys)
    g.SetMarkerStyle(7)
    g.GetXaxis().SetRangeUser(600, 2400)
    g.GetYaxis().SetRangeUser(1e-3, 0.1)
    #g.Draw("ACP")
    return g_excl_xs, g, fs
    

def main():
    from sys import argv
    data = [x for x in [get_limit(filename) for filename in argv[1:]] if x]
    data.sort()
    from pprint import pprint
    pprint(data, width=220)
        
    g2sig, g1sig, exptd = make_uncertaintyband(data)
    
    #exptd.Draw("ACP")
    #wait()
    
    mg = R.TMultiGraph("excl", "excl")
    mg.Add(g2sig, "e3")
    mg.Add(g1sig, "e3")
    
    mg.Add(exptd, "CP")
    
    
    g_excl_xs, g, fs = plot_data_limit("full_log_no_ensemble")
    
    mg.Add(g_excl_xs, "CP")
    
    mg.Draw("A")
    
    for f in fs:
        f.Draw("same")
    #g_excl_xs.Draw("CP")
    
    
    wait()
    
    

if __name__ == "__main__":
    main()
    
