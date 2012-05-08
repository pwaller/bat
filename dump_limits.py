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
    
    mg.Draw("A")
    
    wait()
    
    

if __name__ == "__main__":
    main()
    
