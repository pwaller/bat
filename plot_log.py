#! /usr/bin/env python

import ROOT as R
import numpy as NP

from array import array
from math import sqrt
from rootwait import wait
from sys import argv

def mkarray(x): return array('d', x)

def main():

    input_file = argv[1]

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
    g.Draw("ACP")
    
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
        fnew.Draw("same")
    
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
    g.Draw("ACP")
    
    from IPython import embed; embed()
    wait()
    
if __name__ == "__main__":
    main()
