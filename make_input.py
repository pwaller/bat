#! /usr/bin/env python

import ROOT as R
from ROOT import TFile

def main():
    fdata = TFile("current/data.root")
    firred = TFile("current/sm_diphoton.root")
    fred = TFile("current/Bkg_gg_Template-Quentin.root")
    fgrav = TFile("current/graviton.root")

    foutput_bkg = TFile("current/Bkg_gg_Template.root", "recreate")

    hh_data = fdata.sel_reco_mgg_log.Clone("hh_data")
    hh_data.Write()

    #bkg_reducible_gg = fred.bkg_reducible_gg.Clone("bkg_reducible_gg")
    bkg_reducible_gg = fred.hh_red_400.Clone("bkg_reducible_gg") #hh_red_400")
    bkg_reducible_gg.Write()

    #bkg_irreducible_gg = firred.sel_reco_mgg_log.Clone("bkg_irreducible_gg")
    bkg_irreducible_gg = firred.sel_reco_mgg_log.Clone("bkg_irreducible_gg") #hh_irr_400")
    bkg_irreducible_gg.Write()

    bkg_total_gg = bkg_reducible_gg.Clone("bkg_total_gg")
    bkg_total_gg.Add(bkg_irreducible_gg)
    bkg_total_gg.Write()
    
    for x in "bkg_total_syst_gg bkg_reducible_syst_gg bkg_irreducible_syst_gg bkg_purity_syst_gg".split():
        fred.Get(x).Write()

    foutput_tmpl = TFile("current/Ggg_Templates.root", "recreate")
    
    keys = [(float(k.GetName().split("_")[-1]), k.ReadObj())
            for k in fgrav.Get("limit").GetListOfKeys()]
            
    keys.sort()
    
    arr = R.TObjArray()
    arr.SetName("template")
    for i, (v, k) in enumerate(keys):
        k.SetName("{0}".format(i))
        arr.Add(k)
    
    #arr.Write()
    foutput_tmpl.WriteTObject(arr)
    

if __name__ == "__main__":
    main()
