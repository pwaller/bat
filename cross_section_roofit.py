

import ROOT as R
from ROOT import RooRealVar, RooDataSet, RooArgSet, RooFit, RooGenericPdf, RooArgList, RooLinkedList

massbounds = 600, 4000
xsbounds = 1e-14, 1e10

mass, xs = RooRealVar("mass", "mass", *massbounds), RooRealVar("xs", "xs", *xsbounds)

ds = RooDataSet("data", "data", RooArgSet(mass, xs), RooFit.StoreError(RooArgSet(xs)))

data = [
 #(100, (238656.80456246276, 5021.139940690172)),
 #(200, (4256.158341494928, 212.77361569648355)),
 #(300, (753.6236834753611, 37.5498347216208)),
 #(400, (196.16227118766918, 9.800314801274078)),
 (600, (20.491217905484927, 1.0231064675754569)),
 (800, (4.256970608754619, 0.21145030944619475)),
 (900, (2.1511503582230613, 0.10749890705342333)),
 (1000, (1.2694182781976522, 0.06340890314406443)),
 (1100, (0.643681492794392, 0.03217545331472744)),
 (1500, (0.08287990689532962, 0.004117232695182858)),
 (1750, (0.026150402908862617, 0.0013037026701146663)),
 (2000, (0.008063601507094236, 0.00040306052931164015)),
 (2250, (0.003056637009041868, 0.00015253920943590022)),
 (4000, (1.5147902166034461e-06, 7.553510014249738e-08)),
 #(5000, (6.7752618437184105e-09, 3.3855147194005573e-10)),
 #(6000, (1.673646004156878e-10, 8.333326606880855e-12)),
 #(7000, (4.7172685387069244e-11, 2.356523443723669e-12))
]

for massval, (xsval, xserrval) in data:
    mass.setVal(massval); xs.setVal(xsval); xs.setError(xserrval)
    ds.add(RooArgSet(mass, xs))

#c0pr = RooRealVar("c0pr", "c0pr", 1, 1e-5, 1e5)
c0 = RooRealVar("c0", "c0",  3.13e-06,  1e-8, 1e5)
c1 = RooRealVar("c1", "c1", -1.59e+01, -1e+2, 1e2)
c2 = RooRealVar("c2", "c2",  1.08e+00,  1e-5, 1e1)
c3 = RooRealVar("c3", "c3",  7.66e+01,  0, 1e2)
c4 = RooRealVar("c4", "c4", -4.97e-03, -1e-5, 1)
#c5 = RooRealVar("c5", "c5",  0.01); c5.setConstant()

partial_pdf = RooGenericPdf("xspdf", "xspdf", "(mass^(c1+c2*log(mass))*exp(c3+c4*mass))",
    RooArgList(mass, c1, c2, c3, c4))
    
pdf = R.RooExtendPdf("fullpdf", "fullpdf", partial_pdf, c0)

fr = mass.frame()

ds.plotOnXY(fr, RooFit.YVar(xs))
#pdf.plotOn(fr)

#fr.Draw()
#from IPython import embed; embed()#wait()

#def python_to_root(name, variable):
#    line = '{0}& {1} = *reinterpret_cast<{0}*>((void*)TPython::Eval("{1}"));'
#    R.gInterpreter.ProcessLine(line.format(variable.ClassName(), name))
    
#python_to_root("ds", ds)
#python_to_root("pdf", pdf)
#python_to_root("xs", xs)
#yvar = RooFit.YVar(xs)
#python_to_root("yvar", yvar)

#python_to_root("ll", ll)

#R.gInterpreter.ProcessLine('RooFitResult* result = pdf->chi2FitTo(ds, yvar);')
#R.gInterpreter.ProcessLine('TPython::Bind(result, "result");')

#partial_pdf.fitTo(ds)

ll = RooLinkedList()
ll.Add(RooFit.YVar(xs))
ll.Add(RooFit.Save())
pdfrar = pdf.IsA().DynamicCast(R.RooAbsReal.Class(), pdf)
result = pdfrar.chi2FitTo(ds, RooFit.YVar(xs), RooFit.Save()) #, RooFit.Integrate(R.kTRUE)) #, ll)

#pdf.fitTo(ds)
#partial_pdf.IsA().DynamicCast(R.RooAbsPdf.Class(), partial_pdf).syncNormalization(ds)

pdf.plotOn(fr, RooFit.Normalization(pdf.expectedEvents(pdf.getVariables())))

#partial_pdf.fitTo(ds)
#partial_pdf.plotOn(fr) # "L") #RooFit.YVar(xs))

fr.GetYaxis().SetRangeUser(1e-14, 1e6)
fr.Draw()
R.gPad.SetLogy()
R.gPad.Update()

from IPython import embed; embed()

"""
RooFit v3.50 -- Developed by Wouter Verkerke and David Kirkby 
                Copyright (C) 2000-2011 NIKHEF, University of California & Stanford University
                All rights reserved, please read http://roofit.sourceforge.net/license.txt

[#1] INFO:NumericIntegration -- RooRealIntegral::init(xspdf_Int[mass]) using numeric integrator RooIntegrator1D to calculate Int(mass)
[#1] INFO:Eval -- RooRealVar::getBinning(mass) new range named 'bin' created with default bounds
[#1] INFO:Eval -- RooRealVar::getBinning(xs) new range named 'bin' created with default bounds
[#1] INFO:Minization -- RooMinuit::optimizeConst: activating const optimization
 **********
 **   13 **MIGRAD        2500           1
 **********
 FIRST CALL TO USER FUNCTION AT NEW START POINT, WITH IFLAG=4.
 START MIGRAD MINIMIZATION.  STRATEGY  1.  CONVERGENCE WHEN EDM .LT. 1.00e-03
 FCN=1315.05 FROM MIGRAD    STATUS=INITIATE      151 CALLS         152 TOTAL
                     EDM= unknown      STRATEGY= 1      NO ERROR MATRIX       
  EXT PARAMETER               CURRENT GUESS       STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  c0           2.62986e+03   4.99995e-01   0.00000e+00   6.36705e+02
   2  c1          -1.00000e+01   2.00000e+01   0.00000e+00  -2.37303e+04
   3  c2           1.00000e+00   4.99995e-01   0.00000e+00  -1.05618e+04
   4  c3           7.00000e+01   1.00000e+01   0.00000e+00   0.00000e+00
   5  c4          -1.00000e-02   constant   
 MIGRAD MINIMIZATION HAS CONVERGED.
 MIGRAD WILL VERIFY CONVERGENCE AND ERROR MATRIX.
 EIGENVALUES OF SECOND-DERIVATIVE MATRIX:
        -1.8968e-04  2.8611e-01  1.0000e+00  2.7141e+00
 MINUIT WARNING IN HESSE   
 ============== MATRIX FORCED POS-DEF BY ADDING 0.002904 TO DIAGONAL.
 FCN=182.937 FROM MIGRAD    STATUS=CONVERGED     610 CALLS         611 TOTAL
                     EDM=3.68246e-07    STRATEGY= 1      ERR MATRIX NOT POS-DEF
  EXT PARAMETER                APPROXIMATE        STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  c0           2.84613e+03   7.35639e+01   1.82972e-05   9.93760e-02
   2  c1          -7.30247e+01   2.56213e-01   1.82044e-06   2.55941e+00
   3  c2           5.58166e+00   1.79496e-02   1.76551e-06   2.57367e+00
   4  c3           1.93826e+01   8.85397e+01   5.00000e-01   2.84217e-14
   5  c4          -1.00000e-02   constant   
 EXTERNAL ERROR MATRIX.    NDIM=  25    NPAR=  4    ERR DEF=1
  5.412e+03 -2.148e+00  2.963e-02 -1.850e+04 
 -2.148e+00  6.565e-02 -4.567e-03 -7.665e+03 
  2.963e-02 -4.567e-03  3.222e-04  5.433e+02 
 -1.850e+04 -7.665e+03  5.433e+02  6.852e+15 
ERR MATRIX NOT POS-DEF
 PARAMETER  CORRELATION COEFFICIENTS  
       NO.  GLOBAL      1      2      3      4
        1  0.78164   1.000 -0.114  0.022 -0.000
        2  0.99732  -0.114  1.000 -0.993 -0.000
        3  0.99728   0.022 -0.993  1.000  0.000
        4  0.00037  -0.000 -0.000  0.000  1.000
 ERR MATRIX NOT POS-DEF
 **********
 **   18 **HESSE        2500
 **********
 FIRST CALL TO USER FUNCTION AT NEW START POINT, WITH IFLAG=4.
 COVARIANCE MATRIX CALCULATED SUCCESSFULLY
 FCN=182.937 FROM HESSE     STATUS=OK             23 CALLS         635 TOTAL
                     EDM=9.67132e-07    STRATEGY= 1      ERROR MATRIX ACCURATE 
  EXT PARAMETER                                INTERNAL      INTERNAL  
  NO.   NAME      VALUE            ERROR       STEP SIZE       VALUE   
   1  c0           2.84613e+03   7.67954e+01   3.65944e-06  -1.23177e+00
   2  c1          -7.30247e+01   1.08127e+00   3.64089e-07  -8.18684e-01
   3  c2           5.58166e+00   7.57448e-02   3.53102e-07   1.16596e-01
   4  c3           1.93826e+01   5.28613e+01   5.00000e-01  -6.59027e-01
   5  c4          -1.00000e-02   constant   
 EXTERNAL ERROR MATRIX.    NDIM=  25    NPAR=  4    ERR DEF=1
  5.898e+03 -2.410e+01  1.566e+00 -1.422e+06 
 -2.410e+01  1.169e+00 -8.188e-02  7.298e+04 
  1.566e+00 -8.188e-02  5.738e-03 -5.113e+03 
 -1.422e+06  7.298e+04 -5.113e+03  7.234e+14 
 PARAMETER  CORRELATION COEFFICIENTS  
       NO.  GLOBAL      1      2      3      4
        1  0.80124   1.000 -0.290  0.269 -0.001
        2  0.99985  -0.290  1.000 -1.000  0.003
        3  0.99985   0.269 -1.000  1.000 -0.003
        4  0.00251  -0.001  0.003 -0.003  1.000
[#1] INFO:NumericIntegration -- RooRealIntegral::init(xspdf_Int[mass]) using numeric integrator RooIntegrator1D to calculate Int(mass)

"""
