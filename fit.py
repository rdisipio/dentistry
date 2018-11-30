#!/usr/bin/env python

import os, sys
from ROOT import *

gROOT.LoadMacro("AtlasStyle.C")
gROOT.LoadMacro( "AtlasUtils.C" )
SetAtlasStyle()

gStyle.SetOptStat(0)


# sieve mesh size (a given)
x = [ 15.0, 5.600, 4.000, 2.800, 2.000, 0.850, 0.425, 0.250 ]

infilename = "data.root"
if len(sys.argv) > 1: infilename = sys.argv[1]

ientry = 0
if len(sys.argv) > 2: pentry = int(sys.argv[2])

f = TFile.Open( infilename )

t = f.Get( "data" )
nentries = t.GetEntries()

c = TCanvas( "C", "C", 800, 800 )
c.Draw()

h0 = TH1F( "frame", ";x [mm];Q_{w}^{-} [%]", 20, 0., 20. )
h0.SetMaximum( 130. )
h0.SetMinimum( 0.01 )
h0.SetLineWidth(0)

for ientry in range(nentries):
  if not ientry == (pentry-1): continue

  t.GetEntry( ientry )

  g = TGraphErrors()
  g.SetLineWidth(2)
  
  g.SetName( "g_%i" % ientry )
  
  y = [ t.Sieve0, t.Sieve1, t.Sieve2, t.Sieve3, t.Sieve4, t.Sieve5, t.Sieve6, t.Sieve7] 
  
  g.SetPoint( 0, x[0], y[0] )
  g.SetPoint( 1, x[1], y[1] )
  g.SetPoint( 2, x[2], y[2] )
  g.SetPoint( 4, x[3], y[3] )
  g.SetPoint( 5, x[4], y[5] )
  g.SetPoint( 6, x[5], y[6] )

  # set uncertainties
  dx = 0.
  dy = 0.3
  g.SetPointError( 0, dx, dy*y[0] )
  g.SetPointError( 1, dx, dy*y[1] )
  g.SetPointError( 2, dx, dy*y[2] )
  g.SetPointError( 3, dx, dy*y[3] )
  g.SetPointError( 4, dx, dy*y[4] )
  g.SetPointError( 5, dx, dy*y[5] ) 
  g.SetPointError( 6, dx, dy*y[6] )

  h0.Draw()
  g.Draw( "PE same" )
  
  # Fit 
  # Rosin-Rammler
  # 0.693 = ln(2)
  fit_func = TF1( "fit_func_%i"%ientry, "100. * ( 1 -  TMath::Exp( -0.693*TMath::Power( (x/[0]), [1] ) )  )" )

  fit_func.SetParameter( 0, 0.5 )
  fit_func.SetParameter( 1, 2.0 )
  fit_func.SetLineColor(kBlue)
  fit_func.SetLineStyle(kDashed)

  fit_res = g.Fit( "fit_func_%i"%ientry, "R S  ", "", 0., 20. )
  chi2 = fit_res.Chi2()
  ndf  = fit_res.Ndf()
  
#  if ndf == 0: continue 
  
  chi2xdf = chi2 / ndf
  pvalue  = fit_res.Prob()
  x50  = fit_res.Value(0)
  beta = fit_res.Value(1)
  print "Result: chi2=%.2f : x50=%.1f beta=%.1f" % ( chi2, x50, beta )
  
  txt = TLatex()
  txt.SetNDC()
  txt.SetTextFont(42)
  txt.SetTextSize(0.03)
  txt.DrawLatex( 0.15, 0.90, "Patient id. %i" % pentry )
  txt.DrawLatex( 0.15, 0.85, "#Delta x = %.2f, #Delta y = %.2f" % (dx, dy) )
  txt.DrawLatex( 0.15, 0.80, "#chi2/NDF=%.2f / %i = %.2f, p-value=%.2f" % (chi2,ndf, chi2xdf, pvalue ) )
  txt.DrawLatex( 0.15, 0.75, "x_{50}=%.1f" % x50 )
  txt.DrawLatex( 0.15, 0.70, "#beta=%.1f" % beta )

  c.SaveAs( "img/fit_%i.png" % pentry )
  
  if not fit_res.Status() == 0:
     print "ERROR: fit did not converge"
     exit(1)      
exit(0)
