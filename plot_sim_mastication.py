#!/usr/bin/env python

import os, sys

from ROOT import *

gROOT.LoadMacro("AtlasStyle.C")
gROOT.LoadMacro( "AtlasUtils.C" )
SetAtlasStyle()

gStyle.SetOptStat(0)

def Normalize(h):
    area = h.Integral()
    h.Scale( 100./area )
    h.GetYaxis().SetTitle( "Fraction of particles [%]" )
    h.GetYaxis().SetTitleOffset(1.2)

infilename = "mastication_55_50.root"
if len(sys.argv) > 1: infilename = sys.argv[1]

ntrials        = int(infilename.split(".")[0].split("_")[1])
nparticles_avg = int(infilename.split(".")[0].split("_")[2])

f = TFile.Open( infilename )

c = TCanvas( "C", "C", 2500, 400 )
c.Draw()
c.Divide(5,1)

for itrial in range(1,ntrials+1):
    
    h_radius = f.Get( "all_particle_radius_trial_%i" % itrial )
    h_weight = f.Get( "all_particle_weight_trial_%i" % itrial )
    h_nparticles_per_sieve = f.Get("nparticles_per_sieve_trial_%i" % itrial )
    h_weight_per_sieve     = f.Get("weight_per_sieve_trial_%i" % itrial )
    h_Qw  = f.Get( "Qw_trial_%i" % itrial )
    
    c.cd(1) 
    Normalize(h_radius)
    h_radius.SetMaximum(100.)
    h_radius.Draw()
    
    c.cd(2)
    Normalize( h_weight )
    h_weight.SetMaximum(100.)
    h_weight.Draw()
    
    c.cd(3)
    Normalize(h_nparticles_per_sieve)
    h_nparticles_per_sieve.SetMaximum(100.)
    h_nparticles_per_sieve.Draw()
    
    c.cd(4)
    Normalize(h_weight_per_sieve)
    h_weight_per_sieve.GetYaxis().SetTitle( "Fraction of total weight [%]" )
    h_weight_per_sieve.SetMaximum(100.)
    h_weight_per_sieve.Draw()
    
    c.cd(5)
    frame = TH1F( "frame", "frame", 10, 0., 20 )
    frame.SetMaximum(130.)
    frame.SetMinimum(0.)
    frame.GetXaxis().SetTitle( "x [mm]" )
    frame.GetYaxis().SetTitle("Qw [%]")
    frame.Draw()
    h_Qw.Draw("p same")
    h_Qw.GetYaxis().SetTitleOffset(1.2)
    
    fit_func = TF1( "fit_func_%i"%itrial, "100. * ( 1 -  TMath::Exp( -0.693*TMath::Power( (x/[0]), [1] ) )  )" )
    fit_func.SetParameter( 0, 1.0 )
    fit_func.SetParameter( 1, 2.0 )
    fit_func.SetLineColor(kBlue)
    fit_func.SetLineStyle(kDashed)
    fit_res = h_Qw.Fit( "fit_func_%i"%itrial, "R S  ", "", 0., 20. )
    chi2 = fit_res.Chi2()
    ndf  = fit_res.Ndf()
    chi2xdf = chi2 / ndf
    pvalue  = fit_res.Prob()
    x50  = fit_res.Value(0)
    beta = fit_res.Value(1)
    
    txt = TLatex()
    txt.SetNDC()
    txt.SetTextFont(42)
    txt.SetTextSize(0.03)
    txt.DrawLatex( 0.15, 0.90, "Trial no. %i" % itrial )
    #txt.DrawLatex( 0.15, 0.85, "#Delta x = %.2f, #Delta y = %.2f" % (dx, dy) )
    txt.DrawLatex( 0.15, 0.85, "#chi2/NDF=%.2f / %i = %.2f, p-value=%.2f" % (chi2,ndf, chi2xdf, pvalue ) )
    txt.DrawLatex( 0.15, 0.80, "x_{50}=%.1f" % x50 )
    txt.DrawLatex( 0.15, 0.75, "#beta=%.1f" % beta )

  
  
    c.cd()
    c.SaveAs( "img/sim_mastication_%i_%i_trial_%i.png" % (ntrials, nparticles_avg, itrial) )