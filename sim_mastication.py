#!/usr/bin/env python

import os, sys

import numpy as np
from ROOT import *
from array import array

rho_silicone_rubber = 0.0015 #g/mm3 # Cuttersil by Kulzer
# see data sheet:
# http://msds.kulzer.com/msds/MSDS646_-_CutterSil_Impression_Material_(CDN)_3.pdf

ntrials = 55
nparticles_avg = 50
tau = 3.00
r_max = 10.
w_max = 200.

# Change default parameters from command line input
if len(sys.argv) > 1: ntrials = int(sys.argv[1])
if len(sys.argv) > 2: nparticles_avg = int(sys.argv[2])
if len(sys.argv) > 3: tau = float(sys.argv[3])

def weight( r, rho=rho_silicone_rubber ):
    # assume spherical particles
    V = (4./3.) * 3.1415 * r*r*r
    w = V * rho
    return 1000.*w #mg

# Create random number generator
rng = TRandom3()

sieves_mesh = [ 1000, 5.60, 4.00, 2.80, 2.00, 0.85, 0.425, 0.250 ] #mm
n_sieves = len( sieves_mesh ) - 1

print "INFO: number of sieves:", n_sieves

all_particles = []

ofile = TFile.Open( "mastication.root", "RECREATE" )
_h = {}
_h['all_particle_radius'] = TH1F( "all_particle_radius", ";Radius [mm]", 5, 0., r_max )
_h['all_particle_weight'] = TH1F( "all_particle_weight", ";Weight [mg]", 10, 0., w_max )
_h['weight_per_sieve'] = TH2F( "weight_per_sieve", ";Sieve No.;Weight [mg]", n_sieves, 0.5, n_sieves+0.5, 10, 0., w_max )
_h['nparticles_per_sieve'] = TH1F( "nparticles_per_sieve", ";Sieve", n_sieves, 0.5, n_sieves+0.5 )
_h['Qw'] = TH1F("Qw", ";Sieve", n_sieves, 0.5, n_sieves+0.5 )

for h in _h.values(): h.Sumw2()

for itrial in range(ntrials):
    print "INFO: trial %i/%i" % ( (itrial+1), ntrials )
    
    nparticles = rng.Poisson(nparticles_avg)
    
    # Reset histograms
    for h in _h.values(): h.Reset()
   
    particle_radius = [ rng.Exp(tau) for i in range(nparticles) ]
    particle_weight = [ weight(r)    for r in particle_radius   ]
   
    for x in particle_radius: _h["all_particle_radius"].Fill(x)
    for x in particle_weight: _h["all_particle_weight"].Fill(x)
   
    weight_total = 0.
    weight_per_sieve = np.zeros(n_sieves)
    for i in range(nparticles):
        r = particle_radius[i]
        w = particle_weight[i]
        weight_total += w
       
        #print "INFO: particle r=%.3f mm, w=%.3f mg" % ( r, w )
        for isieve in range(n_sieves):
            r_hi = sieves_mesh[isieve]
            r_lo = sieves_mesh[isieve+1]
            #print isieve, r_hi, " > r > ", r_lo
           
            if (r<r_hi) and (r>r_lo): 
                weight_per_sieve[isieve] = w
                _h['weight_per_sieve'].Fill( isieve-1, w )
                _h['nparticles_per_sieve'].Fill(isieve-1) 
                #print "INFO: w=%.3f mg: falls into sieve %i" % ( w, (isieve) )
                break
    
    # Calculate fraction per sieve (Qw)
    for isieve in range(n_sieves):
        Qw = 100.* weight_per_sieve[isieve] / weight_total
        _h['Qw'].SetBinContent( isieve+1, Qw )
   
    # write out
    for h in _h.values(): 
        hname = "%s_trial_%i" % ( h.GetName(), (itrial+1) )
        h.Write( hname )

#ofile.Write()
ofile.Close()
