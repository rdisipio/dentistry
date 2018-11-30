#!/usr/bin/env python

import sys
from ROOT import *

infilename = "data_percent.csv"
if len( sys.argv ) > 1: infilename = sys.argv[1]

header = "ID/I,Sieve0,Sieve1/f,Sieve2/f,Sieve3/f,Sieve4/f,Sieve5/f,Sieve6/f,Sieve7/f"

f = TFile.Open( "data.root", "RECREATE" )
t = TTree( "data", "data" )
t.ReadFile( infilename, header, ',' )

f.Write()
f.Close()
