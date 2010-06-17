#!/usr/bin/env python
from __future__ import division
from snprankpy_class import DataProperties
import sys, os, string, math, csv

try:
    infilename = sys.argv[1]
    gamma = sys.argv[2]
except:
    print "Usage:", sys.argv[0], "in_all_data", "gamma_value"
    print "Example:", sys.argv[0], "example_data.tab", ".85"
    sys.exit(1)
    
# Create DataProperties object from data file
full_data = DataProperties(infilename)

full_data.pagerank_powermethod(float(gamma))