#!/usr/bin/env python
from __future__ import division
from snprankpy_class import DataProperties
import sys, os, string, math, csv

try:
    infilename = sys.argv[1]
except:
    print "Usage:", sys.argv[0], "in_all_data"
    print "Example:", sys.argv[0], "example_data.tab"
    sys.exit(1)
    
# Create DataProperties object from data file
full_data = DataProperties(infilename)
