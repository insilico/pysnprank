#!/usr/bin/env python
from __future__ import division
from numpy import *
import sys, os, string, math, csv

class DataProperties(object):
    
    def __init__(self, infilename):
        tsv = open(infilename,'r')
        reader = csv.reader(tsv, delimiter="\t")
        
        self.data = [row for row in reader]
        
    def pagerank_powermethod(self):
        matrix = array(self.data[1:])
        print matrix