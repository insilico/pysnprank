#!/usr/bin/env python
from __future__ import division
from numpy import *
import sys, os, string, math, csv

def transpose(a):
    return zip(*a)

class DataProperties(object):
    
    def __init__(self, infilename):
        tsv = open(infilename,'r')
        reader = csv.reader(tsv, delimiter="\t")
        
        data = [row for row in reader]
        
        tsv.close()
        
        data_transpose = transpose(data)
        
        # status symbol
        self.status_key = string.strip(data[0][-1])

        self.num_instances  = len(data) - 1
        
        # Create attribute_name -> data dictionary
        self.attributes = {}
        for row in data_transpose:
            self.attributes[row[0]] = row[1:]
            # remember when constructing data set with sampled attributes
            # to include status_key at the end
            
            
    def pagerank_powermethod(self):
        matrix = array(self.attributes.values())