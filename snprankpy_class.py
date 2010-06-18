#!/usr/bin/env python
from __future__ import division
from numpy import *
import sys, os, string, math, csv

class DataProperties(object):
    
    def __init__(self, infilename):
        tsv = open(infilename,'r')
        reader = csv.reader(tsv, delimiter="\t")
        
        #Store all data
        self.data = [row for row in reader]
        
        #Store the headers
        self.header = self.data[0]
        
    def pagerank_powermethod(self, p):
        #Create a matrix with float values
        matrix = array(self.data[1:],dtype=float32)
        
        #Create a vector of diagonal values
        Gdiag = matrix.diagonal()
        
        #Get sum of the diagonal
        Gtrace = matrix.trace()
        
        #Save dimension of matrix
        [n,n] = matrix.shape
        
        #Get 1xn row vector of column sums
        c = matrix.sum(axis=0)
        
        #Get nx1 column vector of row sums
        r = (matrix.sum(axis=1)).reshape(n,1)
        
        #Get indices of c vector that are not zero
        k = c.nonzero()
        
        #Create matrix of zeros
        D = ones((n,n))
        
        z_t = zeros(n)
        #Create a diagional matrix for formula
        for i in k[0]:
                D[i] = Gtrace/c[i]
                z_t[i] = 1-p
        
        T = (p * matrix * D + Gdiag * z_t) / Gtrace
        
        unit = (ones(n)).reshape(n,1)
        
        x = unit/n;
        
        x = x*((T.sum(axis=0)).reshape(n,1))
        for i in [1,2,3,4,5]:
            lamb = sum(x)
            x = x/lamb
        print x
        print "lambda = ", lamb
            
        
       
                
        