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
        self.SNPs = self.data[0]
        
    def powermethod(self, p):
        #Create a array with float values
        G = array(self.data[1:],dtype=float32)
        
        #Create a vector of diagonal values
        Gdiag = G.diagonal()
        
        #Get sum of the diagonal
        Gtrace = G.trace()
        
        #Save dimension of matrix
        [n,n] = G.shape
        
        #Get 1xn row vector of column sums
        colsum = G.sum(axis=0)
        
        #Get nx1 column vector of row sums
        rowsum = (G.sum(axis=1)).reshape(n,1)
        
        #Get indices of c vector that are not zero
        colsum_nzidx = colsum.nonzero()
        
        #Create matrix of zeros of size nxn
        D = zeros((n,n))
        
        #Create a 1xn row vector of ones
        T_nz = ones(n)
        
        #Create a diagional matrix and 
        #a vector of (1-gamma) of size n for the formula
        for i in colsum_nzidx[0]:
                D[i] = 1/colsum[i]
                T_nz[i] = T_nz[i] - p
        
        T = (p * G * D ) + (Gdiag * T_nz) / Gtrace
        T = T.transpose()
        
        #Reshape row vector into column vector of ones
        unit = (ones(n)).reshape(n,1)
        
        #Initial arbitrary vector
        r = unit/n;
        
        #Cutoff for matrix convergence
        threshold = 10**(-4)
        
        #Initialize Values
        converged = False
        
        while(not converged):
            r_old = r
            r = dot(T,r)
            lamb = sum(r)
            r = r/lamb 
            if all((abs(r-r_old))<threshold):
                converged = True
        print r
        print "lambda = ", lamb
       
                
        
