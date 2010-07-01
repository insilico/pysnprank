#!/usr/bin/env python
from __future__ import division
from numpy import *
import sys, csv, optparse

def normalize(xs):
	return xs/sum(xs)

class SNPrank(object):
	"""Provides functions for running SNPrank algorithm on a GAIN matrix"""
    def __init__(self, infilename):
        #Open file to read
        reader = csv.reader(infilename, delimiter="\t")
        
        #Store the headers
        self.SNPs = reader.next()

        #Store all data
        self.data = [row for row in reader]
        
    def calculate_snprank(self, gamma):
		"""Runs the SNPrank algorithm on the input data, using gamma as the damping factor.
		   Returns the SNPrank scores and diagonal (main effect) of original GAIN matrix."""
        #Create a array with float values
        G = array(self.data,dtype=float64)
        
        #Create a vector of diagonal values
        Gdiag = G.diagonal()
        
        #Save dimension of matrix
        [n,n] = G.shape
        
        #Get 1xn row vector of column sums
        colsum = G.sum(axis=0)
        
        #Get indices of c vector that are not zero
        colsum_nzidx = colsum.nonzero()[0]
        
        #Create matrix of zeros of size nxn
        D = zeros((n,n))
        
        #Create a 1xn row vector of ones
        T_nz = ones(n)
        
        #Create a diagional matrix and 
        #a vector of (1-gamma) of size n for the formula
        for i in colsum_nzidx:
            D[i][i] = 1/colsum[i]
            T_nz[i] -= gamma

        T = (gamma * dot(G,D) ) + (Gdiag.reshape(n,1) * T_nz) / G.trace()
        
        #Initial arbitrary vector
        r = (ones(n)).reshape(n,1)/n;
        
        #Cutoff for matrix convergence
        threshold = 10**(-4)
        
        #Initiate convergence and SNP ranking
        while True:
            r_old, r = r, normalize(dot(T,r))
            if all((abs(r-r_old))<threshold):
                break
                
        #Return SNPrank and diagonal of GAIN matrix
        return r.reshape(1,n)[0],Gdiag
    
    def print_to_file(self, SNPs, snp_rank, ig, output):
        """Output table of SNP names, SNPranks, information gain, and degree to a file."""
        column_header = ('SNPs','SNPrank','IG')
        #Create list of tuples of SNPs, SNPrank and information gain
        temp = zip(SNPs, snp_rank, ig)
        #Sort the list of tuples by snp_rank, i.e. 2nd value in each tuple
        sort_temp = sorted(temp, key=lambda item:item[1], reverse =True)
        #write column headings to file
        output.write("\t".join(head for head in column_header))
        output.write("\n")
        #Write values to file
        for x in range(len(sort_temp)):
            output.write("\t".join(str(value) for value in sort_temp[x]))
            output.write("\n")
        output.close()

#Define main function       
def main():
    #Create option parser
    parser = optparse.OptionParser("usage: %prog -f INFILE -o OUTFILE -g GAMMA")
    #Add options to parser; use defaults if none specified
    parser.add_option("-f", "--file", dest="infile", help="read data from INFILE")
    parser.add_option("-o", "--output", dest="outfile", help="output to OUTFILE")
    parser.add_option("-g", "--gamma", dest="gamma", help="gamma value (default: .85)", default=.85)
    
    (options, args) = parser.parse_args()
    
    #Check to see if filenames are mentioned, open to r/w if not
    infile  = open(options.infile,  'r') if options.infile  else sys.stdin
    outfile = open(options.outfile, 'w') if options.outfile else sys.stdout

    #Create data object from class
    full_data = SNPrank(infile)
    
    #Get SNPrank and information gain from calculate_snprank()
    snprank, IG = full_data.calculate_snprank(float(options.gamma))
    
    #Print to file
    full_data.print_to_file(full_data.SNPs,snprank, IG, outfile)
    
          
if __name__ == "__main__":
    main()
