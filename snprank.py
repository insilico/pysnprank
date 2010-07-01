#!/usr/bin/env python
from __future__ import division
from numpy import *
import sys, csv, optparse

class SNPrank(object):
    
    def __init__(self, infilename):
        #Open file to read
        reader = csv.reader(infilename, delimiter="\t")
        
        #Store all data
        self.data = [row for row in reader]
        
        #Store the headers
        self.SNPs = self.data[0]
        
    def calculate_snprank(self, gamma):
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
            T_nz[i] = T_nz[i] - gamma
        
        T = (gamma * G * D ) + (Gdiag * T_nz) / Gtrace
        T = T.transpose()
        
        #Reshape row vector into column vector of ones
        unit = (ones(n)).reshape(n,1)
        
        #Initial arbitrary vector
        r = unit/n;
        
        #Cutoff for matrix convergence
        threshold = 10**(-4)
        
        #Initialize Values
        converged = False
        
        #Initiate convergence and SNP ranking
        while(not converged):
            r_old = r
            r = dot(T,r)
            lamb = sum(r)
            r = r/lamb 
            if all((abs(r-r_old))<threshold):
                converged = True
                
        #Return SNP_rank and diagonal of GAIN matrix
        return r.reshape(1,n)[0],Gdiag
    
    #Output to file function
    def print_to_file(self, SNPs, snp_rank, ig, output):
        column_header = ('SNPs','SNP_RANK','IG')
        #Create list of tuples of SNPs, snp_rank and InformationGain
        temp = zip(SNPs, snp_rank, ig)
        #Sort the list of tuples by snp_rank, i.e. 2nd value in each tuple
        sort_temp = sorted(temp, key=lambda item:item[1], reverse =True)
        #write column headings to file
        output.write("\t".join(head for head in column_header))
        #Write values to file
        for x in range(len(sort_temp)):
            output.write("\n")
            output.write("\t".join(str(value) for value in sort_temp[x]))
        output.close()

#Define main function       
def main():
    #Create option parser
    parser = optparse.OptionParser("usage: %prog -f FILENAME -o OUTPUT -g GAMMA VALUE")
    #Add options to parser; use defaults if none specified
    parser.add_option("-f", "--file", dest="infile", help="read data from FILENAME", default=sys.stdin, action ='store')
    parser.add_option("-o", "--output", dest="outfile", help="output to FILENAME", default=sys.stdout)
    parser.add_option("-g", "--gamma", dest="gamma", help="gamma value", default=.85)
    
    (options, args) = parser.parse_args()
    
    #Check to see if filenames are mentioned,
    #Open to r/w if not
    if type(options.infile)==type('string'):
        options.infile = open(options.infile, 'r')
    if type(options.outfile)==type('string'):
        options.outfile = open(options.outfile, 'w')
        
    #Create data object from class
    full_data = SNPrank(options.infile)
    
    #Get SNP_rank and InformationGain from powermethod()
    snprank, IG = full_data.calculate_snprank(float(options.gamma))
    
    #Print to file
    full_data.print_to_file(full_data.SNPs,snprank, IG, options.outfile)
    
          
if __name__ == "__main__":
    main()
