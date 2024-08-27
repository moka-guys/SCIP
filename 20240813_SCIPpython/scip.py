# try and copy fhprs code here
# start with the init section - within the class
# then do the if name = main bit - global scope
# should let us input the filename in the commandline at least
# python3 fh.py file.vcf`

import sys 
from samplesheet_parsing import *
from mat_pat import *
import pandas as pd
import math

class SCIP(object):
    """
    Description of class.
    list and description of class attributes. 
    """

    # can add dictionary here - see graemes SCORES dictionary. 
    # dictionary for which allele is which position, which alt and ref?
    
    def __init__(self, sample_sheet):
        self.sample_sheet = sample_sheet
        self.main()

    def main(self):
        """
        Doc string for function
        """
        # convert samplesheet csv into dataframe
        ss = ss_to_df(self.sample_sheet) # if want to print under if name == main, must be self.ss
     
        ss_len = len(ss)

        # for each row in the df, aka, for each patient:
        for row_index in range(0,ss_len,1):

            # assign patient info to variables
            sample,mat_gt,pat_gt,hbb_file,sced_file = row_to_vars(row_index,ss)
            print(sample)
            #print(mat_gt)
            #print(pat_gt)
            #print(sced_file)

            alleles = alleles_of_interest(mat_gt,pat_gt)
            

if __name__ == "__main__":
    sample_sheet = sys.argv[1]
    variable1 = SCIP(sample_sheet) # need to call SCIP class or code will not run
    
    #print("oop doop")
    #print(variable1.ss) # only works if ss is declared as self.ss = blah blah
    #print("oop doop concluded")