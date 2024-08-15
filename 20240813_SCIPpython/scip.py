# try and copy fhprs code here
# start with the init section - within the class
# then do the if name = main bit - global scope
# should let us input the filename in the commandline at least
# python3 fh.py file.vcf`

import sys 
from samplesheet_parsing import *
import pandas as pd

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
        ss = ss_to_df(self.sample_sheet) # if want to print under if name == main, must be self.ss
        
        #print(ss)

        ss_len = len(ss)

        # for each row in the df, aka, for each patient:
        for row_index in range(0,ss_len,1):

            # assign patient info to variables
            s_name,mat_gt,pat_gt,hbb_file,sced_file = row_to_vars(row_index,ss)
            #print(s_name)
            #print(sced_file)

if __name__ == "__main__":
    sample_sheet = sys.argv[1]
    #variable1 = SCIP(sample_sheet)
    
    #print("oop doop")
    #print(variable1.ss) # only works if ss is declared as self.ss = blah blah
    #print("oop doop concluded")