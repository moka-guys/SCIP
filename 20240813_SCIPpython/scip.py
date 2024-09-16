# try and copy fhprs code here
# start with the init section - within the class
# then do the if name = main bit - global scope
# should let us input the filename in the commandline at least
# python3 fh.py file.vcf`

import sys 
from samplesheet_parsing import *
from mat_pat import *
from total_and_alt_count import *
from fetal_frac_calc import *
from fetal_gt_pred import *
import pandas as pd
import math
import os
import matplotlib.pyplot as plt

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
            
            # extract parental alleles
            alleles = alleles_of_interest(mat_gt,pat_gt)
            print(alleles)

            # extract total and alt counts of parental alleles to variables
            S_total, S_alt, C_total, C_alt, E_total, E_alt, D_total, D_alt = total_and_alt_vars(sced_file, alleles)
            total_counts = [S_total, C_total, E_total, D_total]
            alt_counts = [S_alt, C_alt, E_alt, D_alt]
            allele_labels = ["S allele", "C allele", "E allele", "D allele"]
            # TODO remove variables that are empty

            # determine informative snps and calculate fetal fractions, save to txt file
            # TODO softcode depth value
            fetal_frac_output_path = sample + "_fetal_frac_output.txt"
            fetal_frac(390,hbb_file,fetal_frac_output_path)

            # predict genotype using R script conversion

            for total, alt, label in zip(total_counts, alt_counts, allele_labels):
                # selecting only for the relevant alleles according to parental gt
                if total is not None:
                    prediction = gt_prediction(fetal_frac_output_path,total,alt)
                    print("The predicted fetal genotype for the " + label + " allele is: " + prediction)

            

if __name__ == "__main__":
    sample_sheet = sys.argv[1]
    variable1 = SCIP(sample_sheet) # need to call SCIP class or code will not run
    
    #print("oop doop")
    #print(variable1.ss) # only works if ss is declared as self.ss = blah blah
    #print("oop doop concluded")