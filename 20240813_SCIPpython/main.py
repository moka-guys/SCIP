# maternal allele / alleles of interest input
# maternal allele (needed)
# paternal allele (optional)

# using mpileup files as input
# informative snps calcs - output txt file
# alleles of interest calcs - output total and alt count, possibly two.

# paralelisation -
# do we want to run this script once per sample
# or once per bundle of samples?
# if once per bundle, need to put in measures to stop whole thing failing
# if one sample is duff
# if paralell, reading in a csv of mpileup file names + maternal / paternal alleles. 
# if single, command line input? still have input file? 

# consider also, docker container. What will we import into?
from samplesheet_parsing import *
import pandas as pd

# 1: samplesheet_parsing.py

ss = ss_to_df()

# 1.5: loop over rows of samplesheet, pick out patient info
# make function that picks out patient info for a specified row of ss?

ss_len = len(ss)

# for each row in the df, aka, for each patient:
for row_index in range(0,ss_len,1):

    # assign patient info to variables
    s_name,mat_gt,pat_gt,hbb_file,sced_file = row_to_vars(row_index,ss)
    print(s_name)
    print(sced_file)

    # determine alleles of interest for the patient
    # mat_pat.py

    # calculate total_count and alt_count for the alleles of interest
    # convert perl script
    # return values - df?

    # calculate fetal fraction at HBB SNPs
    # convert perl script

    # predict which alleles are present in the fetus
    # feeding data in converted R script.
