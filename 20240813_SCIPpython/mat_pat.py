# function which reads 'AS' 'AC' etc and separates into different alleles
# discards A's
# makes up SNP table of alleles we have to care about
# maternal gt necessary, paternal optional
# what if gts not fully known?

# csv with maternal gt / paternal gt / also columns for file names?

# one function to read in file and assign to variables too. 
# hey thats probably a class

# input is 2 variables - the mat_gt variable and the pat_gt variable. 
# variables should be strings or empty
# need to split strings into component parts. 
# and then ... do stuff. Save strings in list? or set variables to 
# true or false?

import math

def alleles_of_interest(mat_gt,pat_gt):

    mat_alleles = []
    pat_alleles = []

    # TODO check alleles provided are valid aka S C D or E

    # parsing maternal alleles
    if type(mat_gt) == str:
        mat_gt = str(mat_gt) # TODO make sure alleles in upper case
        mat_alleles = list(mat_gt)
        print("The maternal genotype provided is " + mat_gt)
    else:
        print("No maternal genotype provided")

    # parsing paternal alleles
    if type(pat_gt) == str:
        pat_gt = str(pat_gt)
        pat_alleles = list(pat_gt)
        print("The paternal genotype provided is " + pat_gt)
    else:
        print("No paternal genotype provided")
    
    # combine maternal and paternal alleles to produce list of alleles relevant
    # for this patient

    parental_alleles = mat_alleles + list(set(pat_alleles) - set(mat_alleles))

    wt_allele = "A"

    if wt_allele in parental_alleles: parental_alleles.remove(wt_allele)


    # if no alleles / alleles list is empty - error exception here TODO
  
    return(parental_alleles)
#   