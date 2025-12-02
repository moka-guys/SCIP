
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def mat_gt_prediction(fetal_frac,total_count,alt_count):


    # calculate the maternal genotype for this allele, based on total and alt count
    mutant_ratio = alt_count / total_count
    maternal_gt = ""
    if mutant_ratio >= 0.8:
        maternal_gt = "Homozygous mutant"
    elif 0.3<= mutant_ratio <= 0.7:
        maternal_gt = "Heterozygous"
    elif mutant_ratio <= 0.2:
        maternal_gt = "Homozygous wildtype"

    return maternal_gt

def mat_het_gt_prediction(fetal_frac,total_count,alt_count):
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    # read the fetal fraction file into a pandas df
    FL_SNPs = pd.read_csv(fetal_frac, sep = '\t')
    # remove the last five lines which do not contain informative SNP information
    FL_SNPs = FL_SNPs.iloc[:-5]
    
    # Rename the columns
    FL_SNPs.columns = ["Chromosome", "Start", "End", "Num.reads", "A", "A_fraction", 
                    "G", "G_fraction", "C", "C_fraction", "T", "T_fraction"]
    
    # assign number of informative snps to variable
    informative_snps=len(FL_SNPs)

    # Calculate the max value for A_fraction, G_fraction, C_fraction, T_fraction
    FL_SNPs['ColMax'] = FL_SNPs[["A_fraction", "G_fraction", "C_fraction", "T_fraction"]].max(axis=1)
    FL_SNPs['PaternalFraction'] = 100 - FL_SNPs['ColMax']
    FL_SNPs['FetalFraction'] = 2 * FL_SNPs['PaternalFraction']

    # Summary statistics
    mean_pat = FL_SNPs['PaternalFraction'].mean()
    median_pat = FL_SNPs['PaternalFraction'].median()
    IQR_Pat = FL_SNPs['PaternalFraction'].quantile(0.75) - FL_SNPs['PaternalFraction'].quantile(0.25)
    mean_Fet = FL_SNPs['FetalFraction'].mean()
    median_Fet = FL_SNPs['FetalFraction'].median()
    IQR_Fet = FL_SNPs['FetalFraction'].quantile(0.75) - FL_SNPs['FetalFraction'].quantile(0.25)

    # Displaying the table (in Jupyter, pandas automatically outputs it)
    # In a real HTML or PDF generation, use pandas to_html() or to_latex()
    #FL_SNPs

    # SPRT plot setup
    Q0 = 0.5
    Q1 = 0.5 + (median_Fet / 200)
    d = (1 - Q1) / (1 - Q0)
    g = (Q1 * (1 - Q0)) / (Q0 * (1 - Q1))
    Upper_limit = ((np.log(8) / total_count) - np.log(d)) / np.log(g)
    Lower_limit = ((np.log(1/8) / total_count) - np.log(d)) / np.log(g)
    Q3 = 0.5 - (median_Fet / 200)
    Q4 = 0.5
    d_wt = (1 - Q4) / (1 - Q3)
    g_wt = (Q4 * (1 - Q3)) / (Q3 * (1 - Q4))
    Upper_limit_wt = ((np.log(8) / total_count) - np.log(d_wt)) / np.log(g_wt)
    Lower_limit_wt = ((np.log(1/8) / total_count) - np.log(d_wt)) / np.log(g_wt)

    prediction = None

    print(alt_count)

    # Classification based on thresholds + whether there are informative snps
    if (informative_snps == 0):
        prediction = "Prediction not possible, no informative SNPs found"
    elif (alt_count / total_count) > Upper_limit:
        prediction = "Homozygous Mutant"
    elif (alt_count / total_count) > Lower_limit:
        prediction = "Inconclusive between Homozygous Mutant and Heterozygous"
    elif (alt_count / total_count) > Upper_limit_wt:
        prediction = "Heterozygous"
    elif (alt_count / total_count) > Lower_limit_wt:
        prediction = "Inconclusive between Homozygous Wild Type and Heterozygous"
    else:
        prediction = "Homozygous Wild Type"

    return prediction, mean_pat, median_pat, IQR_Pat, mean_Fet, median_Fet, IQR_Fet, FL_SNPs, d, g, d_wt, g_wt


def mat_hom_mut_gt_prediction(fetal_frac,total_count,alt_count):
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    # read the fetal fraction file into a pandas df
    FL_SNPs = pd.read_csv(fetal_frac, sep = '\t')
    # remove the last five lines which do not contain informative SNP information
    FL_SNPs = FL_SNPs.iloc[:-5]
    
    # Rename the columns
    FL_SNPs.columns = ["Chromosome", "Start", "End", "Num.reads", "A", "A_fraction", 
                    "G", "G_fraction", "C", "C_fraction", "T", "T_fraction"]
    
    # assign number of informative snps to variable
    informative_snps=len(FL_SNPs)

    # Calculate the max value for A_fraction, G_fraction, C_fraction, T_fraction
    FL_SNPs['ColMax'] = FL_SNPs[["A_fraction", "G_fraction", "C_fraction", "T_fraction"]].max(axis=1)
    FL_SNPs['PaternalFraction'] = 100 - FL_SNPs['ColMax']
    FL_SNPs['FetalFraction'] = 2 * FL_SNPs['PaternalFraction']

    # Summary statistics
    mean_pat = FL_SNPs['PaternalFraction'].mean()
    median_pat = FL_SNPs['PaternalFraction'].median()
    IQR_Pat = FL_SNPs['PaternalFraction'].quantile(0.75) - FL_SNPs['PaternalFraction'].quantile(0.25)
    mean_Fet = FL_SNPs['FetalFraction'].mean()
    median_Fet = FL_SNPs['FetalFraction'].median()
    IQR_Fet = FL_SNPs['FetalFraction'].quantile(0.75) - FL_SNPs['FetalFraction'].quantile(0.25)

    
    # Displaying the table (in Jupyter, pandas automatically outputs it)
    # In a real HTML or PDF generation, use pandas to_html() or to_latex()
    #FL_SNPs

    # SPRT setup: Mother is BB (mutant)
    # Goal: Test if fetus is BB (null) or AB (heterozygous, alternative)

    f = median_Fet / 100  # Convert fetal fraction from % to proportion

    # Under the null hypothesis (fetus also BB), no A allele should be present
    # this causes a division by 0 later so need to reformat test to be assessing proportion of B
    Q0 = 0.999  # Expected proportion of B allele under H0 - foetus is homozygous mutant
    # with some noise - could calc that somehow using other allele counts?

    # Under the alternative (fetus is AB), fetus contributes half A alleles
    Q1 = 1- (f / 2)

    # SPRT parameters
    d = (1 - Q1) / (1 - Q0)
    g = (Q1 * (1 - Q0)) / (Q0 * (1 - Q1))
    # Upper and lower decision limits
    Upper_limit = ((np.log(8) / total_count) - np.log(d)) / np.log(g)
    Lower_limit = ((np.log(1 / 8) / total_count) - np.log(d)) / np.log(g)


    prediction = None

    # Classification logic for Mother = BB (mutant)
    if informative_snps == 0:
        prediction = "Prediction not possible, no informative SNPs found"
    else:
        observed_ratio = alt_count / total_count  # alt_count = B alleles (observed B allele frequency)

        # If the observed ratio is above the upper limit, it's more likely the fetus is BB
        if observed_ratio > Upper_limit:
            prediction = "Homozygous Mutant"  # Fetus is more likely to be BB

        # If the observed ratio is between the upper and lower limits, it is inconclusive
        elif observed_ratio > Lower_limit:
            prediction = "Inconclusive between Homozygous Mutant and Heterozygous"

        # If the observed ratio is below the lower limit, it's more likely the fetus is AB
        else:
            prediction = "Heterozygous"  # Fetus is more likely to be AB


    return prediction, mean_pat, median_pat, IQR_Pat, mean_Fet, median_Fet, IQR_Fet, FL_SNPs, d, g

def mat_hom_wt_gt_prediction(fetal_frac,total_count,alt_count):
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    # read the fetal fraction file into a pandas df
    FL_SNPs = pd.read_csv(fetal_frac, sep = '\t')
    # remove the last five lines which do not contain informative SNP information
    FL_SNPs = FL_SNPs.iloc[:-5]
    
    # Rename the columns
    FL_SNPs.columns = ["Chromosome", "Start", "End", "Num.reads", "A", "A_fraction", 
                    "G", "G_fraction", "C", "C_fraction", "T", "T_fraction"]
    
    # assign number of informative snps to variable
    informative_snps=len(FL_SNPs)

    # Calculate the max value for A_fraction, G_fraction, C_fraction, T_fraction
    FL_SNPs['ColMax'] = FL_SNPs[["A_fraction", "G_fraction", "C_fraction", "T_fraction"]].max(axis=1)
    FL_SNPs['PaternalFraction'] = 100 - FL_SNPs['ColMax']
    FL_SNPs['FetalFraction'] = 2 * FL_SNPs['PaternalFraction']

    # Summary statistics
    mean_pat = FL_SNPs['PaternalFraction'].mean()
    median_pat = FL_SNPs['PaternalFraction'].median()
    IQR_Pat = FL_SNPs['PaternalFraction'].quantile(0.75) - FL_SNPs['PaternalFraction'].quantile(0.25)
    mean_Fet = FL_SNPs['FetalFraction'].mean()
    median_Fet = FL_SNPs['FetalFraction'].median()
    IQR_Fet = FL_SNPs['FetalFraction'].quantile(0.75) - FL_SNPs['FetalFraction'].quantile(0.25)

    # SPRT setup: Mother is AA (wild-type)
    # Goal: Test if fetus is AA (null) or AB (heterozygous, alternative)

    f = median_Fet / 100  # Convert fetal fraction from % to proportion

    # Under the null hypothesis (fetus also AA), no B allele should be present
    Q0 = 0.001  # Expected proportion of B allele under H0

    # Under the alternative (fetus is AB), fetus contributes half B alleles
    # B alleles will appear in plasma at frequency ≈ fetal fraction / 2
    Q1 = f / 2

    # SPRT parameters
    d = (1 - Q1) / (1 - Q0)
    g = (Q1 * (1 - Q0)) / (Q0 * (1 - Q1))

    # Upper and lower decision limits (log-likelihood thresholds)
    Upper_limit = ((np.log(8) / total_count) - np.log(d)) / np.log(g)
    Lower_limit = ((np.log(1 / 8) / total_count) - np.log(d)) / np.log(g)

    prediction = None

    # Classification logic for Mother = AA
    if informative_snps == 0:
        prediction = "Prediction not possible, no informative SNPs found"
    else:
        observed_ratio = alt_count / total_count  # alt_count = B alleles
        if observed_ratio > Upper_limit:
            prediction = "Heterozygous"
        elif observed_ratio > Lower_limit:
            prediction = "Inconclusive between Homozygous Wild Type and Heterozygous"
        else:
            prediction = "Homozygous Wild Type"


    return prediction, mean_pat, median_pat, IQR_Pat, mean_Fet, median_Fet, IQR_Fet, FL_SNPs, d, g


# # create separate graph functions, save somewhere / pdf report? Ask Karisma re desired output format

#     # Generate SPRT plot

    #     x_vals = np.arange(50, 120001)
    #     Upper_limit_graph_rmd = ((np.log(8) / x_vals) - np.log(d)) / np.log(g)
    #     Lower_limit_graph_rmd = ((np.log(1/8) / x_vals) - np.log(d)) / np.log(g)
    #     Upper_limit_graph_wt = ((np.log(8) / x_vals) - np.log(d_wt)) / np.log(g_wt)
    #     Lower_limit_graph_wt = ((np.log(1/8) / x_vals) - np.log(d_wt)) / np.log(g_wt)
    
    #   g_range = (min(Upper_limit_graph_rmd.min(), Lower_limit_graph_rmd.min(), 
    #             Upper_limit_graph_wt.min(), Lower_limit_graph_wt.min(), alt_count/total_count),
    #         max(Upper_limit_graph_rmd.max(), Lower_limit_graph_rmd.max(), 
    #             Upper_limit_graph_wt.max(), Lower_limit_graph_wt.max(), alt_count/total_count))

    # plt.plot(x_vals, Upper_limit_graph_rmd, color='red', label='Upper Limit RMD')
    # plt.plot(x_vals, Lower_limit_graph_rmd, color='pink', label='Lower Limit RMD')
    # plt.plot(x_vals, Upper_limit_graph_wt, color='green', label='Upper Limit WT')
    # plt.plot(x_vals, Lower_limit_graph_wt, color='blue', label='Lower Limit WT')
    # plt.scatter([total_count], [alt_count / total_count], color='blue', label='Observed')

    # plt.title("Modified SPRT")
    # plt.xlabel("Total number of counts")
    # plt.ylabel("Pr over-represented allele")
    # plt.legend()
    # plt.show()

# #     # Chromosome 11 plot with regions of interest
#     xrange = np.arange(5225264, 5227272)
#     yrange = [0, 1.2]

#     HBB_exon_1 = np.arange(5227071, 5226929, -1)
#     HBB_exon_2 = np.arange(5226799, 5226576, -1)
#     HBB_exon_3 = np.arange(5225726, 5225463, -1)
#     Fetal_Het_alt = [((100 - median_pat) / 100)] * len(xrange)
#     Fetal_Hom_alt = [(0.5 + (median_pat / 100))] * len(xrange)
#     Fetal_Hom_ref = [(0.5 - (median_pat / 100))] * len(xrange)
#     Fetal_Het_ref = [(median_pat / 100)] * len(xrange)
#     Variant_of_interest = [(alt_count / total_count)]

#     plt.plot(xrange, [0.5] * len(xrange), color='black', label='Chr11')
#     plt.plot(HBB_exon_1, [0.5] * len(HBB_exon_1), color='purple', linewidth=3, label='HBB exon 1')
#     plt.plot(HBB_exon_2, [0.5] * len(HBB_exon_2), color='green', linewidth=3, label='HBB exon 2')
#     plt.plot(HBB_exon_3, [0.5] * len(HBB_exon_3), color='orange', linewidth=3, label='HBB exon 3')
#     plt.plot(xrange, Fetal_Het_alt, color='cadetblue', linewidth=4, label='Fetal Het Alt')
#     plt.scatter([5226925], Variant_of_interest, color='red', label='Variant of Interest')

#     plt.title("Chromosome 11 Plot")
#     plt.xlabel("Position")
#     plt.ylabel("Fraction")
#     plt.legend()
#     plt.show()

   

