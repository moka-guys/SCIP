
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def gt_prediction(fetal_frac,total_count,alt_count):
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
    FL_SNPs

    # Print important stats
    # print(f"Interquartile range of the fetal fraction: {IQR_Fet}")
    # print(f"Median of the fetal fraction: {median_Fet}")
    # print(f"Number of SNPs tested: {len(FL_SNPs)}")
    # print(f"Total count: {total_count}")
    # print(f"Alternative count: {alt_count}")
    # print(f"Alternative ratio: {alt_count/total_count}")

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

    # Classification based on thresholds
    if (alt_count / total_count) > Upper_limit:
        prediction = "Homozygous Mutant"
    elif (alt_count / total_count) > Lower_limit:
        prediction = "Inconclusive between Homozygous Mutant and Heterozygous"
    elif (alt_count / total_count) > Upper_limit_wt:
        prediction = "Heterozygous"
    elif (alt_count / total_count) > Lower_limit_wt:
        prediction = "Inconclusive between Homozygous Wild Type and Heterozygous"
    else:
        prediction = "Homozygous Wild Type"

    return prediction



# # create separate graph functions, save somewhere / pdf report? Ask Karisma re desired output format

#     # Generate SPRT plot

    #     x_vals = np.arange(50, 120001)
    #     Upper_limit_graph_rmd = ((np.log(8) / x_vals) - np.log(d)) / np.log(g)
    #     Lower_limit_graph_rmd = ((np.log(1/8) / x_vals) - np.log(d)) / np.log(g)
    #     Upper_limit_graph_wt = ((np.log(8) / x_vals) - np.log(d_wt)) / np.log(g_wt)
    #     Lower_limit_graph_wt = ((np.log(1/8) / x_vals) - np.log(d_wt)) / np.log(g_wt)
#     
#       g_range = (min(Upper_limit_graph_rmd.min(), Lower_limit_graph_rmd.min(), 
#                 Upper_limit_graph_wt.min(), Lower_limit_graph_wt.min(), alt_count/total_count),
#             max(Upper_limit_graph_rmd.max(), Lower_limit_graph_rmd.max(), 
#                 Upper_limit_graph_wt.max(), Lower_limit_graph_wt.max(), alt_count/total_count))

#     plt.plot(x_vals, Upper_limit_graph_rmd, color='red', label='Upper Limit RMD')
#     plt.plot(x_vals, Lower_limit_graph_rmd, color='pink', label='Lower Limit RMD')
#     plt.plot(x_vals, Upper_limit_graph_wt, color='green', label='Upper Limit WT')
#     plt.plot(x_vals, Lower_limit_graph_wt, color='blue', label='Lower Limit WT')
#     plt.scatter([total_count], [alt_count / total_count], color='blue', label='Observed')

#     plt.title("Modified SPRT")
#     plt.xlabel("Total number of counts")
#     plt.ylabel("Pr over-represented allele")
#     plt.legend()
#     plt.show()

#     # Chromosome 11 plot with regions of interest
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

   

