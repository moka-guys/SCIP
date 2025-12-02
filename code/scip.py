import sys 
import re
from total_and_alt_count import *
from fetal_frac_calc import *
from fetal_gt_pred import *
from html_report import *
from resolve_gt_and_clin import *
from reduce_preds import *
import pandas as pd
import math
import os
import matplotlib.pyplot as plt

class SCIP(object):
    """
    Description of class.
    list and description of class attributes. 
    """
    
    def __init__(self):
        #self.sample_sheet = sample_sheet
        self.main()

    def main(self):
        """
        Doc string for function
        """

        # Regular expression to match "SCIP" followed by digits
        match = re.search(r'(SCIP.*?)(?=\.html)', output_path)
        if match:
            # Extract the matched part
            scip_id = match.group(0)
            print("Extracted SCIP ID:", scip_id)
        else:
            print("No SCIP ID found in the string.")
        
        report_name=scip_id + " Report"
               

        alleles = ["S","C","E","D"]

        # extract total and alt counts of SCD alleles to variables
        S_total, S_alt, C_total, C_alt, E_total, E_alt, D_total, D_alt = total_and_alt_vars(sced_file, alleles)
        total_counts = [S_total, C_total, E_total, D_total]

        alt_counts = [S_alt, C_alt, E_alt, D_alt]

        allele_labels = ["S allele", "C allele", "E allele", "D allele"]

        fetal_frac_output_path = "fetal_frac_output.txt"

        x = 100
        
        informative_snp_count = fetal_frac(x,hbb_file,fetal_frac_output_path)
        print("Analysis done with minimum coverage set to " + str(x))
        
        # try:
        #     fetal_frac(500,hbb_file,fetal_frac_output_path)
        # except:
        #     fetal_frac(250,hbb_file,fetal_frac_output_path)

        # Initialise empty html report content variable
        html_content = ""
        html_summary_content = ""

        # generate report title
        report_path = output_path

        # initiate prediction dictionary
        preds = {}
        prediction_possible = True

        # predict genotype using R script conversion
        for total, alt, label in zip(total_counts, alt_counts, allele_labels):
            # selecting only for the relevant alleles according to parental gt
            if total is not None:
                mat_gt_pred = mat_gt_prediction(fetal_frac_output_path,total,alt)
                print("The observed maternal genotype is: " + mat_gt_pred)
                
                if mat_gt_pred == "Heterozygous":
                    pred_and_stats = mat_het_gt_prediction(fetal_frac_output_path,total,alt)
                elif mat_gt_pred == "Homozygous mutant":
                    pred_and_stats = mat_hom_mut_gt_prediction(fetal_frac_output_path,total,alt)
                elif mat_gt_pred == "Homozygous wildtype":
                    pred_and_stats = mat_hom_wt_gt_prediction(fetal_frac_output_path,total,alt)
                else:
                    print("Strange alt:total ratio: " + str(alt/total))


                try:
                    prediction, mean_pat, median_pat, \
                    IQR_pat, mean_Fet, median_Fet, IQR_Fet, FL_SNPs, \
                    d, g, d_wt, g_wt = pred_and_stats[0], pred_and_stats[1], \
                    pred_and_stats[2], pred_and_stats[2], pred_and_stats[4], pred_and_stats[5], pred_and_stats[6], \
                    pred_and_stats[7], pred_and_stats[8], pred_and_stats[9], pred_and_stats[10], pred_and_stats[11]
                except:
                    prediction, mean_pat, median_pat, \
                    IQR_pat, mean_Fet, median_Fet, IQR_Fet, FL_SNPs, \
                    d, g = pred_and_stats[0], pred_and_stats[1], \
                    pred_and_stats[2], pred_and_stats[2], pred_and_stats[4], pred_and_stats[5], pred_and_stats[6], \
                    pred_and_stats[7], pred_and_stats[8], pred_and_stats[9]
                
                if prediction == "Prediction not possible, no informative SNPs found":
                    prediction_possible = False
                
                print("The predicted fetal genotype for the " + label + " allele is: " + prediction)
                label_short = label[0]
                preds[label_short] = prediction

                if mat_gt_pred == "Heterozygous":
                    # generate html content for this allele of interest
                    html_content = html_content + generate_html_content(mean_pat, median_pat, IQR_pat, mean_Fet, \
                                                median_Fet, IQR_Fet, total, alt, FL_SNPs, label, report_name, d, g, d_wt, g_wt)
                else:
                    html_content= html_content + generate_homozygous_html_content(mean_pat, median_pat, IQR_pat, mean_Fet, \
                                                median_Fet, IQR_Fet, total, alt, FL_SNPs, label, report_name, d, g)


                # generate summary html content for this allele of interest
                html_summary_content = html_summary_content + generate_summary_html_content(report_name,label,prediction)


        
        if prediction_possible:
            preds = reduce_preds(preds)
            print(preds)
            gt = resolve_gt(preds)
            print(gt)
            clin_pred = resolve_clin(gt)

            print(clin_pred)

            report_name = report_name + ": " + gt

            # generate html header for report
            html_header = generate_html_header(report_name)

            html_table = generate_html_table(report_name,FL_SNPs,informative_snp_count)

            # combine html contents
            all_html = html_header + html_summary_content + html_table + html_content

            # output html_content to html report
            with open (report_path, 'w') as f:
                f.write(all_html)

        else:
            print("Prediction not possible, no informative SNPs found")            

if __name__ == "__main__":
    output_path = sys.argv[1]
    hbb_file = sys.argv[2]
    sced_file = sys.argv[3]
    hbb_file_155bp = sys.argv[4]
    sced_file_155bp= sys.argv[5]
    SCIP()
