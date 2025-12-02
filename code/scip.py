import sys 
import re
from total_and_alt_count import *
from fetal_frac_calc import *
from fetal_gt_pred import *
from html_report import *
from resolve_gt_and_clin import *
from reduce_preds import *
from prediction import *
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

        os.makedirs(output_dir, exist_ok=True)
        
        report_name=scip_id + " Report"

        report_path = os.path.join(output_dir, f"{scip_id}_Report.html")
        report_path_FE = os.path.join(output_dir, f"{scip_id}_155bp_Report.html")
               

        alleles = ["S","C","E","D"]

        # extract total and alt counts of SCD alleles to variables
        total_counts, alt_counts, allele_labels, fetal_frac_output_path = counts_labels_fetal_frac_path(sced_file, alleles,"fetal_frac_output.txt")

        x = 100
        
        informative_snp_count = fetal_frac(x,hbb_file,fetal_frac_output_path)
        print("Analysis done with minimum coverage set to " + str(x))

        # Initialise empty html report content variable
        html_content = ""
        html_summary_content = ""

        # initiate prediction dictionary
        preds = {}
        prediction_possible = True

        prediction_possible, preds, mat_gt_preds, html_content, html_summary_content, FL_SNPs = scip_pred(report_name,\
            scip_id, fetal_frac_output_path,total_counts,alt_counts,allele_labels)

        
        if prediction_possible:
            preds = reduce_preds(preds)
            gt = resolve_gt(preds)
            clin_pred = resolve_clin(gt)
            
            print("Genotype prediction without foetal enrichment applied: " + gt)
            print("Clinical prediction without foetal enrichment applied: " + clin_pred)

            if clin_pred == "Inconclusive":
                report_name_FE=scip_id + "_155bp Report"
                output_path_FE = output_path + "_155bp"
                total_counts_FE, alt_counts_FE, allele_labels_FE, fetal_frac_output_path_FE = counts_labels_fetal_frac_path(sced_file_155bp, alleles,"fetal_frac_FE.txt")
                informative_snp_count_FE = fetal_frac(x,hbb_file_155bp,fetal_frac_output_path_FE)
                print("FE analysis done with minimum coverage set to " + str(x))

                prediction_possible_FE, preds_FE, html_content_FE, html_summary_content_FE, FL_SNPs_FE = scip_pred_FE(mat_gt_preds, report_name_FE,\
                    output_path_FE, fetal_frac_output_path_FE,total_counts_FE,alt_counts_FE,allele_labels_FE)

                preds_FE = reduce_preds(preds_FE)
                gt_FE = resolve_gt(preds_FE)
                clin_pred_FE = resolve_clin(gt_FE)
                print("Genotype prediction with foetal enrichment applied: " + gt_FE)
                print("Clinical prediction with foetal enrichment applied: " + clin_pred_FE)


                if clin_pred_FE != "Inconclusive":
                    # update report name to include clin prediction
                    report_name_FE=scip_id + "_155bp Report: " + clin_pred_FE
                    # generate html header for report
                    html_header_FE = generate_html_header(report_name_FE)

                    html_table_FE = generate_html_table(report_name_FE,FL_SNPs_FE,informative_snp_count_FE)

                    html_clin_pred_FE = generate_clinical_summary_html(report_name_FE, clin_pred_FE)

                    # combine html contents
                    all_html_FE = html_header_FE + html_clin_pred_FE + html_summary_content_FE + html_table_FE + html_content_FE
                    #report_path_FE = scip_id + "_155bp_report.html"
                    # output html_content to html report
                    with open (report_path_FE, 'w') as f:
                        f.write(all_html_FE)

                else:
                    # if 155bp also makes inconclusive prediction, use the non FE result
                    # update report name to include prediction
                    report_name=scip_id + " Report : " + clin_pred
                    # generate html header for report
                    html_header = generate_html_header(report_name)
                    html_table = generate_html_table(report_name,FL_SNPs,informative_snp_count)
                    html_clin_pred = generate_clinical_summary_html(report_name, clin_pred)

                    # combine html contents
                    all_html = html_header + html_clin_pred + html_summary_content + html_table + html_content

                    # output html_content to html report
                    with open (report_path, 'w') as f:
                        f.write(all_html)

            else: 
                # update report name to include prediction
                report_name=scip_id + " Report: " + clin_pred
                # generate html header for report
                html_header = generate_html_header(report_name)
                html_clin_pred = generate_clinical_summary_html(report_name, clin_pred)

                html_table = generate_html_table(report_name,FL_SNPs,informative_snp_count)

                # combine html contents
                all_html = html_header + html_clin_pred + html_summary_content + html_table + html_content

                # output html_content to html report
                with open (report_path, 'w') as f:
                    f.write(all_html)

        else:
            print("Prediction not possible, no informative SNPs found")            

if __name__ == "__main__":
    output_dir = sys.argv[1]
    scip_id = sys.argv[2]
    hbb_file = sys.argv[3]
    sced_file = sys.argv[4]
    hbb_file_155bp = sys.argv[5]
    sced_file_155bp= sys.argv[6]
    SCIP()
