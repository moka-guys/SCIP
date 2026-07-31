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
import statistics

class SCIP(object):
    """ Runs the SCIP (Sickle Cell in Pregnancy) analysis for predicting foetal sickle cell disease status from 
    allele count data.

    The full workflow includes: extracting allele counts from mpileup files, estimating foetal fraction,
    generating genotype and clinical predicitons and producing an HTML report. Where initial prediction is
    inconclusive, analysis with foetal enrichment is attempted. The final HTML report reflects whichever
    analysis yielded a conclusive clinical prediction.

    Attributes:
        None (analysis is executed immediately on instantiation via main())
    """
    
    def __init__(self):
        """ Initialises the SCIP analysis object and immediately triggers the full analysis workflow.

        No parameters are required at instantiation; all inputs are read from command-line
        arguments parsed in the `__main__` block.
        """
        
        self.main()

    def main(self):
        """
        Coordinates and executes the full SCIP analysis workflow.

        Steps performed:
            1. Extracts total and alternate allele counts from mpileup files for SCD-relevant alleles (S, C, E, D).
            2. Estimates foetal fraction using HBB informative SNPs.
            3. Generates foetal genotype and clinical predictions.
            4. If the initial prediction is inconclusive, re-runs analysis using foetal-enriched (155bp) data.
            5. Writes an HTML report reflecting whichever analysis yielded a conclusive clinical prediction.
            If both analyses are inconclusive, the non-enriched result is reported.

        Raises:
            FileNotFoundError: If any of the required input files (mpileup/HBB) are missing.
            ValueError: If allele count extraction or foetal fraction estimation fails due to malformed input.

        Side effects:
            - Creates the output directory if it does not already exist.
            - Writes one or more HTML report files to `output_dir`.
            - Prints intermediate prediction results and diagnostic messages to stdout.

        Input globals (set from sys.argv in __main__):
            output_dir (str): Path to the directory where reports will be written.
            scip_id (str): Sample identifier used in report naming and file paths.
            hbb_file (str): Path to the HBB mpileup file for foetal fraction estimation.
            sced_file (str): Path to the SCD allele mpileup file (full-length reads).
            hbb_file_155bp (str): Path to the HBB mpileup file for foetal-enriched analysis.
            sced_file_155bp (str): Path to the SCD allele mpileup file (155bp enriched reads).
        """

        def scip_FE_analysis():
            # extract total and alt counts of SCD alleles to variables
            total_counts_FE, alt_counts_FE, allele_labels_FE, fetal_frac_output_path_FE = counts_labels_fetal_frac_path(sced_file_155bp, alleles,"fetal_frac_FE.txt")
            
            # Count informative SNPs in foetally enriched data
            informative_snp_count_FE = fetal_frac(x,hbb_file_155bp,fetal_frac_output_path_FE)
            print("FE analysis done with minimum coverage set to " + str(x))

            prediction_possible_FE, preds_FE, mat_gt_preds, html_content_FE, html_summary_content_FE, FL_SNPs_FE = scip_pred(report_name_FE,\
                scip_id, fetal_frac_output_path_FE,total_counts_FE,alt_counts_FE,allele_labels_FE)

            print(f'prediction_possible_FE: {prediction_possible_FE}')
                    
            # If all values are Homozygous Wild Type, keep only 2 entries
            if all(value == "Prediction not possible, no informative SNPs found" for value in preds_FE.values()): 
                print("Prediction not possible with Foetal Enrichment, no informative SNPs")           
                        
            if prediction_possible_FE:
                preds_FE = reduce_preds(preds_FE)
                print(f'preds_FE after reduce_preds: {preds_FE}')
                gt_FE = resolve_gt(preds_FE)
                print(f'gt_FE: {gt_FE}')
                clin_pred_FE = resolve_clin(gt_FE)
                print("Genotype prediction with foetal enrichment applied: " + gt_FE)
                print("Clinical prediction with foetal enrichment applied: " + clin_pred_FE)
            
                return prediction_possible_FE, clin_pred_FE, preds_FE, FL_SNPs_FE, informative_snp_count_FE, html_content_FE, html_summary_content_FE
        
        def scip_analysis():
            # extract total and alt counts of SCD alleles to variables
            total_counts, alt_counts, allele_labels, fetal_frac_output_path = counts_labels_fetal_frac_path(sced_file, alleles,"fetal_frac.txt")
            
            # Count informative SNPs in non fe data
            informative_snp_count = fetal_frac(x,hbb_file,fetal_frac_output_path)
            print("non-FE analysis done with minimum coverage set to " + str(x))

            prediction_possible, preds, mat_gt_preds, html_content, html_summary_content, FL_SNPs = scip_pred(report_name,\
                scip_id, fetal_frac_output_path,total_counts,alt_counts,allele_labels)

            print(f'prediction_possible: {prediction_possible}')
                    
            # If all values are Homozygous Wild Type, keep only 2 entries
            if all(value == "Prediction not possible, no informative SNPs found" for value in preds.values()): 
                print("Prediction not possible with Foetal Enrichment, no informative SNPs")           
                        
            if prediction_possible:
                preds = reduce_preds(preds)
                print(f'preds after reduce_preds: {preds}')
                gt = resolve_gt(preds)
                print(f'gt: {gt}')
                clin_pred = resolve_clin(gt)
                print("Genotype prediction without foetal enrichment applied: " + gt)
                print("Clinical prediction without foetal enrichment applied: " + clin_pred)
            
                return prediction_possible, clin_pred, preds, FL_SNPs, informative_snp_count, html_content, html_summary_content

        def write_fe_report():
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

        def write_report():
            # update report name to include clin prediction
            report_name=scip_id + " Report: " + clin_pred
            # generate html header for report
            html_header = generate_html_header(report_name)
            html_table = generate_html_table(report_name,FL_SNPs,informative_snp_count)
            html_clin_pred = generate_clinical_summary_html(report_name_FE, clin_pred)
            # combine html contents
            all_html = html_header + html_clin_pred + html_summary_content + html_table + html_content
            # output html_content to html report
            with open (report_path, 'w') as f:
                f.write(all_html)
       



        # make the output directory specified, if it does not already exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Set report name for foetal enrichment analysis
        report_name_FE=scip_id + "_155bp Report"

        # set output report paths for analysis with and without foetal enrichment applied.
        report_path = os.path.join(output_dir, f"{scip_id}_Report.html")
        report_path_FE = os.path.join(output_dir, f"{scip_id}_155bp_Report.html")
               
        # set the names of the SCD alleles being analysed
        alleles = ["S","C","E","D"]

        # set minimum coverage for SNP to be included in foetal fraction calculation
        x = 100
        
        # run analysis using foetally enriched data
        prediction_possible_FE, clin_pred_FE, preds_FE, FL_SNPs_FE, informative_snp_count_FE, html_content_FE, html_summary_content_FE = scip_FE_analysis()
       
        print(f'prediction_possible_FE: {prediction_possible_FE}')

        if prediction_possible_FE:            
            if clin_pred_FE != "Inconclusive":
                # write report if foetal enrichment analysis produces conclusive result
                write_fe_report()
            else:         
                # non FE analysis performed if FE result is inconclusive
                prediction_possible, clin_pred, preds, FL_SNPs, informative_snp_count, html_content, html_summary_content = scip_analysis()
                print(f'prediction_possible: {prediction_possible}')
                if prediction_possible:
                    if clin_pred != "Inconclusive":
                        # if non fe analysis is conclusive, write that report
                        write_report()
                    else:
                        # if non fe analysis is also inconclusive, revert to fe report.
                        write_fe_report()
        else:
            # if no predictions are creating using foetally enriched data, run the analyis on non-fe and report.
            prediction_possible, clin_pred, preds, FL_SNPs, informative_snp_count, html_content, html_summary_content = scip_analysis()
            print(f'prediction_possible: {prediction_possible}')
            if prediction_possible:
                write_report()
            else:
                # if prediction also isn't possible with non fe analysis
                print("prediction not possible")
                   
if __name__ == "__main__":
    output_dir = sys.argv[1]
    scip_id = sys.argv[2]
    hbb_file = sys.argv[3]
    sced_file = sys.argv[4]
    hbb_file_155bp = sys.argv[5]
    sced_file_155bp= sys.argv[6]
    SCIP()
