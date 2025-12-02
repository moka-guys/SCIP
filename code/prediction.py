from fetal_gt_pred import *
from html_report import *

def scip_pred(report_name,output_path, fetal_frac_output_path,total_counts,alt_counts,allele_labels):
    # Initialise empty html report content variable
    html_content = ""
    html_summary_content = ""#
    
    # generate report title
    report_path = output_path

    # initiate prediction dictionary
    preds = {}
    prediction_possible = True
    mat_gt_preds = []

    # predict genotype using R script conversion
    for total, alt, label in zip(total_counts, alt_counts, allele_labels):
        # selecting only for the relevant alleles according to parental gt
        if total is not None:
            mat_gt_pred = mat_gt_prediction(fetal_frac_output_path,total,alt)
            print("The observed maternal genotype is: " + mat_gt_pred)
            mat_gt_preds.append(mat_gt_pred)
            print(mat_gt_preds)
                
            if mat_gt_pred == "Heterozygous":
                print(fetal_frac_output_path)
                print(total)
                print(alt)
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

    print(preds)     

    return prediction_possible, preds, mat_gt_preds, html_content, html_summary_content, FL_SNPs


def scip_pred_FE(mat_gt_preds, report_name,output_path, fetal_frac_output_path,total_counts,alt_counts,allele_labels):
    print(output_path)
    
    # Initialise empty html report content variable
    html_content = ""
    html_summary_content = ""#
    
    # generate report title
    report_path = output_path

    # initiate prediction dictionary
    preds = {}
    prediction_possible = True

    # predict genotype using R script conversion
    for total, alt, label, mat_pred in zip(total_counts, alt_counts, allele_labels, mat_gt_preds):
        # selecting only for the relevant alleles according to parental gt
        if total is not None:
            print("HELLO?: " + mat_pred)
                
            if mat_pred == "Heterozygous":
                pred_and_stats = mat_het_gt_prediction(fetal_frac_output_path,total,alt)
            elif mat_pred == "Homozygous mutant":
                pred_and_stats = mat_hom_mut_gt_prediction(fetal_frac_output_path,total,alt)
            elif mat_pred == "Homozygous wildtype":
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

            if mat_pred == "Heterozygous":
                # generate html content for this allele of interest
                html_content = html_content + generate_html_content(mean_pat, median_pat, IQR_pat, mean_Fet, \
                                            median_Fet, IQR_Fet, total, alt, FL_SNPs, label, report_name, d, g, d_wt, g_wt)
            else:
                html_content= html_content + generate_homozygous_html_content(mean_pat, median_pat, IQR_pat, mean_Fet, \
                                            median_Fet, IQR_Fet, total, alt, FL_SNPs, label, report_name, d, g)


            # generate summary html content for this allele of interest
            html_summary_content = html_summary_content + generate_summary_html_content(report_name,label,prediction)

    print(preds)     

    return prediction_possible, preds, html_content, html_summary_content, FL_SNPs



















# def scip_pred_FE(report_name_FE,output_path, fetal_frac_output_path_FE,total_counts_FE,alt_counts_FE,allele_labels_FE,fetal_frac_output_path, total_counts, alt_counts, allele_labels):
#     # Initialise empty html report content variable
#     html_content = ""
#     html_summary_content = ""#
    
#     # generate report title
#     report_path = output_path

#     # initiate prediction dictionary
#     preds = {}
#     prediction_possible = True

#     # predict genotype using R script conversion
#     for total, total_FE, alt, alt_FE, label, label_FE in zip(total_counts, \
#         total_counts_FE, alt_counts, alt_counts_FE, allele_labels, allele_labels_FE):
#         # selecting only for the relevant alleles according to parental gt
#         if total is not None:
#             # calculate the maternal gt using the non-FE allele counts
#             mat_gt_pred = mat_gt_prediction(fetal_frac_output_path,total,alt)
#             print("The observed maternal genotype is: " + mat_gt_pred)

#             # now make feotal gt prediction and calculate stats using FE allele counts    
#             if mat_gt_pred == "Heterozygous":
#                 print(fetal_frac_output_path_FE)
#                 print(total_FE)
#                 print(alt_FE)
#                 pred_and_stats = mat_het_gt_prediction(fetal_frac_output_path_FE,total_FE,alt_FE)
#             elif mat_gt_pred == "Homozygous mutant":
#                 pred_and_stats = mat_hom_mut_gt_prediction(fetal_frac_output_path_FE,total_FE,alt_FE)
#             elif mat_gt_pred == "Homozygous wildtype":
#                 pred_and_stats = mat_hom_wt_gt_prediction(fetal_frac_output_path_FE,total_FE,alt_FE)
#             else:
#                 print("Strange alt:total ratio: " + str(alt/total))


#             try:
#                 prediction, mean_pat, median_pat, \
#                 IQR_pat, mean_Fet, median_Fet, IQR_Fet, FL_SNPs, \
#                 d, g, d_wt, g_wt = pred_and_stats[0], pred_and_stats[1], \
#                 pred_and_stats[2], pred_and_stats[2], pred_and_stats[4], pred_and_stats[5], pred_and_stats[6], \
#                 pred_and_stats[7], pred_and_stats[8], pred_and_stats[9], pred_and_stats[10], pred_and_stats[11]
#             except:
#                 prediction, mean_pat, median_pat, \
#                 IQR_pat, mean_Fet, median_Fet, IQR_Fet, FL_SNPs, \
#                 d, g = pred_and_stats[0], pred_and_stats[1], \
#                 pred_and_stats[2], pred_and_stats[2], pred_and_stats[4], pred_and_stats[5], pred_and_stats[6], \
#                 pred_and_stats[7], pred_and_stats[8], pred_and_stats[9]
                
#             if prediction == "Prediction not possible, no informative SNPs found":
#                 prediction_possible = False
                
#             print("The predicted fetal genotype for the " + label + " allele is: " + prediction + ", when foetal enrichment is applied.")
#             label_short = label[0]
#             preds[label_short] = prediction

#             if mat_gt_pred == "Heterozygous":
#                 # generate html content for this allele of interest
#                 html_content = html_content + generate_html_content(mean_pat, median_pat, IQR_pat, mean_Fet, \
#                                             median_Fet, IQR_Fet, total, alt, FL_SNPs, label, report_name_FE, d, g, d_wt, g_wt)
#             else:
#                 html_content= html_content + generate_homozygous_html_content(mean_pat, median_pat, IQR_pat, mean_Fet, \
#                                             median_Fet, IQR_Fet, total, alt, FL_SNPs, label, report_name_FE, d, g)


#             # generate summary html content for this allele of interest
#             html_summary_content = html_summary_content + generate_summary_html_content(report_name_FE,label,prediction)

#     print(preds)     

#     return prediction_possible, preds, mat_gt_pred, html_content, html_summary_content, FL_SNPs
