import pandas as pd

#function ss_to_df:
def ss_to_df(sample_sheet):
    ss = pd.read_csv(sample_sheet)
    ss = ss.reset_index() # make sure indexes pair with number of rows
    return ss

#pull out patient details into class when given specific row index
def extract_patient_row(row_index, dataframe):
    patient_row = dataframe.loc[row_index]
    return patient_row

def row_to_vars(row_index, dataframe):
    # extract 1 row of the dataframe
    patient_row = extract_patient_row(row_index, dataframe)

    # pick out the values from that row and assign to variables
    sample= patient_row.loc["Sample"]
    mat_gt = patient_row.loc["Mat_GT"]
    pat_gt = patient_row.loc["Pat_GT"]
    HBB_filename = patient_row.loc["HBB_targets"]
    SCED_filename = patient_row.loc["SCED_alleles"]

    # return variables
    return [sample, mat_gt, pat_gt, HBB_filename, SCED_filename]

############ TESTING ###################
# call ss_to_df and save dataframe to variable, for use in other functions
# ss = ss_to_df()

# patient_row = extract_patient_row(row_index=0, ss_to_df())

# current_patient = patient(patient_row.loc["Sample"])

# print(current_patient.sample_name)


######## CLASS STUFF - CURRENTLY REDUNDANT ###############

# class patient:
#     # currently uncertain as to whether there's any point to this
#     def __init__(self, sample_name, mat_gt,
#                   pat_gt, HBB_filename, SCED_filename):
#         self.sample_name = sample_name
#         self.mat_gt = mat_gt
#         self.pat_gt = pat_gt
#         self.HBB_filename = HBB_filename
#         self.SCED_filename = SCED_filename

# def patient_info_to_object(patient_row):
#     current_patient =  patient(sample_name = patient_row.loc["Sample"],
#                                mat_gt = patient_row.loc["Mat_GT"],
#                                pat_gt = patient_row.loc["Pat_GT"],
#                                HBB_filename = patient_row.loc["HBB_targets"],
#                                SCED_filename = patient_row.loc["SCED_alleles"])
#     return current_patient









