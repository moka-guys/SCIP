import pandas as pd

def identify_allele(row):
    # S allele
    if row.iloc[0] == "chr11" and row.iloc[1] == 5227002:
        return 'S'
    # C allele
    elif row.iloc[0] == "chr11" and row.iloc[1] == 5227003:
        return 'C'
    # E allele
    elif row.iloc[0] == "chr11" and row.iloc[1] == 5226943:
        return 'E'
    # D allele
    elif row.iloc[0] == "chr11" and row.iloc[1] == 5225678:
        return 'D'
    
#print(SCED_pileup_file.apply(identify_allele, axis = 1))



# softcode filename TODO
def total_and_alt_counts(SCED_pileup_file,parental_alleles):
    SCED_pileup_file = pd.read_csv(SCED_pileup_file, delimiter="\t", header = None)

    # add new column to df which contains allele label
    SCED_pileup_file[6] = SCED_pileup_file.apply(identify_allele, axis = 1)

    # subset the df to include only rows for the parental alleles
    filtered_SCED_pileup_file = SCED_pileup_file[SCED_pileup_file[6].isin(parental_alleles)]

    # initiate df rows for total and alt count results
    processed_rows = []

    for index, row in filtered_SCED_pileup_file.iterrows():
        ref_allele = row.iloc[2]
        total_count = int(row.iloc[3])
        calls = row.iloc[4]
        allele = row.iloc[6]
        
        # Replace '.' and ',' with corresponding bases
        if row.iloc[2] == 'A':
            calls = calls.replace('.', 'A').replace(',', 'a')
        elif row.iloc[2] == 'C':
            calls = calls.replace('.', 'C').replace(',', 'c')
        elif row.iloc[2] == 'G':
            calls = calls.replace('.', 'G').replace(',', 'g')
        elif row.iloc[2] == 'T':
            calls = calls.replace('.', 'T').replace(',', 't')
        
        # Count alternate base counts
        ref_count = calls.count(ref_allele.upper()) + calls.count(ref_allele.lower())
        alt_count = total_count - ref_count

        # Append the processed row to the list
        processed_rows.append([allele, total_count, alt_count])

    # convert list of processed rows into a dataframe
    total_alt_count_df = pd.DataFrame(processed_rows, columns=["allele", "total_count", "alt_count"])

    return total_alt_count_df

def counts_to_variables(total_alt_count_df):
    S_total, S_alt, C_total, C_alt, E_total, E_alt, D_total, D_alt = None, None, None, None, None, None, None, None

    for index, row in total_alt_count_df.iterrows():
        if row.iloc[0] == "S":
            S_total = row.iloc[1]
            S_alt = row.iloc[2]
        elif row.iloc[0] == "C":
            C_total = row.iloc[1]
            C_alt = row.iloc[2]
        elif row.iloc[0] == "D":
            D_total = row.iloc[1]
            D_alt = row.iloc[2]
        elif row.iloc[0] == "E":
            E_total = row.iloc[1]
            E_alt = row.iloc[2]
    
    return S_total, S_alt, C_total, C_alt, E_total, E_alt, D_total, D_alt
        
def total_and_alt_vars(SCED_pileup_file, parental_alleles):
    df = total_and_alt_counts(SCED_pileup_file, parental_alleles)
    return counts_to_variables(df)