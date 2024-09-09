import statistics

# function that counts each base in row, calcs fraction
# and total count
#TODO docstring
def process_row(depth, chr, start, end, num_reads, counts, output_file, fetal_fraction_list):
    num_A = counts.count('A') + counts.count('a')
    num_G = counts.count('G') + counts.count('g')
    num_C = counts.count('C') + counts.count('c')
    num_T = counts.count('T') + counts.count('t')

    num_A_Frac = (num_A / num_reads) * 100
    num_G_Frac = (num_G / num_reads) * 100
    num_C_Frac = (num_C / num_reads) * 100
    num_T_Frac = (num_T / num_reads) * 100

    # TODO - currently this does not take into account whether the read depth has reached the minimum set, produces a number I think is misleading
    is_informative = False
    if 1 < num_A_Frac < 20:
        fetal_fraction_list.append(num_A_Frac * 2)
        is_informative = True
    elif 1 < num_G_Frac < 20:
        fetal_fraction_list.append(num_G_Frac * 2)
        is_informative = True
    elif 1 < num_C_Frac < 20:
        fetal_fraction_list.append(num_C_Frac * 2)
        is_informative = True
    elif 1 < num_T_Frac < 20:
        fetal_fraction_list.append(num_T_Frac * 2)
        is_informative = True

    if is_informative:
        output_file.write(f"{chr}\t{start}\t{end}\t{num_reads}\t{num_A}\t{num_A_Frac}\t{num_G}\t{num_G_Frac}\t{num_C}\t{num_C_Frac}\t{num_T}\t{num_T_Frac}\n")
        print(f"{chr}\t{start}\t{end}\t{num_reads}\t{num_A}\t{num_A_Frac}\t{num_G}\t{num_G_Frac}\t{num_C}\t{num_C_Frac}\t{num_T}\t{num_T_Frac}")

    return is_informative

# feed mpileup file and output file path from scip.py / the sample sheet.
def fetal_frac(depth, HBB_mpileup_file, output_file_path):
    # Open files for reading and writing
    with open(HBB_mpileup_file, "r") as r1, open(output_file_path, "w") as fw2:
        fetal_fraction_list = []
        informative_snps = 0

        # Write header
        fw2.write("chr\tstart\tend\tnum_reads\tA\tA_Frac\tG\tG_Frac\tC\tC_Frac\tT\tT_Frac\n")

        for row in r1:
            columns = row.strip().split("\t")
            chr = columns[0]
            start = int(columns[1])
            end = start + 1
            num_reads = columns[3]
            calls = columns[4]

            # Ensure num_reads is numeric
            if not num_reads.isdigit():
                continue

            num_reads = int(num_reads)

            # check whether number of reads meets depth requirement. 
            # If not, move on from row
            if num_reads < depth:
                continue

            # Replace '.' and ',' with corresponding bases
            if columns[2] == 'A':
                calls = calls.replace('.', 'A').replace(',', 'a')
            elif columns[2] == 'C':
                calls = calls.replace('.', 'C').replace(',', 'c')
            elif columns[2] == 'G':
                calls = calls.replace('.', 'G').replace(',', 'g')
            elif columns[2] == 'T':
                calls = calls.replace('.', 'T').replace(',', 't')

            # Process the row
            informative_row = process_row(depth, chr, start, end, num_reads, calls, fw2, fetal_fraction_list)
            if informative_row:
                informative_snps += 1

        # Compute statistics
        if fetal_fraction_list:
            average = statistics.mean(fetal_fraction_list)
            st_dev = statistics.stdev(fetal_fraction_list)
        else:
            average = 0
            st_dev = 0

        fw2.write(f"Min Read Depth Required: {depth}\n")
        fw2.write(f"Fetal Fraction Estimate: {average}\n")
        fw2.write(f"Fetal StdDev: {st_dev}\n")
        fw2.write(f"No. Informative SNPS at Min Read Depth: {informative_snps}\n")
        fw2.write("Informative SNP defined as >2% and <40% Fetal Fraction\n") # this value doesn't take into account whether the SNP reaches the minimum depth.

####### test ###########
#fetal_frac(350,"SCIP259_HBB_targets.mpileup","SCIP259_fetal_frac_output.txt")