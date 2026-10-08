# you cannot have an S mutation and a C mutation on the same HBB gene - bc would result in major structural anomaly and isn't rly seen in practice.
def remove_from_end_until_length(s, target_length, char_to_remove):
    while len(s) > target_length:
        index = s.rfind(char_to_remove)
        if index != -1:
            s = s[:index] + s[index+1:]
        else:
            break
    return s

def remove_duplicates_until_length(s, target_length):
    while len(s) > target_length:
        seen = set()  # Set to track characters we've already added
        result = []  # List to build the resulting string

        for char in s:
            if char not in seen:
                result.append(char)  # Add unique characters to the result
                seen.add(char)  # Mark the character as seen
            # Stop when the result reaches the target length
            if len(result) == target_length:
                break
        s = ''.join(result)
    return s  

def move_char_to_front(s, char_to_move):
    # If the character is in the string, we proceed to move it to the front
    if char_to_move in s:
        # Count how many times the character occurs in the string
        count = s.count(char_to_move)
        # Remove all occurrences of the character and create the result string
        result = char_to_move * count + s.replace(char_to_move, '')
        return result
    return s  # If the character isn't found, return the original string

def genotype_string_processing(s):
    s = remove_from_end_until_length(s,2,"A")
    #s = remove_duplicates_until_length(s,2)
    s = move_char_to_front(s,"A")
    return s

def resolve_gt(prediction_dictionary):
    result = ""

    pred_items = list(prediction_dictionary.items())

    num_of_alleles = len(pred_items)

    allele = ""

    if num_of_alleles == 1:
        allele = pred_items[0][0]
        prediction_to_gt = {"Homozygous Mutant": allele+allele, "Heterozygous": "A" + allele, \
                        "Homozygous Wild Type":"AA","Inconclusive between Homozygous Mutant and Heterozygous" : f"{allele + allele} / A{allele}", \
                        "Inconclusive between Homozygous Wild Type and Heterozygous": f"A{allele} / AA"}  
        prediction = pred_items[0][1]

        result = prediction_to_gt[prediction]

    elif num_of_alleles == 2:
        allele_1 = pred_items[0][0]

        prediction_1 = pred_items[0][1]

        allele_2 = pred_items[1][0]

        prediction_2 = pred_items[1][1]

        prediction_to_gt_1 = {"Homozygous Mutant": [allele_1+allele_1], "Heterozygous": ["A" + allele_1], \
                        "Homozygous Wild Type":["AA"],"Inconclusive between Homozygous Mutant and Heterozygous" : [allele_1 + allele_1, "A" + allele_1], \
                        "Inconclusive between Homozygous Wild Type and Heterozygous": ["A" + allele_1, "AA"]}  
        prediction_to_gt_2 = {"Homozygous Mutant": [allele_2+allele_2], "Heterozygous": ["A" + allele_2], \
                        "Homozygous Wild Type":["AA"],"Inconclusive between Homozygous Mutant and Heterozygous" : [allele_2 + allele_2, "A" + allele_2], \
                        "Inconclusive between Homozygous Wild Type and Heterozygous": ["A" + allele_2, "AA"]}  
        
        result_1 = prediction_to_gt_1[prediction_1]
        result_2 = prediction_to_gt_2[prediction_2]


        allele_1_option_count = len(result_1)
        allele_2_option_count = len(result_2)


        
        if allele_1_option_count == 1 and allele_2_option_count == 1:
            result = str(result_1[0]) + str(result_2[0])
            result = genotype_string_processing(result)
            

        elif allele_1_option_count == 2 and allele_2_option_count == 1:
            option_1 = result_1[0] + result_2[0]
            option_2 = result_1[1] + result_2[0]
            option_1 = genotype_string_processing(option_1)
            option_2 = genotype_string_processing(option_2)
            result = f"Inconclusive between {option_1} and {option_2}"

        elif allele_1_option_count == 1 and allele_2_option_count == 2:
            option_1 = result_1[0] + result_2[0]
            option_2 = result_1[0] + result_2[1]
            option_1 = genotype_string_processing(option_1)
            option_2 = genotype_string_processing(option_2)
            result = f"Inconclusive between {option_1} and {option_2}"

        elif allele_1_option_count == 2 and allele_2_option_count == 2:
            option_1 = result_1[0] + result_2[0]
            option_2 = result_1[0] + result_2[1]
            option_3 = result_1[1] + result_2[0]
            option_4 = result_1[1] + result_2[1]

            option_1 = genotype_string_processing(option_1)
            option_2 = genotype_string_processing(option_2)
            option_3 = genotype_string_processing(option_3)
            option_4 = genotype_string_processing(option_4)

            option_list = [option_1,option_2,option_3, option_4]
            unique_options = []
            seen = set()

            for item in option_list:
                if item not in seen:
                    unique_options.append(item)
                    seen.add(item)


            if len(unique_options) == 4:
                result = f"{unique_options[0]} or {unique_options[1]} or {unique_options[2]} or {unique_options[3]}"
                
            elif len(unique_options) == 3:
                result = f"{unique_options[0]} or {unique_options[1]} or {unique_options[2]}"
            elif len(unique_options) == 2:
                result = f"{unique_options[0]} or {unique_options[1]}"

    return result

def resolve_clin(gt_pred):
    if not 'A' in gt_pred:
        return "Affected"
    else:
        # get rid of spaces
        gt_preds = gt_pred.replace(" ", "").upper()
        # Split composite genotypes like "AS/SS"
        parts = gt_preds.split("/") # produces list of genotypes the prediction is inconclusive between
        for gt in parts:
            if not 'A' in gt:
                return "Inconclusive"
            else:
                continue
        return "Not affected"

        



