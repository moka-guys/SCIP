def reduce_preds(preds):
    #print(preds)

    # String to remove
    target_string = "Homozygous Wild Type"

    # Remove items with value equal to target_string
    filtered_preds = {key: value for key, value in preds.items() if value != target_string}

    #print(filtered_preds)

    return filtered_preds