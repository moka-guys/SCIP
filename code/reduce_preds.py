def reduce_preds(preds):

    # String to remove
    target_string = "Homozygous Wild Type"

    # Remove items with value equal to target_string
    filtered_preds = {key: value for key, value in preds.items() if value != target_string}

    return filtered_preds