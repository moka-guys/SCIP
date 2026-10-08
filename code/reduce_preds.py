def reduce_preds(preds):

    # String to remove
    target_string = "Homozygous Wild Type"

    # If all values are Homozygous Wild Type, keep only 2 entries
    if all(value == target_string for value in preds.values()):
        return dict(list(preds.items())[:2])

    # Otherwise, remove all Homozygous Wild Type entries
    filtered_preds = {key: value for key, value in preds.items() if value != target_string}

    print(filtered_preds)

    return filtered_preds