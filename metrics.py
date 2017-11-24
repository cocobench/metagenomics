#!/usr/bin/env python3


# sensitivity / specificity (can't do specificity without true negatives)
# recall / precision
# f1
# tp / fp / fn (possibly ignore abundance for profiling for now)

# all on different tax levels

# input is bioboxes formats: https://github.com/bioboxes/rfc/tree/master/data-format

# ground truth for each data set in two formats (binning and profiling)

# then compare with what the tool provided


def parse_bioboxes_profiling_file(input_file):
    result_dict = {}
    with open(input_file, 'r') as f:
        for line in f:
            if line[0] == "@" or line[0] == "#" or line[0] == "\n":
                continue
            fields = line.rstrip().split('\t')
            taxid = int(fields[0]) if fields[0].isdigit() else 0
            rank = fields[1]
            name = fields[3].split("|")[-1]
            abundance = float(fields[4])

            try:
                result_dict[rank][taxid] = (name, abundance)
            except KeyError:
                result_dict[rank] = {taxid: (name, abundance)}

    return result_dict


def compare_profiling_dictionaries(rank, ground_truth_dict, tool_dict):
    ground_truth_set = set(ground_truth_dict[rank].keys())
    tool_set = set(tool_dict[rank].keys())

    true_positives = ground_truth_set.union(tool_set)
    false_positves = tool_set.difference(ground_truth_set)
    false_negatives = ground_truth_set.difference(tool_set)

    return len(true_positives), len(false_positves), len(false_negatives)


def calculate_recall(tp, fn):
    return tp / (tp + fn)


def calculate_precision(tp, fp):
    return tp / (tp + fp)


def calculate_f_score(recall, precision):
    return 2 * (precision * recall) / (precision + recall)


def main():
    


if __name__ == '__main__':
    main()
