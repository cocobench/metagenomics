#!/usr/bin/env python3

import argparse
import sys

# todo parser for biobox binning data format
# todo implement taxonomic tree to resolve full lineage for above

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
    try:
        ground_truth_set = set(ground_truth_dict[rank].keys())
    except KeyError:
        ground_truth_set = set()
    try:
        tool_set = set(tool_dict[rank].keys())
    except KeyError:
        tool_set = set()

    true_positives = ground_truth_set.union(tool_set)
    false_positves = tool_set.difference(ground_truth_set)
    false_negatives = ground_truth_set.difference(tool_set)

    return len(true_positives), len(false_positves), len(false_negatives)


def calculate_recall(tp, fn):
    try:
        return tp / (tp + fn)
    except ZeroDivisionError:
        return 0.0


def calculate_precision(tp, fp):
    try:
        return tp / (tp + fp)
    except ZeroDivisionError:
        return 0.0


def calculate_f_score(recall, precision):
    try:
        return 2 * (precision * recall) / (precision + recall)
    except ZeroDivisionError:
        return 0.0


def setup_argument_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('-g', '--ground-truth', required=True, help='Path to the file containing the ground truth')
    parser.add_argument('-t', '--tool-result', required=True, help='Path to the file containing the tool result')
    parser.add_argument('-o', '--output', required=True, help='Path to the output file to write')
    parser.add_argument('-r', '--ranks', type=parse_comma_sep_list, default=['genus', 'species', 'strain'],
                        help='Comma seperated list of ranks to calculate stats on (Default: Genus,Species,Strain)')

    return parser


def parse_comma_sep_list(csl):
    split = csl.split(',')
    return [x.strip().lower() for x in split]


def main():
    parser = setup_argument_parser()
    cli_args = parser.parse_args(sys.argv[1:])

    ground_dict = parse_bioboxes_profiling_file(cli_args.ground_truth)
    tool_dict = parse_bioboxes_profiling_file(cli_args.tool_result)

    with open(cli_args.output, 'w') as outfile:
        outfile.write('RANK\tTP\tFP\tFN\tRECALL\tPRECISION\tF-SCORE\n')
        for rank in cli_args.ranks:
            tp, fp, fn = compare_profiling_dictionaries(rank, ground_dict, tool_dict)
            recall = calculate_recall(tp, fn)
            precision = calculate_precision(tp, fp)
            f_score = calculate_f_score(recall, precision)
            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                rank, tp, fp, fn, recall, precision, f_score))


if __name__ == '__main__':
    main()
