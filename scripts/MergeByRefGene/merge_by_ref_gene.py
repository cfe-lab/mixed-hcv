#! /usr/bin/env python

__author__ = 'rliang'

import HyPhy
import hyphyAlign as align
import csv
import itertools
import os
import logging

import H77_ref
import amino_acids
from collections import defaultdict

logging.basicConfig(level=logging.WARNING)


# Instantiate a single-threaded instance of HyPhy for alignment.
hyphy = HyPhy._THyPhy(os.getcwd(), 1)
# Use the default settings, which are good for HIV amino acid sequences.
align.change_settings(hyphy)


def aa_dict_factory():
    aa_dict = {}
    for aa in amino_acids.aminos:
        aa_dict[aa] = 0
    return aa_dict


def merge_by_ref_gene(input_file, output_file):

    # A table to track the counts by gene.  Each one of these tables has
    # as its keys the gene reference position and as its values another dictionary
    # of amino acid frequencies, e.g.
    # gene_coverage_tables["NS3"][44]["A"] = the total count of A's appearing at
    # position 44 (in reference coordinates) of NS3.
    gene_coverage_tables = {
        "NS3": defaultdict(aa_dict_factory),
        "NS5a": defaultdict(aa_dict_factory),
        "NS5b": defaultdict(aa_dict_factory)
    }

    input_csv = csv.DictReader(input_file)
    # The rows of input_csv look like:
    # rname, gene, aa.pos (relative to the *query*), A, C, D, ... (all amino acids)
    for key, group in itertools.groupby(input_csv, lambda row: (row["rname"], row["gene"])):
        # Now group contains only the entries of the CSV file with the given rname and gene.
        # We'll compute the consensus AA sequence and put it into this dictionary.
        curr_conseq = {}
        curr_gene_ref = ""
        if key[1] == "NS3":
            curr_gene_ref = H77_ref.NS3
        elif key[1] == "NS5a":
            curr_gene_ref = H77_ref.NS5a
        elif key[1] == "NS5b":
            curr_gene_ref = H77_ref.NS5b
        else:
            continue

        logging.info("Processing data with reference: {}, gene: {}".format(key[0], key[1]))

        # We'll need to iterate over this more than once, so:
        rname_gene_rows = list(group)
        for row in rname_gene_rows:
            # row["aa.pos"] is zero-based.
            curr_query_position = int(row["aa.pos"]) + 1

            curr_conseq_aa = ""
            curr_conseq_count = -1
            for amino_acid in amino_acids.aminos:
                curr_aa_count = int(row[amino_acid])
                if curr_aa_count > curr_conseq_count:
                    curr_conseq_count = curr_aa_count
                    curr_conseq_aa = amino_acid

            assert(curr_conseq_count > 0)
            assert(curr_conseq_aa is not "")

            curr_conseq[curr_query_position] = curr_conseq_aa

        conseq_min_position = min(curr_conseq)
        conseq_max_position = max(curr_conseq)
        conseq_list = []

        missing_position = False
        for position in range(conseq_min_position, conseq_max_position+1):
            if position not in curr_conseq:
                conseq_list.append("X")
                missing_position = True

            else:
                conseq_list.append(curr_conseq[position])
        conseq = "".join(conseq_list)

        aligned_conseq, aligned_gene_ref, alignment_score = align.pair_align(hyphy, curr_gene_ref, conseq)

        if missing_position:
            logging.debug("Note: there were missing characters in the consensus sequence.")
            logging.debug("Alignment score: {}".format(alignment_score))
            logging.debug("Query: {}".format(aligned_conseq))
            logging.debug("Reference: {}".format(aligned_gene_ref))

        # We can now create the mapping between reference position and H77 position.
        conseq_to_gene_ref = {}
        conseq_pos = conseq_min_position - 1
        gene_ref_pos = 0
        for zero_based_pos, gene_ref_character in enumerate(aligned_gene_ref):
            one_based_pos = zero_based_pos + 1

            if gene_ref_character is not "-":
                gene_ref_pos += 1

            if aligned_conseq[zero_based_pos] is not "-":
                conseq_pos += 1

            if gene_ref_character is not "-" and aligned_conseq[zero_based_pos] is not "-":
                conseq_to_gene_ref[conseq_pos] = gene_ref_pos

        curr_gene = rname_gene_rows[0]["gene"]

        logging.debug("conseq_to_gene_ref: {}".format(conseq_to_gene_ref))

        for row in rname_gene_rows:
            # Again we have to convert this to 1-based coordinates.  If this consensus coordinate
            # is not in the lookup table, that means that it's an insertion, so we skip it.
            try:
                curr_pos = conseq_to_gene_ref[int(row["aa.pos"]) + 1]
            except KeyError:
                continue

            for aa in amino_acids.aminos:
                gene_coverage_tables[curr_gene][curr_pos][aa] += int(row[aa])


    # We can now create the output, which has columns:
    # gene, ref.pos, A, C, D, ...
    output_csv = csv.DictWriter(
        output_file,
        fieldnames=["gene", "ref.pos"] + amino_acids.aminos,
        lineterminator=os.linesep
    )
    output_csv.writeheader()

    for gene in gene_coverage_tables:
        for position in sorted(gene_coverage_tables[gene]):
            row_to_add = gene_coverage_tables[gene][position]
            row_to_add["gene"] = gene
            row_to_add["ref.pos"] = position

            output_csv.writerow(row_to_add)


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Merge sequences mapping to different reference sequences to the same coordinate system.")

    parser.add_argument("input_csv", type=argparse.FileType("rb", 0),
                        help="Input CSV file to merge")
    parser.add_argument("output_csv", type=argparse.FileType("wb", 0),
                        help="Output CSV file containing merged coverage maps")
    args = parser.parse_args()

    merge_by_ref_gene(args.input_csv, args.output_csv)


if __name__ == "__main__":
    main()