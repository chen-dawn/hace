#!/usr/bin/env python3

"""
Given BAM inputs and the reference FASTA file, count the number of nucleotides in the alignment at each position and plots a mutation heatmap.

python ~/chenlab/Dawn/TRACE/get_pileup_table_trace_edit.py \
    -b TitrationPlate2-A01_S97_aligned_to_MEK1i1_300bp.bam \
    -r ../MEK1i1_300bp.fasta
"""

# from re import A
import numpy as np
from collections import defaultdict
from Bio import SeqIO
import os
import pandas as pd
from collections import Counter
import pysam
import argparse


def get_pileups_table_samfile(samfile, ref, REF_NAME, start_index, end_index):
    """
    Given a samfile, and reference (SeqIO object) we get the pileup table for those reference.
    This is similar to the function below (which has stopped working for some reason), 
    but the difference is that this requires an input read from pysam.Samfile(bam, "rb").
    samfile: SAM file input from PySam.
    ref: Reference FASTA, which is opened as a Seq object.
    REF_NAME: Name of the reference. 
    start_index: Index of the start location.
    end_index: Index of the end location.
    """
    Alist = np.zeros(end_index - start_index)
    Clist = np.zeros(end_index - start_index)
    Glist = np.zeros(end_index - start_index)
    Tlist = np.zeros(end_index - start_index)
    print(start_index, end_index)

    for pos in range(start_index, end_index):
        print("POS:", pos)
        for pileupcolumn in samfile.pileup(REF_NAME, pos, pos + 1):
            for pileupread in pileupcolumn.pileups:
                if pileupcolumn.pos != pos:
                    continue
                if pileupread.query_position == None:
                    # print("PILEUP READ IS NONE", pileupread)
                    continue

                base = pileupread.alignment.query_sequence[pileupread.query_position]
                if base == "A":
                    Alist[pos - start_index] += 1
                elif base == "C":
                    Clist[pos - start_index] += 1
                elif base == "G":
                    Glist[pos - start_index] += 1
                elif base == "T":
                    Tlist[pos - start_index] += 1

    print(Alist, Clist, Glist, Tlist)

    new_ref_list = []
    for i in range(start_index, end_index):
        new_ref_list.append(ref[i])

    pileups_df = pd.DataFrame(
        {
            "pos": list(range(start_index, end_index)),
            "ref": new_ref_list,
            "A": Alist,
            "C": Clist,
            "G": Glist,
            "T": Tlist,
        }
    )

    filtered_df = pileups_df[
        (pileups_df["pos"] >= start_index) & (pileups_df["pos"] <= end_index)
    ]
    return filtered_df.reset_index()


def get_pileups_table(samfile, ref, REF_NAME, start_index, end_index):
    """
    Given a samfile, and reference (SeqIO object) we get the pileup table for those reference.
    samfile: SAM file input from PySam.
    ref: Reference FASTA, which is opened as a Seq object.
    REF_NAME: Name of the reference. 
    start_index: Index of the start location.
    end_index: Index of the end location.
    """
    coverage = samfile.count_coverage(REF_NAME, start_index, end_index)
    Alist = coverage[0]
    Clist = coverage[1]
    Glist = coverage[2]
    Tlist = coverage[3]

    new_ref_list = []
    for i in range(start_index, end_index):
        new_ref_list.append(ref[i])

    # Add 1 to pos so that it is 1 indexed.
    pileups_df = pd.DataFrame(
        {
            "chr": REF_NAME,
            "n_base": list(range(start_index + 1, end_index + 1)),
            "ref_base": new_ref_list,
            "A": Alist,
            "C": Clist,
            "G": Glist,
            "T": Tlist,
        }
    )

    filtered_df = pileups_df[
        (pileups_df["n_base"] >= start_index) & (pileups_df["n_base"] <= end_index)
    ]
    coverage = pileups_df["A"] + pileups_df["C"] + pileups_df["G"] + pileups_df["T"]
    filtered_df = pileups_df[coverage > 0]
    return filtered_df.reset_index()


def subtract_ref_percentage(df):
    """
    Given a data frame (that has the columns pos, ref, A, C, G, T), 
    subtract the ref percentage from the category.
    """
    for i in range(df.shape[0]):
        ref = df["ref"][i]
        df.loc[i, ref.upper()] = 0
    return df


def parse_args():
    parser = argparse.ArgumentParser(description="Create Pileup table from BAM file.")

    parser.add_argument(
        "-b", "--bam", type=str, required=True, help="BAM file path.",
    )

    parser.add_argument(
        "-r",
        "--reference",
        type=str,
        help="Where the reference is located. Default is the MEK1.",
        default="/broad/thechenlab/Dawn/RNAreference/MEK1i1_300bp.fasta",
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        help="Directory where output files will be saved.",
    )

    args = parser.parse_args()
    return args


##################################
# The important part starts here #
##################################

if __name__ == "__main__":
    args = parse_args()
    bam_f = args.bam
    reference_filename = args.reference
    out_dir = args.output_dir
    if out_dir is None:
        out_dir = os.getcwd()

    print("Here are the args:", bam_f, reference_filename, out_dir)

    reference = SeqIO.parse(open(reference_filename), "fasta")
    for fasta in reference:
        print("Reading reference:", reference_filename)
        ref_name, ref = fasta.name, fasta.seq.upper()
        # print(ref_name, ref)

        start_index = 0
        end_index = len(ref)

        # Read the BAM files.
        print("Processing file: ", bam_f)
        basename = os.path.splitext(bam_f)[0]

        samfile = pysam.Samfile(bam_f, "rb")

        pileups_table = get_pileups_table(samfile, ref, ref_name, start_index, end_index)
        if len(pileups_table.index) < 1 or pileups_table['A'].max() < 5:
            print("No pileups found for", ref_name, "Basename is", basename)
            continue

        output_filepath = os.path.join(out_dir, basename + "_with_sequence_" + ref_name+ "_pileup.tsv")
        print("Writing to output: ", output_filepath)
        pileups_table.to_csv(path_or_buf=output_filepath, index=False, sep="\t")
