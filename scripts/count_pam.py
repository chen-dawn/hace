"""
Count number of spacer coverage sites (based on existance of NGG PAM sequences) in the human genome.
"""
import os
import re
import numpy as np
from Bio import SeqIO


def count_bases_covered_by_NGG_PAM_sites(dna_sequence, spacer_length=20):
    # Iterate through every base of the DNA sequence
    dna_coverage_array = np.zeros(len(dna_sequence), dtype=int)
    for i in range(len(dna_sequence)):
        # Check if the current base is part of a PAM site
        if dna_sequence[i:i+2] == "GG":
            # Check if the PAM site is long enough to cover the spacer
            if i - spacer_length >= 0:
                # Mark the spacer region as covered
                dna_coverage_array[i - spacer_length-1:i-1] = 1
        # Also check the reverse complement of the PAM site
        if dna_sequence[i:i+2] == "CC":
            # Check if the PAM site is long enough to cover the spacer
            if i + spacer_length < len(dna_sequence):
                # Mark the spacer region as covered
                dna_coverage_array[i+1:i + spacer_length + 1] = 1

    return np.sum(dna_coverage_array), len(dna_sequence)

def count_bases_covered_by_NGG_PAM_sites_BEwindow(dna_sequence, spacer_length=20):
    # Iterate through every base of the DNA sequence
    dna_coverage_array = np.zeros(len(dna_sequence), dtype=int)
    for i in range(len(dna_sequence)):
        # Check if the current base is part of a PAM site
        if dna_sequence[i:i+2] == "GG":
            # Check if the PAM site is long enough to cover the spacer
            if i - spacer_length >= 0:
                # Mark base 5-8 as covered
                dna_coverage_array[i - 20+4:i-20+9] = 1
        # Also check the reverse complement of the PAM site
        if dna_sequence[i:i+2] == "CC":
            # Check if the PAM site is long enough to cover the spacer
            if i + spacer_length < len(dna_sequence):
                # Mark the spacer region as covered
                dna_coverage_array[i+1+12:i + 17 + 1] = 1

    return np.sum(dna_coverage_array), len(dna_sequence)

# Loop through all the chromosomes of the human genome

genome_filepath = "/Users/dawnxi/Documents/genome/GRCh38.p14.genome.fa"

# output table.
output_table = []
output_table2 = []
output_table3 = []
with open(genome_filepath, "r") as genome_file:
    for record in SeqIO.parse(genome_file, "fasta"):
        # Get the name of the chromosome
        chromosome_name = record.id

        print("Processing chromosome", chromosome_name)
        # Get the sequence of the chromosome
        dna_sequence = str(record.seq)
        # Count the number of NGG PAM sites
        # num_bases_covered, sequence_length = count_bases_covered_by_NGG_PAM_sites(dna_sequence)

        # Also count for if the PAM site is 200 bp spacer length.
        # num_bases_covered2, sequence_length2 = count_bases_covered_by_NGG_PAM_sites(dna_sequence, spacer_length=200)

        num_bases_covered3, sequence_length3 = count_bases_covered_by_NGG_PAM_sites_BEwindow(dna_sequence, spacer_length=20)
        # print("Num of bases covered by PAM sites:", num_bases_covered)
        # print("Sequence length:", sequence_length)
        # Add the results to the output table
        # output_table.append([chromosome_name, num_bases_covered, sequence_length])
        # output_table2.append([chromosome_name, num_bases_covered2, sequence_length2])
        output_table3.append([chromosome_name, num_bases_covered3, sequence_length3])

# Save the output table to a file
# output_filepath = "/Users/dawnxi/Documents/genome/CRISPR_num_PAM_sites_per_chromosome_20bp.tsv"
# with open(output_filepath, "w") as output_file:
#     # Write the header of the table
#     output_file.write("chr\tnum_covered\tseq_len\n")
#     # Write the results for each chromosome
#     for row in output_table:
#         output_file.write("\t".join(map(str, row)) + "\n")

# output_filepath2 = "/Users/dawnxi/Documents/genome/CRISPR_num_PAM_sites_per_chromosome_200bp.tsv"
# with open(output_filepath2, "w") as output_file:
#     # Write the header of the table
#     output_file.write("chr\tnum_covered\tseq_len\n")
#     # Write the results for each chromosome
#     for row in output_table2:
#         output_file.write("\t".join(map(str, row)) + "\n")

output_filepath3 = "/Users/dawnxi/Documents/genome/CRISPR_num_PAM_sites_per_chromosome_BEwindow.tsv"
with open(output_filepath3, "w") as output_file:
    # Write the header of the table
    output_file.write("chr\tnum_covered\tseq_len\n")
    # Write the results for each chromosome
    for row in output_table3:
        output_file.write("\t".join(map(str, row)) + "\n")


# print("Results saved to", output_filepath)
