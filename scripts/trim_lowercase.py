#!/bin/env python
"""
"""
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re


def trim_lowercase(input_filename, output_filename, keep_completely_lowercased_sequences):
    """
    Trim all lowercase letters from the given input FASTA file, remove records
    that were completely lowercase, and write the trimmed FASTA to the given
    output file.
    """
    trimmed_records = []
    for record in SeqIO.parse(input_filename, "fasta"):
        # Remove all lowercase letters in the sequence for this record.
        trimmed_sequence = re.sub("(^[a-z]+|[a-z]+$)", "", str(record.seq))

        # If any uppercase sequence remains, create a new record for it.
        if len(trimmed_sequence) > 0:
            trimmed_records.append(SeqRecord(Seq(trimmed_sequence), id=record.id, description=""))
        elif keep_completely_lowercased_sequences:
            trimmed_records.append(SeqRecord(Seq(str(record.seq)), id=record.id, description=""))

    # Write out all trimmed records to the given output file.
    SeqIO.write(trimmed_records, output_filename, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    parser.add_argument("--keep_completely_lowercased_sequences", action="store_true")
    args = parser.parse_args()

    trim_lowercase(args.input_file, args.output_file, args.keep_completely_lowercased_sequences)
