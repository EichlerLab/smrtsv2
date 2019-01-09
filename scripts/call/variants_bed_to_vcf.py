#!/usr/bin/env python2

import argparse
import datetime
import numpy as np
import pandas as pd
import pysam


def calculate_variant_quality(variant):
    try:
        return int(min(100, round(-10 * np.log10(1 - (variant.contig_support / float(variant.contig_depth))) * np.log(variant.contig_depth), 0)))
    except ZeroDivisionError:
        return 0


def convert_bed_to_vcf(bed_filename, reference_filename, vcf_filename, sample, variant_type):
    # Get variants.
    if variant_type == "sv":
        columns = (0, 1, 2, 3, 4, 5, 8, 11, 12, 15, 16, 20)
        names = ("chr", "start", "end", "sv_call", "event_size", "sv_sequence", "contig", "contig_start", "contig_end", "contig_support", "contig_depth", "repeat_type")
    elif variant_type == "indel":
        # chr1    94824   94827   3       Cttttcttttttttt 1       1       29.04   deletion
        columns = (0, 1, 2, 3, 4, 5, 6, 7, 8)
        names = ("chr", "start", "end", "event_size", "sv_sequence", "contig_support", "contig_depth", "depth", "sv_call")
    elif variant_type == "inversion":
        columns = (0, 1, 2, 3, 4, 5)
        names = ("chr", "start", "end", "sv_call", "contig_support", "contig_depth")
    else:
        raise Exception("Unsupported variant type: %s" % variant_type)

    calls = pd.read_table(bed_filename, header=None, usecols=columns, names=names)
    calls["sample_name"] = sample
    calls["call_id"] = "."
    calls["quality"] = calls.apply(calculate_variant_quality, axis=1)
    calls["filter"] = "PASS"

    # Get the reference base at the position of the variant start.
    reference = pysam.FastaFile(reference_filename)
    calls["reference"] = calls.apply(lambda row: reference.fetch(row.chr, row.start, row.start + 1).upper(), axis=1)

    # Update start position to be 1-based.
    #calls["start"] = calls["start"] + 1  # Start position is the base before the variant

    # Build an INFO field for each call.
    if variant_type == "sv":
        calls["alt"] = calls.apply(lambda row: "<%s>" % row.sv_call[:3].upper(), axis=1)
        calls["info"] = calls.apply(
            lambda row: ";".join(
                ["=".join(map(str, item))
                 for item in (
                        ("END", row.end),
                        ("SVTYPE", row.sv_call[:3].upper()),
                        ("SVLEN", row.event_size),
                        ("CONTIG", row.contig),
                        ("CONTIG_START", row.contig_start),
                        ("CONTIG_END", row.contig_end),
                        ("REPEAT_TYPE", row.repeat_type),
                        ("CONTIG_SUPPORT", row.contig_support),
                        ("CONTIG_DEPTH", row.contig_depth),
                        ("SAMPLES", row.sample_name),
                        ("SEQ", row.sv_sequence)
                 )]
            ),
            axis=1
        )
    elif variant_type == "indel":
        calls["alt"] = calls.apply(lambda row: "<%s>" % row.sv_call[:3].upper(), axis=1)
        calls["info"] = calls.apply(
            lambda row: ";".join(
                ["=".join(map(str, item))
                 for item in (
                        ("END", row.end),
                        ("SVTYPE", row.sv_call[:3].upper()),
                        ("SVLEN", row.event_size),
                        ("CONTIG_SUPPORT", row.contig_support),
                        ("CONTIG_DEPTH", row.contig_depth),
                        ("DP", row.depth),
                        ("SAMPLES", row.sample_name),
                        ("SEQ", row.sv_sequence)
                 )]
            ),
            axis=1
        )
    elif variant_type == "inversion":
        calls["alt"] = "<INV>"
        calls["info"] = calls.apply(
            lambda row: ";".join(
                ["=".join(map(str, item))
                 for item in (
                        ("END", row.end),
                        ("SVTYPE", row.sv_call[:3].upper()),
                        ("SVLEN", row.end - row.start),
                        ("CONTIG_SUPPORT", row.contig_support),
                        ("CONTIG_DEPTH", row.contig_depth),
                        ("SAMPLES", row.sample_name),
                 )]
            ),
            axis=1
        )

    simple_calls = calls[["chr", "start", "call_id", "reference", "alt", "quality", "filter", "info"]].rename_axis({"chr": "#CHROM", "start": "POS", "reference": "REF", "call_id": "ID", "quality": "QUAL", "info": "INFO", "alt": "ALT", "filter": "FILTER"}, axis=1)

    # Save genotypes as tab-delimited file.
    with open(vcf_filename, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##fileDate=%s\n" % datetime.date.strftime(datetime.date.today(), "%Y%m%d"))
        vcf.write("##source=SMRT_SV\n")
        vcf.write('##INFO=<ID=SAMPLES,Number=1,Type=String,Description="Samples with the given variant">' + "\n")
        vcf.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Mean depth of raw reads">' + "\n")
        vcf.write('##INFO=<ID=CONTIG_DEPTH,Number=1,Type=Integer,Description="Total depth of local assemblies">' + "\n")
        vcf.write('##INFO=<ID=CONTIG_SUPPORT,Number=1,Type=Integer,Description="Depth of local assemblies supporting variant">' + "\n")
        vcf.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' + "\n")
        vcf.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">' + "\n")
        vcf.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">' + "\n")
        vcf.write('##INFO=<ID=CONTIG,Number=1,Type=String,Description="Name of alternate assembly contig">' + "\n")
        vcf.write('##INFO=<ID=CONTIG_START,Number=1,Type=Integer,Description="Start coordinate of this variant in the alternate assembly contig">' + "\n")
        vcf.write('##INFO=<ID=CONTIG_END,Number=1,Type=Integer,Description="End coordinate of this variant in the alternate assembly contig">' + "\n")
        vcf.write('##INFO=<ID=REPEAT_TYPE,Number=1,Type=String,Description="Repeat classification of variant content">' + "\n")
        vcf.write('##INFO=<ID=SEQ,Number=1,Type=String,Description="Sequence associated with variant">' + "\n")
        simple_calls.to_csv(vcf, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bed", help="input BED file of variant calls")
    parser.add_argument("reference", help="FASTA file for reference used with SMRT SV")
    parser.add_argument("vcf", help="output VCF file of variant calls")
    parser.add_argument("sample", help="name of sample with variants")
    parser.add_argument("type", help="variant call type", choices=("sv", "indel", "inversion"))
    args = parser.parse_args()

    convert_bed_to_vcf(args.bed, args.reference, args.vcf, args.sample, args.type)
