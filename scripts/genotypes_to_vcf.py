import argparse
import datetime
import numpy as np
import pandas as pd
import pysam


def convert_table_to_vcf(genotypes_filename, calls_filename, reference_filename, vcf_filename):
    # Load genotypes in long format.
    genotypes = pd.read_table(genotypes_filename)

    # Get all SVs.
    vcf = pysam.VariantFile(calls_filename)
    calls = []
    for record in vcf:
        if isinstance(record.info["SVLEN"], tuple):
            sv_length = record.info["SVLEN"][0]
        else:
            sv_length = record.info["SVLEN"]

        calls.append([
            record.chrom,
            record.start,
            record.info["END"],
            record.info["SVTYPE"],
            sv_length,
            record.info["SEQ"],
            record.info["CONTIG"],
            record.info["CONTIG_START"],
            record.info["CONTIG_END"],
            record.info["CONTIG_SUPPORT"],
            record.info["CONTIG_DEPTH"],
            record.info["REPEAT_TYPE"]
        ])

    calls = pd.DataFrame(calls, columns=("chr", "start", "end", "sv_call", "event_size", "sv_sequence", "contig", "contig_start", "contig_end", "support", "depth", "repeat_type"))
    calls["call_id"] = calls.apply(lambda row: "-".join(map(str, (row.chr, row.start, row.end, row.sv_call))), axis=1)
    calls["quality"] = calls.apply(lambda row: int(min(100, round(-10 * np.log10(1 - (row.support / float(row.depth))) * np.log(row.depth), 0))), axis=1)

    # Get the reference base at the position of the variant start.
    reference = pysam.FastaFile(reference_filename)
    calls["reference"] = calls.apply(lambda row: reference.fetch(row.chr, row.start, row.start + 1).upper(), axis=1)

    # Update start position to be 1-based.
    calls["start"] = calls["start"] + 1

    # Assign an id to each call.
    genotypes["call_id"] = genotypes.apply(lambda row: "-".join(map(str, (row.chr, row.start, row.end, row.sv_call))), axis=1)
    sample_specific_fields = ["sample", "genotype", "genotype_quality", "genotype_likelihoods", "concordant", "discordant"]
    genotypes = genotypes[sample_specific_fields + ["call_id"]]

    # Build a genotype field for each call/sample combination.
    format = "GT:GQ:PL:DPR:DPA"
    genotypes["genotype_data"] = genotypes.apply(lambda row: ":".join(map(str, [row.genotype, row.genotype_quality, row.genotype_likelihoods, int(round(row.discordant, 0)), int(round(row.concordant, 0))])), axis=1)

    # Pivot genotypes to wide format to match VCF style.
    wide_genotypes = genotypes.pivot("call_id", "sample", "genotype_data")
    wide_genotypes = wide_genotypes.sort_index(axis=0)
    wide_genotypes.insert(0, "FORMAT", format)

    # Build an INFO field for each call.
    calls["info"] = calls.apply(lambda row: ";".join(["=".join(map(str, item)) for item in (("END", row.end), ("SVTYPE", row.sv_call), ("SVLEN", row.event_size), ("CONTIG", row.contig), ("CONTIG_START", row.contig_start), ("CONTIG_END", row.contig_end), ("REPEAT_TYPE", row.repeat_type), ("DPA", row.support), ("DP", row.depth), ("SEQ", row.sv_sequence))]), axis=1)

    calls["sv_sequence"] = calls.apply(lambda row: "<DEL>" if row["sv_call"] == "deletion" else "<INS>", axis=1)
    simple_calls = calls[["chr", "start", "call_id", "reference", "sv_sequence", "quality", "info"]].rename_axis({"chr": "#CHROM", "start": "POS", "reference": "REF", "call_id": "ID", "quality": "QUAL", "info": "INFO", "sv_sequence": "ALT"}, axis=1)

    # Annotate distinct calls with wide genotypes.
    calls_with_genotypes = simple_calls.merge(wide_genotypes, how="inner", left_on="ID", right_index=True)

    # Add reference, alternate allele, quality, and filter fields.
    calls_with_genotypes.insert(6, "FILTER", "PASS")

    # Save genotypes as tab-delimited file.
    with open(vcf_filename, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##fileDate=%s\n" % datetime.date.strftime(datetime.date.today(), "%Y%m%d"))
        vcf.write("##source=SMRT_SV_genotyper\n")
        vcf.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth of local assemblies">' + "\n")
        vcf.write('##INFO=<ID=DPA,Number=1,Type=Integer,Description="Depth of local assemblies supporting variant">' + "\n")
        vcf.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' + "\n")
        vcf.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">' + "\n")
        vcf.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of this variant">' + "\n")
        vcf.write('##INFO=<ID=CONTIG,Number=1,Type=String,Description="Name of alternate assembly contig">' + "\n")
        vcf.write('##INFO=<ID=CONTIG_START,Number=1,Type=Integer,Description="Start coordinate of this variant in the alternate assembly contig">' + "\n")
        vcf.write('##INFO=<ID=CONTIG_END,Number=1,Type=Integer,Description="End coordinate of this variant in the alternate assembly contig">' + "\n")
        vcf.write('##INFO=<ID=REPEAT_TYPE,Number=1,Type=String,Description="Repeat classification of variant content">' + "\n")
        vcf.write('##INFO=<ID=SEQ,Number=1,Type=String,Description="Sequence associated with variant">' + "\n")
        vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + "\n")
        vcf.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">' + "\n")
        vcf.write('##FORMAT=<ID=PL,Number=G,Type=String,Description="Genotype likelihood">' + "\n")
        vcf.write('##FORMAT=<ID=DPR,Number=1,Type=String,Description="Read depth supporting reference allele">' + "\n")
        vcf.write('##FORMAT=<ID=DPA,Number=1,Type=String,Description="Read depth supporting alternate allele">' + "\n")
        calls_with_genotypes.to_csv(vcf, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genotypes", help="tab-delimited file of sample, call coordinates in query sequence, call type, and genotypes in standard format (e.g., 0/0)")
    parser.add_argument("calls", help="VCF of original SV calls with repeat type from SMRT SV")
    parser.add_argument("reference", help="FASTA file for reference used with SMRT SV")
    parser.add_argument("output", help="genotypes in VCF format")
    args = parser.parse_args()

    convert_table_to_vcf(args.genotypes, args.calls, args.reference, args.output)
