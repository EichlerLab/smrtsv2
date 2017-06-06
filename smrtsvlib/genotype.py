"""
Functions to support genotyping.
"""

import gzip
import math
import pysam
import time

import pandas as pd

import smrtsvlib


# Convert genotype names to VCF genotype annotation
GENOTYPE_TO_GT = {
    'HOM_REF': '0/0',
    'HET': '0/1',
    'HOM_ALT': '1/1',
    'NO_CALL': './.'
}


def get_sample_column(table_file_name, sample_name):
    """
    Get a VCF column as a Pandas Series for one sample.

    :param table_file_name: Name of genotyped features table file output by the genotyper after applying the genotyping
        model and annotating the genotype for no-call variants.

    :param sample_name: Name of the sample. This is saved to the name of the returned Series object.

    :return: Series of a column to add to the VCF for one genotyped sample.
    """

    df_gt = pd.read_table(
        table_file_name, header=0,
        usecols=('GENOTYPE', 'REF_COUNT', 'REF_COUNT', 'DEN_HOM_REF', 'DEN_HET', 'DEN_HOM_ALT'),
        names=('GENOTYPE', 'REF_COUNT', 'REF_COUNT', 'HOM_REF', 'HET', 'HOM_ALT')
    )

    # Set genotype (GT), genotype quality (GQ), and genotype likelihood (GL)
    df_gt['GT'] = df_gt['GENOTYPE'].apply(lambda genotype: GENOTYPE_TO_GT[genotype])

    df_gt['GQ'] = df_gt.apply(
        lambda row: int(10 * -math.log10(1 - row[row['GENOTYPE']])) if df_gt['GENOTYPE'] != 'NO_CALL' else '.',
        axis=1
    )

    df_gt['GL'] = df_gt.apply(lambda row: '{HOM_REF},{HET},{HOM_ALT}'.format(**row), axis=1)

    # Get a series representing the column to be added to the VCF
    sample_column = df_gt.apply(lambda row: '{GT}:{GQ}:{GS}:{REF_COUNT}:{ALT_COUNT}'.format(**row), axis=1)
    sample_column.name = sample_name

    # Return
    return sample_column


def vcf_table(vcf_file_name):
    """
    Get a table of the VCF file as a DataFrame. Contains each VCF record and the correct column names. Each record
    is complete up to the "FORMAT" column, and there are no samples. Genotype calls are appended as columns to
    make the output VCF file for the genotype pipeline.

    DataFrame does not include any header content other than the column names (as the column names of the frame).

    :param vcf_file_name: Name of the gzipped VCF file to load.

    :return: A DataFrame of VCF variant records.
    """

    # Read
    with gzip.open(vcf_file_name) as in_file:
        df_var = pd.read_table(in_file, header=None, comment='#')

    # Get only variant information (truncate FORMAT and any samples)
    df_var = df_var.loc[:, 0:7]

    # Add columns and FORMAT as the genotyper will write them
    df_var.columns = ('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')

    df_var['FORMAT'] = 'GT:GQ:GL:DPR:DPA'

    # Return base VCF records
    return df_var


def vcf_header_lines(vcf_file_name):
    """
    Get header lines for the genotype output VCF from the input variant-call VCF file.

    :param vcf_file_name: Name of a variant VCF file.

    :return: A list of VCF headers as strings with newlines. Does not include the column heading line at the end
        of the headers.
    """

    header_list = list()

    with pysam.VariantFile(vcf_file_name) as vcf_file:
        # Add VCF version if missing
        header_list.append(vcf_version(vcf_file))

        # Set file date
        header_list.append('##fileDate={}\n'.format(time.strftime("%Y%m%d")))

        # Set source
        header_list.extend(vcf_get_source_list(vcf_file))
        header_list.append('##source=SMRTSV_Genotyper_{}\n'.format(smrtsvlib.__version__))

        # Get header elements excluding FORMAT tags
        for header_element in vcf_file.header.records:

            # Replace FORMAT tags
            if header_element.type == 'FORMAT':
                continue

            # Source and date handled
            if header_element.type == 'GENERIC' and header_element.key.lower() in {'fileformat', 'source', 'filedate'}:
                continue

            # Write record
            header_list.append(str(header_element))

        # Add FORMAT tags written by the genotyper
        header_list.extend(vcf_get_format_tags())

    # Return header lines
    return header_list


def vcf_version(vcf_file):
    """
    Get VCF version header line or create one.

    :param vcf_file: Open pysam VCF file.

    :return: VCF version header line.
    """
    for header_element in vcf_file.header.records:
        if header_element.type == 'GENERIC' and header_element.key.lower() == 'fileformat':
            return str(header_element)

    return '##fileformat=VCFv4.2\n'


def vcf_get_source_list(vcf_file):
    """
    Get a list of generic source annotations in the VCF headers.

    :param vcf_file: Open pysam VCF file.

    :return: A list of source annotations.
    """

    source_list = list()

    for header_element in vcf_file.header.records:
        if header_element.type == 'GENERIC' and header_element.key.lower() == 'source':
            source_list.append(str(header_element))

    return source_list


def vcf_get_format_tags():
    """
    Get a list of VCF format tags records for the VCF header.

    :return: List of VCF format tags for the VCF header.
    """
    return [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n',
        '##FORMAT=<ID=GL,Number=G,Type=String,Description="Genotype likelihood">\n',
        '##FORMAT=<ID=DPR,Number=1,Type=Float,Description="Average read depth over reference allele breakpoints">\n',
        '##FORMAT=<ID=DPA,Number=1,Type=Float,Description="Average read depth over alternate allele breakpoints">\n'
    ]
