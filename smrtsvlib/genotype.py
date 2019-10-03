#!/usr/bin/env python3

"""
Functions to support genotyping.
"""

import gzip
import math
import pysam
import time

import numpy as np
import pandas as pd

import smrtsvlib


# Convert genotype names to VCF genotype annotation
GENOTYPE_TO_GT = {
    'HOM_REF': '0/0',
    'HET': '0/1',
    'HOM_ALT': '1/1',
    'NO_CALL': './.'
}


def get_sample_column(table_file_name, sample_name, sex='U'):
    """
    Get a VCF column as a Pandas Series for one sample.

    :param table_file_name: Name of genotyped features table file output by the genotyper after applying the genotyping
        model and annotating the genotype for no-call variants.
    :param sex: "M", "F", or "U" depending on if the sample is male, female, or unknown. For males, chrX is never het.
        for females, chrY is absent. For males or unknown, chrY is never het.

    :param sample_name: Name of the sample. This is saved to the name of the returned Series object.

    :return: Series of a column to add to the VCF for one genotyped sample.
    """

    df_gt = pd.read_csv(
        table_file_name, sep='\t', header=0,
        usecols=('#CHROM', 'CALLABLE', 'BP_REF_COUNT', 'BP_ALT_COUNT', 'HOM_REF', 'HET', 'HOM_ALT')
    )

    df_gt = df_gt.loc[:, ('#CHROM', 'CALLABLE', 'BP_REF_COUNT', 'BP_ALT_COUNT', 'HOM_REF', 'HET', 'HOM_ALT')]

    # Adjust density estimates on sex
    if sex == 'M':
        adjust_chrx_for_males(df_gt)
        adjust_chry(df_gt, False)

    elif sex == 'F':
        adjust_chry(df_gt, True)

    elif sex == 'U':
        adjust_chry(df_gt, False)

    # Set genotype (GT), genotype quality (GQ), and genotype likelihood (GL)
    df_gt['CLASS'] = df_gt.apply(
        lambda row: str(np.argmax(row[['HOM_REF', 'HET', 'HOM_ALT']])) if row['CALLABLE'] else 'NO_CALL',
        axis=1
    )

    df_gt['GT'] = df_gt['CLASS'].apply(lambda gt_class: GENOTYPE_TO_GT[gt_class])

    df_gt['GQ'] = df_gt.apply(
        lambda row: (
            int(10 * -math.log10(1 - row[row['CLASS']])) if row[row['CLASS']] < 1 else 255
        ) if row['CALLABLE'] else '.',
        axis=1
    )

    df_gt['GL'] = df_gt.apply(lambda row: '{HOM_REF:.4f},{HET:.4f},{HOM_ALT:.4f}'.format(**row), axis=1)

    # Get a series representing the column to be added to the VCF
    sample_column = df_gt.apply(lambda row: '{GT}:{GQ}:{GL}:{BP_REF_COUNT:.1f}:{BP_ALT_COUNT:.1f}'.format(**row), axis=1)
    sample_column.name = sample_name

    # Return
    return sample_column


def adjust_chrx_for_males(df_gt):
    """
    Adjust class densities on chrX variants to only hom-ref or hom-alt calls.

    :param df_gt: Genotype table.
    """

    for index in df_gt.index:
        if df_gt.loc[index, '#CHROM'] in {'chrX', 'X'}:
            df_gt.loc[index, 'HET'] = 0.0

            density_norm = sum(df_gt.loc[index, ['HOM_REF', 'HOM_ALT']])
            df_gt.loc[index, 'HOM_REF'] /= density_norm
            df_gt.loc[index, 'HOM_ALT'] /= density_norm


def adjust_chry(df_gt, is_female):
    """
    Adjust the chrY densities.

    :param df_gt: Genotype table.
    :param is_female: `True` if sample is female and chrY calls are always hom-ref, and `False` if samples are male
        or unknown and chrY calls are hom-ref or hom-alt (never het).
    """

    if is_female:
        for index in df_gt.index:
            if df_gt.loc[index, '#CHROM'] in {'chrY', 'Y'}:
                df_gt.loc[index, 'HOM_REF'] = 1.0
                df_gt.loc[index, 'HET'] = 0.0
                df_gt.loc[index, 'HOM_ALT'] = 0.0

    else:
        for index in df_gt.index:
            if df_gt.loc[index, '#CHROM'] in {'chrY', 'Y'}:
                df_gt.loc[index, 'HET'] = 0.0

                density_norm = sum(df_gt.loc[index, ['HOM_REF', 'HOM_ALT']])
                df_gt.loc[index, 'HOM_REF'] /= density_norm
                df_gt.loc[index, 'HOM_ALT'] /= density_norm


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
        df_var = pd.read_csv(in_file, sep='\t', header=None, comment='#')

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


def preprocess_manifest(manifest_file_name):
    """
    Pre-process the sample manifest and normalize all sex columns.

    :param manifest_file_name: File name to read.

    :return: A table with columns "SAMPLE", "SEX", and "DATA".
    """

    # Read manifest
    df = pd.read_csv(manifest_file_name, sep='\t', header=0)

    df.columns = [colname.upper() for colname in df.columns]

    for colname in ['SAMPLE', 'SEX', 'DATA']:
        if colname not in df.columns:
            raise RuntimeError(
                'Sample manifest "{}" does not contain a column with the label "{}"'.format(manifest_file_name, colname)
            )

    # Normalize sex
    df['SEX'] = df['SEX'].apply(lambda sex: 'U' if sex is None else sex.upper())

    if any(df['SEX'].apply(lambda sex: sex not in {'M', 'F', 'U'})):
        raise RuntimeError(
            'Sample manifest "{}" contains sexes that are not "M", "F", or "U"'.format(manifest_file_name)
        )

    # Check reference column
    if 'REF' not in df.columns:
        df['REF'] = 'NA'

    df['REF'] = df['REF'].fillna('NA')

    df['REF'] = df['REF'].apply(lambda val: val.upper() if val.upper() == 'NA' else val)

    # Order and set index
    df = df.loc[:, ['SAMPLE', 'SEX', 'DATA', 'REF']]
    df.set_index('SAMPLE', inplace=True)

    # Return
    return df
