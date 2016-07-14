"""
Genotype a set of SV calls in a given set of BAMs.
"""
import argparse
import numpy as np
import pandas as pd
from scipy.stats import binom
np.random.seed(1)

# Input file looks like this with SV coordinates in local assembly at columns 6-8.
# chr1    350793  350794  insertion       40      chr1-350756-362784|utg7180000000000|merged      37      77
CHROMOSOME=0
START=1
END=2
EVENT_TYPE=3
EVENT_LENGTH=4
CONTIG_NAME=5
CONTIG_START=6
CONTIG_END=7


def calculate_genotype_likelihoods(reference_support, alternate_support, chromosomes, sample_sexes, homozygous_probability, heterozygous_probability):
    total_reads = reference_support + alternate_support
    gl_hom_ref = binom.pmf(n=total_reads, k=reference_support, p=homozygous_probability)
    gl_het = binom.pmf(n=total_reads, k=reference_support, p=heterozygous_probability)
    gl_hom_alt = binom.pmf(n=total_reads, k=reference_support, p=1 - homozygous_probability)

    # Set probabilities for heterozygous and homozygous calls in hemizygous
    # chromosomes to 0 depending on the sex of the sample.
    gl_het[np.array((chromosomes == "chrX") & (sample_sexes == "M"), dtype=bool)] = 0
    gl_het[np.array((chromosomes == "chrY"), dtype=bool)] = 0
    gl_hom_ref[np.array((chromosomes == "chrY") & (sample_sexes == "F"), dtype=bool)] = 0
    gl_hom_alt[np.array((chromosomes == "chrY") & (sample_sexes == "F"), dtype=bool)] = 0

    return np.array((gl_hom_ref, gl_het, gl_hom_alt)).T


def genotype(df, output, samples, min_depth, homozygous_binomial_probability, heterozygous_binomial_probability):
    # Define genotype classes.
    genotype_classes = ("0/0", "1/0", "1/1")

    # Calculate binomial genotype likelihoods.
    binomial_genotype_probabilities = calculate_genotype_likelihoods(
        df["discordant"].astype(int),
        df["concordant"].astype(int),
        df["chr"],
        df["sex"],
        homozygous_binomial_probability,
        heterozygous_binomial_probability
    )
    df["genotype"] = np.array(genotype_classes)[np.argmax(binomial_genotype_probabilities, axis=1)]

    # Scale binomial probabilities and convert to Phred score space.
    binomial_genotype_probabilities_sum = binomial_genotype_probabilities.sum(axis=1)
    binomial_genotype_likelihoods = -10 * np.log10(1 - (binomial_genotype_probabilities / binomial_genotype_probabilities_sum[:, None]))

    # Set missing and infinite likelihoods to reasonable values.
    binomial_genotype_likelihoods[np.isnan(binomial_genotype_likelihoods)] = 0
    max_non_infinite_likelihood = np.ceil(binomial_genotype_likelihoods[~np.isinf(binomial_genotype_likelihoods)].max())
    binomial_genotype_likelihoods[np.isinf(binomial_genotype_likelihoods)] = max_non_infinite_likelihood

    # Convert likelihoods to integers.
    binomial_genotype_likelihoods = binomial_genotype_likelihoods.astype(int)

    # Convert tuple of likelihoods per site to a comma-delimited list for VCF.
    df["genotype_likelihoods"] = [",".join(row) for row in binomial_genotype_likelihoods.astype(str)]

    # Save the likelihood of the final genotype as the genotype quality.
    df["genotype_quality"] = np.max(binomial_genotype_likelihoods, axis=1)

    # Convert genotypes for low coverage regions to unknowns and set likelihoods
    # to 0.
    df.loc[(df["concordant"] < min_depth) & (df["discordant"] < min_depth), "genotype"] = "./."
    df.loc[(df["concordant"] < min_depth) & (df["discordant"] < min_depth), "genotype_quality"] = 0
    df.loc[(df["concordant"] < min_depth) & (df["discordant"] < min_depth), "genotype_likelihoods"] = ".,.,."

    # Save genotypes.
    df.to_csv(output, sep="\t", index=False, header=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("genotypes", help="BED file of SV calls with reference coordinates in columns 1-3, SV type in 4, length in 5, and contig coordinates in 6-8")
    parser.add_argument("samples", help="tab-delimited list of samples with column names 'sample' and 'sex' where sex is either 'M' or 'F'")
    parser.add_argument("output", help="file to emit VCF output to; use /dev/stdout for piping")
    parser.add_argument("--min_depth", type=int, default=5, help="minimum depth required across at least one haplotype to genotype a region")
    parser.add_argument("--homozygous_binomial_probability", type=float, default=0.95, help="binomial probability for a homozygous genotype")
    parser.add_argument("--heterozygous_binomial_probability", type=float, default=0.4, help="binomial probability for a heterozygous genotype")
    args = parser.parse_args()

    complete_df = pd.read_table(args.genotypes)
    samples_df = pd.read_table(args.samples)

    # Annotate sex from sample manifest to genotype input if necessary.
    if "sex" not in complete_df.columns:
        complete_df = complete_df.merge(samples_df, on="sample", how="left")

    with open(args.output, "w") as output:
        genotype(
            complete_df,
            output,
            samples_df,
            args.min_depth,
            args.homozygous_binomial_probability,
            args.heterozygous_binomial_probability
        )
