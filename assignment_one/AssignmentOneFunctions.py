# Functions required for assignment one
import numpy as np
import pandas as pd
from scipy.stats import chi2
from collections import Counter
from itertools import combinations_with_replacement, combinations


def get_allozyme_data(file='data/MyAllozymeDataset.csv'):
    data = pd.read_csv(file)
    return data


def get_sequence_data(file='data/MySequenceDiversityData.csv'):
    data = pd.read_csv(file, skiprows=3)
    print(data)
    return data


def calc_allele_frequencies(est):
    # Give me only 1 column
    gene_pool = ''

    # Solve Gene Pool!
    for dude in est:
        gene_pool += str(dude)

    # Find the number of unique alleles
    alleles = list(set(gene_pool))
    print('Alleles in Population: {}'.format(alleles))
    print('Number of Alleles per Locus: {}'.format(len(alleles)))

    # Count the frequency of each unique allele
    freq_table = pd.DataFrame()
    for p in alleles:
        p = str(p)
        print('Freq p{} is {}'.format(
                p,
                str(gene_pool.count(p) / len(gene_pool))
        ))
        freq_table[p] = [gene_pool.count(p) / len(gene_pool)]

    return freq_table, alleles, len(alleles)  # Allele Frequency Table, Alleles, Alleles per Locus


def calc_observed_heterozygosity(est):
    # Observed Hetero = Heteros / N
    num_hetero = 0
    hetero_genotypes = []

    for genotype in est:
        # Listen to this logic bro. If the character in the first position is not the same as the character in the
        # second position then it most be hetero!
        i = str(genotype)
        if i[0] != i[1]:
            num_hetero += 1  # Hetero
            hetero_genotypes.append(i)

    observed_heterozygosity = num_hetero / len(est)
    print('Observed Heterozygosity: {}'.format(observed_heterozygosity))

    return observed_heterozygosity, hetero_genotypes


def calc_expected_heterozygosity(freq_table):
    # Expected Hetero = 1 - Sum(proportion of allele ^ 2)
    expected_heterozygosity = 1 - (freq_table**2).sum(axis=1).values[0]
    print('Expected Heterozygosity: {}'.format(expected_heterozygosity))

    return expected_heterozygosity


def calc_f_coef(h_observed, h_expected):
    f_is = abs(h_expected - h_observed) / h_expected
    print('F_is: {}'.format(f_is))
    return f_is


def chi_test(alleles, freq_table, est, n=1):
    # Combinations: https://docs.python.org/2/library/itertools.html
    # Calculate Observed and the Expected then calc the chi stat
    allele_pool = ''.join(alleles)
    individuals = Counter(est)
    observed_table = pd.DataFrame()
    expected_table = pd.DataFrame()

    for genotype in combinations_with_replacement(allele_pool, 2):
        genotype_string = genotype[0] + genotype[1]

        # EXPECTED
        if genotype[0] != genotype[1]:
            # These are the Heteros, multiple their expected frequencies by 2
            exp_hetero_count = 2 * freq_table[genotype[0]].values * freq_table[genotype[1]].values * n
            expected_table[genotype_string] = exp_hetero_count

        else:
            # These are the Homos, square their expected frequencies
            exp_homo_count = (freq_table[genotype[0]].values ** 2) * n
            expected_table[genotype_string] = exp_homo_count

        # OBSERVED
        if genotype[0] == genotype[1]:
            # Homo therefore order doesn't matter
            obs_homo_count = individuals[int(genotype_string)]
            observed_table.set_value(0, genotype_string, obs_homo_count)
        else:
            # Hetero therefore order does matter
            # I'm using extended slice syntax to reverse the string [::-1] This makes synonymous genotypes the same
            obs_hetero_count = individuals[int(genotype_string)]
            obs_hetero_count += individuals[int(genotype_string[::-1])]
            observed_table.set_value(0, genotype_string, obs_hetero_count)

    # Calculate Chi sq
    chi_table = ((observed_table - expected_table) ** 2) / expected_table
    chi_sq_statistic = chi_table.sum(axis=1).values[0]
    df = len(expected_table.columns) - 2
    p_value = chi2.sf(chi_sq_statistic, df)

    print('\n Expected Numbers')
    print(expected_table)

    print('\n Observed Numbers')
    print(observed_table)

    print('\n CHI^2: {}'.format(chi_sq_statistic))
    print('df: {}, p: {}'.format(df, p_value))

    return chi_sq_statistic, df, p_value


def calc_segregating_sites(sequence_data, sequence_length=1000):
    s = len(sequence_data) / sequence_length
    return s


def pairwise_mismatches(site):
    # Segregating site pairwise mismatches
    mismatches = 0
    combos = 0

    for pair in combinations(site, 2):
        combos += 1
        if pair[0] != pair[1]:
            mismatches += 1

    return mismatches, combos


def calc_nucleotide_diversity(seg_site, genome_size):
    # Calculate Nucleotide Diversity.
    # pi = sum(number of mismatches) / sum(number of combos) OR sum(num of mismatches) / (num of combos * num of sites)

    site_mismatches = []
    total_mismatch = 0
    site_combos = 0

    for pos in seg_site:
        site = ''.join(np.array(seg_site[pos]))
        i_mismatches, i_combos = pairwise_mismatches(site)
        total_mismatch += i_mismatches
        site_combos = i_combos
        site_mismatches.append(i_mismatches)

    print('\nSite Mismatches \n{}'.format(site_mismatches))

    pi = total_mismatch / (site_combos * genome_size)
    print('Pi: {}'.format(pi))

    return pi, site_mismatches


def run_analysis_part_1(file='data/MyAllozymeDataset.csv'):
    allozyme_data = get_allozyme_data(file)
    print(allozyme_data.columns)

    for est in allozyme_data.columns:
        print('-------------------------------------------')
        print(est)

        # Allele Frequency and Number of Alleles / Locus
        allele_freq_table, alleles, _ = calc_allele_frequencies(allozyme_data[est])

        # Heterozygosity / Locus
        h_i, _ = calc_observed_heterozygosity(allozyme_data[est])
        h_s = calc_expected_heterozygosity(allele_freq_table)

        # Inbreed Coef
        calc_f_coef(h_i, h_s)

        # Count things
        chi_test(alleles, allele_freq_table, allozyme_data[est], len(allozyme_data))

        # Just a space
        print()


def run_analysis_part_2(genome_size=1000, file='data/MySequenceDiversityData.csv'):
    sequence_data = get_sequence_data(file)

    # Calc Segregating Site 'S'
    s = calc_segregating_sites(sequence_data)

    # Calc Mismatches at each seg site in the same order as the original dataset
    # Nucleotide Diversity 'pi'
    pi, site_mismatches = calc_nucleotide_diversity(sequence_data, genome_size)


if __name__ == '__main__':
    # Test Functions to make sure nothing was broken.
    run_analysis_part_1()
    run_analysis_part_2()
