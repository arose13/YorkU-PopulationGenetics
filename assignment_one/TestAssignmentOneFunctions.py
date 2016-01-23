import unittest
from assignment_one.AssignmentOneFunctions import *


# Since this is for marks it is imperative that we test the functions produce the correct results
class TestFunctions(unittest.TestCase):

    def test_correct_data(self):
        assert len(get_allozyme_data()) == 20

    def test_allele_freq_calculator(self):
        test_pop = get_allozyme_data()['EST1']
        output_freq_table, output_alleles, output_a = calc_allele_frequencies(test_pop)

        self.assertEqual(output_a, 2)
        self.assertEqual(output_freq_table['1'].values, 0.05)
        self.assertEqual(output_freq_table['2'].values, 0.95)

    def test_observed_heterozygosity_measures(self):
        test_pop = get_allozyme_data()['EST2']
        output_obs_hetero, output_hetero_genotypes = calc_observed_heterozygosity(test_pop)

        self.assertEqual(len(output_hetero_genotypes), 5)
        self.assertEqual(output_obs_hetero, 0.25)

    def test_expected_heterozygosity_measures(self):
        test_pop = get_allozyme_data()['EST1']
        output_freq_table, _, _ = calc_allele_frequencies(test_pop)
        output_expected_hetero = calc_expected_heterozygosity(output_freq_table)

        self.assertAlmostEqual(output_expected_hetero, 0.095, places=4)

    def test_f_is(self):
        test_pop = get_allozyme_data()['EST1']
        output_freq_table, _, _ = calc_allele_frequencies(test_pop)
        output_obs_hetero, _ = calc_observed_heterozygosity(test_pop)
        output_expected_hetero = calc_expected_heterozygosity(output_freq_table)
        output_f_is = calc_f_coef(output_obs_hetero, output_expected_hetero)

        self.assertAlmostEqual(output_f_is, 0.052631578, places=4)

    def test_chi_sq_calc(self):
        test_pop = get_allozyme_data()['EST1']
        output_freq_table, output_alleles, _ = calc_allele_frequencies(test_pop)
        output_chi_stat, output_df, output_p_value = chi_test(output_alleles, output_freq_table, test_pop, n=len(test_pop))

        self.assertAlmostEqual(output_chi_stat, 0.0554, places=4)

    def test_pairwise_mismatch(self):
        self.assertEqual(pairwise_mismatches('AACC')[0], 4)
        self.assertEqual(pairwise_mismatches('TTTA')[0], 3)
        self.assertEqual(pairwise_mismatches('TGGA')[0], 5)

    def test_pairwise_mismatch_num_combos(self):
        self.assertEqual(pairwise_mismatches('ABCD')[1], 6)

    def test_calc_nucleotide_diversity(self):
        seg_site_data = get_sequence_data()
        pi, mismatch_list = calc_nucleotide_diversity(seg_site=seg_site_data, genome_size=1000)

        self.assertEqual(pi, 0.0068)
