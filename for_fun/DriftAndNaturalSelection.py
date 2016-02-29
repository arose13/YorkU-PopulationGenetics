# Methods for the Beneficial Allele Fixation Notebook
import warnings
import numpy as np
import matplotlib.pyplot as graph
import seaborn as sns

warnings.simplefilter('ignore')
sns.set(font_scale=1.5, rc={'figure.figsize': (14, 9)})
sns.set_style(style='whitegrid')

DOMINANT = 'Dominant'
RECESSIVE = 'Recessive'
ADDITIVE = 'Additive'


class NSWithDriftSimulation(object):
    """
    Monte Carlo Natural Selection Genetic Drift Simulation to figure out the odds of a beneficial allele being fixed.

    NOTE:
        - Drift occurs before selection occurs.
        - Population size is 1/2 of the allele pool size
    """

    def __init__(self, size_p=20, size_q=20, p_rf=1.0, q_rf=0.9, type_of_trait=DOMINANT):
        """
        :param type_of_trait: of allele p
        :return:
        """
        # Hyperparameters
        self.max_generations = 2000

        self.type_of_trait = type_of_trait
        self.size_p = size_p
        self.size_q = size_q
        self.allele_pool_size = size_p + size_q

        self.pp_rf, self.pq_rf, self.qq_rf = self.calc_relative_fitness(p_rf, q_rf, type_of_trait)

    def run(self, number_of_runs=1, max_generations=None):
        # Some Stats
        number_lost = 0
        num_of_generations = []

        if max_generations is not None:
            self.max_generations = max_generations

        graph.ylabel('Proportion of p')
        graph.xlabel('Number of Generations')

        # Selection With Drift
        for i in range(number_of_runs):
            results, lost = self.single_run()
            num_of_generations.append(len(results))
            number_lost += lost
            graph.plot(np.array(results) / self.allele_pool_size, linewidth=1.0, alpha=0.4)

        num_of_generations = np.array(num_of_generations)

        graph.plot(
            [0, num_of_generations.max()],
            [1, 1],
            linewidth=1.0,
            linestyle='--',
            label='Fixation Point',
            color='k'
        )

        # Selection Alone (Benchmark)
        graph.plot(
            np.array(self.selection_only_benchmark(break_point=num_of_generations.max())) / self.allele_pool_size,
            linewidth=1.5,
            color='k',
            label='Selection Only (Non-Discrete)'
        )

        graph.title(
            '{} Simulation ({}% of runs lost allele p)'.format(
                self.type_of_trait,
                round(number_lost * 100 / number_of_runs, 2)
            )
        )
        graph.legend(loc=4)
        graph.xlim(0, num_of_generations.max())
        graph.ylim(0, 1.05)
        graph.show()

        return num_of_generations, number_lost

    def single_run(self):
        p, q = self.size_p, self.size_q

        i = 0
        lost = 0
        run_results = [p]

        while 0 < p < self.allele_pool_size:
            # Continue as long as the population has not become fixed.
            p, q = self.i_drift(p, q)
            p, q = self.i_selection(p, q)

            run_results.append(p)

            if p == 0:
                lost = 1

            i += 1  # Loop counter for potential infinite termination
            if i == self.max_generations:
                break

        return run_results, lost

    def selection_only_benchmark(self, break_point):
        i = 0
        p, q = self.size_p, self.size_q
        benchmark_results = [p]

        while 0 < p < self.allele_pool_size:
            p, q = self.i_selection(p, q, discretization=False)
            benchmark_results.append(p)

            i += 1
            if i == break_point:
                break

        return benchmark_results

    def calc_genotypic_frequencies(self, size_p, size_q):
        pp = (size_p / self.allele_pool_size) ** 2
        pq = (size_p / self.allele_pool_size) * (size_q / self.allele_pool_size) * 2
        qq = (size_q / self.allele_pool_size) ** 2

        return pp, pq, qq

    @staticmethod
    def calc_relative_fitness(p_rf, q_rf, type_of_trait):
        pp_rf = p_rf
        qq_rf = q_rf

        if type_of_trait is DOMINANT:
            pq_rf = p_rf
        elif type_of_trait is RECESSIVE:
            pq_rf = q_rf
        elif type_of_trait is ADDITIVE:
            pq_rf = (p_rf + q_rf) / 2.0
        else:
            # Assume Dominant
            pq_rf = p_rf

        return pp_rf, pq_rf, qq_rf

    @staticmethod
    def i_drift(p, q):
        n = p + q
        p_prime = np.random.binomial(n, p / n, 1)[0]
        q_prime = n - p_prime

        return p_prime, q_prime

    def i_selection(self, p, q, discretization=True):
        # Calculate the genotypic frequencies (HWE)
        pp, pq, qq = self.calc_genotypic_frequencies(p, q)
        genotypic_freq = np.array([pp, pq, qq])

        # Calculate the frequencies of survivors
        genotypic_freq = genotypic_freq * np.array([self.pp_rf, self.pq_rf, self.qq_rf])
        mean_relative_fitness = np.sum(genotypic_freq)

        # Calculate new frequencies after selection
        genotypic_freq = genotypic_freq / mean_relative_fitness

        # Convert frequencies to numbers and resolve allele numbers
        population_genotypes = genotypic_freq * (self.allele_pool_size / 2)

        if discretization:
            population_genotypes = np.around(population_genotypes)

        [pp_prime, pq_prime, qq_prime] = population_genotypes.tolist()

        p = (pp_prime * 2) + pq_prime
        q = (qq_prime * 2) + pq_prime
        return p, q


if __name__ == '__main__':
    simulation = NSWithDriftSimulation(size_p=1000, size_q=100000, type_of_trait=DOMINANT)
    simulation.run(number_of_runs=1000)
