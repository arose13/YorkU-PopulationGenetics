import pandas as pd
from pandas.stats.api import ols

# Trait Keys
[parent1, parent2, offspring] = ['parent1', 'parent2', 'offspring']


def get_parent_offspring_data(file: str = 'data/anthony_parent_offspring.csv') -> pd.DataFrame:
    data = pd.read_csv(file)
    return data


def get_qtl_data(file: str = 'data/anthony_qtl_31.csv') -> pd.DataFrame:
    data = pd.read_csv(file)
    return data


def init_trait_data(trait: str, data=get_parent_offspring_data()) -> pd.DataFrame:
    trait_data = pd.DataFrame(columns=[parent1, parent2, offspring])

    trait_data[parent1] = data['p1_' + trait]
    trait_data[parent2] = data['p2_' + trait]
    trait_data[offspring] = data['o_' + trait]

    trait_data['parent_mid'] = (trait_data[parent1] + trait_data[parent2]) / 2

    return trait_data


def init_qtl_marker_data(chromosome: str, data=get_qtl_data()) -> pd.DataFrame:
    summary_stats = pd.DataFrame(columns=['Position', 'Slope', 'p-value'])

    all_markers = data.columns.values.tolist()

    i = 0
    for i_marker in all_markers:
        if chromosome in i_marker:
            # Calculate Stats
            position_data = pd.DataFrame(columns=['genotype', 'phenotype'])
            position_data['genotype'] = data[i_marker]
            position_data['phenotype'] = data['Trait_1']

            i_results = ols(x=position_data['genotype'], y=position_data['phenotype'])
            summary_stats.loc[i] = [i*10, i_results.beta['x'], i_results.p_value['x']]

            # Iterate
            i += 1

    return summary_stats
