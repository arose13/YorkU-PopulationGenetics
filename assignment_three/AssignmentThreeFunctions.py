# Part 1
import pandas as pd

# Trait Keys
[parent1, parent2, offspring] = ['parent1', 'parent2', 'offspring']


def get_data(file: str = 'data/anthony_parent_offspring.csv') -> pd.DataFrame:
    data = pd.read_csv(file)
    return data


def init_trait_data(trait: str, data=get_data()):
    trait_data = pd.DataFrame(columns=[parent1, parent2, offspring])

    trait_data[parent1] = data['p1_' + trait]
    trait_data[parent2] = data['p2_' + trait]
    trait_data[offspring] = data['o_' + trait]

    trait_data['parent_mid'] = (trait_data[parent1] + trait_data[parent2]) / 2

    return trait_data
