import numpy as np
import pandas as pd
from pandas.stats.api import ols
from sklearn.ensemble import RandomForestRegressor

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


def select_by_markers(markers: list, qtl_data=get_qtl_data()):
    data = pd.DataFrame(columns=markers)
    data['phenotype'] = qtl_data['Trait_1']

    for marker in markers:
        data[marker] = qtl_data[marker]

    return data


def assessing_marker_importance(x_train, y_train, labels):
    rf = RandomForestRegressor(
        n_estimators=5000,
        criterion='mse',
        random_state=1992,
        n_jobs=-1
    )
    rf.fit(x_train, y_train)


    importance_matrix = rf.feature_importances_
    indices = np.argsort(importance_matrix)[::-1]

    output_array = []
    data = pd.DataFrame(columns=['Marker', 'Importance'])
    for feature in range(x_train.shape[1]):
        output_array.append('{}) **{}** {}%'.format(
            feature + 1,
            labels[feature],
            round(importance_matrix[indices[feature]] * 100, 3)
        ))
        data.loc[feature] = [labels[feature], importance_matrix[indices[feature]]]

    return output_array, data
