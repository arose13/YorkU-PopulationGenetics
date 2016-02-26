from itertools import combinations


site = 'ANTHONY'
mismatches = 0
combos = 0

for comparison in combinations(site, 2):
    print(comparison)
    combos += 1
    if comparison[0] != comparison[1]:
        mismatches += 1

print('Number of mismatches: {}'.format(mismatches))
print('Number of combos: {}'.format(combos))
