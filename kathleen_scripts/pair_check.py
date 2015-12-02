# check whether all disease pairs are represented

# import modules
import pandas as pd
import os

# import original disease pairs dataset
path_to_orig = '/Users/kathleen/Box Sync/UCSF/Fall Classes 2015/BMI_206/Project/DiseaseDiseaseBMI206/data/DataS4_disease_pairs.tsv'
orig_pairs_df = pd.read_csv(path_to_orig, sep='\t', header = 0, comment='#')

# import new disease pairs dataset
path_to_new = '/Users/kathleen/Box Sync/UCSF/Fall Classes 2015/BMI_206/Project/DiseaseDiseaseBMI206/bmi/output_distances.tsv'
new_pairs_df = pd.read_csv(path_to_new, sep='\t', header = None, comment='#')

# check whether each value from the disease A column of the
# original pairs set is in the first column of the new set

# convert cols of interest into lists
new_col_1 = list(new_pairs_df[0])
new_col_1_rev = [w.replace('.txt', '') for w in new_col_1]
old_col_1 = list(orig_pairs_df['disease_A'])

# # compare
# for val in old_col_1:
#     if val in new_col_1_rev:
#         continue
#     else:
#         print val

# check whether each value from the disease B column of the
# original pairs set is in the second column of the new set

# convert cols of interest into lists
new_col_2 = list(new_pairs_df[1])
new_col_2_rev = [w.replace('.txt', '') for w in new_col_1]
old_col_2 = list(orig_pairs_df['disease_B'])

# compare
for val in old_col_2:
    if val in new_col_2_rev:
        continue
    else:
        print val

# for some reason this has a hard time with "intestinal diseases" even though it's definitely there


# this doesn't work because disease names are different

