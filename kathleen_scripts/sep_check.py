# check whether <sAB> are equal for each pair

# import modules
import pandas as pd
import os

# import original disease pairs dataset
path_to_orig = '/Users/kathleen/Box Sync/UCSF/Fall Classes 2015/BMI_206/Project/DiseaseDiseaseBMI206/data/DataS4_disease_pairs.tsv'
orig_pairs_df = pd.read_csv(path_to_orig, sep='\t', header = 0, comment='#')

# import new disease pairs dataset
path_to_new = '/Users/kathleen/Box Sync/UCSF/Fall Classes 2015/BMI_206/Project/DiseaseDiseaseBMI206/reanalysis_output/output_distances.tsv'
new_pairs_df = pd.read_csv(path_to_new, sep='\t', header = None, comment='#')
# replace '.txt' from each disease name
new = new_pairs_df[0].str.replace('.txt','')
del new_pairs_df[0]
new_pairs_df[0] = new
new1 = new_pairs_df[1].str.replace('.txt','')
del new_pairs_df[1]
new_pairs_df[1] = new1
new_pairs_df.columns = ['d_AB', 's_AB', 'dis_1', 'dis_2']
# new_pairs_df.to_csv('new_s_output.csv', sep = '\t', header=True, index = False)

# find pairs and compare their sABs
# df1 = orig_pairs_df.loc[orig_pairs_df['disease_A'] == 'abnormalities, multiple']
# df2 = df1.loc[df1['disease_B'] == 'actinomycetales infections']
# print df2

for i, r in orig_pairs_df.iterrows():
    df1 = new_pairs_df.loc[new_pairs_df['dis_1'] == (r['disease_A'] or r['disease_B'])]
    df2 = df1.loc[df1['dis_2'] == (r['disease_B'] or r['disease_A'])]
    s_AB_dif = df2['s_AB'] - r['s_AB (observed)']
    with open('output_s_AB.tsv', 'ab') as f: 
        row = str(r['disease_A']) + '\t' + str(r['disease_B']) + '\t' + str(s_AB_dif) + '\n'
        f.write((row))



