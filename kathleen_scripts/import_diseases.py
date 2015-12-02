# import disease genes list and create .txt files with list of genes for each disease

# import modules
import sys, getopt, csv
import argparse
import pandas as pd
from itertools import chain
import operator

# set up parsing input and output files
 
parser = argparse.ArgumentParser(description='Parse input and output files')
parser.add_argument('-i', help='Input (should be genes list file from paper',required=True)
parser.add_argument('-od', help='Output directory', required=True, action='store')
args = parser.parse_args()

# show values #

print ("Genes file: %s" % args.i)
print ("Output directory: %s" % args.od)

# read original genes file into Pandas dataframe

genes_df = pd.read_csv(args.i, sep='\t', header=0, comment='#')
genes_df['GWAS_genes'].fillna(str(0), inplace=True)

# separate diseases into their own dataframes and read these to 
# their own .csv files (with 1 column each)

result = genes_df.groupby('disease')[['OMIM_genes','GWAS_genes']].apply(lambda x: x.values.flatten().tolist()).to_dict()
# print result
for key, values in result.iteritems():
    filename = str(key)
    # print values
    values_list = reduce(operator.add, values).split(';')
    # print values_list
    with open(args.od + filename + '.txt', 'w') as infile:
        infile.write('\n'.join(values_list))



    # print gwas_df.dtype
    # new_df = pd.concat([omim_df, gwas_df], axis = 1)
    # new_df.to_csv(i['disease'])
    # print i['OMIM_genes']




#     omim_df = genes_df['OMIM_genes']
# omim_df.tolist()
# chain = list(itertools.chain.from_iterable(omim_df)) 
# gwas_df = genes_df['GWAS_genes']

# print omim_df