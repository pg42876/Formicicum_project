import pandas as pd

kbase_genes = pd.read_csv('../Results/kbase_uniprotinfo.tsv', sep='\t')

kbase_genes['Entry'].tolist()




