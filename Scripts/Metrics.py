import pandas as pd
from Scripts.utils import ReconstructionTool, get_reactions, get_genes, get_metabolites, calculate_quality_metrics

formicicum_metabolites = pd.read_csv('../Results/formicicum_metabolites.txt', sep='\n', header=None)[0].tolist()
func_ids_df = pd.read_csv('../Results/metanetx_functional_ids.tsv', sep='\t', index_col=0)

conversion = pd.read_csv('../Scripts/Xrefs files/compounds-conversion.csv')

conversion['External ID'].isin(formicicum_metabolites).sum() #Verdadeiros

reactions = conversion[conversion['External ID'].isin(formicicum_reactions)]['Internal ID'].tolist()
#print(reactions)

metrics = []
for tool in ['aureme', 'merlin_blast', 'merlin_bit', 'carveme', 'kbase']:
    func_reactions = func_ids_df['metabolites'].tolist()
    lista = func_reactions[0].split(',')
    print(lista)
    TPs = len([ide for ide in lista if ide in reactions])
    FPs = len([ide for ide in lista if ide not in reactions])
    FNs = len([ide for ide in reactions if ide not in lista])
    #print(TPs, FNs, FPs)
    metrics.append([tool, TPs, FPs, FNs] + list(calculate_quality_metrics(TPs, FPs, FNs)))
pd.DataFrame(metrics, columns=['tool', 'TPs', 'FPs', 'FNs', 'Precision', 'Recall', 'F1 score', 'Jaccard distance']).to_excel(
    '../Results/quality_metrics_compounds.xlsx', index=False)