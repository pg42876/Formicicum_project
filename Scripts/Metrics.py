import pandas as pd
from Scripts.utils import ReconstructionTool, get_reactions, get_genes, get_metabolites, calculate_quality_metrics

with open(f'../Results/formicicum_reactions.txt', 'r') as f:

    formicicum = f.readlines()

    print(formicicum)
    print(formicicum[0])

    for reaction in formicicum:
        print(reaction)

print(formicicum)

metrics = []
for tool in ['merlin_blast', 'merlin_bit', 'carveme', 'kbase']:
    tool_reactions = open(f'../Results/kegg_functional_ids.tsv', 'r')
    #reactions = tool_reactions['reactions']
    TPs = len([ide for ide in tool_reactions if ide in formicicum])
    FPs = len([ide for ide in tool_reactions if ide not in formicicum])
    FNs = len([ide for ide in formicicum if ide not in tool_reactions])
    print(TPs, FNs, FPs)
    metrics.append([tool, TPs, FPs, FNs] + list(calculate_quality_metrics(TPs, FPs, FNs)))
pd.DataFrame(metrics, columns=['tool', 'TPs', 'FPs', 'FNs', 'Precision', 'Recall', 'F1 score', 'Jaccard distance']).to_excel('../Results/quality_metrics_reactions.xlsx', index=False)
