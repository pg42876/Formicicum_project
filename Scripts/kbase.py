import pandas as pd

import Scripts
from Scripts import Metrics, merlin_blast, merlin_bit

## KBase

kbase = pd.read_csv('../Results/kbase_uniprotinfo.tsv', sep ='\t')
kbase_genes = kbase['Entry'].tolist()
print(kbase_genes)

formicicum = pd.read_csv('../Results/formicicum_uniprotinfo.tsv', sep ='\t')
formicicum_genes = formicicum['Entry'].tolist()
print(formicicum_genes)

TPs = [ide for ide in kbase_genes if ide in formicicum_genes]
FPs = [ide for ide in kbase_genes if ide not in formicicum_genes]
FNs = [ide for ide in formicicum_genes if ide not in kbase_genes]

print(len(TPs))
print(len(FPs))
print(len(FNs))

def calculate_quality_metrics(tp, fp, fn):
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1_score = 2 * (precision * recall) / (precision + recall)
    jaccard_distance = 1 - ( tp / fp + fp + fn)
    return precision, recall, f1_score, jaccard_distance

print(calculate_quality_metrics(len(TPs), len(FPs), len(FNs)))


## Create DataFrame

KB = calculate_quality_metrics(len(TPs), len(FPs), len(FNs))
CV = Scripts.Metrics.CV
BIT = Scripts.merlin_bit.BIT
MRL = Scripts.merlin_blast.MRL

precision, recall, f1_score, jaccard_distance =

DataFrame = pd.DataFrame({'KBase' : [KB],
                          'CarveMe' : [CV],
                          'merlin' : [MRL],
                          'merlin(bit)' : [BIT]},
                         columns = ['KBase', 'CarveMe', 'merlin', 'merlin(bit)'])
print(DataFrame)


