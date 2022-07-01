import pandas as pd

merlin_blast = pd.read_csv('../Results/merlin_blast_uniprotinfo.tsv', sep ='\t')
merlin_genes = merlin_blast['Entry'].tolist()
print(merlin_genes)

formicicum = pd.read_csv('../Results/formicicum_uniprotinfo.tsv', sep ='\t')
formicicum_genes = formicicum['Entry'].tolist()
print(formicicum_genes)

TPs = [ide for ide in merlin_genes if ide in formicicum_genes]
FPs = [ide for ide in merlin_genes if ide not in formicicum_genes]
FNs = [ide for ide in formicicum_genes if ide not in merlin_genes]

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

MRL = calculate_quality_metrics(len(TPs), len(FPs), len(FNs))