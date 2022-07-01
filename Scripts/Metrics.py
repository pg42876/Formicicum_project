import pandas as pd

carveme = pd.read_csv('../Results/carveme_uniprotinfo.tsv', sep ='\t')
carveme_genes = carveme['Entry'].tolist()
print(carveme_genes)

formicicum = pd.read_csv('../Results/formicicum_uniprotinfo.tsv', sep ='\t')
formicicum_genes = formicicum['Entry'].tolist()
print(formicicum_genes)

TPs = [ide for ide in carveme_genes if ide in formicicum_genes]
FPs = [ide for ide in carveme_genes if ide not in formicicum_genes]
FNs = [ide for ide in formicicum_genes if ide not in carveme_genes]

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

CV = calculate_quality_metrics(len(TPs), len(FPs), len(FNs))