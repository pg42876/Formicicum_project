import pandas as pd

## Genes

with open('../Results/aureme_genes.txt', 'r') as a:
    aureme = a.readlines()

with open("../Results/carveme_genes.txt", "r") as b:
    carveme = b.readlines()

with open("../Results/kbase_genes.txt", "r") as c:
    kbase = c.readlines()

with open("../Results/merlin_blast_genes.txt", "r") as d:
    merlin_blast = d.readlines()

with open("../Results/merlin_bit_genes.txt", "r") as e:
    merlin_bit = e.readlines()

with open("../Results/functional_ids.tsv", "r") as tsv:
    total_genes = [line for line in tsv]

genes = aureme + carveme + merlin_bit + merlin_blast + kbase

output_genes = []

for line in total_genes:
    for gene in genes:
        if gene in line:
            if line.endswith('\n'):
                output_genes.append(line)
            else:
                output_genes.append('{0}\n'.format(line))
            print('Found line {0} that matched word {1}'.format(line, gene))

with open('../Results/output_genes.tsv', 'w') as output_file:
    output_file.writelines(output_genes)

output_genes = pd.read_csv('../Results/output_genes.tsv', sep = '\t').T
print(output_genes)
#