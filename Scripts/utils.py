import re
from enum import Enum
import cobra.io.sbml

class Utils:

    @staticmethod
    def get_metabolite_ids(model_path):
        model = cobra.io.sbml.read_sbml_model(model_path)
        metabolite_list = []

        for metabolite in model.metabolites:
            metabolite_id = re.sub("_[a-z]*$", "", metabolite.id)
            metabolite_list.append(metabolite_id)

        metabolite_set = set(metabolite_list)

        with open("metabolites_list_%s.txt" % model.name, "w") as file:
            for met_id in metabolite_set:
                file.write(met_id + "\n")

    @staticmethod
    def jaccard_distance(generated_elements, reference_elements):

        jd = 0

        if type(generated_elements) != set:
            raise TypeError("generated_elements should be a set instead of " + str(type(generated_elements)))

        elif type(reference_elements) != set:
            raise TypeError("reference_elements should be a set instead of " + str(type(reference_elements)))

        else:
            intersection = generated_elements & reference_elements
            union = generated_elements | reference_elements
            jd = 1 - (abs(len(intersection))/abs(len(union)))

        return jd

    @staticmethod
    def ratio(generated_elements, reference_elements):

        ratio = 0

        if type(generated_elements) != set:
            raise TypeError("generated_elements should be a set instead of " + str(type(generated_elements)))

        elif type(reference_elements) != set:
            raise TypeError("reference_elements should be a set instead of " + str(type(reference_elements)))

        else:

            intersection = generated_elements & reference_elements
            difference = generated_elements - reference_elements
            ratio = abs(len(intersection)) / abs(len(difference))

        return ratio


class Type(Enum):
    GENES = 1
    REACTIONS = 2
    METABOLITES = 3

class ReconstructionTool(Enum):
    AUTOKEGGREC = 1
    AUREME = 2
    CARVEME = 3
    MERLIN = 4
    MERLIN_BIT = 5
    MODELSEED = 6
    PATHWAYTOOLS = 7
    RAVEN = 8
    T_GONDII_CURATED = 9


def get_genes(xml, tool=None):
    if tool == 'kbase':
        return [gene.id.replace('_RS', '_') for gene in xml.genes]
    if tool == 'carveme':
        parts = [gene.id.split('_') for gene in xml.genes]
        return [f'{part[0]}|{part[1]}|{part[2]}_{part[3]}' for part in parts]
    return [gene.id for gene in xml.genes]


def get_reactions(xml, tool=None):
    if tool == 'merlin_blast':
        return [reaction.id.split('__')[0] for reaction in xml.reactions]
    if tool == 'kbase':
        return [reaction.annotation['seed.reaction'] for reaction in xml.reactions
                if 'seed.reaction' in reaction.annotation.keys()]
    if tool == 'aureme':
        return [reaction.id.split('__')[0] for reaction in xml.reactions]
    return [reaction.id for reaction in xml.reactions]


def get_cross_reference_reactions(xml, conversion_df, tool=None):
    reactions = get_reactions(xml, tool=tool)
    metanetx_reactions = conversion_df[(conversion_df['External ID'].isin(reactions))]['Internal ID'].tolist()
    kegg_reactions = conversion_df[(conversion_df['Internal ID'].isin(metanetx_reactions)) &
                                   (conversion_df['Source'] == 'KEGG')]['External ID'].tolist()
    print(f'Found {len(metanetx_reactions)} reactions.')
    return metanetx_reactions, kegg_reactions


def get_metabolites(xml, tool=None):
    if tool == 'carveme':
        return [metabolite.id.replace('__', '#').split('_')[0].replace('#', '__') for metabolite in xml.metabolites]
    if tool == 'kbase':
        return [metabolite.id.split('_')[0] for metabolite in xml.metabolites]
    if tool in ['aureme', 'merlin_blast']:
        return [metabolite.id.split('__')[0] for metabolite in xml.metabolites]
    return [metabolite.id for metabolite in xml.metabolites]


def get_cross_reference_metabolites(xml, conversion_df, tool=None):
    metabolites = get_metabolites(xml, tool=tool)
    metanetx_metabolites = conversion_df[(conversion_df['External ID'].isin(metabolites))]['Internal ID'].tolist()
    kegg_metabolites = conversion_df[(conversion_df['Internal ID'].isin(metanetx_metabolites)) &
                                     (conversion_df['Source'] == 'KEGG')]['External ID'].tolist()
    print(f'Found {len(metanetx_metabolites)} reactions.')
    return metanetx_metabolites, kegg_metabolites


def calculate_quality_metrics(tp, fp, fn):
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1_score = 2 * (precision * recall) / (precision + recall)
    jaccard_distance = 1 - (tp / (fp + fp + fn))
    return precision, recall, f1_score, jaccard_distance
