import re

from Bio import SeqIO


class ProteinSequence:
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        self._sequence = value

    @property
    def header(self):
        return self._sequence

    @header.setter
    def header(self, value):
        self._header = value

    @property
    def locus_tag(self):
        return self._locus_tag

    @locus_tag.setter
    def locus_tag(self, value):
        self._locus_tag = value


class FileUtils:

    @staticmethod
    def read_fasta_file(input_file):
        sequences = []
        fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
        for sequence in fasta_sequences:
            protein_sequence = str(sequence.seq)
            description = sequence.description
            new_sequence = ProteinSequence(description, protein_sequence)

            match = re.search("locus_tag=[0-9A-Za-z_]+", description)
            if match:
                locus_tag = match.group(0)
                locus_tag = locus_tag.split("=")[1]
                new_sequence.locus_tag = locus_tag

            sequences.append(new_sequence)

        return sequences

    @staticmethod
    def generate_gene_map(gene_file,output_file):
        sequences = FileUtils.read_fasta_file(gene_file)
        with open(output_file, 'w') as handle:
            for sequence in sequences:
                sbml_id = sequence.locus_tag
                handle.write(sbml_id)
                handle.write("\t")
                handle.write(sequence.locus_tag)
                handle.write("\n")


    @staticmethod
    def change_gene_ids_to_locus_tag(input_file, output_file):
        sequences = FileUtils.read_fasta_file(input_file)
        with open(output_file, 'w') as handle:

            for sequence in sequences:

                handle.write(">")
                handle.write(sequence.locus_tag)
                handle.write("\n")
                handle.write(sequence.sequence)
                handle.write("\n")

    @staticmethod
    def convert_genebank_file_to_fasta(file_path, output_file):
        sequences = []
        with open(file_path) as input_handle:
            for record in SeqIO.parse(input_handle, "genbank"):
                features = record.features
                for feature in features:
                    if "translation" in feature.qualifiers and "locus_tag" in feature.qualifiers:
                        sequence = ProteinSequence(feature.qualifiers["locus_tag"][0],
                                                   feature.qualifiers["translation"][0])
                        sequence.locus_tag = feature.qualifiers["locus_tag"][0]
                        sequences.append(sequence)

        with open(output_file, 'w') as handle:

            for sequence in sequences:

                handle.write(">")
                handle.write(sequence.locus_tag)
                handle.write("\n")
                handle.write(sequence.sequence)
                handle.write("\n")

