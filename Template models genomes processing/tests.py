from unittest import TestCase


from utils import FileUtils


class TestFileUtils(TestCase):

    def test_read_fasta_file(self):
        FileUtils.read_fasta_file("Llactis/Llactis.faa")

    def test_change_genes_name(self):
        FileUtils.change_gene_ids_to_locus_tag("Llactis/Llactis.faa", "Llactis/Llactis_locus_tag.faa")

    def test_generate_gene_map(self):
        FileUtils.generate_gene_map("Llactis/Llactis.faa", "dict_genes.txt")

    def test_read_genebank_file(self):
        FileUtils.convert_genebank_file_to_fasta("Pputida/Pputida.gb", "Pputida/Pputida_converted.faa")


