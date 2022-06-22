from analysis_scripts.utils import FileUtils


def convert_ecoli_genome_to_have_only_the_locus_tag():
    FileUtils.change_gene_ids_to_locus_tag("Ecoli/Ecoli.faa", "Ecoli/Ecoli_with_locust_tag.faa")

def convert_pputida_genome_to_have_only_the_locus_tag():
    FileUtils.change_gene_ids_to_locus_tag("Pputida/Pputida.faa", "Pputida/Pputida_with_locust_tag.faa")

def convert_pfalciparum_genome_to_have_only_the_locus_tag():

    with open("Pfalciparum/Pfalciparum.faa", "r") as file:
        lines = file.readlines()

        with open("Pfalciparum/Pfalciparum_locus_tag.faa", "w") as file_locus:
            for line in lines:
                if ">" in line:
                    locus_tag = line.split(" ")[0]
                    file_locus.write(locus_tag + "\n")
                else:
                    file_locus.write(line)


convert_pfalciparum_genome_to_have_only_the_locus_tag()