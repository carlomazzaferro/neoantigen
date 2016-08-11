import os
import re
import glob
import pandas as pd
from varcode import load_vcf


class HandleReaders(object):

    def __init__(self, reader_names):
        """
        Attribute of class
        :param reader_names: desired names of the reader objects
        """
        self.readers = reader_names

    def create_collection_from_readers(self, path_to_dir):
        """
        Given a directory as input, the method will return a list of reader objects that matches the given keys
        :param path_to_dir: directory where the fiktered vcf files reside
        :return: list of reader objects
        """
        os.chdir(path_to_dir)
        reader_list = []
        files = glob.glob("*.vcf")
        files.sort()
        self.readers.sort()

        for i in range(0, len(files)):
            if len(files) != len(self.readers):
                print files
                raise ValueError("Number of reader names and number of vcf files doesn't match")

            self.readers[i] = load_vcf(files[i])
            reader_list.append(self.readers[i])

        return reader_list

    def return_list_coding_effects(self, variant_collection):
        """
        Obtain a list of variant.Effect objects (defined in Varcode's documentation), given a reader object.
        Method removes duplicate and non-coding effects.
        :param variant_collection: collection of variants from vcf reader
        :return: list of coding effects
        """

        list_of_effects = []
        for i in variant_collection:
            list_of_effects.append(i.effects())  #list of effects

        list_of_coding_effects = []
        for i in list_of_effects:
            list_of_coding_effects.append(i.drop_silent_and_noncoding())   #list of coding effects

        effect_clean = self.remove_empty_variants(list_of_coding_effects)  #rem duplicates and noncoding effects
        return effect_clean

    def return_protein_list(self, effects):
        """
        From grch_37, varcode can extract the effect on protein sequences that a mutation can cause.
        :param effects: list of effects
        :return: list of protein sequences
        """
        prot_list = self.get_protein_list(effects)
        return prot_list

    def return_dataframe(self, prot_list, effects):
        """
        Create a nicely displayable dataframe that contains, for each variant description, an associated protein
        sequences that is mutated by the variant.
        :param prot_list: list of mutated proteins
        :param effects: list of variant effects
        :return: dataframe of lists with removed duplicates
        """
        rem_dup = self.rem_duplicates(prot_list, effects)
        return rem_dup

    def generate_fasta_file(self, df, file_path_and_name):

        unique_prot_list = list(df.prot.values)
        var_names = df.variants.apply(lambda x: str(x))
        unique_var_name = list(var_names)
        unique_var_name = self.clean_list_coding_effects(unique_var_name)

        with open(file_path_and_name, 'w') as outfile:
            for i in range(0, len(unique_var_name)):
                outfile.write(unique_var_name[i] + "\n" + unique_prot_list[i] + "\n")
        outfile.close()

    @staticmethod
    def clean_list_coding_effects(effects_clean):
        long_string = []
        for i in effects_clean:
            l = re.split(r"[\s,]+", i)
            long_string.append('>' + '-'.join(l))

        return long_string

    @staticmethod
    def remove_empty_variants(variant_list):
        # Remove empty
        non_empty = [x for x in variant_list if x != []]
        # Flatten list
        flattened_list = [item for sublist in non_empty for item in sublist]
        return flattened_list

    @staticmethod
    def get_protein_list(variant_list):
        # Varcode's method variant.mutant_protein_sequence will output the proteins affected by it.
        list_of_proteins = []
        for i in range(0, len(variant_list)):
            list_of_proteins.append(str(variant_list[i].mutant_protein_sequence))

        return list_of_proteins

    @staticmethod
    def rem_duplicates(protein_list, variant_list):
        # Remove duplicate variants/proteins
        df = pd.DataFrame({'variants': variant_list, 'prot': protein_list})
        df = df.drop_duplicates(subset='variants')
        df = df.drop_duplicates(subset='prot')
        return df
