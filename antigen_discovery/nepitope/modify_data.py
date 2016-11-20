from nepitope import pep_utils
from nepitope import merge_workflows
import re
import pandas
import random
import glob
import os


class ModifyData(object):

    reference = 'StreptococcusPyogenes_reference'

    def __init__(self, top_swaps_df, fasta_ref, fasta_out_dir, list_dfs):

        self.top_swaps = top_swaps_df
        self.fasta_ref = fasta_ref
        self.fasta_out_dir = fasta_out_dir
        self.df_list = list_dfs
        self.swaps_dic = self._hash_swaps()
        self.fasta_dic = self._hash_fasta()

    def get_modified_fasta(self):
        """
        In progress. Will be used to get a fasta file modified at more than one location
        :return:
        """

        zipped = self.fasta_dic
        orig_and_swaps = self.swaps_dic
        swapped_fasta = []

        for i in zipped:
            if i[0][1::] == self.reference:
                i[1] = self.replace_all(i[1], orig_and_swaps)
            swapped_fasta.append(i)

        self._write_fasta_out(swapped_fasta)

    def get_modified_fastas_single_swap(self):
        """
        Create a modified fasta file for each swap 'entry' (row in top_swaps_df). These files, alongside with
        a modified lsss_1 (list of protein dataframes) can then go into the PairwiseComp class and return the csv
        files containing the orthogonality matrices (one per each peptide swap).
        :return: None
        """

        zipped = self.fasta_dic
        orig_and_swaps = self.swaps_dic

        for orig, swap in orig_and_swaps.items():
            swapped_fasta = []
            for i in zipped:

                if i[0][1::] == self.reference:
                    if not orig in i[1]:
                        print(orig + '\n', i[1] + '\n')
                    new = i[1].replace(orig, swap)
                    swapped_fasta.insert(0, [i[0], new])
                swapped_fasta.append(i)
            del swapped_fasta[1]
            self._write_fasta_out(swapped_fasta, orig, swap)

    def _write_fasta_out(self, swapped_fasta, orig, swap):

        with open(self.fasta_out_dir + orig + '-' + swap + '.fasta', 'w') as out:
            for i in swapped_fasta:
                out.write(i[0] + '\n')
                out.write(i[1] + '\n')

    def _hash_swaps(self):

        peps_to_swap = self._get_swap_peps()
        orig_peps = self._get_orig_peps()

        return dict(zip(orig_peps, peps_to_swap))

    def _hash_fasta(self):

        idx, seq = pep_utils.create_separate_lists(self.fasta_ref)
        zipped = list(zip(idx, seq))

        return [list(x) for x in zipped]

    @staticmethod
    def replace_all(text, dic):
        for i, j in dic.items():
            text = text.replace(i, j)
        return text

    def get_singly_modified_df_list(self, original, swap):
        """
        Modify single entry in the list of protein dataframes.
        :param original: Original peptide
        :param swap: Peptide to be swapped
        :return: Modified list of dataframes
        """

        list_df = self.df_list

        for idx, df in enumerate(list_df):
            if df.ID.unique()[0] == self.reference:
                new_df = df
                new_df.Peptide = new_df.Peptide.str.replace(original, swap)
                ix = new_df.loc[new_df.Peptide == swap].index
                new_df.set_value(ix, 'Affinity Level', 'No')
                list_df[idx] = new_df

        return list_df

    def _get_swap_peps(self):

        #TODO: something better than random choice for the peptide swap

        swaps = []
        list_swaps = list(self.top_swaps['top scoring peptides'].values)
        for i in list_swaps:
            swap = list(filter(None, re.sub('[^A-Za-z0-9]+', ',', i).split(',')))
            swaps.append(random.choice(swap))

        return swaps

    def _get_orig_peps(self):
        return self.top_swaps['original peptide'].values.tolist()

    def create_signly_modified_csvs(self, csv_dir):
        """
        Returns a csv file for each entry in the top_swaps_df. Must be run after getting the swapped fastas with
        get_modified_fastas_single_swap.

        :param csv_dir: Directory where the csv files will be created
        :return: List of dataframes (each dataframe is also written to csv).
        """

        fastas = glob.glob(self.fasta_out_dir + '*')
        exchange_pairs = self._exchange_pairs()
        list_dfs = []

        for pair in exchange_pairs:
            for fasta in fastas:
                if os.path.splitext(os.path.basename(fasta))[0].split('-')[0] == pair[0]:
                    my_fasta = fasta
                    ls_mod = self.get_singly_modified_df_list(pair[0], pair[1])
                    pwcomp_pair_swap = merge_workflows.PairwiseComp(ls_mod, 5, my_fasta)
                    df_swap = pwcomp_pair_swap.pipe_run()
                    df_swap = df_swap.drop('Matches Loc', 1)
                    df_swap.to_csv(csv_dir + pair[0] + '-' + pair[1] + '.csv')
                    list_dfs.append(df_swap)

        return list_dfs

    def _exchange_pairs(self):

        exchange_pairs = []
        for i, j in self.swaps_dic.items():
            exchange_pairs.append([i, j])
        return exchange_pairs
