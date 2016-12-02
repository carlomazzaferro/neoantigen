import csv
import pandas
import sys
sys.path.append('/Users/carlomazzaferro/Documents/Code/neoantigen/')
from mhc_parser import utilities, methods
from tabulate import tabulate
import os
from shutil import rmtree
import importlib
importlib.reload(methods)


class PredictionCollection(object):

    threshold_levels = {'High': [0,50],
                        'Intermediate': [50,500],
                        'Low': [500, 5000],
                        'No': [5000, 100000]}

    cols = ['Pos', 'Peptide', 'Protein', 'Affinity_level', 'rank', 'core_pep',
            'h_avg_ranks', 'n_binders', 'allele', 'nmer']

    def __init__(self, files, fasta_file):

        self.files = files
        self.fasta = fasta_file
        self.protein_list = self.protein_from_fasta()
        self.dictionary_collection = self.dic_initiate()

    def return_protein_df_list(self):
        list_dfs = []

        for prot in self.protein_list:
            prot_df = self.protein_dataframe(prot)
            list_dfs.append(prot_df)

        return list_dfs

    def dic_initiate(self):
        prediction_collection = {}
        for protein in self.protein_list:
            prediction_collection[protein] = []

        return prediction_collection

    def protein_from_fasta(self):

        names, seqs = utilities.create_separate_lists(self.fasta)
        names = [name[0:15] for name in names]
        return names

    def sequence_from_fasta(self):
        names, seqs = utilities.create_separate_lists(self.fasta)
        return seqs

    def protein_dataframe(self, protein_id):

        protein_data = self.dictionary_collection[protein_id]
        pred_data = []

        for dat in protein_data:
            elemnts = []
            for elem in dat.data_list:
                elemnts.append(elem)
            elemnts.append(dat.nmer)
            pred_data.append(elemnts)

        return pandas.DataFrame(pred_data, columns=self.cols)

    def digest_multiple(self):

        pred_collection = []
        i = 0
        for file in self.files:
            print(file, i)
            i = i + 1
            pred_list = self.digest(file)
            pred_collection.extend(pred_list)

        return pred_collection

    def digest(self, file):

        with open(file) as csv_file:

            iterator = csv.reader(csv_file, delimiter='\t')
            allele = self.get_allele(iterator)
            _ = next(iterator)

            single_preds = []

            for line in iterator:
                pred = self.pass_to_pred_class(line, allele)
                single_preds.append(pred)

        return single_preds

    def pass_to_pred_class(self, line, allele):

        data = line
        data.append(allele)
        pred_object = Prediction(data)
        self.dictionary_collection[data[2]].append(pred_object)

        return Prediction(data)

    def predict_swaps(self, protein_id, threshold, net_mhc_path, dirty_mode=True):
        """

        Create predictions for each swap in a specific protein's PredictionCollection.
        Only peptides that score above a specific affinity level are considered and predicted on.
        This works by grouping swaps into a fasta file and automatically running netMHC on them.
        This creates an excel file with the predictions for each swap which is the parsed with the
        Prediction class.

        :param protein_id: the protein of interest
        :param threshold: in affinity (nM)
        :param net_mhc_path: path to local netMHC install
        :param dirty_mode: if True, keep temporary files
        :return:
        """

        threshold_range = self._return_threshold_level(threshold)
        prot_dic = self.dictionary_collection[protein_id]
        tmp_dir = "".join([os.path.dirname(self.fasta), '/temp/'])
        os.mkdir(tmp_dir)
        filtered_pred = []

        for pred in prot_dic:
            if (pred.affinity_level > threshold_range[0]) & (pred.affinity_level < threshold_range[1]):
                filtered_pred.append(pred)

        for filt in filtered_pred:
            tmp_fasta_file = methods.write_to_fasta(filt, protein_id, tmp_dir)
            methods.run_mhc(filt, tmp_fasta_file, net_mhc_path)

        for filt in filtered_pred:
            pred_file = methods.retrive_xls(filt, tmp_dir)
            self.parse_swap_preds(pred_file, filt)

        self.cleanup(dirty_mode, tmp_dir)

    def parse_swap_preds(self, pred_file, pred):
        single_preds = self.digest(pred_file)
        swaps = pred.Swap.swaps
        pred.Swap.swaps = [{i: j} for i, j in zip(swaps, single_preds)]

    def filter_low_aff(self, protein_id, threshold):
        return methods.filter_low_affinity(self.dictionary_collection, protein_id, threshold)

    @staticmethod
    def cleanup(mode, tmp_dir):
        if mode:
            rmtree(tmp_dir)


    @staticmethod
    def get_allele(iterator):
        return list(filter(None, next(iterator)))[0]

    def _return_threshold_level(self, threshold):

        if not threshold:
            return [0,50]
        if isinstance(threshold, str):
            return self.threshold_levels[threshold]
        if isinstance(threshold, int):
            return [0, threshold]


class Prediction(object):

    threshold_levels = {'High': [0,50],
                        'Intermediate': [50,500],
                        'Low': [500, 5000],
                        'No': [5000, 100000]}

    def __init__(self, data):
        """

        :param data: [peptide, allele, nmer, pos, protein, core_pep, h_avg_ranks, N_binders, affinity_lev]
        """
        self.data_list = self.check_data(data)
        self.original_position = self.data_list[0]
        self.peptide = self.data_list[1]
        self.protein = self.data_list[2]
        self.affinity_level = float(self.data_list[3])
        self.rank = self.data_list[4]
        self.core_pep = self.data_list[5]
        self.h_avg_ranks = self.data_list[6]
        self.n_binders = self.data_list[7]
        self.allele = self.data_list[8]
        self.nmer = len(self.peptide)

        self.lists_to_print = [["Peptide Seq", self.peptide],
                               ["Allele", self.allele],
                               ["Peptide Length", self.nmer],
                               ["Original Pos", self.original_position],
                               ["Affinity Level", self.affinity_level],
                               ["Core Pep", self.core_pep],
                               ["Average Rank", self.h_avg_ranks],
                               ["Number of High Binders", self.n_binders],
                               ["Protein of Origin", self.protein]]

        self.Swap = self.get_swaps()

    def __str__(self):
        return 'Basic Prediction Info: \n' + tabulate(self.lists_to_print)

    def get_swaps(self):
        if self.affinity_level < 500:
            return Swap(self.peptide, self.nmer, self.allele, self.protein, self.original_position)
        else:
            return []

    def check_data(self, data):
        #TODO: check data
        return data


class Swap(object):

    #TODO: test this.

    list_AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    def __init__(self, peptide, nmer, allele, protein, pos):

        self.original_peptide = peptide
        self.original_position = pos
        self.protein = protein
        self.allele = allele
        self.nmer = nmer
        self.swaps = self.gen_swaps()

    def swaps_func(self):
        for index, _ in enumerate(self.original_peptide):
            for letter in self.list_AA:
                temp = list(self.original_peptide)
                temp[index] = letter
                yield ''.join(temp)

    def generate_all_variants(self):
        for i in range(len(self.original_peptide)):
            head = self.original_peptide[:i]
            tail = self.original_peptide[i + 1:]
            for letter in self.list_AA:
                yield head + letter + tail

    def gen_swaps(self):
        return [v for v in self.generate_all_variants() if v != self.original_peptide]

