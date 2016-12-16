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

    cols = ['Pos', 'Peptide', 'ID', 'nM', 'Rank', 'Core',
            'H_Avg_Ranks', 'N_binders', 'Allele', 'Nmer']

    def __init__(self, files, fasta_file, threshold=None):

        self.files = files
        self.fasta = fasta_file
        self.threshold = self._return_threshold_level(threshold)
        self.protein_list = self.protein_from_fasta()
        self.seqs_list = self.sequence_from_fasta()
        self.dictionary_collection = self.dic_initiate()

    def dic_initiate(self):
        prediction_collection = {}
        prot_seq_dic = dict(zip(self.protein_list, self.seqs_list))

        for protein in self.protein_list:
            for key in prot_seq_dic.keys():
                if key.startswith(protein):

                    prediction_collection[protein] = {}
                    prediction_collection[protein]['Predictions'] = []
                    prediction_collection[protein]['ProtSeq'] = prot_seq_dic[key]
                    prediction_collection[protein]['High Affinity Ranges'] = []

        return prediction_collection

    def protein_from_fasta(self):

        names, seqs = utilities.create_separate_lists(self.fasta)
        names = [name[0:15] for name in names]
        return names

    def sequence_from_fasta(self):
        names, seqs = utilities.create_separate_lists(self.fasta)
        return seqs

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
        self.dictionary_collection[data[2]]['Predictions'].append(pred_object)
        if pred_object.nM < self.threshold[1]:
            self.dictionary_collection[data[2]]['High Affinity Ranges'].append(self.get_ranges(pred_object))

        return Prediction(data)

    @staticmethod
    def get_ranges(pred_object):

        final_pos = pred_object.Pos + pred_object.Nmer - 1
        ranges = [list(range(pred_object.Pos, final_pos + 1))]

        return ranges

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
        prot_dic = self.dictionary_collection[protein_id]['Predictions']
        tmp_dir = "".join([os.path.dirname(self.fasta), '/temp/'])
        os.mkdir(tmp_dir)
        filtered_pred = []

        for pred in prot_dic:
            if (pred.nM > threshold_range[0]) & (pred.nM < threshold_range[1]):
                filtered_pred.append(pred)

        for filt in filtered_pred:
            #create fasta file containging swpas and run netMHC on the fly on the newly created fasta
            tmp_fasta_file = methods.write_to_fasta(filt, protein_id, tmp_dir)
            methods.run_mhc(filt, tmp_fasta_file, net_mhc_path)

        for filt in filtered_pred:
            #parse the prediction to the Swap class
            pred_file = methods.retrive_xls(filt, tmp_dir)
            self.parse_swap_preds(pred_file, filt)

        self.cleanup(dirty_mode, tmp_dir)

    def parse_swap_preds(self, pred_file, pred):
        single_preds = self.digest(pred_file)
        swaps = pred.Swap.swaps
        pred.Swap.swaps = [{i: j} for i, j in zip(swaps, single_preds)]

    @staticmethod
    def cleanup(mode, tmp_dir):
        if mode:
            rmtree(tmp_dir)

    @staticmethod
    def get_allele(iterator):
        return list(filter(None, next(iterator)))[0]

    def _return_threshold_level(self, threshold):

        if not threshold:
            return [0, 50]
        if isinstance(threshold, str):
            return self.threshold_levels[threshold]
        if isinstance(threshold, int):
            return [0, threshold]

    def filter_low_aff(self, protein_id, threshold):
        return methods.filter_low_affinity(self.dictionary_collection, protein_id, threshold)

    def filter_all(self, threshold):
        return methods.filter_all(self.dictionary_collection, threshold)

    def protein_dataframe(self, protein_id):
        return methods.to_df(self.dictionary_collection, self.cols, protein_id)

    def return_protein_df_list(self):
        list_dfs = []
        for prot in self.protein_list:
            prot_df = self.protein_dataframe(prot)
            list_dfs.append(prot_df)

        return list_dfs


class Prediction(object):

    threshold_levels = {'High': [0, 50],
                        'Intermediate': [50, 500],
                        'Low': [500, 5000],
                        'No': [5000, 100000]}

    def __init__(self, data):
        """

        :param data: [peptide, allele, nmer, pos, protein, core_pep, h_avg_ranks, N_binders, affinity_lev]
        """
        self.data_list = self.check_data(data)
        self.Pos = int(self.data_list[0])
        self.Peptide = self.data_list[1]
        self.ID = self.data_list[2]
        self.nM = float(self.data_list[3])
        self.Rank = int(round(float(self.data_list[4])))
        self.Core = self.data_list[5]
        self.H_Avg_Ranks = int(round(float(self.data_list[6])))
        self.N_binders = int(self.data_list[7])
        self.Allele = self.data_list[8]
        self.Nmer = len(self.Peptide)
        self.lists_to_print = [["Peptide Seq", self.Peptide],
                               ["Allele", self.Allele],
                               ["Peptide Length", self.Nmer],
                               ["Original Pos", self.Pos],
                               ["Affinity Level", self.nM],
                               ["Core Pep", self.Core],
                               ["Average Rank", self.H_Avg_Ranks],
                               ["Number of High Binders", self.N_binders],
                               ["Protein of Origin", self.ID]]

        self.Swap = self.get_swaps()

    def __str__(self):
        return 'Basic Prediction Info: \n' + tabulate(self.lists_to_print)

    def get_swaps(self):
        if self.nM < 500:
            return Swap(self.Peptide, self.Nmer, self.Allele, self.ID, self.Pos)
        else:
            return []

    def check_data(self, data):
        #TODO: check data
        return data


class Swap(object):

    list_AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    def __init__(self, peptide, nmer, allele, protein, pos):

        self.original_peptide = peptide
        self.original_position = pos
        self.protein = protein
        self.allele = allele
        self.nmer = nmer
        self.swaps = self.gen_swaps()

    def generate_all_variants(self):
        for i in range(len(self.original_peptide)):
            #Slice original peptdie into head, tail; place letter form list of amino acids in between.
            head = self.original_peptide[:i]
            tail = self.original_peptide[i + 1:]
            for letter in self.list_AA:
                yield head + letter + tail    #Proves to be faster than "".join due to small length.
                                              #Returns a generator

    def gen_swaps(self):
        #Generator to list, exclude parent peptide
        return [v for v in self.generate_all_variants() if v != self.original_peptide]

