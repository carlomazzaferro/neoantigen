import csv
from tabulate import tabulate


class PredictionCollection(object):

    def __init__(self, files):

        self.files = files
        self.collection = self.digest_multiple()

    def digest_multiple(self):

        pred_collection = []

        for file in self.files:
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

    @staticmethod
    def pass_to_pred_class(line, allele):
        data = line
        data.append(allele)
        return Prediction(data)

    @staticmethod
    def get_allele(iterator):
        return list(filter(None, next(iterator)))[0]


class Prediction(object):

    def __init__(self, data):
        """

        :param data: [peptide, allele, nmer, pos, protein, core_pep, h_avg_ranks, N_binders, affinity_lev]
        """
        self.data_list = self.check_data(data)

        self.original_position = self.data_list[0]
        self.peptide = self.data_list[1]
        self.protein = self.data_list[2]
        self.affinity_level = self.data_list[3]
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

        self.Swap = Swap(self.peptide)

    def __str__(self):
        return 'Basic Prediction Info: \n' + tabulate(self.lists_to_print)



    def check_data(self, data):
        #TODO: check data
        return data


class Swap(object):

    list_AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    def __init__(self, peptide):
        self.original_peptide = peptide
        self.swaps = self._create_swaps()

    def _create_swaps(self):
        list_peps = []
        for i in range(len(self.original_peptide)):
            for k in range(len(self.list_AA)):
                list_peps.append(self._insert_aa(self.original_peptide, i, self.list_AA[k]))

        return [i for i in list_peps if i != self.original_peptide]

    @staticmethod
    def _insert_aa(string, index, aa):
        hash_string = list(string)
        del hash_string[index]
        hash_string.insert(index, aa)
        return "".join(hash_string)