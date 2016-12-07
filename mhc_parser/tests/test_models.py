import unittest
import csv
import os
import sys
import glob
sys.path.append('/Users/carlomazzaferro/Documents/Code/neoantigen/')

from mhc_parser import utilities, methods, models

class TestModels(unittest.TestCase):

    def setUp(self):

        self.peptide = 'DCQEGHILY'
        self.amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        self.allele = 'HLA-A0102'
        self.nmer = len(self.peptide)
        self.protein_name = 'Streptococcus_Pyogenes'
        self.affinity_level = 740
        self.original_pos = 1200
        self.fasta_test = '/Users/carlomazzaferro/Documents/Code/neoantigen/mhc_parser/tests/' \
                             'test_files/test_fasta.fasta'
        self.reduced_fasta = '/Users/carlomazzaferro/Documents/Code/neoantigen/mhc_parser/tests/' \
                             'test_files/test_fasta_reduced.fasta'
        self.prot_names_test_fasta = ['S__pyogenes_Cas9', 'Staphylococcus_aureus_Cas9', 'S_CRISPR_1_thermophilus_Cas9', 'N__meningitidis_Cas9']
        self.prots_initals_fasta = ['MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIGALLFDSGE',
                                    'MKRNYILGLDIGITSVGYGIIDYETRDVIDAGVRLFKEANVENNEGRRSK',
                                    'MSDLVLGLDIGIGSVGVGILNKVTGEI',
                                    'MAAFKPNPINYILGLDIGIASVGWAMVEIDEDENPICLIDLGVRVFERAEVPKTGDSLAMARRLARSVRR']

        self.net_mhc_file_test = glob.glob('/Users/carlomazzaferro/Documents/Code/neoantigen/mhc_parser/' \
                                 'tests/test_files/mhc_preds_test_fasta_reduced/*B1501_9.xls')

        self.model_methods = models.PredictionCollection(self.net_mhc_file_test, self.reduced_fasta)
        self.init_dic = self.model_methods.dic_initiate()
        self.top_excel_lines = [['', '', '', 'HLA-B1501', '', '', '', ''],
                                ['Pos', 'Peptide', 'ID','nM', 'Rank','Core', 'H_Avg_Ranks','N_binders'],
                                [0, 'MDKKYSIGL', 'Test1', 9694.8, 14, 'MDKKYSIGL', 14, 0],
                                [1, 'DKKYSIGLD', 'Test1', 42058.3, 95, 'DKKYSIGLD', 95, 0],
                                [2, 'KKYSIGLDI', 'Test1', 22207.5, 31, 'KKYSIGLDI', 31, 0],
                                [3, 'KYSIGLDIG', 'Test1', 29965.4, 48, 'KYSIGLDIG', 48, 0]]

        self.swap_letters = ['A', 'B', 'C', 'D', 'E']
        self.swap_data = {'ACC': ['BCC', 'CCC', 'DCC', 'ECC',
                                  'AAC', 'ABC', 'ADC', 'AEC',
                                  'ACA', 'ACB', 'ACD', 'ACE'],

                          'BDE': ['ADE', 'CDE', 'DDE', 'EDE',
                                  'BAE', 'BBE', 'BCE', 'BEE',
                                  'BDA', 'BDB', 'BDC', 'BDD']}

    def test_swapping(self):
        self.assertEqual(len(self.peptide), self.nmer)

    def test_names_seq_retrieval_from_fasta(self):
        names, seqs = utilities.create_separate_lists(self.fasta_test)
        self.assertEqual(names, self.prot_names_test_fasta)
        for idx, seq in enumerate(seqs):
            self.assertIn(self.prots_initals_fasta[idx], seq)

    def test_dic_initiate(self):
        init_dic = self.model_methods.dic_initiate()
        self.assertEqual(list(init_dic.keys()).sort(), ['Test4', 'Test3', 'Test2', 'Test1'].sort())

    def test_digest(self):

        #Reproduction of parsing function on trial data
        iterator = iter(self.top_excel_lines)
        allele = self.model_methods.get_allele(iterator)
        _ = next(iterator)
        single_preds = []

        for line in iterator:
            pred = self.model_methods.pass_to_pred_class(line, allele)
            single_preds.append(pred)

        for i in single_preds:
            print(i.nM)
        print(self.model_methods.files)

        #parsing using the actual module
        multi_preds = self.model_methods.digest_multiple()[0:4]

        for i, _ in enumerate(single_preds):
            self.assertEqual(single_preds[i].Allele, multi_preds[i].Allele)
            self.assertEqual(single_preds[i].Nmer, multi_preds[i].Nmer)
            self.assertEqual(single_preds[i].nM, multi_preds[i].nM)
            self.assertEqual(single_preds[i].Peptide, multi_preds[i].Peptide)
            self.assertEqual(single_preds[i].ID, multi_preds[i].ID)
            self.assertEqual(single_preds[i].lists_to_print, multi_preds[i].lists_to_print)

    @staticmethod
    def get_allele(iterator):
        return list(filter(None, next(iterator)))[0]

    def test_swaps_generation(self):
        original_pep1 = list(self.swap_data.keys())[0]
        original_pep2 = list(self.swap_data.keys())[1]
        swaps_1 = self.swap_data[original_pep1]
        swaps_2 = self.swap_data[original_pep2]

        self.assertEqual(swaps_1, self.gen_swaps(original_pep1))
        self.assertEqual(swaps_2, self.gen_swaps(original_pep2))

    #Functions taken verbatim from models class
    def generate_all_variants(self, pep):
        for i in range(len(pep)):
            head = pep[:i]
            tail = pep[i + 1:]
            for letter in self.swap_letters:
                yield head + letter + tail

    def gen_swaps(self, pep):
        return [v for v in self.generate_all_variants(pep) if v != pep]

if __name__ == '__main__':
    unittest.main()