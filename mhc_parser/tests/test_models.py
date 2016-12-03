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

    def test_swapping(self):
        self.assertEqual(len(self.peptide), self.nmer)

    def test_names_seq_retrieval_from_fasta(self):
        names, seqs = utilities.create_separate_lists(self.fasta_test)
        self.assertEqual(names, self.prot_names_test_fasta)
        for idx, seq in enumerate(seqs):
            self.assertIn(self.prots_initals_fasta[idx], seq)

    def test_dic_initiate(self):
        init_dic = self.model_methods.dic_initiate()
        self.assertEqual(list(init_dic.keys()).sort(), ['Test4', 'Test3','Test2', 'Test1'].sort())

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
            print(i.affinity_level)
        print(self.model_methods.files)

        #parsing using the actual module
        multi_preds = self.model_methods.digest_multiple()[0:4]

        for i, _ in enumerate(single_preds):
            self.assertEqual(single_preds[i].allele, multi_preds[i].allele)
            self.assertEqual(single_preds[i].nmer, multi_preds[i].nmer)
            self.assertEqual(single_preds[i].affinity_level, multi_preds[i].affinity_level)
            self.assertEqual(single_preds[i].peptide, multi_preds[i].peptide)
            self.assertEqual(single_preds[i].protein, multi_preds[i].protein)
            self.assertEqual(single_preds[i].lists_to_print, multi_preds[i].lists_to_print)

    @staticmethod
    def get_allele(iterator):
        return list(filter(None, next(iterator)))[0]

if __name__ == '__main__':
    unittest.main()