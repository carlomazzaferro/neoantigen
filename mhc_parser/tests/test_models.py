import unittest
import csv
import os
import sys
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
        self.fasta_test = 'test_fasta.fasta'
        self.prot_names_test_fasta = ['S__pyogenes_Cas9', 'Staphylococcus_aureus_Cas9', 'S_CRISPR_1_thermophilus_Cas9', 'N__meningitidis_Cas9']
        self.prots_initals_fasta = ['MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIGALLFDSGE',
                                    'MKRNYILGLDIGITSVGYGIIDYETRDVIDAGVRLFKEANVENNEGRRSK',
                                    'MSDLVLGLDIGIGSVGVGILNKVTGEI',
                                    'MAAFKPNPINYILGLDIGIASVGWAMVEIDEDENPICLIDLGVRVFERAEVPKTGDSLAMARRLARSVRR']

    def test_swapping(self):
        self.assertEqual(len(self.peptide), self.nmer)

    def test_names_seq_retrieval_from_fasta(self):
        names, seqs = utilities.create_separate_lists(self.fasta_test)
        self.assertEqual(names, self.prot_names_test_fasta)
        for idx, seq in enumerate(seqs):
            self.assertIn(self.prots_initals_fasta[idx], seq)

    def test_dic_initiate(self):


if __name__ == '__main__':
    unittest.main()

    """
    fl = '/Users/carlomazzaferro/Documents/Code/neoantigen/antigen_discovery/tests/fasta_base_new_prots_HLA-A0101_9.xls'
    with open(fl) as csvf:
        data = []
        rd = csv.reader(csvf, delimiter='\t')
        allele = list(filter(None, next(rd)))
        ids = next(rd)
        for i in rd:
            data.append(i)

        for row in rd:
            print(row)
    """