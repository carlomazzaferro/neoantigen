import unittest


class TestModels(unittest.TestCase):

    def setUp(self):

        self.peptide = 'DCQEGHILY'
        self.amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        self.allele = 'HLA-A0102'
        self.nmer = len(self.peptide)
        self.protein_name = 'Streptococcus_Pyogenes'
        self.affinity_level = 740
        self.original_pos = 1200

    def test_swapping(self):
        self.assertEqual(len(self.peptide), self.nmer)
        self.assertEqual()



if __name__ == '__main__':
    unittest.main()


import csv
fl = '/Users/carlomazzaferro/Documents/Code/neoantigen/antigen_discovery/tests/fasta_base_new_prots_HLA-A0101_9.xls'
with open(fl) as csvf:
    data = []
    rd = csv.reader(csvf, delimiter='\t')
    allele = list(filter(None, next(rd)))
    ids = next(rd)
    for i in rd:
        data.append(i)

 #   for row in rd:
#        print(row)