import unittest
import sys
sys.path.append('/Users/carlomazzaferro/Documents/CCBB/neoantigen/VarP-master/')

import os
import glob
from VarP import downstream_analysis
from shutil import copyfile
import pandas as pd


class TestHandleReaders(unittest.TestCase):

    def setUp(self):
        self.test_directory = os.path.dirname(os.path.realpath('__file__'))
        self.test_pvac_file = os.path.join(self.test_directory,'Test_Tsv_File.tsv')
        self.test_gene_data_file = os.path.join(self.test_directory,'Test_Uveal_Genes.result')
        self.collection_name = 'Somatic_Filtered_Variants'
        self.db_name = 'Variant_Prioritization_Workflow'
        self.enricher = downstream_analysis.DataEnrichment(self.db_name, self.collection_name,
                                                           self.test_pvac_file, self.test_gene_data_file)
        self.pvac_headers = ['Chromosome', 'Start', 'Stop', 'Ensembl Gene ID']
        self.gene_data_headers = ['gene_id', 'FPKM']
        self.col_headers = ['Normal Ref Count', 'Normal Var Count', 'Tumor DNA Ref Count',
                            'Tumor DNA Var Count', 'Tumor RNA Ref Count', 'Tumor RNA Var Count',
                            'Gene Name', 'HLA Allele', 'Peptide Length', 'MT Epitope Seq', 'WT Epitope Seq']

    def test_get_dfs(self):
        df1, df2 = self.enricher.get_dfs()
        gene_data_df = pd.read_csv(self.test_gene_data_file, sep='\t')

        if set(self.pvac_headers).issubset(set(df1.columns)) is False:
            raise ValueError("Dataframe headers are not present in pVAC-seq file")

        if set(self.gene_data_headers).issubset(set(gene_data_df.columns)) is False:
            raise ValueError("Dataframe headers are not present in gene expression file")

    def test_integrate_data(self):
        out_df_normal = self.enricher.integrate_data(type_data="Normal")
        out_df_rna = self.enricher.integrate_data(type_data="Tumor RNA")
        out_df_dna = self.enricher.integrate_data(type_data="Tumor DNA")

        joint = list(out_df_dna.columns) + list(out_df_normal.columns) + list(out_df_rna.columns)

        if set(self.col_headers).issubset(set(joint)) is False:
            raise ValueError("Dataframe headers are not present in out file")

    def test_retrieve_multiple(self):
        list_out = self.enricher.retrieve_multiple()
        print list_out[0]
        if set(['Start', 'End', 'Chr']).issubset(set(list(list_out[0][0].keys()))) is False:
            raise ValueError("Dataframe headers are not present in gene expression file")


if __name__ == '__main__':
    unittest.main()
