"""
import sys
sys.path.append('/Users/carlomazzaferro/Documents/CCBB/neoantigen/VarP-master/')

import glob
from VarP import utils

reader_names = ['Tumor_RNA_Reader', 'Tumor_Targeted_Reader',
           'Normal_DNA_Reader', 'Normal_Blood_Reader',
           'Somatic_Mutect_Reader']


#path_to_files = '/Volumes/Seagate Backup Plus Drive 1/vcf_files/varcode_to_test/'


myhandler = utils.HandleReaders(reader_names)
list_collections = myhandler.create_collection_from_readers(path_to_files)

type(list_collections[4])

my_collection = list_collections[4]


list_coding_effects = myhandler.return_list_coding_effects(my_collection)
protein_list = myhandler.return_protein_list(list_coding_effects)

dataframe = myhandler.return_dataframe(protein_list, list_coding_effects)

#Fasta file
myhandler.generate_fasta_file(dataframe, path_to_files+'peps.txt')

"""
import sys
sys.path.append('/Users/carlomazzaferro/Documents/CCBB/neoantigen/VarP-master/')
import VarP
from VarP import downstream_analysis
import pandas
from pymongo import MongoClient
pvac_seq_out_file = '/Volumes/Seagate Backup Plus Drive/vcf_files/not_annotated/OUT_PVAQ_SEQ/MutecT_filtered.tsv'

#This one contains gene expression data.
gene_expression_file = '/Volumes/Seagate Backup Plus Drive/bam_files_extras/TCGA_normalized_deseq_RNAseq_counts.csv'

df = pandas.read_csv(pvac_seq_out_file, sep='\t')
gene_expression_file = pandas.read_csv(gene_expression_file)

db_name = 'Variant_Prioritization_Workflow'
collection_name = 'Somatic_Filtered_Variants'

enricher = downstream_analysis.DataEnrichment(self.collection_name, self.db_name,
                                                           self.test_pvac_file, self.test_gene_data_file)

to_query = downstream_analysis.retrieve_chr_start_end(df)
mylist = enricher.retrieve_multiple(to_query)


