
import sys
sys.path.append('/Users/carlomazzaferro/Documents/CCBB/neoantigen/VarP-master/')

import glob
from VarP import utils

reader_names = ['Tumor_RNA_Reader', 'Tumor_Targeted_Reader',
           'Normal_DNA_Reader', 'Normal_Blood_Reader',
           'Somatic_Mutect_Reader']


path_to_files = '/Volumes/Seagate Backup Plus Drive 1/vcf_files/varcode_to_test/'


myhandler = utils.HandleReaders(reader_names)
list_collections = myhandler.create_collection_from_readers(path_to_files)

type(list_collections[4])

my_collection = list_collections[4]


list_coding_effects = myhandler.return_list_coding_effects(my_collection)
protein_list = myhandler.return_protein_list(list_coding_effects)

dataframe = myhandler.return_dataframe(protein_list, list_coding_effects)

#Fasta file
myhandler.generate_fasta_file(dataframe, path_to_files+'peps.txt')

