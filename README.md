

```python
from VarP import utils

"""
Set reader names. Since these and the files will be sorted automatically
by name, make sure the names of readers and files are similar
"""

reader_names = ['Tumor_RNA_Reader', 'Tumor_Targeted_Reader',
           'Normal_DNA_Reader', 'Normal_Blood_Reader',
           'Somatic_Mutect_Reader']

path_to_files = '/Volumes/Seagate Backup Plus Drive/vcf_files/varcode_to_test/'

#Initialize class HandleReaders.
myhandler = utils.HandleReaders(reader_names)

"""
Create a list of variant collections. See varcode's documentation for
more info on this type of data
"""
list_collections = myhandler.create_collection_from_readers(path_to_files)

type(list_collections[4])  #variant.Collection

#Select one to work on
my_collection = list_collections[4]

#Obtain codon effects from variant collection
list_coding_effects = myhandler.return_list_coding_effects(my_collection)

#Obtain protein list from variant collection
protein_list = myhandler.return_protein_list(list_coding_effects)

#Return a dataframe for easy viz
dataframe = myhandler.return_dataframe(protein_list, list_coding_effects)

#Generate fasta file for post-processing
myhandler.generate_fasta_file(dataframe)
```