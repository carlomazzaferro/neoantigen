from variantannotation import MongoDB_querying
from pymongo import MongoClient


collection_name = ['ANNOVAR_MyVariant_chunks', 'ANNOVAR_MyVariant_full',
                   'My_Variant_Info_Collection_Full', 'My_Variant_Info_Collection_Chunks']
db_name = 'My_Variant_Database'

l_of_variant_lists = []
for i in range(0, len(collection_name)):
    l_of_variant_lists.append(MongoDB_querying.query(collection_name[i], db_name))
    print l_of_variant_lists




"""
filter1 = df.query('ThousandGenomeAll < 0.05 or ThousandGenomeAll != ThousandGenomeAll', engine='python')

filter2 = filter1.query('esp6500siv2_all != esp6500siv2_all or esp6500siv2_all < 0.05 or cosmic70 == cosmic70', engine='python')

filter3 = filter2.query('Func_knownGene == "exonic" or Func_knownGene == "exonic;splicing" or Func_knownGene == "splicing"', engine='python')

filter4 = filter3.query('ExonicFunc_knownGene != "synonymous SNV"', engine='python')

filtered_DP = df.query('Genotype_call.DP > 10', engine='python')
"""

##########TEST VARCODE##########

from varcode import Variant
from pyensembl import ensembl_grch37


def extract_info(list_rare_variants):
    for i in l_of_variant_lists[0]:


l_of_variant_lists[0][0]["hgvs_key"]
