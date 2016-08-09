import os
import sys
#sys.path.append('/Users/carlomazzaferro/Documents/Code/variantannotation-master')

from variantannotation import annotate_batch
from variantannotation import myvariant_parsing_utils
from variantannotation import mongo_DB_export
from variantannotation import annovar_subprocess

#set paths
filepath = "/Users/carlomazzaferro/Desktop/CSV to be tested"
csv_file = "Tumor_targeted_processed.csv"
vcf_file = "Tumor_targeted_seq.vqsr.vcf"
os.chdir(filepath)

"""
ANNOVAR_PATH = '/database/annovar/'
IN_PATH = '/data/Nof1/file.vcf'
OUT_PATH = '/data/ccbb_internal/interns/Carlo/annovar_results'

"""

#1. Get csv file: run annovar
"""
annovar_subprocess.run_annovar(ANNOVAR_PATH, IN_PATH, OUT_PATH)
"""

#METHOD 1: by chunks, iteratively.
chunksize = 1000
step = 0
collection_name = 'ANNOVAR_MyVariant_chunks'
db_name = 'My_Variant_Database'

#Get variant list. Should always be the first step after running ANNOVAR
open_file = myvariant_parsing_utils.VariantParsing()
list_file = open_file.get_variants_from_vcf(vcf_file)

#Run process, export to MongoDB in-built
as_batch = annotate_batch.AnnotationMethods()
as_batch.by_chunks(list_file, chunksize, step, csv_file, collection_name, db_name)


#---------------#--------------#---------------#--------------#---------------#--------------#---------------#

#METHOD 2: usign full file, and holding it in memory (OK for smaller files)   ##TEST THIS##

#get variant list. Should always be the first step after running ANNOVAR
open_file = myvariant_parsing_utils.VariantParsing()
list_file = open_file.get_variants_from_vcf(vcf_file)

#Run process, data saved to joint_list
as_one_file = annotate_batch.AnnotationMethods()
joint_list = as_one_file.full_file(list_file, csv_file)

#Name Collection & DB
collection_name = 'ANNOVAR_MyVariant_full'
db_name = 'My_Variant_Database'

#Export
exporting_function = mongo_DB_export.export
exporting_function(joint_list, collection_name, db_name)


#---------------#--------------#---------------#--------------#---------------#--------------#---------------#
#METHOD 3: ignore annovar, get data solely from myvariant (much faster, requires nothing but a VCF file.
#will however be incomplete (some variants will have no information).

#Get variant list form vcf file
open_file = myvariant_parsing_utils.VariantParsing()
list_file = open_file.get_variants_from_vcf(vcf_file)

#Run process
my_variants = annotate_batch.AnnotationMethods()
myvariant_data = my_variants.my_variant_at_once(list_file)

#Name Collection & DB
collection_name = 'My_Variant_Info_Collection_Full'
db_name = 'My_Variant_Database'

#Export
exporting_function = mongo_DB_export.export
exporting_function(myvariant_data, collection_name, db_name)

#---------------#--------------#---------------#--------------#---------------#--------------#---------------#
#METHOD 4: ignore annovar, Get data solely from myvariant (much faster, requires nothing but a VCF file.
#will however be incomplete (some variants will have no information).
#Do so BY CHUNKS. Export function is built in the methods myvariant_chunks

chunksize = 1000
step = 0

#Get variant list from vcf file
open_file = myvariant_parsing_utils.VariantParsing()
list_file = open_file.get_variants_from_vcf(vcf_file)

#Name Collection & DB
collection_name = 'My_Variant_Info_Collection_Chunks'
db_name = 'My_Variant_Database'

#Run process, export to MongoDB in-built
my_variants = annotate_batch.AnnotationMethods()
myvariant_data = my_variants.myvariant_chunks(list_file, chunksize, step, collection_name, db_name)


