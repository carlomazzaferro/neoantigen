import os
import sys
#sys.path.append('/Users/carlomazzaferro/Documents/Code/variantannotation-master')

from variantannotation import annotate_batch
from variantannotation import myvariant_parsing_utils
from variantannotation import mongo_DB_export
from variantannotation import create_output_files
from variantannotation import utilities
from variantannotation import MongoDB_querying


#set paths
collection_name = 'Test_Normal_Targeted'
db_name = 'My_Variant_Database'

#set paths
filepath = "/Volumes/Seagate Backup Plus Drive/vcf_files/"

VCF_FILE_NAMES = ['Tumor_RNAseq_variants.vcf', 'Tumor_targeted_seq.vcf',
                  'normal_targeted_seq.vcf', 'somatic_mutect_old.vcf',
                  'normal_blood_WGS.vqsr.vcf']

CSV_FILE_NAME = ['Tumor_RNAseq_variants.hg19_multianno.csv', 'Tumor_targeted_seq.hg19_multianno.csv',
                 'normal_targeted_seq.hg19_multianno.csv','somatic_mutect_old.hg19_multianno.csv',
                 'normal_blood_WGS.hg19_multianno.csv', ]
os.chdir(filepath)

#ANNOVAR_PATH = '/database/annovar/'
#IN_PATH = '/data/Nof1/file.vcf'
#OUT_PATH = '/data/ccbb_internal/interns/Carlo/annovar_results'

#1. Get csv file: run annovar

#num_lines = sum(1 for line in open(csv_file))
#num_lines
#utilities.run_annovar(ANNOVAR_PATH, IN_PATH, OUT_PATH)


#METHOD 1: by chunks, iteratively.
chunksize = 10000
step = 0
db_name = 'LOCAL_Variant_Prioritization_Workflow'

collection_name = ['Tumor_RNAseq_Variants', 'Tumor_Targeted_Variants',
                  'Normal_Targeted_Variants', 'Somatic_Mutect_Variants',
                  'Normal_Blood_WGS_Variants']

#Get variant list. Should always be the first step after running ANNOVAR
#open_file = myvariant_parsing_utils.VariantParsing()
#list_file = open_file.get_variants_from_vcf(vcf_file)


#Run process, export to MongoDB in-built
as_batch = annotate_batch.AnnotationMethods()
as_batch.by_chunks(VCF_FILE_NAMES[4], CSV_FILE_NAME[4], collection_name[4], db_name, step, chunksize)

#Apply filter(s).
filter_collection = MongoDB_querying.Filters(db_name, collection_name[4])

rare_cancer_variants = filter_collection.rare_cancer_variant()
rare_disease_variants = filter_collection.rare_disease_variant()
cadd_phred_high_impact_variants = filter_collection.rare_high_impact_variants()

#Create 4 output files: annotated vcf, annotated csv, filtered vcf, filtered csv
#Annotated vcf and csv, unfiltered. Will contain all info coming from annovar and myvariant

out_unfiltered_vcf_file = filepath + "/varcode_to_test/" + collection_name[4] + "unfilterd_vcf_annotated.vcf"
out_unfiltered_csv_file = filepath + "/varcode_to_test/" + collection_name[4] + "normal_targ_unfiltered_csv_annotated.csv"

rare_cancer_variants_csv = filepath + "/varcode_to_test/" + collection_name[4] + "normal_targ_rare_cancer_vars.csv"
rare_cancer_variants_vcf = filepath + "/varcode_to_test/" + collection_name[4] + "normal_targ_rare_cancer_vars.vcf"

rare_disease_variants_csv = filepath + "/varcode_to_test/" + collection_name[4] + "normal_targ_rare_disease_vars.csv"
rare_diseasw_variants_vcf = filepath + "/varcode_to_test/" + collection_name[4] + "normal_targ_rare_disease_vars.vcf"

cadd_phred_high_impact_variants_csv = filepath + "/varcode_to_test/" + collection_name[4] + "normal_targ_cadd_phred_high_impact_variants.csv"
cadd_phred_high_impact_variants_vcf = filepath + "/varcode_to_test/" + collection_name[4] + "normal_targ_cadd_phred_high_impact_variants.vcf"

in_vcf_file = filepath + "normal_blood_WGS.vqsr.vcf.gz"

#Create writer object
my_writer_1 = create_output_files.FileWriter(db_name, collection_name[4])
#Write collection to csv and vcf
my_writer_1.generate_unfiltered_annotated_csv(out_unfiltered_csv_file)
my_writer_1.generate_unfiltered_annotated_vcf(in_vcf_file, out_unfiltered_vcf_file)


#Crete writer object for filtered lists:
my_writer_2 = create_output_files.FileWriter(db_name, collection_name[4])

#cancer variants filtered files
my_writer_2.generate_annotated_csv(rare_cancer_variants, rare_cancer_variants_csv)
my_writer_2.generate_annotated_vcf(rare_cancer_variants, in_vcf_file, rare_cancer_variants_vcf)

#disease variants filtered files
my_writer_2.generate_annotated_csv(rare_disease_variants, rare_disease_variants)
my_writer_2.generate_annotated_vcf(rare_disease_variants, rare_disease_variants)

#high impact cadd_phredd filtered files
my_writer_2.generate_annotated_csv(cadd_phred_high_impact_variants, cadd_phred_high_impact_variants_csv)
my_writer_2.generate_annotated_vcf(cadd_phred_high_impact_variants, cadd_phred_high_impact_variants_vcf)


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

#Generate output files
out_vcf_file = filepath + "/Tumor_RNAseq_rare_variants_ANNOTATED_FULL.vcf"
out_csv_file = filepath + "/Tumor_RNAseq_rare_variants_ANNOTATED_FULL.csv"
in_vcf_file = filepath + "/Tumor_RNAseq_rare_variants_VCF.vcf"
create_output_files.generate_annotated_vcf(joint_list, in_vcf_file, out_vcf_file)
create_output_files.generate_annotated_csv(joint_list, out_csv_file)

#Filtering


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
list_ids = list(myvariant.get_hgvs_from_vcf(vcf_file))

list_file = open_file.get_variants_from_vcf(vcf_file)

#Name Collection & DB
collection_name = 'My_Variant_Info_Collection_Chunks'
db_name = 'My_Variant_Database'

#Run process, export to MongoDB in-built
my_variants = annotate_batch.AnnotationMethods()
myvariant_data = my_variants.myvariant_chunks(list_file, chunksize, step, collection_name, db_name)


out_vcf_file = filepath + "/Tumor_RNAseq_rare_variants_ANNOTATED_MYV_FULL.vcf"
out_csv_file = filepath + "/Tumor_RNAseq_rare_variants_ANNOTATED_MyV_FULL.csv"
in_vcf_file = filepath + "/Tumor_RNAseq_rare_variants_VCF.vcf"
create_output_files.generate_annotated_vcf(myvariant_data, in_vcf_file, out_vcf_file)
create_output_files.generate_annotated_csv(myvariant_data, out_csv_file)




########DEBUG#########
import os
collection_name = 'Test_Normal_Targeted'
db_name = 'My_Variant_Database'

#set paths
filepath = "/Volumes/Seagate Backup Plus Drive/vcf_files"
csv_file = "normal_targeted_seq.hg19_multianno.csv"
vcf_file = "normal_targeted_seq.vcf"
os.chdir(filepath)


from variantannotation import myvariant_parsing_utils
from variantannotation import csv_to_df
from variantannotation import annovar_processing
from variantannotation import utilities

open_file = myvariant_parsing_utils.VariantParsing()
list_file = open_file.get_variants_from_vcf(vcf_file)


df = csv_to_df.parse_to_df(csv_to_df.open_and_parse(csv_file))
list1 = annovar_processing.get_list_from_annovar_csv(df, list_file[0:5000])
open_file = myvariant_parsing_utils.VariantParsing()
from_myvariant = open_file.get_dict_myvariant(list_file[0:5000])
utilities.final_joint(list1, from_myvariant)
joined_list = list1

from pymongo import MongoClient
client = MongoClient()
db = client.My_Variant_Database
collection = db.Test_Normal_Targeted

all_my_data = list(collection.find({}))

chr_vars = []
location_vars_ant = []
location_vars_pos = []

for i in range(0, len(all_my_data)):
    if all_my_data[i]['Chr'] == 'chrMT':
        chr_vars.append('chrM')
    else:
        chr_vars.append(all_my_data[i]['Chr'].encode('ascii','ignore'))
    location_vars_ant.append(all_my_data[i]['Start'] + 1)
    location_vars_pos.append(all_my_data[i]['Start'] - 1)


import vcf
in_vcf_file = filepath + "/somatic_mutect_old.vcf"
vcf_output_path = "/Users/carlomazzaferro/Desktop/test.vcf"

vcf_reader = vcf.Reader(filename=in_vcf_file)
vcf_writer = vcf.Writer(open(vcf_output_path, 'w'), vcf_reader)



import itertools
import myvariant

list_ids = []
reading = vcf.Reader(open('/Volumes/Seagate Backup Plus Drive/vcf_files/normal_blood_WGS.vqsr.vcf', 'r'))

for record in itertools.islice(reading, 0, 10):
    values = ','.join(str(v) for v in record.ALT)
    list_ids.append(myvariant.format_hgvs(record.CHROM, record.POS, record.REF, record.ALT))


sl = ['AT', 'TSTEF']

values = ','.join(str(v) for v in sl)
values