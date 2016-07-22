import os
import sys

#quick and dirty way of importing functions
sys.path.append("/Users/carlomazzaferro/Documents/Bioinformatics Internship/Python Codes")

from utilities import convert
from utilities import final_joint
import myvariant_parsing_utils
import annotate_batch

import mongo_DB_export
import annovar_subprocess


#set by user
filepath = "/Users/carlomazzaferro/Desktop/CSV to be tested"
csv_file = "Tumor_targeted_processed.csv"
vcf_file = "Tumor_targeted_seq.vqsr.vcf"

os.chdir(filepath)



#METHOD 1: by chunks, iteratively
#TODO: fix joint_list format
chunksize = 1000
step = 0


open_file = myvariant_parsing_utils.VariantParsing()
list_file = open_file.get_variants_from_vcf(vcf_file)


as_batch = annotate_batch.AnnotationMethods()
joint_list, length = as_batch.by_chunks(list_file, chunksize, step, csv_file)

#METHOD 2: parallelize everything. Add method to class AnnotationMethods to do so.

#get variant list. Should always be the first step after running ANNOVAR
open_file = myvariant_parsing_utils.VariantParsing()
list_file = open_file.get_variants_from_vcf(vcf_file)

as_one_file = annotate_batch.AnnotationMethods()
joint_list = as_one_file.full_file(list_file, csv_file)








#METHOD A: ANNOVAR + MYVARIANT.INFO

#1. Get csv file: run annovar
annovar_subprocess.run_annovar()


"""
#ANNOVAR
filtered_annovar = get_list_from_annovar_csv(file_name, vcf_file)

#test on smaller subsample, a huge one will take a few hours to query all the data from myvariant.info
from_annovar = filtered_annovar[0:200]
"""

#MYVARIANT.INFO()
variant_list = myvariant_parsing_utils.get_variants_from_vcf(vcf_file)
from_myvariant = myvariant_parsing_utils.get_dict_myvariant(variant_list, chunksize)

#Join data
final_joint(from_annovar,from_myvariant)
joined_list = from_annovar

#From unicode to string
joined_list = convert(joined_list)



#---------------#--------------#---------------#--------------#---------------#--------------#---------------#
#Ignore
#METHOD B: Just MYVARIANT.INFO
import os
import sys

#quick and dirty way of importing functions
sys.path.append("/Users/carlomazzaferro/Documents/Bioinformatics Internship/Python Codes/")

from myvariant_methodB import annoatation_myvariant

#set by user
filepath = "/Users/carlomazzaferro/Desktop/CSV to be tested"
vcf_file = "Tumor_targeted_seq.vqsr.vcf"
os.chdir(filepath)

from_myvariant = annoatation_myvariant(vcf_file, first_variant = 1000, last_variant= 2000)


#---------------#--------------#---------------#--------------#---------------#--------------#---------------#


"MongoDB export"

from pymongo import MongoClient

client = MongoClient()

db = client.test_database

collection = db.test_variant_annotation
collection.insert_many(from_annovar)

collection.find_one({'HGVS_id': 'chrMT:g.146T>C'})


""""
DEBUG

import myvariant

list_ids = list(myvariant.get_hgvs_from_vcf(vcf_file))
mv = myvariant.MyVariantInfo()
myvariant_dot_info = mv.getvariants(list_ids)

for i in range(0, len(list_ids)):
    if 'M' in list_ids[i]:
        list_ids[i] = list_ids[i][0:4] + 'T' + list_ids[i][4::]

list_of_wrong_IDS =[]
for i in range(0, len(myvariant_dot_info)):
    if myvariant_dot_info[i].keys()[1] == 'notfound':
        if myvariant_dot_info[i].values()[0][0] != 'c':
            list_of_wrong_IDS.append(i)

        if 'M' in myvariant_dot_info[i].values()[0]:
            #myvariant_dot_info[i].values()[0] = myvariant_dot_info[i].values()[0][0:4] + 'T' + myvariant_dot_info[i].values()[0][4::]
            print myvariant_dot_info[i].values()[0]


for i in from_myvariant:
    print i.values()
from_myvariant[22]
"""