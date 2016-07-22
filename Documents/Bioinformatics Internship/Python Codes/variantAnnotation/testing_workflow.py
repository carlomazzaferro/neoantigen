import pandas
import re
import os
import sys
import myvariant
from itertools import chain
#quick and dirty way of importing functions
sys.path.append("/Users/carlomazzaferro/Documents/Bioinformatics Internship/Python Codes/")
from variantAnnotation import genotype_calling
from variantAnnotation import utilities
from variantAnnotation import csv_to_df
#set by user
filepath = "/Users/carlomazzaferro/Desktop/CSV to be tested"
file_name = "Tumor_targeted_processed.csv"
vcf_file = "Tumor_targeted_seq.vqsr.vcf"

os.chdir(filepath)


as_list = csv_to_df.open_and_parse(file_name)
df = csv_to_df.parse_to_df(as_list)

df.Chr = df.Chr.replace(to_replace='chrM', value='chrMT')
print 'Processing knownGene info ...'
utilities.split_string(df, "Func.knownGene")
utilities.split_string(df, "ExonicFunc.knownGene")

print 'Processing nci60 info ...'
df["nci60"] = utilities.to_float(df, "nci60")


df["tfbsConsSites"] = df["tfbsConsSites"].dropna().apply(utilities.cell_to_dict)


df["targetScanS"] = df["targetScanS"].dropna().apply(utilities.cell_to_dict)

print 'Processing genomicSuperDups info ...'
df["genomicSuperDups"] = df["genomicSuperDups"].dropna().apply(utilities.cell_to_dict)


print 'Processing cytoBand info ...'
df["cytoBand"] = df["cytoBand"].dropna().apply(utilities.split_cytoband)
df["cytoBand"] = df["cytoBand"].dropna().apply(utilities.lists_to_dict)


print 'Creating hgvs key ...'
# Myvariant.info offers a light-weight, conveninent method to create HGVS ID's from a vcf file. This will be of primary importance
# in order to be able to then retrieve the data from myvariant.info's website.\
list_ids = list(myvariant.get_hgvs_from_vcf(vcf_file))
expanded_list = utilities.expand_list(list_ids)
df['hgvs_key'] = pandas.Series(expanded_list)

print 'Processing genotype call info ...'
my_sample_id = df["Otherinfo"].dropna()

print 'Processing genotype call info ...'
my_sample_id = df["Otherinfo"].dropna().apply(genotype_calling.split_sample_ID)
#my_sample_split = my_sample_id
genotype_call =my_sample_id.apply(lambda x: x[-2::])


dict_split = genotype_call.apply(genotype_calling.return_dict)
df['Otherinfo'] = dict_split
df = df.rename(columns={'Otherinfo': 'Genotype_Call'})



myList2 = list(genotype_call)

genotype_call2 = genotype_call.apply(genotype_calling.split_gen_call)
dict_split = genotype_call2.apply(genotype_calling.parse_genotype)
df['Otherinfo'] = dict_split
df = df.rename(columns={'Otherinfo': 'Genotype_Call'})


df = utilities.modify_df(df)

print 'Transforming to JSON from dataFrame'
# Clean up dataframe
df_final = df.where((pandas.notnull(df)), None)