
import os
import sys
import math
#quick and dirty way of importing functions
sys.path.append("/Users/carlomazzaferro/Documents/Bioinformatics Internship/Python Codes")

from utilities import convert
from utilities import final_joint
import myvariant_parsing_utils
from annovar_processing import get_list_from_annovar_csv

from mongo_DB_export import export

import annovar_subprocess
import csv_to_df
#set by user


filepath = "/Users/carlomazzaferro/Desktop/CSV to be tested"
csv_file = "Tumor_targeted_processed.csv"
vcf_file = "Tumor_targeted_seq.vqsr.vcf"

os.chdir(filepath)

chunksize = 1000
step = 0


open_file = myvariant_parsing_utils.VariantParsing()
variant_list = open_file.get_variants_from_vcf(vcf_file)

joint_list = []

while step * chunksize < len(variant_list):
    chunk_ids = variant_list[chunksize * step:chunksize * (step + 1)]

    df = csv_to_df.parse_to_df(csv_to_df.open_and_parse_chunks(csv_file, chunksize, step))

    from_annovar = get_list_from_annovar_csv(df, chunk_ids)
    open_file = myvariant_parsing_utils.VariantParsing()

    from_myvariant = open_file.get_dict_myvariant(chunk_ids)

    final_joint(from_annovar, from_myvariant)

    joined_list = from_annovar

    # From unicode to string
    joined_list = convert(joined_list)

    joint_list.append(joined_list)
    print len(joint_list)

    print 'Parsing to MongoDB ...'
    export(joined_list)
    step = step + 1
    print 'Step: {} of {}'.format(step, (len(variant_list) / 1000) - 1)


import pandas
import csv
from itertools import islice


def open_and_parse_chunks(file_name, chunksize, step):

    listoflists = []
    with open(file_name, 'rb') as csvfile:

        reader = csv.reader(csvfile, delimiter=',')
        header = next(reader)
        lines_of_interest = list(reader)[step*chunksize:(step+1)*chunksize]

        for i in lines_of_interest:
            if len(i) == len(header):
                listoflists.append(i)
            else:
                del i[4]
                listoflists.append(i)

        listoflists.insert(0, header)
    return listoflists

lol = open_and_parse_chunks(csv_file, 1000, 10)




def parse_to_df(listoflists):

    df = pandas.DataFrame(listoflists[1::], columns=listoflists[0])
    # Keep only required columns
    df = df[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.knownGene', 'Gene.knownGene', 'GeneDetail.knownGene',
             'ExonicFunc.knownGene', 'tfbsConsSites', 'cytoBand', 'targetScanS',
             'genomicSuperDups', 'nci60', 'Otherinfo']]
    df = df.replace({'.': None})  # None values are easier to deal with

    return df

