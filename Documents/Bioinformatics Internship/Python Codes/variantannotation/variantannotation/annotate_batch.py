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

class AnnotationMethods(object):

    def __init___(self, variant_list, chunksize, step, file_name):

        self.file_name = file_name
        self.variant_list = variant_list
        self.chunksize = chunksize
        self.step = step

    def by_chunks(self, variant_list, chunksize, step, file_name):

        while step*chunksize < len(variant_list):

            chunk_ids = variant_list[chunksize*step:chunksize*(step+1)]

            df = csv_to_df.parse_to_df(csv_to_df.open_and_parse_chunks(file_name, chunksize, step))

            from_annovar = get_list_from_annovar_csv(df, chunk_ids)
            open_file = myvariant_parsing_utils.VariantParsing()

            from_myvariant = open_file.get_dict_myvariant(chunk_ids)

            final_joint(from_annovar, from_myvariant)
            joined_list = from_annovar

            # From unicode to string
            joined_list = convert(joined_list)

            print 'Parsing to MongoDB ...'
            export(joined_list)
            step = step + 1
            print 'Step: {} of {}'.format(step, (len(variant_list)/1000)-1)

        return joined_list, chunksize*step

    def full_file(self, variant_list, file_name):

        df = csv_to_df.parse_to_df(csv_to_df.open_and_parse(file_name))
        from_annovar = get_list_from_annovar_csv(df, variant_list)
        open_file = myvariant_parsing_utils.VariantParsing()
        from_myvariant = open_file.get_dict_myvariant(variant_list)
        final_joint(from_annovar, from_myvariant)
        joined_list = from_annovar
        joined_list = convert(joined_list)

        return joined_list

    def myvariant_chunks(self, variant_list, chunksize, step,):

        while step * chunksize < len(variant_list):

            chunk_ids = variant_list[chunksize * step:chunksize * (step + 1)]
            open_file = myvariant_parsing_utils.VariantParsing()
            from_myvariant = open_file.get_dict_myvariant(chunk_ids)
            export(from_myvariant)

        return 'Done'

    def my_variant_at_once(self, variant_list):

        open_file = myvariant_parsing_utils.VariantParsing()
        from_myvariant = open_file.get_dict_myvariant(variant_list)

        return from_myvariant




  #  def parallelized(self, variant_list, file_name):
