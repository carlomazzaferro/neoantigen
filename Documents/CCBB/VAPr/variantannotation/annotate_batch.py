from utilities import convert
from utilities import final_joint
import myvariant_parsing_utils
from annovar_processing import get_list_from_annovar_csv
from mongo_DB_export import export
import csv_to_df
import logging


class AnnotationMethods(object):


    @staticmethod
    def by_chunks(vcf_file, csv_file, collection_name, db_name, step, chunksize):

        """
        Export data to MongoDB by chunks, iteratively. It retrives annovar data from a csv file and myvariant data from
        their web-server, joins it, and exports it.
        This method is well-fitted for large files. Only the 1000 documents are held in memory and processed at a time,
        instead of attempting to parse and process an entire csv file at once.
        As soon as the method is run, the data will automatically be stored to it. Database and collection name should
        be specified, and there must be a running MongoDB connection. The script will set up a client to communicate
        between python (through pymongo) and the the database.

        :param vcf_file: from where the variant ids will be retrieved
        :param csv_file: name of the collection to which store the data
        :param chunksize: set at 1000 usually for simplicity, since the queries to myvariant.info also have a cap of
        1000 variants per query
        :param step: just a counter object

        :param collection_name:
        :param db_name: name of the database to which store the collection
        :return: nothing, simply exports data to MongoDB
        """

        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(message)s',
                            datefmt='%a, %d %b %Y %H:%M:%S')



        logging.info('calculating length of file...')

        variant_list = [0, 1] * chunksize

        num_lines = sum(1 for line in open(csv_file))
        seconds = 67*(num_lines/10000)
        m, s = divmod(seconds, 60)
        h, m = divmod(m, 60)

        time_message = 'File is {} lines long. This will take approximately {} hours, {} ' \
                       'minutes to run.'.format(num_lines, h, m)
        logging.info(time_message)

        while len(variant_list) >= chunksize:

            logging.getLogger("urllib3").setLevel(logging.WARNING)
            logging.getLogger("requests").setLevel(logging.WARNING)

            open_file = myvariant_parsing_utils.VariantParsing()
            variant_list = open_file.get_variants_from_vcf_chunk(vcf_file, chunksize, step)

            step_message = 'Step: {} of {}'.format(step+1, (num_lines / chunksize)+1)
            logging.info(step_message)

            df = csv_to_df.parse_to_df(csv_to_df.open_and_parse_chunks(csv_file, len(variant_list), step))
            from_annovar = get_list_from_annovar_csv(df, variant_list)
            from_myvariant = open_file.get_dict_myvariant(variant_list)

            final_joint(from_annovar, from_myvariant)
            joined_list = from_annovar

            # From unicode to string
            joined_list = convert(joined_list)

            logging.info('Parsing to MongoDB ...')
            export(joined_list, collection_name, db_name)
            step += 1

        return 'Finished!'

    @staticmethod
    def full_file(variant_list, file_name):
        """
        Retrives annovar data from a csv file and myvariant data from their web-server and return a list of documents
        containing info from both services. All done at once.
        This method is well-fitted for smaller files.

        :param variant_list: list of HGVS variant ID's. Usually retrived beforehand using the method get_variants_from_vcf
        from the class VariantParsing.
        :param file_name: name of csv_file to be used (coming from annovar's output)
        :return:
        """

        df = csv_to_df.parse_to_df(csv_to_df.open_and_parse(file_name))
        from_annovar = get_list_from_annovar_csv(df, variant_list)

        open_file = myvariant_parsing_utils.VariantParsing()
        from_myvariant = open_file.get_dict_myvariant(variant_list)

        final_joint(from_annovar, from_myvariant)
        joined_list = from_annovar
        joined_list = convert(joined_list)

        print 'Finished!'
        return joined_list

    @staticmethod
    def myvariant_chunks(variant_list, chunksize, step, collection_name, db_name):
        """
        Export data to MongoDB by chunks, iteratively. It retrives myvariant data from their web-server solely.
        Annovar, and a processed csv file are not required.
        This method is well-fitted for large files. Only the 1000 documents are held in memory and processed at a time,
        instead of attempting to hold in memory data from all the variants at the same time.
        As soon as the method is run, the data will automatically be stored to it. Database and collection name should
        be specified, and there must be a running MongoDB connection. The script will set up a client to communicate
        between python (through pymongo) and the the database.

        :param variant_list: list of HGVS variant ID's. Usually retrived beforehand using the method get_variants_from_vcf
        from the class VariantParsing.
        :param chunksize: set at 1000 usually for simplicity, since the queries to myvariant.info also have a cap of
        1000 variants per query
        :param step: just a counter object
        :param collection_name: name of the collection to which store the data
        :param db_name: name of the database to which store the collection
        :return: nothing, simply exports data to MongoDB
        """

        while step * chunksize < len(variant_list):

            chunk_ids = variant_list[chunksize * step:chunksize * (step + 1)]

            open_file = myvariant_parsing_utils.VariantParsing()
            from_myvariant = open_file.get_dict_myvariant(chunk_ids)

            export(from_myvariant, collection_name, db_name)
            step = step + 1

            print 'Step: {} of {}'.format(step, (len(variant_list) / chunksize) + 1)

        return 'Finished!'

    @staticmethod
    def my_variant_at_once(variant_list):
        """
        Retrives myvariant data from their web-server and return a list of documents. All done at once.
        This method is well-fitted for smaller files.

        :param variant_list: list of HGVS variant ID's. Usually retrived beforehand using the method get_variants_from_vcf
        from the class VariantParsing.
        :return: list of documents containing the joint data, can be exporter easily later on
        """
        open_file = myvariant_parsing_utils.VariantParsing()
        from_myvariant = open_file.get_dict_myvariant(variant_list)

        print 'Finished!'
        return from_myvariant
