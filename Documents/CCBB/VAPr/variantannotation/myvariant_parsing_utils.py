import utilities
import itertools
import vcf
import myvariant


class VariantParsing(object):

    def __init___(self, vcf_file, variant_list):
        self.file = vcf_file
        self.variant_list = variant_list

    def get_variants_from_vcf(self, vcf_file):
        """
        Retrieves variant names from a vcf file.

        :param vcf_file: initial file containing genomic data.
        :return: a list of variants formatted according to HGVS standards
        """

        list_ids = list(myvariant.get_hgvs_from_vcf(vcf_file))
        expanded_list = utilities.expand_list(list_ids)
        variant_list = self.complete_chromosome(expanded_list)

        return variant_list

    def get_variants_from_vcf_chunk(self, vcf_file, chunksize, step):
        """
        Retrieves variant names from a LARGE vcf file.

        :param vcf_file: initial file containing genomic data.
        :return: a list of variants formatted according to HGVS standards
        """
        list_ids = []
        reading = vcf.Reader(open(vcf_file, 'r'))

        for record in itertools.islice(reading, step*chunksize, (step+1)*chunksize):
            values = ','.join(str(v) for v in record.ALT)
            list_ids.append(myvariant.format_hgvs(record.CHROM, record.POS, record.REF, values))

        expanded_list = utilities.expand_list(list_ids)
        variant_list = self.complete_chromosome(expanded_list)

        return variant_list


    def get_dict_myvariant(self, variant_list):
        """
        Function designated to place the queries on myvariant.info servers.

        :param variant_list: list of HGVS variant ID's. Usually retrived beforehand using the method get_variants_from_vcf
        from the class VariantParsing.
        :return: list of dictionaries. Each dictionary contains data about a single variant.
        """

        mv = myvariant.MyVariantInfo()
        # This will retrieve a list of dictionaries
        variant_data = mv.getvariants(variant_list, as_dataframe=False)
        variant_data = self.remove_id_key(variant_data)
        return variant_data

    @staticmethod
    def complete_chromosome(expanded_list):
        for i in range(0, len(expanded_list)):
            if 'M' in expanded_list[i]:
                one = expanded_list[i].split(':')[0]
                two = expanded_list[i].split(':')[1]
                if 'MT' not in one:
                    one = 'chrMT'
                expanded_list[i] = one + ':' + two
        return expanded_list

    @staticmethod
    def remove_id_key(variant_data):
        for i in variant_data:
            if "_id" in i:
                del i["_id"]
            if "query" in i:
                del i["query"]
        return variant_data







#-------------#-------------#-------------#-------------#-------------#-------------#-------------#-------------


"""





def get_variants_from_vcf(vcf_file):
    list_ids = list(myvariant.get_hgvs_from_vcf(vcf_file))
    expanded_list = utilities.expand_list(list_ids)
    final_list = complete_chromosome(expanded_list)

    return final_list


def complete_chromosome(list_variants):
    for i in range(0, len(list_variants)):
        if 'M' in list_variants[i]:
            one = list_variants[i].split(':')[0]
            two = list_variants[i].split(':')[1]
            if 'MT' not in one:
                one = 'chrMT'
            list_variants[i] = one + ':' + two
    return list_variants


def get_dict_myvariant(my_variants):
    # This will retrieve a list of dictionaries
    data1 = mv.getvariants(my_variants, as_dataframe=False)
    data1 = remove_id_key(data1)
    return data1

def remove_id_key(from_myvariant):
    for i in from_myvariant:
        if "_id" in i:
            del i["_id"]
        if "query" in i:
            del i["query"]
    return from_myvariant


def dict_to_list(dictionary):
    list_filtered = []
    for key in dictionary.keys():
        list_filtered.append({key: dictionary[key]})

    return list_filtered

def clean_from_annovar(list_variants):
    for i in list_variants:
        if i[0] != 'c':
            del i

    return list_variants

"""