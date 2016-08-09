import re


#split by looking for separators: commas, colons.
def split_sample_ID(x):
    p = re.compile('[\t]+')
    x = p.split(x)
    return x


def split_gen_call(x):
    p = re.compile('[:,]+')
    x = p.split(x)
    return x


def get_dict_keys(x):
    keys_dict = x[0].split(':')
    return keys_dict


def get_info(x):
    info_dict = x[-1].split(':')
    for i in range(0, len(info_dict)):

        if ',' in info_dict[i]:
            info_dict[i] = list(info_dict[i].split(','))
            info_dict[i] = map(float, info_dict[i])
    return info_dict


def return_dict(x):
    """
    This list of functions is aimed at recreating a dictionary from the extremely messy column 'Otherinfo' that appears
    in annovar-processed csv files. The columns specifically contains genotype call data, which is in genral organized
    as follows:

    {'AD': [0.0, 267.0],
     'DP': 267,
     'GQ': 99,
     'GT': '1/1',
     'PL': [11155.0, 803.0, 0.0]},

    Where the AD, DP, GQ, GT, PL values are described at the top of the vcf file input, i.e:
    ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed"
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)"
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF

    :param x: entry of a pandas.Series as a string
    :return: dictionary
    """

    if x[0][0:2] == 'GT':
        del x[1]
    else:
        del x[0]

    my_dict = dict(zip(get_dict_keys(x), get_info(x)))
    if 'DP' in my_dict.keys():
        my_dict['DP'] = int(my_dict['DP'])

    if 'GQ' in my_dict.keys():
        my_dict['GQ'] = int(my_dict['GQ'])

    return my_dict
