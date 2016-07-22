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
    info_dict = x[1].split(':')
    for i in range(0, len(info_dict)):

        if ',' in info_dict[i]:
            info_dict[i] = list(info_dict[i].split(','))
            info_dict[i] = map(float, info_dict[i])
    return info_dict


def return_dict(x):

    my_dict = dict(zip(get_dict_keys(x), get_info(x)))
    if 'DP' in my_dict.keys():
        my_dict['DP'] = int(my_dict['DP'])

    if 'GQ' in my_dict.keys():
        my_dict['GQ'] = int(my_dict['GQ'])

    return my_dict