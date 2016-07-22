
import pandas
import re
from itertools import chain
import collections


def split_string(dataframe, column):
    """ General Split String Function"""
    dataframe[column] = dataframe[column].str.split(pat=';', expand=False)


def to_float(dataframe, column):
    """General Function to return floats"""
    dataframe[column] = dataframe[column].dropna().astype(float)
    dataframe[column] = dataframe[column].where(pandas.notnull(dataframe[column]), None)
    return dataframe[column]


def cell_to_dict(s):
    """When separators are '=' and ';'"""
    as_dict =  dict(item.split("=") for item in s.split(";"))
    as_dict["Score"] = float(as_dict["Score"])
    return as_dict

def split_cytoband(x):
    """Specific for splitting cytoband"""

    letters = ['X', 'Y', 'p', 'q']     #Possible letters: Chrm X,Y, or arm regions p (short), or q (long)
    spliced = re.split('(\D+)', x)     #Split letters and numbers
    spliced = filter(None, spliced)
    if any(letter in spliced[0] for letter in letters): spliced[0] = map(None, spliced[0])

    if type(spliced[0]) == type(letters):
        first = spliced[0][0]
        second = spliced[0][1]
        spliced[0] = first
        spliced.insert(1,second)
    return spliced

def lists_to_dict(x):
    """Specific for cytoband. See tests for clearer visualization of the modifications being done"""

    cyto_dict = {'Chromosome': None, 'Arm': None, 'Region': None, 'Band': None, 'Sub_Band': None}
    cyto_dict['Chromosome'] = x[0]
    del x[0]
    cyto_dict['Arm'] = x[0]
    del x[0]

    if '.' in x:
        cyto_dict['Sub_Band'] = int(x[-1])
        del x[-2::]
        if len(x[0]) == 2:
            cyto_dict['Region'] = int(x[0][0])
            cyto_dict['Band'] = int(x[0][1])
        else:
            cyto_dict['Region'] = int(x[0])
        clean_dict = {k: v for k, v in cyto_dict.iteritems() if v is not None}
        return clean_dict

    else:
        if len(x[0]) == 2:
            cyto_dict['Region'] = int(x[0][0])
            cyto_dict['Band'] = int(x[0][1])
        else:
            cyto_dict['Region'] = int(x[0])
        clean_dict = {k: v for k, v in cyto_dict.iteritems() if v is not None}
        return clean_dict


def scrub_dict(d):
    """Clean up dict - after converting everything to dict"""
    if type(d) is dict:
        return dict((k, scrub_dict(v)) for k, v in d.iteritems() if v and scrub_dict(v))
    else:
        return d


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

def modify_df(df):
    """Replace dot notation by underscore notation in columns. This will allow proper and easy parsing
    to the MongoDB database (dot notation is not allowed)"""
    cols = df.columns
    cols = cols.map(lambda df: df.replace('.', '_') if isinstance(df, (str, unicode)) else df)
    df.columns = cols
    return df


def expand_list(list_ids):
    """Expand list of HGVS IDS when there are duplicate ids as the same list item"""
    for i in range(0, len(list_ids)):
        if ',' in list_ids[i]:
            old = list_ids[i].split(',')[0]
            new = old.split('s')[0] + list_ids[i].split(',')[1]
            list_ids[i] = old
            list_ids.insert(i + 1, new)

    return list_ids

def convert(data):
    """Given a data structure of unicode-type variables, convert to string each piece of data, recursively"""
    try:
        if isinstance(data, basestring):
            return str(data)
        elif isinstance(data, collections.Mapping):
            return dict(map(convert, data.iteritems()))
        elif isinstance(data, collections.Iterable):
            return type(data)(map(convert, data))
        else:
            return data
    except:
        pass

"""
def clean_dict(list_filtered):
    #Function for removing bad keys and introducing the MT notation for chromosomes.
    for i in range(0, len(list_filtered)):
        if 'M' in list_filtered[i]['hgvs_key']:
            one = list_filtered[i]['hgvs_key'].split(':')[0]
            two = list_filtered[i]['hgvs_key'].split(':')[1]
            if 'MT' not in one:
                one = 'chrMT'

            list_filtered[i]['hgvs_key'] = one + ':' + two

    return list_filtered
"""



def final_joint(list1, list2):
    print 'Joining lists ...'
    for i in range(0, len(list2)):
        list1[i].update(list2[i])
    print 'DONE'
    return list2

def get_paths(PATHS_file):
    with open('/Users/carlomazzaferro/Documents/Bioinformatics Internship/Python Codes/variantAnnotation/PATHS.txt',
              'r') as f:
        content = f.readlines()
        annovar_path = content[1].split('=')[1].strip()
        input_vcf_path = content[2].split('=')[1].strip()
        output_csv_path = content[3].split('=')[1].strip()


        return annovar_path, input_vcf_path, output_csv_path

#def get_chunks(chunksize, variant_list):
 #   while (len(variant_list) > 1000):

