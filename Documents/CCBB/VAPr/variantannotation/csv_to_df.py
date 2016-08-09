import pandas
import csv
import itertools


def open_and_parse_chunks(file_name, chunksize, step):
    """
    Self-explanatory. The list of lists output was chosen for simplicity and easiness to parse to pandas df.

    :param file_name: name of csv file
    :param chunksize: size of chunk to be processed at a time
    :param step: counter object
    :return: list of lists coming from csv file
    """
    listoflists = []
    with open(file_name, 'rb') as csvfile:

        reader = csv.reader(csvfile, delimiter=',')
        header = next(reader)

        for i in itertools.islice(reader,step*chunksize,(step+1)*chunksize):
            if len(i) == len(header):
                listoflists.append(i)
            else:
                del i[4]
                listoflists.append(i)

        listoflists.insert(0, header)
    return listoflists


def open_and_parse(file_name):
    """Same as before with the differentail of parsing the entire file at once"""

    listoflists = []
    with open(file_name, 'rb') as csvfile:
        # tp = read_csv('large_dataset.csv', iterator=True, chunksize=1000)
        reader = csv.reader(csvfile, delimiter=',')
        header = next(reader)
        for row in reader:
            if len(row) == len(header):
                listoflists.append(row)
            else:
                del row[4]
                listoflists.append(row)

        listoflists.insert(0, header)
    return listoflists

def parse_to_df(listoflists):
    """
    :param listoflists: list of lists coming from previous step
    :return: a pandas dataframe with selected columns
    """

    df = pandas.DataFrame(listoflists[1::], columns=listoflists[0])
    # Keep only required columns

    df = df[['Chr',
             'Start',
             'End',
             'Ref',
             'Alt',
             'Func.knownGene',
             'Gene.knownGene',
             'GeneDetail.knownGene',
             'ExonicFunc.knownGene',
             'tfbsConsSites',
             'cytoBand',
             'genomicSuperDups',
             '1000g2015aug_all',
             'ESP6500si_ALL',
             'cosmic70',
             'nci60',
             'Otherinfo']]

    df = df.replace({'.': None})  # None values are easier to deal with

    return df

