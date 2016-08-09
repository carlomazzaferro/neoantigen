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


def open_and_parse(file_name):

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

    df = pandas.DataFrame(listoflists[1::], columns=listoflists[0])
    # Keep only required columns
    df = df[['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.knownGene', 'Gene.knownGene', 'GeneDetail.knownGene',
             'ExonicFunc.knownGene', 'tfbsConsSites', 'cytoBand', 'targetScanS',
             'genomicSuperDups', 'nci60', 'Otherinfo']]
    df = df.replace({'.': None})  # None values are easier to deal with

    return df

