import pandas
from variantannotation import genotype_calling
from variantannotation import utilities
import logging


def get_list_from_annovar_csv(df, chunk_ids):

    length = len(chunk_ids)
    #message_transform = 'Transforming data from %d variants from csv to JSON...' % length

    #logging.INFO(message_transform)

    df = df.rename(columns={'1000g2015aug_all': 'ThousandGenomeAll'})
    df.Chr = df.Chr.replace(to_replace='chrM', value='chrMT')
    df['Start'] = pandas.to_numeric(df['Start'])
    df['End'] = pandas.to_numeric(df['End'])
    df["nci60"] = utilities.to_float(df, "nci60")
    df["ThousandGenomeAll"] = utilities.to_float(df, "ThousandGenomeAll")
    df["ESP6500si_ALL"] = utilities.to_float(df, "ESP6500si_ALL")
    df["tfbsConsSites"] = df["tfbsConsSites"].dropna().apply(utilities.cell_to_dict)

    utilities.split_string(df, "Func.knownGene")
    utilities.split_string(df, "ExonicFunc.knownGene")

    df["genomicSuperDups"] = df["genomicSuperDups"].dropna().apply(utilities.cell_to_dict)
    df["cytoBand"] = df["cytoBand"].dropna().apply(utilities.split_cytoband)
    df["cytoBand"] = df["cytoBand"].dropna().apply(utilities.lists_to_dict)
    df['hgvs_key'] = pandas.Series(chunk_ids)

    my_sample_id = df["Otherinfo"].dropna().apply(genotype_calling.split_sample_ID)
    genotype_call = my_sample_id.apply(lambda x: x[-3::])
    dict_split = genotype_call.apply(genotype_calling.return_dict)
    df['Otherinfo'] = dict_split
    df = df.rename(columns={'Otherinfo': 'Genotype_Call'})
    df = utilities.modify_df(df)


    #Clean up dataframe
    df_final = df.where((pandas.notnull(df)), None)
    list_dict = df_final.T.to_dict().values()

    #Attempt to transform dataframe to dictionary
    #Set the ID to be the HGVS_ID

    for i in range(0, len(list_dict)):
        list_dict[i] = utilities.scrub_dict(list_dict[i])

    #list_filtered = []
    #for key in filtered.keys():
    #    list_filtered.append({key: filtered[key]})
    return list_dict


#REFACTOR THIS CODEEEEE AASAP


def get_df_from_annovar_csv(df, chunk_ids):

    df = df.rename(columns={'1000g2015aug_all': 'ThousandGenomeAll'})
    df.Chr = df.Chr.replace(to_replace='chrM', value='chrMT')
    df['Start'] = pandas.to_numeric(df['Start'])
    df['End'] = pandas.to_numeric(df['End'])
    df["nci60"] = utilities.to_float(df, "nci60")
    df["ThousandGenomeAll"] = utilities.to_float(df, "ThousandGenomeAll")
    df["ESP6500si_ALL"] = utilities.to_float(df, "ESP6500si_ALL")
    df["tfbsConsSites"] = df["tfbsConsSites"].dropna().apply(utilities.cell_to_dict)

    utilities.split_string(df, "Func.knownGene")
    utilities.split_string(df, "ExonicFunc.knownGene")

    #df["targetScanS"] = df["targetScanS"].dropna().apply(utilities.cell_to_dict)

    df["genomicSuperDups"] = df["genomicSuperDups"].dropna().apply(utilities.cell_to_dict)
    df["cytoBand"] = df["cytoBand"].dropna().apply(utilities.split_cytoband)
    df["cytoBand"] = df["cytoBand"].dropna().apply(utilities.lists_to_dict)
    df['hgvs_key'] = pandas.Series(chunk_ids)

    my_sample_id = df["Otherinfo"].dropna().apply(genotype_calling.split_sample_ID)
    genotype_call = my_sample_id.apply(lambda x: x[-2::])
    dict_split = genotype_call.apply(genotype_calling.return_dict)
    df['Otherinfo'] = dict_split
    df = df.rename(columns={'Otherinfo': 'Genotype_Call'})


    #Clean up dataframe
    df = utilities.modify_df(df)
    df_final = df.where((pandas.notnull(df)), None)

    return df_final

#after myvariant data has been obtained
#def join_df_chunks(df):
