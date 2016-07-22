import pandas
import myvariant
import logging
from variantAnnotation import genotype_calling
from variantAnnotation import utilities


def get_list_from_annovar_csv(df, chunk_ids):

    df.Chr = df.Chr.replace(to_replace='chrM', value='chrMT')
    df['Start'] = pandas.to_numeric(df['Start'])
    df['End'] = pandas.to_numeric(df['End'])


    print 'Processing knownGene info ...'
    utilities.split_string(df, "Func.knownGene")
    utilities.split_string(df, "ExonicFunc.knownGene")


    print 'Processing nci60 info ...'
    df["nci60"] = utilities.to_float(df, "nci60")


    print 'Processing tfbsConsSites info ...'
    df["tfbsConsSites"] = df["tfbsConsSites"].dropna().apply(utilities.cell_to_dict)


    print 'Processing targetScanS info ...'
    df["targetScanS"] = df["targetScanS"].dropna().apply(utilities.cell_to_dict)


    print 'Processing genomicSuperDups info ...'
    df["genomicSuperDups"] = df["genomicSuperDups"].dropna().apply(utilities.cell_to_dict)


    print 'Processing cytoBand info ...'
    df["cytoBand"] = df["cytoBand"].dropna().apply(utilities.split_cytoband)
    df["cytoBand"] = df["cytoBand"].dropna().apply(utilities.lists_to_dict)


    print 'Creating hgvs key ...'
    df['hgvs_key'] = pandas.Series(chunk_ids)

    print 'Processing genotype call info ...'
    my_sample_id = df["Otherinfo"].dropna().apply(genotype_calling.split_sample_ID)
    genotype_call = my_sample_id.apply(lambda x: x[-2::])
    dict_split = genotype_call.apply(genotype_calling.return_dict)
    df['Otherinfo'] = dict_split
    df = df.rename(columns={'Otherinfo': 'Genotype_Call'})


    df = utilities.modify_df(df)

    print 'Transforming to JSON from dataFrame'
    #Clean up dataframe
    df_final = df.where((pandas.notnull(df)), None)
    list_dict = df_final.T.to_dict().values()

    #Attempt to transform dataframe to dictionary
    #Set the ID to be the HGVS_ID

    print 'cleaning up...'
    for i in range(0, len(list_dict)):
        list_dict[i] = utilities.scrub_dict(list_dict[i])

    #list_filtered = []
    #for key in filtered.keys():
    #    list_filtered.append({key: filtered[key]})

    print 'Done'
    return list_dict