import myvariant
import pandas


#Myvariant.info offers a light-weight, conveninent method to create HGVS ID's. This will be of primary importance
#in order to be able to then retrieve the data from myvariant.info's website.


def annovar_to_hgvs(dataframe):
    chrom, pos, ref, alt = dataframe[1], dataframe[2],dataframe[3], dataframe[4]
    return myvariant.format_hgvs(chrom, pos, ref, alt)


#Needed because myvariants' function only works for non None values: a deletion is detected only if ref[0] = alt
#or alt[0] = ref. For instance: ref:ACG, alt:A. In our df, such deletion is represented instead as: ref: CG, alt: None
#Hence the padding with an arbitrary letter, in this case 'c'.
def padding(dataframe):
    to_pad =  dataframe[['Ref','Alt']]
    format = lambda x: 'c' if x == None else 'c' + x
    to_pad = to_pad[pandas.isnull(to_pad).any(axis=1)].applymap(format)
    #print to_pad
    return to_pad

