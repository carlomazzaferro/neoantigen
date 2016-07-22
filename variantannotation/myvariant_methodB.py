import myvariant


def annoatation_myvariant(vcf_file, first_variant = 0, last_variant = 1):
    """Produces a list of dictionaries containing data from myvariant. It gets the HGVS id straight from the VCF
    file and queries the remaining information. It also cleans up the dictionaries from spurious entries"""

    if first_variant != None:

        list_ids = list(myvariant.get_hgvs_from_vcf(vcf_file))
        for i in range(0, len(list_ids)):
            if 'M' in list_ids[i]:
                list_ids[i] = list_ids[i][0:4] + 'T' + list_ids[i][4::]

        mv = myvariant.MyVariantInfo()
        myvariant_dot_info = mv.getvariants(list_ids[first_variant:last_variant])

        #clean up list.
        #Some wrong queries that appear usually
        # { u'query': u'AGTGT', u'notfound': True}
        #  {u'query': u'GCACACACACA', u'notfound': True}

        list_of_wrong_IDS = []
        for i in range(0, len(myvariant_dot_info)):
            if myvariant_dot_info[i].keys()[1] == 'notfound':
                if myvariant_dot_info[i].values()[0][0] != 'c':
                    list_of_wrong_IDS.append(i)

            if myvariant_dot_info[i].keys()[1] == 'query':
                if myvariant_dot_info[i]['query'][0] != 'c':
                    list_of_wrong_IDS.append(i)



        myvariant_dot_info = [x for i,x in enumerate(myvariant_dot_info) if i not in list_of_wrong_IDS]

    return myvariant_dot_info
