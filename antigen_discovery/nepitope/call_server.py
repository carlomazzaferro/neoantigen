import sys
import requests
import os


def post_request(input_file, nmer_lenght, output_file):


    allele_supertypes = ['HLA-A*01:01',
                         'HLA-A*02:01',
                         'HLA-A*03:01',
                         'HLA-A*24:02',
                         'HLA-A*26:01',
                         'HLA-B*07:02',
                         'HLA-B*08:01',
                         'HLA-B*27:05',
                         'HLA-B*39:01',
                         'HLA-B*40:01',
                         'HLA-B*58:01',
                         'HLA-B*15:01']

    fo = open(input_file)
    data = {
        'sequence_text': fo.read(),
        'method':        'smm',
        'allele':        'HLA-A*01:01',
        'length':        nmer_lenght,
    }

    response = requests.post('http://tools-api.iedb.org/tools_api/mhci/', data=data)
    if response.status_code != 200:
        sys.exit("Error posting request to IEDB.\n%s" % response.text)

    #if response.status_code == 500:

        """
                #Retry once
        response = requests.post('http://tools-api.iedb.org/tools_api/mhci/', data=data)
        if response.status_code == 500:
            #Retry a second time
            response = requests.post('http://tools-api.iedb.org/tools_api/mhci/', data=data)
        """



    tmp_output_file = output_file + '.tmp'
    tmp_output_filehandle = open(tmp_output_file, 'w')
    tmp_output_filehandle.write(response.text)
    tmp_output_filehandle.close()
    os.replace(tmp_output_file, output_file)

    input_file.close()

