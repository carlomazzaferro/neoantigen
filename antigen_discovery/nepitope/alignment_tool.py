from urllib.request import urlopen
from urllib.error import HTTPError
from urllib.error import URLError

from Bio.Blast import NCBIXML
import errno
import os

filepath = '/Users/carlomazzaferro/Desktop/BINF_rand_files/CAS9_stuff/test.fasta'
ID = 'Q99ZW2'
website = 'uniprot'




def fetch_protein(filepath, website=None, ID=None, fetch_fasta=False):
    dict_websites = {'uniprot': 'http://www.uniprot.org/uniprot/',
                     'refseq': 'http://www.ncbi.nlm.nih.gov/protein/'
                     # TODO: add other possibilities
                     }

    if website not in dict_websites:
        raise ValueError("Website invalid or not supported")

    if fetch_fasta is True:
        url = dict_websites[website] + ID + '.fasta'

    fasta = handle_request(url)
    write_file(fasta, filepath)

    return 'File written to ' + filepath


def handle_request(url):
    try:
        urlopen(url)
    except HTTPError as e:
        print(e.code)
    except URLError as e:
        print(e.args)

    # if no errors, proceed:
    response = urlopen(url)
    html = response.read()
    return html


def write_file(html_input, filepath):
    if isinstance(html_input, bytes):
        html_input = html_input.decode("utf-8")
    if isinstance(html_input, str):
        with open(filepath, 'w') as outfile:
            outfile.write(html_input)
    else:
        print(type(html_input))
        raise ValueError("HTML returned bad format")


def viz_blast(xml_filepath):
    result_handle = open(xml_filepath)
    blast_records = NCBIXML.parse(result_handle)
    for record in blast_records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 0.1:
                    print('****Alignment****')
                    print('sequence:', alignment.title)
                    print('length:', alignment.length)
                    print(hsp.query[0:75] + '...')
                    print(hsp.match[0:75] + '...')
                    print(hsp.sbjct[0:75] + '...')
                    print('identities: {}'.format(hsp.identities))
                    print('PERCENTAGE: {}'.format(float(hsp.identities)/alignment.length))


def viz_blast_2(xml_filepath):
    result_handle = open(xml_filepath)
    blast_records = NCBIXML.parse(result_handle)
    i = 0
    for record in blast_records:
        for dat in record.descriptions:

            print (dat.title + '\n')
            i += 1
            print (i)
            """
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < 0.1:
                        print('****Alignment****')
                        print('sequence:', record.posted_date)
                        print('length:', alignment.length)
                        print(hsp.query[0:75] + '...')
                        print(hsp.match[0:75] + '...')
                        print(hsp.sbjct[0:75] + '...')
                        print('identities: {}'.format(hsp.identities))
                        print('PERCENTAGE: {}'.format(float(hsp.identities)/alignment.length))
            """


def write_blast_filtered_output(xml_in_filepath, fasta_out, treshold, extra_data_file):
    result_handle = open(xml_in_filepath)
    blast_records = NCBIXML.parse(result_handle)
    with open(fasta_out, 'w') as outfile, open(extra_data_file, 'w') as outfile_2:
        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    PERCENTAGE = float(hsp.identities)/alignment.length*100
                    if PERCENTAGE > treshold:
                        outfile.write('>' + alignment.accession + '\n')
                        outfile.write(hsp.sbjct + '\n')
                        outfile_2.write('>' + alignment.accession + '\n')
                        outfile_2.write(alignment.title + '\n')
                        outfile_2.write(str(alignment.length) + '\n')
                        outfile_2.write(str(hsp.identities) + '\n')


def create_lists(fasta_file):
    with open(fasta_file) as infile:
        all_list = []
        peptide = ""
        for line in infile:

            if line.startswith('>'):
                all_list.append(line.rstrip())
            else:
                peptide += line.rstrip()
            if next(infile).startswith('>'):
                all_list.append(peptide)
                peptide = ""

    return all_list


def create_single_fasta(fasta_file, dir_name):
    with open(fasta_file) as infile:
        all_list = []
        peptide = ""
        for line in infile:

            if line.startswith('>'):
                all_list.append(line.rstrip())
            else:
                peptide += line.rstrip()
                all_list.append(peptide)
                peptide = ""

    mkdir_p(os.path.dirname(fasta_file) + '/' + dir_name)
    for i in range(len(all_list)):
        with open('fasta_out_single_prot_%i.fasta' %i, 'w') as fasta:
            fasta.write(all_list[i] + '/n')
            fasta.write(all_list[i+1])


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def make_windows(fasta_file, out_fasta_path, nmers=None):
    lines = create_lists(fasta_file)
    length = len(lines[1])
    for j in nmers:
        for i in range(0, (length / j) + 1):
            with open(out_fasta_path + 'nmerized_%i_%i' % (j, i), 'w') as outfile:
                for line in lines:
                    if '>' in line:
                        outfile.write(line + '\n')
                    else:
                        outfile.write(line[i:i + j] + '\n')


