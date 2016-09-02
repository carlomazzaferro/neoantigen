import urllib2
from Bio.Blast import NCBIXML


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
        urllib2.urlopen(url)
    except urllib2.HTTPError, e:
        print(e.code)
    except urllib2.URLError, e:
        print(e.args)

    # if no errors, proceed:
    response = urllib2.urlopen(url)
    html = response.read()
    return html


def write_file(html_input, filepath):
    if isinstance(html_input, basestring):
        with open(filepath, 'w') as outfile:
            outfile.write(html_input)
    else:
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


def write_blast_filtered_output(xml_in_filepath, fasta_out, treshold):
    result_handle = open(xml_in_filepath)
    blast_records = NCBIXML.parse(result_handle)
    with open(fasta_out, 'w') as outfile:
        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    PERCENTAGE = float(hsp.identities)/alignment.length*100
                    if PERCENTAGE > treshold:
                        outfile.write('>' + alignment.accession + '\n')
                        outfile.write(hsp.sbjct + '\n')


def create_lists(fasta_file):
    with open(fasta_file) as infile:
        all_list = []
        peptide = ""
        for line in infile:

            if line.startswith('>'):
                all_list.append(line.rstrip())
            else:
                peptide = peptide + line.rstrip()
            if next(infile).startswith('>'):
                all_list.append(peptide)
                peptide = ""

    return all_list


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


