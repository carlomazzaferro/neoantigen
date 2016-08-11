import subprocess
import re
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class NetMHCIIPanPredictor(Predictor):
    """netMHCIIpan predictor"""

    def __init__(self, data=None):
        Predictor.__init__(self, data=data)
        self.name = 'netmhciipan'
        self.colnames = ['pos','HLA','peptide','Identity','Pos','Core',
                         '1-log50k(aff)','Affinity','Rank']
        self.scorekey = '1-log50k(aff)'
        self.cutoff = .426
        self.operator = '>'
        self.rankascending = 0

    def readResult(self, res):
        """Read raw results from netMHCIIpan output"""

        data=[]
        res = res.split('\n')[19:]
        ignore=['Protein','pos','']
        for r in res:
            if r.startswith('-'): continue
            row = re.split('\s*',r.strip())[:9]
            if len(row)!=9 or row[0] in ignore:
                continue
            data.append(dict(zip(self.colnames,row)))
        return data

    def prepareData(self, df, name):
        """Prepare netmhciipan results as a dataframe"""

        df = df.convert_objects(convert_numeric=True)
        #df = df.apply(pd.to_numeric)#, errors='ignore')
        df['name'] = name
        df.rename(columns={'Core': 'core','HLA':'allele'}, inplace=True)
        df = df.drop(['Pos','Identity','Rank'],1)
        df = df.dropna()
        df['allele'] = df.allele.apply( lambda x: self.convert_allele_name(x) )
        self.getRanking(df)
        self.data = df
        return

    def runSequence(self, seq, length, allele, overlap=1):
        """Run netmhciipan for a single sequence"""

        seqfile = createTempSeqfile(seq)
        cmd = 'netMHCIIpan -s -length %s -a %s -f %s' %(length, allele, seqfile)
        #print cmd
        temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        rows = self.readResult(temp)
        df = pd.DataFrame(rows)
        return df

    def predict(self, sequence=None, peptides=None, length=11, overlap=1,
                    allele='HLA-DRB1*0101', name='',
                    pseudosequence=None):
        """Call netMHCIIpan command line"""

        #assume allele names are in standard format HLA-DRB1*0101
        try:
            allele = allele.split('-')[1].replace('*','_')
        except:
            print('invalid allele')
            return
        allele = allele.replace(':','')
        if peptides != None:
            res = pd.DataFrame()
            for p in peptides:
                temp = self.runSequence(p, len(p), allele)
                res = res.append(temp,ignore_index=True)
        else:
            res = self.runSequence(sequence, length, allele, overlap)
        if len(res)==0:
            return res
        self.prepareData(res, name)
        #print self.data[self.data.columns[:7]][:5]
        return self.data

    def getAlleles(self):
        """Get available alleles"""

        cmd = 'netMHCIIpan -list'
        try:
            temp = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
        except:
            print('netmhciipan not installed?')
            return []
        alleles=temp.split('\n')[34:]
        #print sorted(list(set([getStandardmhc1Name(i) for i in alleles])))
        return alleles

    def convert_allele_name(self, r):
        """Convert allele names to internally used form"""

        if not r.startswith('HLA'):
            return 'HLA-'+r.replace(':','')
        else:
            return r.replace(':','')


def createTempSeqfile(sequences, seqfile='tempseq.fa'):

    if isinstance(sequences, str):
        sequences=[sequences]
    out = open(seqfile, 'w')
    i=1
    for seq in sequences:
        SeqIO.write(SeqRecord(Seq(seq),id='temp%s'%i,
                    description='temp'), out, 'fasta')
        i+=1
    out.close()
    return seqfile


