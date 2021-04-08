import subprocess
import urllib.request
import xml.etree.ElementTree as ET
from sys import argv
import os

class accession:
    pid = ''
    paired = False
    runinfo = {}
    biosample = {}
    fastq = ''
    
    def __init__(self,pid):
        self.pid = pid
        
    def getRunInfo(self):
        url='https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&amp;db=sra&amp;rettype=runinfo&amp;term=' + self.pid
        u = urllib.request.urlopen(url)
        x = u.readlines()
        x0 = x[0].decode().strip().split(',')
        x1 = x[1].decode().strip().split(',')
        for i in range(0, len(x0)):
            self.runinfo[x0[i]] = x1[i]

    def getBioSample(self):
        if not self.runinfo:
            getRunInfo()
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=' + self.runinfo['BioSample']
        u = urllib.request.urlopen(url)
        tree = ET.parse(u)
        root = tree.getroot()
        for attr in root[0][5]:
            self.biosample[attr.attrib['attribute_name']] = attr.text
        bkeys = ['ecotype', 'tissue', 'dev_stage', 'sample_type', 'treatment', 'geo_loc_name',
                 'genotype', 'age']
        for k in bkeys:
            if not k in self.biosample:
                self.biosample[k]=''

    def getData(self):
        if not self.runinfo:
            getRunInfo()
        cmd = ['fastq-dump']
        if self.runinfo['LibraryLayout'] == 'PAIRED':
            cmd.append('--split-spot')
        cmd.append(self.pid)
        subprocess.run(cmd)
        self.fastq = self.pid + '.fastq'

    def sortmerna(self):
        if not self.fastq:
            getData()
        cmd = ['sortmerna','--ref',os.environ['SORTMERNADB'],'-a','12','--fastx','--log',
               '--aligned',self.pid + '_aligned', '--other', self.pid + '_other',
               '--reads',self.fastq]
        if self.runinfo['LibraryLayout'] == 'PAIRED':
            cmd.append('--paired_in')
        print(cmd)
        subprocess.run(cmd)
        if self.runinfo['LibraryLayout'] == 'PAIRED':
            cmd = ['unmerge-paired-reads.sh', self.pid + '_other.fastq',
                   self.pid + '_other_1.fq',  self.pid + '_other_2.fq']
        else:
            cmd = ['mv', self.pid + '_other.fastq', self.pid + '_other.fq']
        subprocess.run(cmd)

    def trimmomatic(self):
        PE3=os.environ['TRIMMOMATIC_HOME'] + '/adapters' + '/TruSeq3-PE.fa'
        PE32=os.environ['TRIMMOMATIC_HOME'] + '/adapters' + '/TruSeq3-PE-2.fa'
        SE3=os.environ['TRIMMOMATIC_HOME'] + '/adapters' + '/TruSeq3-SE.fa'
        PE2=os.environ['TRIMMOMATIC_HOME'] + '/adapters' + '/TruSeq2-PE.fa'
        SE2=os.environ['TRIMMOMATIC_HOME'] + '/adapters' + '/TruSeq2-SE.fa'
        EXE=os.environ['TRIMMOMATIC_HOME'] + '/trimmomatic-0.36.jar'
        if self.runinfo['LibraryLayout'] == 'PAIRED':
            cmd = ['java','-jar',EXE, 'PE','-threads','12', self.pid + '_other_1.fq', self.pid + '_other_2.fq',
                   self.pid + '_other_trimmed_1.fq', self.pid + '_other_orphans_1.fq',
                   self.pid + '_other_trimmed_2.fq', self.pid + '_other_orphans_2.fq',
                   'ILLUMINACLIP:'+PE3+':2:30:10', 'ILLUMINACLIP:'+PE32+':2:30:10',
                   'ILLUMINACLIP:'+PE2+':2:30:10', 'LEADING:3','TRAILING:3','SLIDINGWINDOW:4:20',
                   'MINLEN:36']
        else:
            cmd = ['java','-jar',EXE, 'SE','-threads','12', self.pid + '_other.fq', self.pid + '_other_trimmed.fq',
                   'ILLUMINACLIP:'+SE3+':2:30:10', 'ILLUMINACLIP:'+SE2+':2:30:10',
                   'LEADING:3','TRAILING:3','SLIDINGWINDOW:4:20', 'MINLEN:36']
        subprocess.run(cmd)

    def salmon(self):
        if self.runinfo['LibraryLayout'] == 'PAIRED':
            cmd = ['salmon','quant','-l','IU','-1',self.pid + '_other_trimmed_1.fq',
                   '-2',self.pid + '_other_trimmed_2.fq','-i',
                   '/mnt/picea/storage/reference/Arabidopsis-thaliana/TAIR10/indices/salmon/TAIR10_mrna_20101214_updated_new.inx',
                   '--output', self.pid]
        else:
            cmd = ['salmon','quant','-l','U','-r',self.pid + '_other_trimmed.fq','-i',
                   '/mnt/picea/storage/reference/Arabidopsis-thaliana/TAIR10/indices/salmon/TAIR10_mrna_20101214_updated_new.inx',
                   '--output', self.pid]
        subprocess.run(cmd)
        
    def printTabular(self):
        r = self.runinfo
        b = self.biosample
        print(self.pid,
              r['Experiment'],r['Submission'],r['Sample'],r['BioSample'],r['BioProject'],r['LibraryLayout'],r['Platform'],
              r['spots'],r['avgLength'],r['bases'],r['InsertSize'],r['InsertDev'],r['spots_with_mates'],r['LibrarySource'],
              r['LibraryStrategy'],r['LibrarySelection'],r['TaxID'],r['size_MB'],r['ReleaseDate'],b['ecotype'],b['tissue'],
              b['dev_stage'],b['sample_type'],b['treatment'],b['geo_loc_name'],b['genotype'],b['age'],sep='\t')

    def cleanup(self):
        if self.runinfo['LibraryLayout'] == 'PAIRED':
            torm = [self.pid + '_other_trimmed_1.fq', self.pid + '_other_trimmed_2.fq',
                    self.pid + '_other_orphans_1.fq', self.pid + '_other_orphans_2.fq',
                    self.pid + '_other.fastq', self.pid + '_aligned.fastq', self.pid + '.fastq',
                    self.pid + '_other_1.fq', self.pid + '_other_2.fq']
        else:
             torm = [self.pid + '_other_trimmed.fq',
                     self.pid + '_other.fq', self.pid + '_aligned.fastq', self.pid + '.fastq']
        for f in torm:
            os.remove(f)
        subprocess.run(['mv',self.pid + '_aligned.log',self.pid])
        
if __name__ == '__main__':
    
    x = accession(argv[1])
    x.getRunInfo()
    print(x.runinfo)
    x.getBioSample()
    print(x.biosample)
    #x.printTabular()
    x.getData()
    x.sortmerna()
    x.trimmomatic()
    x.salmon()
    x.cleanup()
