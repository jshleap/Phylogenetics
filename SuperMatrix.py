#!/usr/bin/python
"""
SuperMatrix Copyright (C) 2012 Jose Sergio Hleap

This  script will go over the fasta files created by SeqFetcher (or at least in the same format),
will create new fasta files with consensus sequences per species and transform them to nexus
using Joseph Hughes (University of Glasgow) Consensus.pl and Fasta2Nexus.pl, respectively.
by , concatenate all the sequences using the 'NEXUS' module of biopython ('NEXUS:
An extensible file format for systematic information' Maddison, Swofford, Maddison. 1997. Syst. Biol.
46(4):590-621, and will launch a RAxML (Stamatakis, A. (2006). RAxML-VI-HPC: maximum likelihood-based 
phylogenetic analyses with thousands of taxa and mixed models. Bioinformatics, 22(21), 2688-2690.)
either in fester (a private cluster) or locally.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: jshleap@squalus.org

Requires:
1) Consensus.pl, which in turn requires BioPerl libraries. Available at http://www.vcu.edu/csbc/bnfo601/Scenarios/Blast/consensus.pl
2) Biopython
3) Fasta2Nexus.pl, which in turn requires BioPerl libraries
4) RAxML, available at http://sco.h-its.org/exelixis/web/software/raxml/index.html
5) perl available at http://raven.iab.alaska.edu/~ntakebay/teaching/programming/perl-scripts/fasta2nexus.pl
6) A utils folder available at https://github.com/AlexSafatli/LabBlouinTools

Points 1,3 should be in the same folder of the selected fasta files. The requirement 6 should be either in the same folder
or in the python path

"""

import sys,os
from utils.FASTAnet import FASTAstructure as F
from properAlignment import if_prot
from subprocess import Popen, PIPE
from glob import glob
from Bio.Nexus import Nexus
import cPickle as pickle


class consensus:
    def __init__(self,prefix):
        self.F= F(prefix+'.fasta')
        self.prefix=prefix
        self.species=[]
        self.fsp=self.F.sequences.keys()
        self.sps={}
        self.typ=None
        self.main()
    
    def get_sps(self):
        for k in self.fsp:
            bl= k.split('_')
            if 'cf' in bl:
                sp='_'.join(bl[:3])
            else:
                sp='_'.join(bl[:2])
            if sp not in self.species:
                self.species.append(sp)

    def consense_by_sp(self):
        for s in self.species:
            self.sps[s]=[]
            for k in self.fsp:
                if s in k:
                    self.sps[s].append(k)
    def fst2fas(self):
        out=open(self.prefix+'.fas','w')
        fst=glob('*.fst')
        for fs in fst:
            out.write('>'+fs[:-4]+'\n')
            fo=open(fs).read().split('\n')
            if fo[0] == '':
                out.write('\n'.join(fo[2:]))
            else:
                out.write('\n'.join(fo[1:]))
        os.system('rm *.fst')
        
    def Consensus(self):
        typ=None
        for s in self.sps:
            if 'temp.fasta' in os.listdir(os.getcwd()):
                print 'Temporary file not removed'
            self.F.writeSequences('temp.fasta',self.sps[s])
            p=if_prot('temp.fasta')
            if p:
                print '\t is protein:'
                pr=Popen('perl Consensus.pl -in %s -out %s '%('temp.fasta',s+'.fst'),shell=True)
                pr.wait()
                typ='prot'
            else:
                print '\t is either nucleotide or RNA:'
                pr=Popen('perl Consensus.pl -in %s -out %s -iupac'%('temp.fasta',s+'.fst'),shell=True)
                pr.wait()
                typ='nuc'
            r=Popen('rm temp.fasta',shell=True)
            r.wait()
        self.fst2fas()
        n=Popen('perl Fasta2Nexus.pl %s.fas %s.nex'%(self.prefix,self.prefix),shell=True)
        n.wait()
        self.typ=typ
    
    def main(self):
        self.get_sps()
        self.consense_by_sp()
        self.Consensus()
    #    convert2nexus()
        
def Concatenate(prefix):
    file_list = glob('*.nex')   
    nexi =  [(fname, Nexus.Nexus(fname)) for fname in file_list]
    combined = Nexus.combine(nexi)
    combined.write_nexus_data(filename=open('btCOMBINED.nex', 'w'))
    combined.export_phylip(prefix+'.phy')
        
        
def run_prottest(aln):
    cwd= os.getcwd()
    cp = Popen('cp %s $PROTTEST_HOME'%(aln),shell=True)
    cp.wait()
    pthome=os.environ.get('PROTTEST_HOME')
    os.chdir(pthome)
    print 'Running prottest in %s file'%(aln)
    pt = Popen('java -jar prottest-3.2.jar -i %s -Dayhoff -DCMut -JTT -MtREV -WAG -RtREV -CpREV '\
               '-VT -Blosum62 -MtMam -all-distributions -AICC'%(aln),shell=True, stdout=PIPE,
               stderr=PIPE)
    o,e = pt.communicate()
    print e
    model, gamma, inv, frec = parse_prottest(o)
    os.remove(aln)
    os.chdir(cwd)
    return model, gamma, inv, frec

def parse_prottest(stdout):
    gamma=False
    inv=False
    frec=False
    query='Best model according to AICc:'
    m=stdout[stdout.find(query)+29:stdout.find('Sample Size:')].strip()
    m=m.split('+')
    model=m[0].upper()
    if 'G' in m[1:]:
        gamma=True
    elif 'I' in m[1:]:
        inv=True
    elif 'F' in m[1:]:
        frec=True
    return model, gamma, inv, frec

def parse_nexus(nexusfile):
    parts=[]
    f=open(nexusfile).read().split('begin sets;')
    for e in f[1].split('\n'):
        if 'charset' not in e:
            continue
        else:
            be=e.split()
            gene=be[1][:-4].replace("'","")
            rang=be[-1][:-1]
            parts.append(gene+' '+'= '+rang)
    return parts

def partition_file(prefix,types,g):
    fout=open(prefix+'.part','w')
    parts=parse_nexus('btCOMBINED.nex')
    pickle.dump((parts,types),open('types.pckl','wb'))
    for k,v in types.iteritems():
        for el in parts:
            if k in el:
                if v[0] == 'nuc':
                    s='DNA'
                else:
                    s=v[1]
                    if v[2]:
                        if g:
                            continue
                        else:
                            g=True
                    elif v[3]:
                        s+='I'
                    elif v[4]:
                        s+='F'
                fout.write('%s, %s\n'%(s,el))
    fout.close()

                    
    
    
def run_raxml(fileprefix,cluster=False):
    print 'Running raxml on the combined dataset'
    if not cluster:
        rax=Popen('raxmlHPC -m GTRGAMMA -n %s -s %s -f s -q %s -k -f a -x 12345 -N 100'
                  %(fileprefix,fileprefix+'.phy', fileprefix+'.part'),shell=True,stdout=PIPE,
                  stderr=PIPE)
        o,e = rax.communicate()
        print o,e
    else:
        cmd='raxmlHPC-HYBRID-SSE3 -T 4 -m GTRGAMMA -n %s -s %s -f s -q %s -k -f a -x 12345 -N '\
            '100'%(fileprefix,fileprefix+'.phy', fileprefix+'.part')
        qsc=Popen('qscript -now -name %s -np 300 -cmd %s'%(fileprefix,cmd),stdout=PIPE, stderr=PIPE,
                  shell=True)
        o,e=qsc.communicate()
        print o,e
    
prefix = sys.argv[1]
try:
    c= sys.argv[2]
    cluster=True
except:
    cluster=False
    
files = glob('*.fasta')
types={}
g=False
for fi in files:
    print fi
    if not os.path.isfile('types.pkcl'):
        model, gamma, inv, frec = None,None,None,None
        t=consensus(fi[:-6])   
        if t.typ =='prot':
            model, gamma, inv, frec = run_prottest(fi)
        types[fi[:-6]]=(t.typ,model, gamma, inv, frec)
        pickle.dump(types,open('Dtps.pkcl','wb'))
    else:
        tup = pickle.load(open('types.pckl'))
        parts,types = tup[0],tup[1]
    

Concatenate(prefix)
partition_file(prefix,types,g)
run_raxml(prefix,cluster)
print 'Done'


