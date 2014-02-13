#!/usr/bin/python
"""
mergeFasta Copyright (C) 2013 Jose Sergio Hleap

This script will get two or more <prefix>.fasta files and its corresponding 
<prefix>SpList.txt and merge them. This is useful when SeqFetcher.py creates 
two files with different names, but that refers to the same gene (i.e. cytb and
cytochrome b). It requires a <prefix>SpList.txt which is a file linking a species
name with all accession numbers downloaded. It is of the form:
   Species1: Accesion1, Accesion2, ... , Accession_n
   Species2: Accesion1, Accesion2, ... , Accession_n
   .
   .
   .
   Speciesn: Accesion1, Accesion2, ... , Accession_n
   

The fasta file should have the .fas extension and should contain the species name 
and accession as sequence header like:
>Species1_accesion1
<SEQUENCE BIT>
>Species1_accession2
<SEQUENCE BIT>
>Species2_accession3
<SEQUENCE BIT>
.
.
.

This files are also created by SeqFetcher.py, freely available at 
https://github.com/jshleap/Phylogenetics

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

"""
import sys,os

def merge_fastas(l):
    lists=[]
    for e in l:
        if os.path.isfile(e+'.fas'):
            f = open(e+'.fas').read().split('\n>')
        else:
            f = open(e+'.fasta').read().split('\n>')
        lists.extend(f)
    s = set(lists)
    fout=open(l[0]+'_merged.fas','w')
    for g in s:
        if g == '' or g == ' ':
            continue
        fout.write('>'+g+'\n')
    fout.close()

print 'Usage: merger.py <prefix1> <prefix2> .... <prefix_n>'
print 'The prefixes are the common part between the fasta file and the list file. Fasta files should have the .fas extension' 
prefixes= sys.argv[1:]
sp={}
for p in prefixes:
    f = open(p+'SpList.txt')
    for l in f:
        if l == '':
            continue
        bl=l.split(':')
        acc=bl[1][:-1].split(',')
        if bl[0] not in sp:
            sp[bl[0]]=[]
            for e in acc:
                if e == '':
                    continue
                if e not in sp[bl[0]]:
                    sp[bl[0]].append(e.strip())
        else:        
            for e in acc:
                if e == '':
                    continue                
                if e not in sp[bl[0]]:
                    sp[bl[0]].append(e.strip())

fil = open(prefixes[0]+'SpList_merged.txt','w')
for k,v in sp.iteritems():
    fil.write(k+': '+','.join(v)+',\n')#the last comma is to mimic the normal file
merge_fastas(prefixes)
#clean up
if not os.path.isdir('ori_merged_fastas'):
    os.system('mkdir ./ori_merged_fastas')
for pr in prefixes:
    os.system('mv %s.fas ./%s'%(pr,'ori_merged_fastas'))
    os.system('mv %sSpList.txt ./%s'%(pr,'ori_merged_fastas'))



