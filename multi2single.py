#!/usr/bin/python
"""
multi2single Copyright (C) 2012 Jose Sergio Hleap
 This Script allows you to collapse multiple species into a single leaf. It will
 return a file with a list of unique names, a csv file with Species,gene,type,
 Accession, and Source (the last one if the format is apropriate), and a new tree
 with "collapsed" leafs. This script assumes that you have a bunch of tree files
 in newick format with the extension .tree, and the name must follow <gene>_<type>
 where gene is the label for the given marker and type is the label for the type of
 data (either prot for protein or nuc for nucleotide)

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
1) Dendropy available at http://pythonhosted.org/DendroPy/
"""

import dendropy as D
from glob import glob as G

files = G('*.tree')

#equival=[]
spcount={}
uniquenames=[]
csvou= open('rajoidei.csv','w')
csvou.write('Species \t gene \t type \t Accession \t Source \n')
for f in files:
	print 'Processing %s'%(f)
	gene = f.split('_')
	gen = gene[0]
	typ = gene[1].strip('.tree')
	newtree=open('./%s_%s.newick'%(gen,typ),'w')
	if 'prot' in typ:
		typ = 'Protein'
	else:
		typ = 'Nucleotide'
	done = []
	T=D.Tree.get_from_path(f,'newick')
	for t in T.taxon_set:
		bl = t.label.split()
		if 'BOLD' in bl:
			source = 'BOLD'
		else:
			source = 'NCBI'
		if not ' '.join(bl[:2]) in done:
			sp = ' '.join(bl[:2])
			done.append(sp)
			t.label = sp
			spcount[sp]=0
			if sp not in uniquenames:
				uniquenames.append(sp)
		else:
			sp = ' '.join(bl[:2])
			t.label = sp+str(spcount[sp])
			spcount[sp]+=1
			
		if source == 'BOLD':
			acc= '-'.join(bl[3:])
		else:
			acc = bl[2]
			
		csvou.write(sp+'\t'+gen+'\t'+typ+'\t'+acc+'\t'+source+'\n')
	newtree.write(T.as_newick_string()+';')
	newtree.close()
csvou.close()
un=open('species.uniquenames','w')
un.write('\n'.join([x.replace(' ','_') for x in uniquenames]))
un.close()