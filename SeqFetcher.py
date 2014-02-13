#!/usr/bin/python
'''
This script will get the sequences from NCBI. Read the NCBI code of conduct before.
This script will return a set of files of unaligned ([].fas) and aligned ([].fasta,
if specified) sequences for a given search. It also return a list of species and a
pickled file (GenesnSpecies.pckl), if the program crash in the writing part.

SeqFetcher Version 1.2 Copyright (C) 2013  Jose Sergio Hleap

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
1) Biopython

'''
#importing bit########################################################################
import sys,glob,os,time,unicodedata,string,types,htmlentitydefs,urllib,httplib,\
       fnmatch,optparse
import cPickle as pickle
from Bio import Entrez
#from collections import defaultdict
from urllib import urlencode
from subprocess import Popen, PIPE
# End importing#######################################################################
# Start Definitions ##################################################################
class SPcheck:
	'''
	Conversion of mike's access to CAS
	'''
	def __init__(self, SEARCH='Urotrygon aspidurus',TYPE='Species'):
		if len(SEARCH.split()) == 2:
			self.SEARCH = SEARCH
		else:
			newname = SEARCH.split()[0] + ' ' + SEARCH.split()[1]
			self.SEARCH = newname
		self.TYPE   = TYPE
		self.data   = None
		self.URL = None
		self.CGI = None
		self.HEADERS = None
		self.connection = None
		self.f = None
		self.validname = ''
		self.SetData()
		self.SetGlobals()
		self.SetConnection()
		self.GetConnectionResult()
		self.GetValidName()
		print '\t\tTesting the validity of the scientific name "%s"'%(self.SEARCH)
		
	def SetData(self):
		# The POST data includes the search parameters and the submit button 
		#(which probably isn't even used)
		self.data = {'tbl':self.TYPE,'contains':self.SEARCH,'Submit':'Search'}
		# Encode POST data
		self.data = urllib.urlencode(self.data)
		return self.data

	def SetGlobals(self):
		self.URL    = "researcharchive.calacademy.org"
		self.CGI = "/research/Ichthyology/catalog/fishcatmain.asp"
		self.HEADERS = {
		    "Accept"          : "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
		    "Accept-Encoding" : "gzip, deflate",
		    "Accept-Language" : "en-us,en;q=0.5",
		    "Cache-Control"   : "max-age=0",
		    "Content-Length"  : str(len(self.data)),
		    "Content-Type"    : "application/x-www-form-urlencoded",
		    "Host"            : "researcharchive.calacademy.org",
		    "Referer"         : "http://researcharchive.calacademy.org/research/"\
		    "Ichthyology/catalog/fishcatmain.asp",
		    "User-Agent"      : "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.6; rv:11.0)"\
		    " Gecko/20100101 Firefox/11.0"
		}
		return

	def SetConnection(self):
		# Send POST request and get resulting data
		self.connection = httplib.HTTP(self.URL,80)
		self.connection.putrequest("POST", self.CGI)
		for header,value in self.HEADERS.iteritems(): 
			self.connection.putheader(header,value)
		self.connection.endheaders()
		self.connection.send(self.data)
		reply, msg, hdrs = self.connection.getreply()
		if reply != 200:
			print "ERROR: Website returned error",reply,msg
			sys.exit(0)

		return self.connection

	def GetConnectionResult(self):
		self.f = self.connection.getfile()
		return self.f
	
	def CurrentStatus(self, HTML, query):
		#  when not available
		if ('No matches found') in HTML:
			name = False
		elif 'Current status' not in HTML:
			name = False
		else:
			current = HTML.split('<b>Current status:</b>')
			if current[1].find('Valid as <i>') == -1:
				name = current[1][current[1].find('Synonym of <i>')+12:current[1].find('</i>')]
			else:
				name = current[1][current[1].find('Valid as <i>')+12:current[1].find('</i>')]
		return name	
	
	def cleanDuplicates(self,names):
		b = False
		count={}
		s = set(names)
		for i in s:
			count[i]=names.count(i)
		if self.SEARCH in s:
			b=True
		return b, count, s
	
	def GetValidName(self):
		names=[]
		for l in self.f:
			if '<p class="result"><b>' in l:
				if not self.CurrentStatus(l,self.SEARCH):
					pass
				else:
					name = self.CurrentStatus(l,self.SEARCH)
					names.append(name)
		valid, count, s = self.cleanDuplicates(names)
		if valid == True:
			self.validname = self.SEARCH
		elif  len(s) == 0:
			pass		
		
		elif len(s) != 1:
			self.validname = max(count)

		else:
			self.validname = names[0]
			
class RecordsFetcher:
	def __init__(self, db, term, filename):
		self.genes = {}#defaultdict(list)
		self.sps = {}
		self.spl = []
		self.record = None
		self.db=db
		self.listfname=filename
		#check if a pikled file is in path
		if os.path.isfile('GenesnSpecies.pckl'):
			self.genes,self.sps = pickle.load(open('GenesnSpecies.pckl'))
		else:
			if self.listfname:
				self.term=None
				self.splist()
				self.loop_over_list()
			else:
				self.term = term
				self.get_records(self.term)
	def splist(self):
		'''
		take the file with the species name and parse it
		'''
		f = open(self.listfname)
		for line in f:
			bline = line.replace(':',' ').strip().split()
			if len(bline) >= 2:
				if 'cf.' not in bline:
					self.spl.append(str(bline[0])+' '+str(bline[1]))
				else:
					self.spl.append(str(bline[0])+' '+str(bline[1])+' '+str(bline[2]))
			elif len(bline) == 1:
				self.spl.append(bline[0])
				print 'Only genus will be queried'
			else:
				print 'List file not in apropriate format. Exiting the program.'
				sys.exit()
				
	def get_records(self,term):
		'''get record'''
		c=0
		print 'Getting records for the search "%s"...'%(term)
		while True:
			try:
				c+=1
				h = Entrez.esearch(db=self.db, term=term, retmode='xml', retmax='1000')
				rec=Entrez.read(h)
				break
			except:
				if c == 20:
					print 'There is a problem with either the query or NCBI API.'\
					      'Check it out and re-run the program'
					sys.exit()
				else:
					print 'Awaiting for NCBI response'
					time.sleep(2)
		idschecked=[]
		# Iterate over the records
		for I in rec['IdList']:
			#check if has not been assign
			if I in idschecked:
				continue
			idschecked.append(I)
			print '\tGI: '+I
			seqs = {}
			#get individual record
			while True:
				try:
					handle=Entrez.efetch(db=options.db, id=I, retmode='xml')
					record=Entrez.read(handle)
					break
				except:
					print '\t\tAwaiting NCBI response'
					time.sleep(10)
			if not record:
				'No record found for %s'%(term)
				continue
			# deal with the whole genome sequencing
			definition=record[0]['GBSeq_definition']
			if  'genome' in definition:
				if 'shotgun' in definition or 'WGS' in definition:
					print '\t\tSkipping this entry. Is a shotgun genome sequence'
					continue
				else:
					self.solve_genome(record)
			else:
				#Get the accession number
				ID = record[0]['GBSeq_primary-accession']
				# Get the gene name
				gene = self.get_gene(record)
				# deal with unasigned DNA
				if not gene:
					continue
				if not gene in self.genes:
					self.genes[gene]={}
								
				
				#Get the sequence and organize them acording accession
				try:
					seq  = record[0]['GBSeq_sequence']
				except:
					continue
				
				if ID not in self.genes[gene]:
					self.genes[gene][ID]=seq
				#	seqs[ID]=seq
				self.checknfill_species(record)
		pickle.dump([self.genes, self.sps], open('GenesnSpecies.pckl','w'))
				
	def loop_over_list(self):
		for e in self.spl:
			if options.delay:
				print 'Waiting %f seconds to send new request'%(round(secs,2))
				time.sleep(secs)
			self.get_records(e+'[organism]')
	
	def get_gene(self,record):
		if len(record[0]['GBSeq_feature-table']) > 1:
			table = record[0]['GBSeq_feature-table'][1]
		else:
			table = record[0]['GBSeq_feature-table'][0]
			
		try:
			qualsname=table['GBFeature_quals'][0]['GBQualifier_name']
			if qualsname == 'product' or qualsname == 'gene':
				gene=table['GBFeature_quals'][0]['GBQualifier_value'].lower()
			else:
				#just break
				record[999]
		except:
			try:
				if table['GBFeature_key'] == 'CDS' and (table['GBFeature_quals'][2]['GBQualifier_name'] == 'product'):
					gene=record[0]['GBSeq_feature-table'][1]['GBFeature_quals'][2]['GBQualifier_value'].lower()
				else:
					record[999]
			except:
				try:
					gene=record[0]['GBSeq_keywords'][0]
				except:
					gene = False
		return gene

	def	checknfill_species(self,record):
		'''Get and check the species'''
		if options.check:
			species = SPcheck(SEARCH=record[0]['GBSeq_organism'])
			if species.validname == '' or record[0]['GBSeq_organism'].split >= 3:
				species = record[0]['GBSeq_organism']
			else:
				species = species.validname
		else:
			species = record[0]['GBSeq_organism']
			
		for k in self.genes.iterkeys():
			for l in self.genes[k].keys():
				if not l in self.sps:
						self.sps[l]= species
	def solve_genome(self,record):
		genome = record[0]['GBSeq_sequence']
		quals = record[0]['GBSeq_feature-table']
		for e in quals:
			iquals= e['GBFeature_quals'][0]['GBQualifier_name']
			if iquals == 'product' or iquals == 'gene':
				gene = e['GBFeature_quals'][0]['GBQualifier_value']
				fro = int(e['GBFeature_intervals'][0]['GBInterval_from'])-1
				to = int(e['GBFeature_intervals'][0]['GBInterval_to'])
				ID = e['GBFeature_intervals'][0]['GBInterval_accession']
				if to > len(genome):
					seq = genome[fro:]
				else:					
					seq = genome[fro:to]
				if not gene in self.genes:
					self.genes[gene]={}
				if ID not in self.genes[gene]:
					self.genes[gene][ID]=seq
				self.checknfill_species(record)


	def Filter(self):
		'''
		will create an OK folder where the genes with more than 4 sps are found. 
		Will check for that number and store the fasta files in the folder
		'''
		print 'Using filter'
		print 'WARNING! the files are explored by gene name which can have synonyms.'\
		      'Please check the folder for those synonyms and merge them. There is a'\
		      'sccript called merger.py to do so.'
		os.mkdir('OK')
		files = fnmatch.filter(os.listdir(os.getcwd()),'*SpList.txt')
		for f in files:
			if '(' in f:
				f = f.replace('(','').replace(')','_')
			s=[]
			#'''
			#fi=open(f).read().split('\n>')
			#for el in fi:
				#s.append('_'.join(el.split('_')[:2]))'''
			countlines = Popen('cat %s | wc -l'%(f),shell=True, stdout=PIPE)
			o,e = countlines.communicate()
			#s=len(set(s))
			if int(o) >= int(options.filter) :
				mv = Popen('mv %s %s.fas ./OK'%(f,f[:-10]),shell=True)
				
	def writeFastas(self):
		'''
		write a fasta file for each gene and a species and accession number lists
		'''
		for g in self.genes.keys():
			try:
				outf=open(g.replace(' ','_')+'.fas', 'w')
				fout = open(g.replace(' ','_')+'SpList.txt','w')
			except:
				continue
			newsp={}
			for e, v in self.genes[g].iteritems():
				acc= e
				sp = self.sps[acc]
				if not sp in newsp:
					newsp[sp]=[acc]
				else:
					newsp[sp].append(acc)
				sp = sp.split()
				sp = sp[0][0] + '._' + sp[1]
				se = v
				outf.write('>'+sp+'_'+acc+'\n'+se+'\n')
			outf.close()
			for s, a in newsp.iteritems():
				fout.write(s+': ')
				for e in a:
					fout.write(e+',')
				fout.write('\n')
			fout.close()

def executeMuscle():
	'''
	get the fasta files and align them using muscle
	'''
	
	fastafiles = glob.glob('*.fas')
	for f in fastafiles:
		os.system('muscle -in %s -out %s'%(f,f[:-4]+'.fasta'))

# End of definitions##################################################################
# Aplication of the code #############################################################
if __name__ == "__main__":
	# Command line input #############################################################
	opts = optparse.OptionParser(usage='%prog [options] email')
	opts.add_option("-t", "--term", dest="term", action="store",
		              help="Any term that NCBI accepts giving the database provided."\
	                  "If more than one word is used, use quotation (see http://www."\
	                  "ncbi.nlm.nih.gov/books/NBK3837/).", 
	                  default = 'Rhinobatos[organism]')
	opts.add_option("-L", "--list", dest="filename", default=None,
		              help = "Use this option if a list of scientific names to be"\
		              " downloaded is provided ")
	opts.add_option("-d", "--database", dest="db", default='protein', action="store",
	                  help= "Database to be queried. The term used must be suitable"\
	                  " for the database used (see http://www.ncbi.nlm.nih.gov/book"\
	                  "s/NBK3837/).")
	opts.add_option("-a", "--align", dest="muscle",action="store_true", default=False,
	                  help="Use Muscle (Edgar 2004) to aling the fasta files in the"\
	                  " execution folder")
	opts.add_option("-y", "--delay", dest="delay", action="store", type="int",
	                help="Delay each retrieval for some given seconds."\
	                "Very important for large queries.",default=None)
	opts.add_option("-F","--filter", dest="filter", type="int", action="store",default=0,
	                help="Use this option if you want to check which genes have more than"\
	                " a given number of species. An OK folder will be created and the fas"\
	                "ta files of the genes with more than the specied species number will"\
	                " be stored.")
	opts.add_option("-c","--checksp", dest='check', action="store_false", default=True,
	                help="Use this option if you want to turn off the species checking."\
	                "The species check in this program is only intended for fishes, se"\
	                "arching the Eschmeyer catalog of fishes. If other class of organi"\
	                "sms is required or if the Eschmeyer classification is not what you"\
	                " are looking for, use this option.")
	options, args = opts.parse_args()
	# End command line input ##########################################################
	
	# Tell NCBI who you are
	Entrez.email = args[0]
	#Create the instance of SeqFetcher
	S = RecordsFetcher(options.db, options.term, options.filename)
	#write the fasta files
	S.writeFastas()
	if options.muscle:
		executeMuscle()
	if options.filter > 0:
		S.Filter()