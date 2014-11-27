#!/usr/bin/python
"""
GarliRunner Copyright (C) 2012 Jose Sergio Hleap

This script will iterate over fasta files in a given directory, create the garli.conf and run garli for each gene.
Each Fasta file should be an alignment for each gene.

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

It requires:
1) Garli
2) JModelTest or ProtTest
3) if in a cluster passTosub.py available from Alex Safatli at https://github.com/AlexSafatli/LabBlouinTools
"""
#importing bit####################################################################################################
import os, sys
from labblouin.passToqsub
from subprocess import PIPE, Popen
# End importing####################################################################################################
#Some definitions##################################################################################################

def CheckRequiredExecutables(PATH=os.getcwd(), datatype='nuc',mpi=False):
	'''
	will check if you have the mdoeltest and garli installed and in your PATH
	'''
	print 'Checking required executables and paths . . .'
	if datatype == 'nuc':
		jm = os.getenv('JMODELTEST_HOME')
		if not jm:
			print 'jModelTest not in your PATH. Please set the PATH to $JMODELTEST_HOME with export'
			print 'or set it in your bash file (.bashrc in Unix) e.g. export JMODELTEST_HOME=<PATH_TO_JMODELTEST>'
			sys.exit(-1)
		print jm
	else:
		pt = os.getenv('PROTTEST_HOME')
		if not pt:
			print 'ProtTest not in your PATH. Please set the PATH to $PROTTEST_HOME with export'
			print 'or set it in your bash file (.bashrc in Unix)e.g. export PROTTEST_HOME=<PATH_TO_PROTTEST>'
			sys.exit(-1)
		print pt
	
	p,f = Popen('which Garli-2.0',shell=True, stdout=PIPE).communicate()
	if not p.strip():
		print'No executable of GARLI in your PATH (or older version of Garli). Please copy the executable to the bin folder or assign the PATH with export'
		print 'or set it in your bash file (.bashrc in Unix)'
		sys.exit(-1)
	print p
	
	if mpi:
		garlimpi=os.getenv('GARLIMPI_HOME')
		q,g = Popen('which mpirun',shell=True, stdout=PIPE).communicate()
		if not q.strip() or not garlimpi:
			print 'Trying to use MPI with MPI applications. Install mpirun and related resources. Set a GARLIMPI_HOME to the executable with export'
			print 'e.g. export GARLIMPIHOME=[PATH_TO_GARLI], or set it in your bash file (.bashrc in Unix)'
			print 'MAKE SURE YOUR GARLI IS COMPILED FOR MPI'
			sys.exit(-1)
	
def get_models(prefix, datatype, bayesian):
	'''
	This function will execute jModelTest and get the substitution
	'''
	print 'Running model testing . . .'
	stderr=open(prefix+'.stderr','w')
	if datatype == 'nuc':
		jmthome=Popen('echo $JMODELTEST_HOME',shell=True, stdout=PIPE)
		jmthome.wait()
		o=jmthome.communicate()[0]
		
		if bayesian:
			jmodel=Popen('java -jar %s/jModelTest.jar -d %s.fasta  -f -i -g 4 -s 11 -S BEST -BIC -uLnL > %s.model'%(o.strip(), prefix,prefix),stderr=PIPE, shell=True)
			jmodel.wait()
		else:
			jmodel=Popen('java -jar %s/jModelTest.jar -d %s.fasta  -f -i -g 4 -s 11 -S BEST -AICc -uLnL > %s.model'%(o.strip(),prefix,prefix),stderr=PIPE, shell=True)
			jmodel.wait()		
		stderror = jmodel.communicate()[1]
	elif datatype == 'prot':
		cwd = os.getcwd()
		ph = os.getenv('PROTTEST_HOME')
		print ph
		#ph = Popen('echo $PROTTEST_HOME',shell=True, stdout=PIPE)
		#prottesthome = ph.communicate()[0]
		#os.chdir(prottesthome.strip())
		os.chdir(ph)
		if bayesian:
			ptmodel = Popen('java -jar prottest-3.2.jar -i %s -Dayhoff -JTT -WAG -MtREV -MtMam -all-distributions -BIC > %s.model'%(cwd+'/%s.fasta'%(prefix),prefix), shell=True,stderr=PIPE)
			ptmodel.wait()
		else:
			ptmodel = Popen('java -jar prottest-3.2.jar -i %s -Dayhoff -JTT -WAG -MtREV -MtMam -all-distributions -AICC > %s.model'%(cwd+'/%s.fasta'%(prefix),prefix), shell=True,stderr=PIPE)
			ptmodel.wait()
		stderror = ptmodel.communicate()[1]
		mv = Popen('mv %s.model %s.model'%(prefix,cwd+'/'+prefix),shell=True)
		mv.wait()
		os.chdir(cwd)   
	stderr.write(stderror)
	
def create_n_parse_modelfile(prefix, datatype, bayesian):
	'''
	will parse the model file created by jModelTest and return a string with the ratematrix bit for the garli.conf. It will also return a file with the
	jmodeltest parameters
	'''
	print 'Creating the configuration file for %s . . .'%(prefix)
	if not os.path.isfile('%s.model'%(prefix)):
		get_models(prefix, datatype, bayesian)
		
	modout=open('models.txt','a')
	if bayesian:
		modline='Dataset \t MODEL SELECTION \t Model \t f(a) \t f(c) \t f(g) \t f(t) \t Kappa \t titv \t Ra \t Rb \t Rc \t Rd \t Re \t Rf \t pInv \t gamma \t -lnL \t K \t BIC \t AICcw\n'
	else:
		modline='Dataset \t MODEL SELECTION \t Model \t f(a) \t f(c) \t f(g) \t f(t) \t Kappa \t titv \t Ra \t Rb \t Rc \t Rd \t Re \t Rf \t pInv \t gamma \t -lnL \t K \t AICc \t AICcw\n'
	garliline=''	
	if datatype == 'nuc':
		bline = open('%s.model'%(prefix)).read().split('::Best Models::\n \n\tModel \t\tf(a) \tf(c) \tf(g) \tf(t) \tkappa \ttitv \tRa\tRb\tRc\tRd\tRe\tRf\tpInv \tgamma\n----------------------------------------------------------------------------------------------------------------------------------------\n')
		bmodel = bline[1].replace('\t',' ').replace(' ', '\t')
		bbmodel=bmodel.split()
		if bayesian:
			aic=bline[0][bline[0].find('* BIC MODEL SELECTION : Selection uncertainty')+222: bline[0].find('* AICc MODEL SELECTION : Selection uncertainty')+222+86].split()
		else:
			aic=bline[0][bline[0].find('* AICc MODEL SELECTION : Selection uncertainty')+222: bline[0].find('* AICc MODEL SELECTION : Selection uncertainty')+222+86].split()
		aic=aic[1]+'\t'+aic[2]+'\t'+aic[3]+'\t'+aic[5]+'\n'
		if modline in open('models.txt','r').read():
			modline= '%s\t%s'%(prefix,bmodel+aic)
		else:
			modline+='%s\t%s'%(prefix,bmodel+aic)
		f = open(prefix+'.model').read().split('-----\n \n Model selected: \n')
		sli = f[1][:f[1].find('\n \nML tree')+23]
		model= sli.split('\n')
		garliline+='datatype = nucleotide\nratematrix = (%s)\n'%(model[1][model[1].find('= ')+2:].replace('',' '))
		if bbmodel[2] == bbmodel[3] == bbmodel[4] == bbmodel[5]:
			garliline+='statefrequencies = equal\n'
		else:
			garliline+='statefrequencies = estimate\n'
			
		if bbmodel[-1].strip() == 'N\A':
			garliline+='ratehetmodel = none\n'
		else:
			garliline+='ratehetmodel = gamma\nnumratecats = 4\n'
			
		if bbmodel[-2].strip() == 'N\A':
			garliline+='invariantsites = none\n\n'
		else:
			garliline+='invariantsites = estimate\n\n'
			
	if datatype == 'prot':
		if bayesian:
			bline = open('%s.model'%(prefix)).read().split('Best model according to BIC:')
			modelsel='BIC'
		else:
			bline = open('%s.model'%(prefix)).read().split('Best model according to AICc:')
			modelsel='AICc'
		model = bline[1][:bline[1].find('\n')]
		modelsel+='\t%s\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\t'%(model.strip())
		
		if modline in open('models.txt','r').read():
			modline= '%s\t%s'%(prefix,modelsel)
		else:
			modline+='%s\t%s'%(prefix,modelsel)
		
		if model.strip().split('+')[0].lower() == 'jtt':
			usemodel='jones'
		else:
			usemodel = model.strip().split('+')[0].lower()
		garliline+='datatype = aminoacid\nratematrix = %s\n'%(usemodel)
		if '+F' in model:
			garliline+='statefrequencies = empirical\n'
		else:
			garliline+='statefrequencies = equal\n'
		
		if not '+G' in model:
			garliline+='ratehetmodel = none\nnumratecats = 1\n'
			gamma='N/A'
		else:
			garliline+='ratehetmodel = gamma\nnumratecats = 4\n'
			gamma=bline[1].split('alpha (+G):')
			gamma = gamma[1][:gamma[1].find('\n')].strip()
		if not '+I' in model:
			garliline+='invariantsites = none\n\n'
			inv='N/A'
		else:
			garliline+='invariantsites = estimate\n\n'
			inv=bline[1].split('p-inv (+I):')
			inv = inv[1][:inv[1].find('\n')].strip()
			
	modout.write(modline)
	modout.close()		
	return garliline

def process_outgroups(spslistfile,outgroups):
	outgroupslist=[]
	if outgroups:
		print 'Getting the outgroups ...'
		f=open(spslistfile)
		for line in f:
			if not line.split()[0] in outgroups:
				continue
			else:
				genusini=line.split()[0][0]+'._'
				if line.split()[1] == 'cf.':
					species=line.split()[1]+'_'+line.split()[2]
				else:
					species=line.split()[1]
				name = genusini+species
				outgroupslist.append(name)
	return outgroupslist

def process_fasta(prefix):
	'''
	this will process the fasta and return the numer of taxa and a dictionary with the indexes of each sps
	'''
	print 'Processing the FASTA . . .'
	spsindex={}
	usedsps=[]
	f=open(prefix+'.fasta').read().split('>')
	taxa=0
	for e in f:
		if e == '':
			continue
		else:
			taxa+=1
			temp = e[:e.find('\n')].split('_')
			if len(temp) > 2:
				name = temp[0]+'_'+temp[1]+'_'+temp[2]
			else:
				name = temp[0]+'_'+temp[1]
			if not name in usedsps:
				usedsps.append(name)
				spsindex[name]= str(taxa)
			else:
				spsindex[name]+=' %s'%(str(taxa))
	return taxa, spsindex
			
def general_line(prefix,datatype, bayesian, bootstrapreps,outgroups,spslistfile,mpi):
	'''
	This will create the first block for the Garli configuration file
	'''
	taxa, spsindex = process_fasta(prefix)
	general='[general]\ndatafname = %s.fasta\nconstraintfile = none\nstreefname = stepwise\nattachmentspertaxon = %d\n'%(prefix,taxa*2)
	general+='ofprefix = %s\nrandseed = -1\navailablememory = 512\nlogevery = 10\nsaveevery = 100\nrefinestart = 1\n'%(prefix+'.'+datatype)
	general+='outputeachbettertopology = 0\noutputcurrentbesttopology = 0\nenforcetermconditions = 1\n'
	if bootstrapreps == 0:
		general+='genthreshfortopoterm = %d\nscorethreshforterm = 0.05\nsignificanttopochange = 0.01\noutputphyliptree = 0\n'%(20000)
	else:
		general+='genthreshfortopoterm = %d\nscorethreshforterm = 0.05\nsignificanttopochange = 0.01\noutputphyliptree = 0\n'%(1000)
	if not outgroups:
		general+='outputmostlyuselessfiles = 0\nwritecheckpoints = 0\nrestart = 0\noutgroup = 1'
	else:
		general+='outputmostlyuselessfiles = 0\nwritecheckpoints = 0\nrestart = 0\noutgroup ='
		outgroupsl = process_outgroups(spslistfile,outgroups)
		if not outgroupsl:
			general += str(1)
		else:
			for e in outgroupsl:
				for k in spsindex.iterkeys():
					if not e == k[:k.rfind('_')]:
						continue
					else:
						general+= ' %s'%(spsindex[k])
	if mpi:
		general+='\noutputsitelikelihoods = 0\ncollapsebranches = 1\nsearchreps = 1	\n\n'
	else:
		general+='\noutputsitelikelihoods = 0\ncollapsebranches = 1\nsearchreps = 3	\n\n'
		
	general += create_n_parse_modelfile(prefix, datatype, bayesian)
	general += master_line(bootstrapreps,mpi)
	return general

def master_line(bootstrapreps,mpi):
	'''
	will create the masterline of the configfile
	'''
	internal= 1
	if bootstrapreps > 0:
		internal = 0
	master ='[master]\nnindivs = 4\nholdover = 1\nselectionintensity = 0.5\nholdoverpenalty = 0\nstopgen = 5000000\n'
	master+='stoptime = 5000000\nstartoptprec = 0.5\nminoptprec = 0.01\nnumberofprecreductions = 10\n'
	master+='treerejectionthreshold = 50.0\ntopoweight = 1.0\nmodweight = 0.05\nbrlenweight = 0.2\n'
	master+='randnniweight = 0.1\nrandsprweight = 0.3\nlimsprweight =  0.6\nintervallength = 100\nintervalstostore = 5\n\n'
	master+='limsprrange = 6\nmeanbrlenmuts = 5\ngammashapebrlen = 1000\ngammashapemodel = 1000\nuniqueswapbias = 1.0\n'
	master+='distanceswapbias = 1.0\n\n'
	if mpi:
		boots='bootstrapreps = %s\nresampleproportion = 1.0\ninferinternalstateprobs = %d'%(str(round(bootstrapreps/runs)),internal)
	else:
		boots='bootstrapreps = %s\nresampleproportion = 1.0\ninferinternalstateprobs = %d'%(str(bootstrapreps),internal)
	
	return master+boots

def execute_all(prefix,datatype='nuc',bayesian=False,bootstrapreps=0,outgroups=[],mpi=False,spslistfile=None,processors=4,runs=20,cluster=False):
	'''
	Operative funtion to run all
	'''
	#CheckRequiredExecutables(os.getcwd(), datatype, mpi)
	
	configfile=open('garli.conf','w')
	configfilecontent = general_line(prefix,datatype, bayesian, bootstrapreps,outgroups,spslistfile,mpi)
	configfile.write(configfilecontent)
	configfile.close()
	os.mkdir(prefix)
	mv= Popen('mv %s.* garli.conf ./%s/'%(prefix,prefix),shell=True)
	mv.wait()
	cwd=os.getcwd()
	os.chdir(prefix)-MPI=
	if mpi and cluster:
		print 'Using MPI version in cluster (needed) . . .'
		cmd = passToqsub.returnScript('Garli-2.0 -%d'%(runs),'%s'%(prefix),processors,processors)
		launch = Popen(cmd,shell=True)
		launch.wait()
	elif mpi and not cluster:
		print 'Using MPI version . . .'
		#os.rename('%s.conf'%(os.getcwd()+'/'+prefix),os.getcwd()+'/'+'garli.conf')
		ou =os.getenv('GARLIMPI_HOME')
		print ou
		cl='mpirun -num-cores %d -np %d %sGarli-2.0 -20'%(nc,nc,ou.strip())
		print 'Executing Garli . . .'
		garli= Popen('%s'%(cl),shell=True, stderr=PIPE, stdout=PIPE)
		#garli.wait()
		out,err=garli.communicate()
		Err=open('garli.stderr','w')
		Err.write(err)
		Err.write('command line passed:%s'%(cl))
		Err.close()
	else:
		garli= Popen('Garli-2.0',shell=True, stderr=PIPE, stdout=PIPE)
		garli.wait()
	
	if not cluster:
		print 'Running Sumtrees . . .'
		sumtrees= Popen('sumtrees.py --to-newick --rooted -m* --extract-edges=%s -f0.5 -p -d0 %s.%s.run*.boot.tre > %s.70consensus.sumtrees'%(prefix+'.edges',prefix,datatype,prefix),shell=True,stderr=PIPE, stdout=PIPE)
		sumtrees.wait()
		o,e=sumtrees.communicate()
		binsumtrees= Popen('sumtrees.py --to-newick --no-meta-comments --no-taxa-block --no-support  --no-node-ages --no-summary-metadata -m* -f0.5 -p -d0 %s.%s.run*.boot.tre > %s.70consensus.bin.sumtrees'%(prefix,datatype,prefix),shell=True,stderr=PIPE, stdout=PIPE)
		p,f=binsumtrees.communicate()
		O=open('sumtrees.stdout','w')
		O.write(o+'\n'+p)
		O.close()
		E=open('sumtrees.stderr','w')
		E.write(e+'\n'+f)
		E.close()
	
# End of definitions###############################################################################################
# Aplication of the code ##########################################################################################
if __name__ == "__main__":
	if len(sys.argv) == 1 or '-help' in sys.argv:
		print 'usage GarliRunner.py [prefix] [option]'
		print 'Options:'
		print '\t-datatype=XXX : will run the scripts for nucleotides (nuc) or aminoacids (prot) ( Default: nuc )'
		print '\t-bayesian : Boolean. Will use BIC when choosing the subtitution model( Default: False )'
		print '\t-bootstrapreps=XXX : Will perform boostrap and with XXX replicates ( Default : No )'
		print '\t-outgroups=XXX,YYY,...,ZZZ: Genus of the outgroups. This option has to be used with the Species list ( Default: No )'
		print '\t-splist=XXX : XXX is the filename with the species names to be used. Outgroups should be in this file ( Default: No )'
		print '\t-MPI=XXX,YYY : will run MPI in multiprocessors. XXX is the numer of cores and process to use. YYY is the number of',		
		print '\times to run the config file. You will need an mpirun application ( Default: no )'
		print '\t-cluster : If is going to use a cluster( Default : no )'
		
	# Default Parameters ###################################################
	datatype = 'nuc'
	bayesian = False
	bootstrapreps=0
	outgroups = False
	splist = False
	mpi=False
	cluster = False
	prefix = sys.argv[1]
	# Command line input ###################################################
	for arg in sys.argv[1:]:
		if arg.startswith('-datatype='):
			datatype=arg[10:]
		if arg.startswith('-outgroups='):
			outgroups = str(arg[11:].split(','))[1:-1]
		if arg == '-bayesian':
			bayesian = True
		if arg.startswith('-bootstrapreps='):
			bootstrapreps = int(arg[15:])
		if arg.startswith('-splist='):
			splist= True
			spslistfile=arg[8:]
		if arg.startswith('-MPI='):
			mpi=True
			nc=int(arg[5:].split(',')[0])
			runs=int(arg[5:].split(',')[1])
		if arg == '-cluster':
			cluster = True
	execute_all(prefix,datatype,bayesian,bootstrapreps,outgroups,mpi,spslistfile,int(nc),int(runs),cluster)
