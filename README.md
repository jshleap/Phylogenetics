Phylogenetics
=============

Python scripts for phylogenetics. This is mainly work on creation of trees and super trees, and some utilities related with them. Phylogenetic building need appropriate software, so this is not intended for tree buiding, but more to glue techniques. It ONLY WORKS RELIABLY IN LINUX SYSTEMS!!!!! Feel free to change all paths, end of line characters, etc.. to make it work on others (Yes, i did not use os.path modules to make it happen).

It contains:


GarliRunner.py:
--------------
  This script will iterate over fasta files in a given directory, create the garli.conf and run garli for each gene.
  Each Fasta file should be an alignment for each gene.
  
  It requires:
  1) Garli
  2) JModelTest or ProtTest
  3) if in a cluster passTosub.py available from Alex Safatli at https://github.com/AlexSafatli/LabBlouinTools


mergeFasta.py:
-------------
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
  SEQUENCE BIT
  >Species1_accession2
  SEQUENCE BIT
  >Species2_accession3
  SEQUENCE BIT
  .
  .
  .
  
  
  This files are also created by SeqFetcher.py, freely available in this repo.


multi2single.py:
---------------
  This Script allows you to collapse multiple species into a single leaf. It will
  return a file with a list of unique names, a csv file with Species,gene,type,
  Accession, and Source (the last one if the format is apropriate), and a new tree
  with "collapsed" leafs. This script assumes that you have a bunch of tree files
  in newick format with the extension .tree, and the name must follow <gene>_<type>
  where gene is the label for the given marker and type is the label for the type of
  data (either prot for protein or nuc for nucleotide).
  Requires:
  1) Dendropy available at http://pythonhosted.org/DendroPy/
  

properAlignment.py:
------------------
  This script will go over the files created by SeqFetcher.py and will execute 
  translatorX on the coding sequences, giving the apropriate genetic code. For 
  non-coding or RNA, will use muscle.
  Requires:
  1) Biopython
  2) translatorx.pl available at http://pc16141.mncn.csic.es/cgi-bin/translatorx_vLocal.pl
  3) perl
  4) muscle available at http://www.drive5.com/muscle/


SeqFetcher.py:
-------------
  This script will get the sequences from NCBI. Read the NCBI code of conduct before.
  This script will return a set of files of unaligned ([].fas) and aligned ([].fasta,
  if specified) sequences for a given search. It also return a list of species and a
  pickled file (GenesnSpecies.pckl), if the program crash in the writing part.
  Requires:
  1) Biopython
  
  
SuperMatrix.py:
--------------
  This  script will go over the fasta files created by SeqFetcher (or at least in the same format),
  will create new fasta files with consensus sequences per species and transform them to nexus
  using Joseph Hughes (University of Glasgow) Consensus.pl and Fasta2Nexus.pl, respectively.
  by , concatenate all the sequences using the 'NEXUS' module of biopython ('NEXUS:
  An extensible file format for systematic information' Maddison, Swofford, Maddison. 1997. Syst. Biol.
  46(4):590-621, and will launch a RAxML (Stamatakis, A. (2006). RAxML-VI-HPC: maximum likelihood-based 
  phylogenetic analyses with thousands of taxa and mixed models. Bioinformatics, 22(21), 2688-2690.)
  either in fester (a private cluster) or locally.
  
  Requires:
  1) Consensus.pl, which in turn requires BioPerl libraries. Available at           
     http://www.vcu.edu/csbc/bnfo601/Scenarios/Blast/consensus.pl
  2) Biopython
  3) Fasta2Nexus.pl, which in turn requires BioPerl libraries
  4) RAxML, available at http://sco.h-its.org/exelixis/web/software/raxml/index.html
  5) perl available at http://raven.iab.alaska.edu/~ntakebay/teaching/programming/perl-scripts/fasta2nexus.pl
  6) A utils folder available at https://github.com/AlexSafatli/LabBlouinTools
  
  Points 1,3 should be in the same folder of the selected fasta files. The requirement 6 should be either in the same     folder or in the python path
