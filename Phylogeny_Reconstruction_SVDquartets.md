#Phylogeny Reconstruction using SVDquartets
24 January 2019

Mao Li, Laura L. Klein, Keith Duncan, Ni Jiang, Daniel H. Chitwood, Jason Londo, Allison J. Miller, Christopher N. Topp. Characterizing grapevine (Vitis spp.) inflorescence architecture using X-ray imaging: implications for understanding bunch density.

###Contents
These analyses were performed using a Mac computer and a high performance computing cluster. Individuals used for phylogeny reconstruction (N=99; Table 1) were extracted from a the full species dataset in Klein et al. 2018 (doi:10.1002/ajb2.1033). The full species dataset file is available upon request. The dataset used in this study (rachisphy\_renamed\_taxap.nex) is available in this github folder.

###Required files:
rachisphy\_renamed\_taxap.nex

_Note:_
This file was prepared for SVDquartets by adding the following to the end of the nexus file:

	begin sets;
	taxpartition species = 
	V_acerifolia: 1-9,
	V_aestivalis: 10,
	V_amurensis: 11-12,
	V_cinerea: 13-25,
	V_coignetiae: 26,
	V_labrusca: 27-38,
	V_palmata: 39,
	V_riparia: 40-86,
	V_rupestris: 87-96,
	V_vulpina: 97-98
	  ;
	END;
This is a standard nexus file format that indicates to the program which species individuals should be considered (i.e., the flag `taxpartition species`).

###Programs used:
PAUP\* version 4.0a (build 163). _Note that PAUP\* can only be used with the most current version. Future analyses may be in a different version._

##SVDquartets program commands (executed in PAUP\*):
	
	# generate a log file for the run
	log file=rachisphy.log;
	# load the executable file 
	exe rachisphy_renamed_taxap.nex
	# run SVDquartets with the desired parameters
	svdq evalq=all taxpartition=species bootstrap;
	# save the tree file
	SaveTrees file = Vitis_rachises_SVDquartet_MR_indv.tre format = Newick brLens = yes supportValues = Both trees = all

At the end of the run, SVDquartets generates a tree file, a format that captures the phylogenetic tree topology as a coded sequence.