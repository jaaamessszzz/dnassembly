"""
GGpart.py

Modified from code originally written by Will DeLoache.
"""

import sys, os
import copy
import pdb
from Bio import SearchIO, SeqIO
#from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Restriction import BsmBI, BbsI, BsaI, Restriction, FormattedSeq
import more_itertools

from .GGfrag import *
from .MergeSegments import *
from .OHfinder import *
from .AssemblyInstructions import *
from .removeRS import *
from .findTemplates import *

from ...dna.part import Part, DNA
from ...utils.benchlingAPI import *
from ...design.partdesigner import *
from ..cloning import StickyEndAssembly

#from CodonOptimize.codonOptimize import codonOptimize
#from makeBLASTdb import makeBLASTdb

# BASE_PATH is the absolute path of ../../.. relative to this script location
#BASE_PATH = reduce (lambda l,r: l + os.path.sep + r, os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep )[:-3] )
#sys.path.append( os.path.join( BASE_PATH, "main" ) )
#sys.path.append('../../../main')
#from siteinfo import *

#Define static variables
gBlockMaxSize = 3000
PCAMaxSize = 800
gBlockMinSize = 125
oligoAssemblySize = 100
annealingLength = 20
oligoTM = 48

#Overhangs of the entry vector
entry_ohs = ('TTCT', 'AACG')

#Sequence between entry vector OH and part ohs
entry_extensions = ('GAAGACTC', 'GAGTCTTC') #Changed the right extension to include a BbsI cut site. Not sure why this wasn't included originally?

#part assembly overhangs
partOH_5 = {'1': 'TCAA',
'2a': 'ACGA',
'2b': 'AACG',
'3a': 'CACC',
'3b': 'TTCT',
'3c': 'AGCA',
'3d': 'AGGC',
'3e': 'TCCA',
'4a': 'ATCC',
'4b': 'GTAA',
'5': 'GCTA',
'6': 'ACTC',
'7': 'ATTG'
}

partOH_3 = {'1': 'ACGA',
'2a': 'AACG',
'2b': 'CACC',
'3a': 'TTCT',
'3b': 'AGCA',
'3c': 'AGGC',
'3d': 'TCCA',
'3e': 'ATCC',
'4a': 'GTAA',
'4b': 'GCTA',
'5': 'ACTC',
'6': 'ATTG',
'7': 'TCAA'
}

#Sequence between part and part overhangs
partextension_5 = {'1': '',
'2a': '',
'2b': '',
'3a': '',
'3b': '',
'3c': 'GT',
'3d': '',
'3e': 'GT',
'4a': '',
'4b': 'GGTACC', #Add KpnI
'5': '',
'6': '',
'7': '',
'Custom': '',
}

partextension_3 = {'1': '',
'2a': '',
'2b': 'ACTCTTCTGGTCCCCACAGACTCAGAGAGAACCCGCTAGC', # This is where the standardized 5' UTR and 3' UTR can go
'3a': 'GG',
'3b': '',
'3c': 'TC',
'3d': '',
'3e': 'GG',
'4a': 'TA',
'4b': '',
'5': '',
'6': '',
'7': '',
'Custom': '',
}

multigene_ohs = ["CTGA", "CCAA", "GATG", "GTTC", "GGTA", "AAGT", "AGCA"] # - NEED TO FIX, is this even used?

# Sequence to add to ends of a PCR in order to digest with BsmBI for entry into Stage 1 vector
tails = {"BsmBI":("GCATCGTCTCA", "TGAGACGGCAT"),
		 "BbsI":("GCATGAAGACTC", "AGGTCTTCGCAT")}

# part entry vector is required for GGpart instantiation
# entryVectors = {"OPL4772":"gcatcaaatgaaactgcaatttattcatatcaggattatcaataccatatttttgaaaaagccgtttctgtaatgaaggagaaaactcaccgaggcagttccataggatggcaagatcctggtatcggtctgcgattccgactcgtccaacatcaatacaacctattaatttcccctcgtcaaaaataaggttatcaagtgagaaatcaccatgagtgacgactgaatccggtgagaatggcaaaagtttatgcatttctttccagacttgttcaacaggccagccattacgctcgtcatcaaaatcactcgcatcaaccaaaccgttattcattcgtgattgcgcctgagccagacgaaatacgcgatcgctgttaaaaggacaattacaaacaggaatcgaatgcaaccggcgcaggaacactgccagcgcatcaacaatattttcacctgaatcaggatattcttctaatacctggaatgctgtttttccggggatcgcagtggtgagtaaccatgcatcatcaggagtacggataaaatgcttgatggtcggaagaggcataaattccgtcagccagtttagtctgaccatctcatctgtaacatcattggcaacgctacctttgccatgtttcagaaacaactctggcgcatcgggcttcccatacaagcgatagattgtcgcacctgattgcccgacattatcgcgagcccatttatacccatataaatcagcatccatgttggaatttaatcgcggcctcgacgtttcccgttgaatatggctcatactcttcctttttcaatattattgaagcatttatcagggttattgtctcatgagcggatacatatttgaatgtatttagaaaaataaacaaataggggttccgcgcacatttccccgaaaagtgccagatacctgaaacaaaacccatcgtacggccaaggaagtctccaataactgtgatccaccacaagcgccagggttttcccagtcacgacgttgtaaaacgacggccagtcatgcataatccgcacgcatctggaataaggaagtgccattccgcctgacctttctagagacgttctagagcacagctaacaccacgtcgtccctatctgctgccctaggtctatgagtggttgctggataactttacgggcatgcataaggctcgtataatatattcagggagaccacaacggtttccctctacaaataattttgtttaacttttactagatcacacaggaaagtactagatgcgtaaaggagaagaacttttcactggagttgtcccaattcttgttgaattagatggtgatgttaatgggcacaaattttctgtcagtggagagggtgaaggtgatgcaacatacggaaaacttacccttaaatttatttgcactactggaaaactacctgttccatggccaacacttgtcactactttcggttatggtgttcaatgctttgcgagatacccagatcatatgaaacagcatgactttttcaagagtgccatgcccgaaggttatgtacaggaaagaactatatttttcaaagatgacgggaactacaagacacgtgctgaagtcaagtttgaaggtgatacccttgttaatagaatcgagttaaaaggtattgattttaaagaagatggaaacattcttggacacaaattggaatacaactataactcacacaatgtatacatcatggcagacaaacaaaagaatggaatcaaagttaacttcaaaattagacacaacattgaagatggaagcgttcaactagcagaccattatcaacaaaatactccaattggcgatggccctgtccttttaccagacaaccattacctgtccacacaatctgccctttcgaaagatcccaacgaaaagagagaccacatggtccttcttgagtttgtaacagctgctgggattacacatggcatggatgaactatacaaatagtagtactagagccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttatacgtctcgaacgaggctaggtggaggctcagtgatgataagtctgcgatggtggatgcatgtgtcatggtcatagctgtttcctgtgtgaaattgttatccgctcagagggcacaatcctattccgcgctatccgacaatctccaagacattaggtggagttcagttcggcgtatggcatatgtcgctggaaagaacatgtgagcaaaaggccagcaaaaggccaggaaccgtaaaaaggccgcgttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaagaacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctctattcaacaaagccgccgtcccgtcaagtcagcgtaaatgggtagggggcttcaaatcgtcctcgtgataccaattcggagcctgcttttttgtacaaacttgttgataatggcaattcaaggatcttcacctagatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatatgagtaaacttggtctgacagttagaaaaactcatcga"}
# entryVectorDigests = {"OPL4772":"aggctaggtggaggctcagtgatgataagtctgcgatggtggatgcatgtgtcatggtcatagctgtttcctgtgtgaaattgttatccgctcagagggcacaatcctattccgcgctatccgacaatctccaagacattaggtggagttcagttcggcgtatggcatatgtcgctggaaagaacatgtgagcaaaaggccagcaaaaggccaggaaccgtaaaaaggccgcgttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaagaacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctctattcaacaaagccgccgtcccgtcaagtcagcgtaaatgggtagggggcttcaaatcgtcctcgtgataccaattcggagcctgcttttttgtacaaacttgttgataatggcaattcaaggatcttcacctagatccttttaaattaaaaatgaagttttaaatcaatctaaagtatatatgagtaaacttggtctgacagttagaaaaactcatcgagcatcaaatgaaactgcaatttattcatatcaggattatcaataccatatttttgaaaaagccgtttctgtaatgaaggagaaaactcaccgaggcagttccataggatggcaagatcctggtatcggtctgcgattccgactcgtccaacatcaatacaacctattaatttcccctcgtcaaaaataaggttatcaagtgagaaatcaccatgagtgacgactgaatccggtgagaatggcaaaagtttatgcatttctttccagacttgttcaacaggccagccattacgctcgtcatcaaaatcactcgcatcaaccaaaccgttattcattcgtgattgcgcctgagccagacgaaatacgcgatcgctgttaaaaggacaattacaaacaggaatcgaatgcaaccggcgcaggaacactgccagcgcatcaacaatattttcacctgaatcaggatattcttctaatacctggaatgctgtttttccggggatcgcagtggtgagtaaccatgcatcatcaggagtacggataaaatgcttgatggtcggaagaggcataaattccgtcagccagtttagtctgaccatctcatctgtaacatcattggcaacgctacctttgccatgtttcagaaacaactctggcgcatcgggcttcccatacaagcgatagattgtcgcacctgattgcccgacattatcgcgagcccatttatacccatataaatcagcatccatgttggaatttaatcgcggcctcgacgtttcccgttgaatatggctcatactcttcctttttcaatattattgaagcatttatcagggttattgtctcatgagcggatacatatttgaatgtatttagaaaaataaacaaataggggttccgcgcacatttccccgaaaagtgccagatacctgaaacaaaacccatcgtacggccaaggaagtctccaataactgtgatccaccacaagcgccagggttttcccagtcacgacgttgtaaaacgacggccagtcatgcataatccgcacgcatctggaataaggaagtgccattccgcctgacct"}


class GGpart():
	def __init__(self, partName, leftPartType, rightPartType, seq, part_entry_vector, method=None, fiveprime=None, threeprime=None,
				 possibleTemplates={}, plasmidName="", inputSeqType="DNA", rsToRemove=["BbsI","BsmBI"], removeHomology=False,
				 skipPartPlasmid=False, maxPrimerLength=60, databaseTemplates=[], addOHs=[]):
		"""
		:param part_entry_vector: DNAssembly Plasmid
		"""
		self.partName = partName
		self.leftPartType = leftPartType
		self.rightPartType = rightPartType
		self.inputSeq = seq

		self.method = "None" if method is None else method

		self.plasmidName = plasmidName

		self.inputSeqType = inputSeqType
		self.rsToRemove = rsToRemove
		self.removeHomology = removeHomology

		self.maxPrimerLength = maxPrimerLength
		self.skipPartPlasmid = skipPartPlasmid

		self.possibleTemplates = possibleTemplates
		for each in self.possibleTemplates:
			self.possibleTemplates[each] = self.possibleTemplates[each].replace("\r","").replace("\n","").replace(" ","")

		self.additionalOHs = addOHs

		self.databaseTemplates = databaseTemplates
		self.warnings = []
		self.errors = []

		## These will all be initialized later unless the an error is thrown in which case they will be left blank
		self.vectorName = ""
		self.vectorSeq = ""
		self.vectorDigest = ""
		self.maxOverhangConstraints = 1
		self.plasmidSeq = ""
		self.enzyme = ""
		self.subclone = ""
		self.partSeq = ""
		self.GGfrags = []
		self.fragments = []
		self.removeRS_tuples = None

		if self.checkInput_beforeOptimization():
			#Generate partSeq by removing RS and codon optimizing
			self.partSeq = ""
			self.removeInternalRS()
			self.partSeq = self.partSeq.lower() #change all bases to lowercase

			if self.checkInput_afterOptimization():
				#Initialize GG fragment by making it a single piece
				self._initGGfrag(fiveprime, threeprime)
				#Find the best overhangs
				self.pickGGfrags()

				#Combine adjacent Oligo Assemblies if possible
				if self.method in ["None", "Oligo Assembly"]:
					self.mergeOligoAssemblies()

				self.assembled_plasmid = self.performPartAssembly(part_entry_vector)

				#List of AssemblyInstruction objects
				self.fragments = self.getAssemblyInstructions()


	################ Other functions #################

	# REDUNDANT METHOD! - not used anywhere and can be removed
    # Check to make sure input is in the form of DNA letters
	def _verify_alphabet(sequence):  # todo: what.
		letters = sequence.alphabet.letters
		if not letters:
			raise ValueError("Alphabet does not define letters.")
		for letter in sequence:
			if letter not in letters:
				return False
		return True

	def checkInput_beforeOptimization(self):
		"""Check to ensure the input sequence is DNA and does not break part 1 and 5 construction
		Specifically checks that BsmBI sites are not removed from connector parts.
		"""
		#seq = Seq(self.inputSeq.upper())
		#if not (_verify_alphabet(seq)):
		#    self.errors.append("The DNA sequence you entered in invalid.")
		if self.leftPartType in ["1", "5"] and "BsmBI" in self.rsToRemove:
			self.warnings.append("You chose to remove BsmBI sites from a type connector part. Doing so will break multigene assembly.")
		if self.rightPartType in ["1", "5"] and "BsmBI" in self.rsToRemove:
			self.warnings.append("You chose to remove BsmBI sites from a type connector part. Doing so will break multigene assembly.")
		if len(self.errors) != 0:
			return False
		else:
			return True

	def removeInternalRS(self): #goto removeRS.py
		"""Use custom script to remove defined RS. By default these are defined as BbsI and BsmBI"""
		self.partSeq, self.removeRS_tuples = removeRS(self.inputSeq, self.rsToRemove)

	def checkInput_afterOptimization(self):
		f_seq = FormattedSeq(Seq(self.partSeq), True)

		#Check to make sure BbsI/BsmBI sites aren't present
		if BbsI.search(f_seq) and BsmBI.search(f_seq):
			#comment the line below to let through parts with BbsI and BsmBI
			#be careful, though, these assemblies may be problematic
			#self.errors.append("Your part contains both BbsI and BsmBI sites and cannot be assembled using golden gate.")
			pass
		elif self.leftPartType in ['1', '2a', '2b', '3a', '3b', '3c', '3d', '3e', '4a', '4b', '5', '6', '7'] and BbsI.search(f_seq):
			self.errors.append("Your part contains a BbsI site which must be removed prior to assembly.")

		elif self.rightPartType in ['1', '2a', '2b', '3a', '3b', '3c', '3d', '3e', '4a', '4b', '5', '6', '7'] and BbsI.search(f_seq):
			self.errors.append("Your part contains a BbsI site which must be removed prior to assembly.")

		elif self.leftPartType in ['2a', '2b', '3a', '3b', '3c', '3d', '3e', '4a', '4b', '6', '7'] and BsmBI.search(f_seq):
			self.errors.append("Your part contains a BsmBI site which must be removed prior to assembly.")

		elif self.rightPartType in ['2a', '2b', '3a', '3b', '3c', '3d', '3e', '4a', '4b', '6', '7'] and BsmBI.search(f_seq):
			self.errors.append("Your part contains a BsmBI site which must be removed prior to assembly.")

		#Check to make sure connector parts have BsmBI sites
		bsmBIFor = self.partSeq.upper().count(BsmBI.site)
		bsmBIRev = self.partSeq.upper().count(str(Seq(BsmBI.site).reverse_complement()))
		if str(self.leftPartType) == "1" or str(self.rightPartType) == "1":
			if bsmBIFor + bsmBIRev != 1:
				self.warnings.append("Your type 1 part should have exactly 1 BsmBI site. Consider modifying for multigene assembly.")
		if str(self.leftPartType) == "5" or str(self.rightPartType) == "5":
			if bsmBIFor + bsmBIRev != 1:
				self.warnings.append("Your type 5 part should have exactly 1 BsmBI site. Consider modifying for multigene assembly.")

		#Warn if ORF has a start codon
		if self.leftPartType in ["3a"]:
			if str(self.partSeq.upper())[:3] != "ATG":
				self.warnings.append("Your part is missing a start codon.")
		#Warn if ORF has a stop codon
		if str(self.rightPartType) in ['3a', '3b', '3c', '3d', '3e', '4a']:
			if len(self.partSeq) % 3 != 0:
				self.warnings.append("Your part appears to be out of frame (length is not a multiple of 3). If this is a coding sequence, check to make sure it is correct.")
			if Seq(self.partSeq).translate().find("*") > -1:
				self.warnings.append("Your part has a stop codon. If it is not removed, the part cannot be used for making N-terminal fusions.")

		if len(self.errors) != 0:
			return False
		else:
			return True

	def _initGGfrag(self, fiveprime, threeprime):
		#fiveprime and threeprime are variables that will get passed in the case of custom part design
		#only set to true if you have to gel purify the fragments
		self.subclone = False

		frag_seq = self.partSeq
		#These are the strings that will get appended onto the ends of the fratment
		fiveprimeOH = ""
		threeprimeOH = ""
		fiveprimeExt = ""
		threeprimeExt = ""

		#Custom Part Design
		if self.leftPartType == "Custom" or self.rightPartType == "Custom":
			self.enzyme = "BsmBI"

			fiveprimeOH = entry_ohs[0].upper()
			threeprimeOH = entry_ohs[1].upper()
			fiveprimeExt = entry_extensions[0].upper() + fiveprime[:4].upper() + partextension_5[self.leftPartType]
			threeprimeExt = partextension_3[self.rightPartType] + threeprime[-4:].upper() + entry_extensions[1].upper()

		#Add a section here to catch 4a-4b parts and add a stop codon
		elif self.leftPartType == "4a" and self.rightPartType == "4b":
			self.enzyme = "BsmBI"
			fiveprimeOH = entry_ohs[0].upper()
			threeprimeOH = entry_ohs[1].upper()
			fiveprimeExt = entry_extensions[0].upper() + partOH_5[self.leftPartType] + 'TAA' #Insert additional cut site to enable subcloning here
			threeprimeExt = partextension_3[self.rightPartType] + partOH_3[self.rightPartType] + entry_extensions[1].upper()

		else:
			self.enzyme = "BsmBI"
			fiveprimeOH = entry_ohs[0].upper()
			threeprimeOH = entry_ohs[1].upper()
			fiveprimeExt = entry_extensions[0].upper() + partOH_5[self.leftPartType] + partextension_5[self.leftPartType]
			threeprimeExt = partextension_3[self.rightPartType] + partOH_3[self.rightPartType] + entry_extensions[1].upper()

		#Set vector
		# todo: get rid of this dependency
		self.vectorName = "OPL4772" #Remember to hardcode in this vector

		newFrag = Part.GGfrag(fiveprimeOH, fiveprimeExt, frag_seq, threeprimeExt, threeprimeOH)
		self.GGfrags = [newFrag]

	def pickGGfrags(self):
		allowableSize = 0
		forced_method = ""
		if self.method == "gBlocks":
			allowableSize = gBlockMaxSize - 24
		elif self.method == "Oligo Assembly":
			allowableSize = gBlockMinSize - 3
		elif self.method == "PCA":
			allowableSize = PCAMaxSize - 24

		# Divide the initialized GGfrag into chunks based on assembly method
		if self.method in ["gBlocks", "Oligo Assembly", "PCA"]:
			self.GGfrags[0].forced_method = self.method
			self.GGfrags = divideBySize(self.GGfrags[0], allowableSize)
		else:
			#pdb.set_trace()
			self.GGfrags = divideByIndexTuples(self.GGfrags[0], self.removeRS_tuples) #Uppercase should come from the removeRS
			self._mergeFragments()
			### Split sequence without template into gBlock sized chunks###
			if self.method == "None":
				tempFrags = []
				for index, each in enumerate(self.GGfrags):
					if len(each.seq) > 8:
						results = searchSeqBenchling(each.seq)
						if results:
							for i in range(len(results)):
								if results[i]['entityRegistryId']:
									alias = results[i]['entityRegistryId']
								else:
									continue

								seq = results[i]['bases']
								#linear = not(results[i]['isCircular'])
								self.possibleTemplates[alias] = seq

					if findTemplate(each.seq, self.possibleTemplates):
						each.forced_method = "PCR"
						tempFrags.append(each)
					else:
						if len(each) < gBlockMinSize:
							each.forced_method = "Oligo Assembly"
							tempFrags += [each]
						#If length is just over gBlock length, make the assembly include one oligo assembly to save on a gBlock
						elif 0 < len(each) % (gBlockMaxSize - 24) < gBlockMinSize - 4 * len(each)/(gBlockMaxSize-24):
							oligoAssemSize = max(20, len(each) % (gBlockMaxSize - 24) + 4 * len(each)/(gBlockMaxSize-24))
							if index == 0:
								oligoFrag = Part.GGfrag(each.sequence[:4], "", each.sequence[4:oligoAssemSize], "", "", forced_method="Oligo Assembly")
							else:
								oligoFrag = Part.GGfrag("", "", each.sequence[:oligoAssemSize], "", "", forced_method="Oligo Assembly")
							if index == len(self.GGfrags)-1:
								gBlockFrag = Part.GGfrag("", "", each.sequence[oligoAssemSize:-4], "", each.sequence[-4:], forced_method="gBlocks")
							else:
								gBlockFrag = Part.GGfrag("", "", each.sequence[oligoAssemSize:], "", "", forced_method="gBlocks")
							tempFrags += [oligoFrag]
							tempFrags += divideBySize(gBlockFrag, gBlockMaxSize - 24)
						else:
							each.forced_method = "gBlocks"
							tempFrags += divideBySize(each, gBlockMaxSize - 24)

				self.GGfrags = tempFrags

		# Merge adjacent gBlocks
		def gblock_merge():
			core_seq = ''
			for is_first, is_last, item in more_itertools.mark_ends(merged_seq):
				if is_first:
					fiveprimeOH = item.fiveprimeOH
					extension_5 = item.fiveprimeExt
				if is_last:
					extension_3 = item.threeprimeExt
					threeprimeOH = item.threeprimeOH
				if not is_first:
					core_seq += item.fiveprimeExt
				core_seq += item.seq
				if not is_last:
					core_seq += item.threeprimeExt

			merged_frag = Part.GGfrag(fiveprimeOH=fiveprimeOH, fiveprimeExt=extension_5, seq=core_seq,
									  threeprimeExt=extension_3, threeprimeOH=threeprimeOH, forced_method="gBlocks")
			tempFrags.append(merged_frag)


		if len(self.GGfrags) > 1:
			print('Attempting to merge gBlocks')
			tempFrags = []
			merged_seq = []
			for frag in self.GGfrags:
				if frag.forced_method == "gBlocks":
					merged_seq.append(frag)
				else:
					if len(merged_seq) > 0:
						gblock_merge()
						merged_seq = []
					tempFrags.append(frag)
			if len(merged_seq) > 0:
				gblock_merge()
			self.GGfrags = tempFrags

		self.pickGGfragOHs()
		if self.maxOverhangConstraints > 12:
			self.warnings.append("Due to its complexity, this assembly uses overhangs that are somewhat similar and may be prone to misannealing.")
		for each in self.GGfrags:
			each.assignTails(tails[self.enzyme], self.enzyme)

	#combines items from self.GGfrags where possible based on maxPrimer length
	#works by converting the GGfrags into segments of form [oh+ext, seq, ext+oh]
	#calls mergeSegments to combine them and then turns them back into GGfrags.
	#If I have time, this should be written since it is stupid to convert into segments
	def _mergeFragments(self):
		#Maximum length of primer
		maxPrimerLength = self.maxPrimerLength - 11 #11 = length of type IIs tail
		if len(self.GGfrags) <= 1: #no need to merge
			return
		#a segment is just a ist of 3 sequences, the middle piece is the template and the flanking pieces are extra sequence
		segments = []
		for fragment in self.GGfrags:
			newSeg = [fragment.fiveprimeOH + fragment.fiveprimeExt, fragment.seq, fragment.threeprimeExt + fragment.threeprimeOH]
			segments.append(newSeg)
		mergedSegs = mergeSegments(segments, maxPrimerLength, annealingLength)

		#convert segments into GGfrags
		newFrags = []
		for i in range(len(mergedSegs)):
			fiveprimeOH = ""
			fiveprimeExt = mergedSegs[i][0]
			seq = mergedSegs[i][1]
			threeprimeExt = mergedSegs[i][2]
			threeprimeOH = ""
			if i == 0:
				fiveprimeOH = self.GGfrags[i].fiveprimeOH
				fiveprimeExt = mergedSegs[i][0][4:]
			if i == len(mergedSegs) - 1:
				threeprimeOH = self.GGfrags[-1].threeprimeOH
				threeprimeExt = mergedSegs[i][2][:-4]
			newFrag = Part.GGfrag(fiveprimeOH, fiveprimeExt, seq, threeprimeExt, threeprimeOH)  # todo: do these need to be annotated with the assembly method?
			newFrags.append(newFrag)
		self.GGfrags = newFrags

	#Find the best overhangs
	def pickGGfragOHs(self):
		allowableSize = 0
		if self.method == "gBlocks":
			allowableSize = gBlockMaxSize - 22
		elif self.method == "Oligo Assembly":
			allowableSize = gBlockMinSize

		#Find possible overhangs at each junction
		possOHs = []
		if self.method in ["gBlocks", "Oligo Assembly"]:
			possOHs = findPossOH_byFragLength(self.GGfrags, allowableSize, annealingLength)
		else:
			possOHs = findPossOH_byPrimerLength(self.GGfrags, self.maxPrimerLength - 11, annealingLength, gBlockMaxSize, self.enzyme)
		if possOHs == False:
			self.errors.append("No valid overhangs could be found. Consider increasing the maximum primer length or specifying a different assembly strategy.")
			return False

		#Compile all additional overhangs to avoid conflicts with
		#AHN: I think this is only necessary if you are trying to skip the part assembly step? May comment out
		additionalOHs = list(self.additionalOHs)

		'''
		if self.skipPartPlasmid and self.GGtype != "Circular":
			if self.GGtype in ["mV", "mVa", "mVb", "V", "Va", "Vb"]:
				additionalOHs += list(multigene_ohs)
			else:
				for each in part_ohs:
					if part_ohs[each][0] not in additionalOHs:
						additionalOHs.append(part_ohs[each][0])
			if self.GGfrags[0].fiveprimeOH not in additionalOHs:
				additionalOHs.append(self.GGfrags[0].fiveprimeOH)
			if self.GGfrags[-1].threeprimeOH not in additionalOHs:
				additionalOHs.append(self.GGfrags[-1].threeprimeOH)
		'''

		#Add existing overhangs at beginning and end to the list
		existingOHs = []
		existingOHs.append((self.GGfrags[0].fiveprimeOH, 0))
		existingOHs.append((self.GGfrags[-1].threeprimeOH, 0))

		#Find the best set of overhangs
		constraints = 1
		bestOHs = False
		while not bestOHs:
			if constraints > 16:
				self.errors.append("No valid overhangs could be found")
				return False
			bestOHs = findBestOverhangs(existingOHs, possOHs, additionalOHs, constraints)
			constraints += 1
		#reorder the overhangs so that the 3' OH is last
		bestOHs.append(bestOHs.pop(1))
		if constraints - 1 > self.maxOverhangConstraints:
			self.maxOverhangConstraints = constraints - 1
		self.GGfrags = reindexGGfrags(self.GGfrags, bestOHs)

	def mergeOligoAssemblies(self):
		fragsMerged = True
		while fragsMerged:
			fragsMerged = False
			for i in range(len(self.GGfrags)-1):
				currFrag = self.GGfrags[i]
				nextFrag = self.GGfrags[i+1]
				if (len(currFrag) + len(nextFrag)) < oligoAssemblySize:
					fiveprimeOH = currFrag.fiveprimeOH
					fiveprimeExt = currFrag.fiveprimeExt + currFrag.seq + currFrag.threeprimeExt + currFrag.threeprimeOH + nextFrag.fiveprimeExt
					seq = nextFrag.seq
					threeprimeExt = nextFrag.threeprimeExt
					threeprimeOH = nextFrag.threeprimeOH
					mergedFrag = Part.GGfrag(fiveprimeOH, fiveprimeExt, seq, threeprimeExt, threeprimeOH)
					self.GGfrags.pop(i)
					self.GGfrags[i] = mergedFrag
					fragsMerged = True
					break

	def performPartAssembly(self, part_entry_vector):
		"""
		Use StickyEndAssembly to perform assembly
		Assumes BsmBI for part assembly

		:param part_entry_vector: Plasmid instance for part entry vector
		"""
		input_fragments = copy.deepcopy(self.GGfrags)
		part_assembly = StickyEndAssembly(input_fragments + [part_entry_vector], [BsmBI])
		part_assembly.digest()
		assembled_part = part_assembly.perform_assembly()
		return assembled_part

	def getAssemblyInstructions(self):
		assemIns = []
		primers, methods = self.getPrimersAndMethods()
		# pdb.set_trace()

		for primer_set, method, fragment in zip(primers, methods, self.GGfrags):
			assemIns.append(AssemblyInstructions(method, primer_set, fragment, self.possibleTemplates))
		return assemIns

	def getPrimersAndMethods(self): #Can a single fragment have multiple assembly methods/primers?
		primers = []
		methods = []
		for each in self.GGfrags:
			if each.forced_method == "gBlocks":
				primers.append(['',''])
				methods.append("gBlocks")
			elif each.forced_method == "Oligo Assembly":
				primers.append(getOligoAssemPrimers(each, self.maxPrimerLength))
				methods.append("Oligo Assembly")
			elif each.forced_method == "PCR":
				primers.append(getPCRprimers(each, self.maxPrimerLength, annealingLength, oligoTM=oligoTM))
				methods.append("PCR")
				if not findTemplate(each.seq, self.possibleTemplates):
					message = "PCR templates could not be found for this assembly. If you have a template, add it to the templates text box. If you don't have a template, consider using an assembly method other than PCR."
					if message not in self.errors:
						self.errors.append(message)
			else:
				if not self.possibleTemplates: #Change - if no template, then
					if len(each) <= gBlockMinSize:
						primers.append(getOligoAssemPrimers(each, self.maxPrimerLength))
						methods.append("Oligo Assembly")
					#elif each.getLength() + len(tails[self.enzyme][0]) + len(tails[self.enzyme][1]) >  gBlockMaxSize:
						#primers.append(each.getPCAprimers(self.maxPrimerLength))
						#methods.append("PCA")
					else:
						primers.append(['',''])
						methods.append("gBlocks")
				else:
					if len(each) <= oligoAssemblySize:
						primers.append(getOligoAssemPrimers(each, self.maxPrimerLength))
						methods.append("Oligo Assembly")
					else:
						primers.append(getPCRprimers(each, self.maxPrimerLength, annealingLength, oligoTM=oligoTM))
						methods.append("PCR")
		return primers, methods
