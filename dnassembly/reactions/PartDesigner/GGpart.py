"""
GGpart.py

Modified from code originally written by Will DeLoache.
"""

import sys, os
from GGfrag import GGfrag
from Bio.Alphabet import IUPAC
from MergeSegments import mergeSegments
from OHfinder import *
from Bio.Seq import Seq
from Bio.Restriction import BsmBI, BsaI, Restriction, FormattedSeq
from AssemblyInstructions import AssemblyInstructions, findTemplate
from RemoveRS import removeRS
from CodonOptimize.codonOptimize import codonOptimize
from Bio import SearchIO, SeqIO
from findTemplates import findTemplates, reindexCircular
from makeBLASTdb import makeBLASTdb
import MySQLdb

# BASE_PATH is the absolute path of ../../.. relative to this script location
BASE_PATH = reduce (lambda l,r: l + os.path.sep + r, os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep )[:-3] )
sys.path.append( os.path.join( BASE_PATH, "main" ) )
sys.path.append('../../../main')
from siteinfo import *

#from XMLoutput import getXML

#Overhangs of the entry vector - NEED TO FIX
entry_ohs = {"1": ("TCGG", "GACC"),
			"234": ("TCGG", "GACC"),
			"2": ("TCGG", "GACC"),
			"3": ("TCGG", "GACC"),
			"3alpha": ("TCGG", "GACC"),
			"3beta": ("TCGG", "GACC"),
			"3a": ("TCGG", "GACC"), #new TTK
			"3b": ("TCGG", "GACC"), #new TTK
			"3c": ("TCGG", "GACC"), #new TTK
			"3d": ("TCGG", "GACC"), #new TTK
			"3d_1": ("TCGG", "GACC"), #expanded ICD
			"3d_2": ("TCGG", "GACC"), #expanded ICD
			"4": ("TCGG", "GACC"),
			"4a": ("TCGG", "GACC"),
			"4b": ("TCGG", "GACC"),
			"5": ("TCGG", "GACC"),
			"67": ("TCGG", "GACC"),
			"6": ("TCGG", "GACC"),
			"6a": ("TCGG", "GACC"),
			"6b": ("TCGG", "GACC"),
			"7": ("TCGG", "GACC"),
			"6b7": ("TCGG", "GACC"),
			"7a": ("TCGG", "GACC"),
			"7b": ("TCGG", "GACC"),
			"234r": ("TCGG", "GACC"),
			"12a": ("TCGG", "GACC"),
			"45": ("TCGG", "GACC")}

#Sequence between entry vector OH and part ohs - NEED TO FIX
entry_extensions = {"1": ("TCTCA", "TGA"),
					"234": ("TCTCA", "TGA"),
					"2": ("TCTCA", "TGA"),
					"3": ("TCTCA", "TGA"),
					"3alpha": ("TCTCA", "TGA"),
					"3beta": ("TCTCA", "TGA"),
					"3a": ("TCTCA", "TGA"), #new TTK
					"3b": ("TCTCA", "TGA"), #new TTK
					"3c": ("TCTCA", "TGA"), #new TTK
					"3d": ("TCTCA", "TGA"), #new TTK
					"3d_1": ("TCTCA", "TGA"), #expanded ICD
					"3d_2": ("TCTCA", "TGA"), #expanded ICD
					"4": ("TCTCA", "TGA"),
					"4a": ("TCTCA", "TGA"),
					"4b": ("TCTCA", "TGA"),
					"5": ("TCTCA", "TGA"),
					"67": ("TCTCA", "TGA"),
					"6": ("TCTCA", "TGA"),
					"6a": ("TCTCA", "TGA"),
					"6b": ("TCTCA", "TGA"),
				    "6b7": ("TCTCA", "TGA"),
					"7": ("TCTCA", "TGA"),
					"7a": ("TCTCA", "TGA"),
					"7b": ("TCTCA", "TGA"),
					"234r": ("", ""),
				    "12a": ("TCTCA", "TGA"),
			        "45": ("TCTCA", "TGA")}

#part assembly overhangs - NEED TO FIX
part_ohs = {"1": ("CCCT", "AACG"),
			"234": ("AACG", "GCTG"),
			"2": ("AACG", "CATG"), #updated to create concensus Kozak
			"3": ("CATG", "ATCC"),
			"3alpha": ("CATG", "TTCT"), #original YTK overhangs
			"3beta": ("TTCT", "ATCC"), #original YTK overhangs
			"3a": ("CATG", "TTCT"), #new TTK
			"3b": ("TTCT", "GGAA"), #new TTK
			"3c": ("GGAA", "ATCA"), #new TTK
			"3d": ("ATCA", "ATCC"), #new TTK
			"3d_1": ("ATCA", "CTCC"), #expanded ICD
			"3d_2": ("CTCC", "ATCC"), #expanded ICD
			"4": ("ATCC", "GCTG"),
			"4a": ("ATCC", "TGGC"),
			"4b": ("TGGC", "GCTG"),
			"5": ("GCTG", "TACA"),
			"67": ("TACA", "CCCT"),
			"6": ("TACA", "CCGA"),
			"6a": ("TACA", "GAGT"),
			"6b": ("GAGT", "CCGA"),
			"6b7": ("GAGT", "CCCT"),
			"7": ("CCGA", "CCCT"),
			"7a": ("CCGA", "CAAT"),
			"7b": ("CAAT", "CCCT"),
			"234r": ("AACG", "GCTG"),
		    "12a": ("CCCT", "CGAC"),
			"45": ("ATCC", "TACA")}

#Sequence between part and part overhangs - NEED TO FIX
part_extensions = {"1": ("GAATTCGCATCTAGA", ""),
					"234": ("", ""),
					"2": ("", "GCCGCCAC"), #Removed BglII cut site, add Kozak concensus
					"3": ("", "GG"),
					"3alpha": ("", "GG"),
					"3beta": ("", "GG"),
					"3a": ("", "GG"), #new TTK
					"3b": ("", ""), #new TTK
					"3c": ("GT", "GG"), #new TTK
					"3d": ("", "GG"), #new TTK
					"3d_1": ("", "GG"), #expanded ICD
					"3d_2": ("", "GG"), #expanded ICD
					"4": ("TAA", ""),
					"4a": ("", "TAA"),
					"4b": ("", ""),
					"5": ("", "ACTAGTGCACTGCAG"),
					"67": ("GCGGCCGC", "GCGGCCGC"), #add NotI sites
					"6": ("", ""),
					"6a": ("", ""),
					"6b": ("", ""),
					"6b7": ("", ""),
					"7": ("GCGGCCGC", "GCGGCCGC"),#add NotI sites
					"7a": ("GCGGCCGC", "GCGGCCGC"),#add NotI sites
					"7b": ("", ""),
					"234r": ("TGAGACC", "GGTCTCA"),
		                        "12a": ("GAATTCGCATCTAGA", ""),
		                 	"45": ("TAA", "ACTAGTGCACTGCAG")}

multigene_ohs = ["CTGA", "CCAA", "GATG", "GTTC", "GGTA", "AAGT", "AGCA"] # - NEED TO FIX

#Sequence to add to ends of a PCR in order to digest with a type IIs enzyme - NEED TO FIX
tails = {"BsmBI":("GCATCGTCTCA", "TGAGACGGCAT"),
		"BsaI":("GCATGGTCTCA", "TGAGACCGCAT")}

#NEED TO CHANGE TO OPL
entryVectors = {"pML967":"cgtctctgaccagaccaataaaaaacgcccggcggcaaccgagcgttctgaacaaatccagatggagttctgaggtcattactggatctatcaacaggagtccaagcgagctcgatatcaaattacgccccgccctgccactcatcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatggtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgaaacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagcgatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccgtctttcattgccatacgaaattccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccattttagcttccttagctcctgaaaatctcgataactcaaaaaatacgcccggtagtgatcttatttcattatggtgaaagttggaacctcttacgtgcccgatcaatcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtagtcggcgagacggaaagtgaaacgtgatttcatgcgtcattttgaacattttgtaaatcttatttaataatgtgtgcggcaattcacatttaatttatgaatgttttcttaacatcgcggcaactcaagaaacggcaggttcggatcttagctactagagaaagaggagaaatactagatgcgtaaaggcgaagagctgttcactggtgtcgtccctattctggtggaactggatggtgatgtcaacggtcataagttttccgtgcgtggcgagggtgaaggtgacgcaactaatggtaaactgacgctgaagttcatctgtactactggtaaactgccggttccttggccgactctggtaacgacgctgacttatggtgttcagtgctttgctcgttatccggaccatatgaagcagcatgacttcttcaagtccgccatgccggaaggctatgtgcaggaacgcacgatttcctttaaggatgacggcacgtacaaaacgcgtgcggaagtgaaatttgaaggcgataccctggtaaaccgcattgagctgaaaggcattgactttaaagaagacggcaatatcctgggccataagctggaatacaattttaacagccacaatgtttacatcaccgccgataaacaaaaaaatggcattaaagcgaattttaaaattcgccacaacgtggaggatggcagcgtgcagctggctgatcactaccagcaaaacactccaatcggtgatggtcctgttctgctgccagacaatcactatctgagcacgcaaagcgttctgtctaaagatccgaacgagaaacgcgatcatatggttctgctggagttcgtaaccgcagcgggcatcacgcatggtatggatgaactgtacaaatgaccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata",
				"pWCD1669":"cgtctctgaccagtccctatcagtgatagagattgacatccctatcagtgatagagatactgagcacggatctgaaagaggagaaaggatctatggcgagtagcgaagacgttatcaaagagttcatgcgtttcaaagttcgtatggaaggttccgttaacggtcacgagttcgaaatcgaaggtgaaggtgaaggtcgtccgtacgaaggtactcagaccgctaaactgaaagttaccaaaggtggtccgctgccgttcgcttgggacatcctgtccccgcagttccagtacggttccaaagcttacgttaaacacccggctgacatcccggactacctgaaactgtccttcccggaaggtttcaaatgggaacgtgttatgaacttcgaagacggtggtgttgttaccgttacccaggactcctccctgcaagacggtgagttcatctacaaagttaaactgcgtggtactaacttcccgtccgacggtccggttatgcagaaaaaaaccatgggttgggaagcttccaccgaacgtatgtacccggaagacggtgctctgaaaggtgaaatcaaaatgcgtctgaaactgaaagacggtggtcactacgacgctgaagttaaaaccacctacatggctaaaaaaccggttcagctgccgggtgcttacaaaaccgacatcaaactggacatcacctcccacaacgaagactacaccatcgttgaacagtacgaacgtgctgaaggtcgtcactccaccggtgcttaataaggatctccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttataagtcggcgagacgaccaataaaaaacgcccggcggcaaccgagcgttctgaacaaatccagatggagttctgaggtcattactggatctatcaacaggagtccaagcgagctcgatatcaaattacgccccgccctgccactcatcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatggtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgaaacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagcgatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccgtctttcattgccatacgaaattccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccattttagcttccttagctcctgaaaatctcgataactcaaaaaatacgcccggtagtgatcttatttcattatggtgaaagttggaacctcttacgtgcccgatcaatcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgt",
				"pWCD2515":"cgtctctgaccatgaaagtgaaacgtgatttcatgcgtcattttgaacattttgtaaatcttatttaataatgtgtgcggcaattcacatttaatttatgaatgttttcttaacatcgcggcaactcaagaaacggcaggttcggatcttagctactagagaaagaggagaaatactagatgcgtaaaggcgaagagctgttcactggtgtcgtccctattctggtggaactggatggtgatgtcaacggtcataagttttccgtgcgtggcgagggtgaaggtgacgcaactaatggtaaactgacgctgaagttcatctgtactactggtaaactgccggttccttggccgactctggtaacgacgctgacttatggtgttcagtgctttgctcgttatccggaccatatgaagcagcatgacttcttcaagtccgccatgccggaaggctatgtgcaggaacgcacgatttcctttaaggatgacggcacgtacaaaacgcgtgcggaagtgaaatttgaaggcgataccctggtaaaccgcattgagctgaaaggcattgactttaaagaagacggcaatatcctgggccataagctggaatacaattttaacagccacaatgtttacatcaccgccgataaacaaaaaaatggcattaaagcgaattttaaaattcgccacaacgtggaggatggcagcgtgcagctggctgatcactaccagcaaaacactccaatcggtgatggtcctgttctgctgccagacaatcactatctgagcacgcaaagcgttctgtctaaagatccgaacgagaaacgcgatcatatggttctgctggagttcgtaaccgcagcgggcatcacgcatggtatggatgaactgtacaaatgaccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttataagtcggcgagacgaccaataaaaaacgcccggcggcaaccgagcgttctgaacaaatccagatggagttctgaggtcattactggatctatcaacaggagtccaagcgagctcgatatcaaattacgccccgccctgccactcatcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatggtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgaaacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagcgatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccgtctttcattgccatacgaaattccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccattttagcttccttagctcctgaaaatctcgataactcaaaaaatacgcccggtagtgatcttatttcattatggtgaaagttggaacctcttacgtgcccgatcaatcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgt"
				}

entryVectorDigests = {"pML967":"agaccaataaaaaacgcccggcggcaaccgagcgttctgaacaaatccagatggagttctgaggtcattactggatctatcaacaggagtccaagcgagctcgatatcaaattacgccccgccctgccactcatcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatggtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgaaacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagcgatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccgtctttcattgccatacgaaattccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccattttagcttccttagctcctgaaaatctcgataactcaaaaaatacgcccggtagtgatcttatttcattatggtgaaagttggaacctcttacgtgcccgatcaatcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtag",
					"pWCD1669":"agtccctatcagtgatagagattgacatccctatcagtgatagagatactgagcacggatctgaaagaggagaaaggatctatggcgagtagcgaagacgttatcaaagagttcatgcgtttcaaagttcgtatggaaggttccgttaacggtcacgagttcgaaatcgaaggtgaaggtgaaggtcgtccgtacgaaggtactcagaccgctaaactgaaagttaccaaaggtggtccgctgccgttcgcttgggacatcctgtccccgcagttccagtacggttccaaagcttacgttaaacacccggctgacatcccggactacctgaaactgtccttcccggaaggtttcaaatgggaacgtgttatgaacttcgaagacggtggtgttgttaccgttacccaggactcctccctgcaagacggtgagttcatctacaaagttaaactgcgtggtactaacttcccgtccgacggtccggttatgcagaaaaaaaccatgggttgggaagcttccaccgaacgtatgtacccggaagacggtgctctgaaaggtgaaatcaaaatgcgtctgaaactgaaagacggtggtcactacgacgctgaagttaaaaccacctacatggctaaaaaaccggttcagctgccgggtgcttacaaaaccgacatcaaactggacatcacctcccacaacgaagactacaccatcgttgaacagtacgaacgtgctgaaggtcgtcactccaccggtgcttaataaggatctccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttataag",
					"pWCD2515":"atgaaagtgaaacgtgatttcatgcgtcattttgaacattttgtaaatcttatttaataatgtgtgcggcaattcacatttaatttatgaatgttttcttaacatcgcggcaactcaagaaacggcaggttcggatcttagctactagagaaagaggagaaatactagatgcgtaaaggcgaagagctgttcactggtgtcgtccctattctggtggaactggatggtgatgtcaacggtcataagttttccgtgcgtggcgagggtgaaggtgacgcaactaatggtaaactgacgctgaagttcatctgtactactggtaaactgccggttccttggccgactctggtaacgacgctgacttatggtgttcagtgctttgctcgttatccggaccatatgaagcagcatgacttcttcaagtccgccatgccggaaggctatgtgcaggaacgcacgatttcctttaaggatgacggcacgtacaaaacgcgtgcggaagtgaaatttgaaggcgataccctggtaaaccgcattgagctgaaaggcattgactttaaagaagacggcaatatcctgggccataagctggaatacaattttaacagccacaatgtttacatcaccgccgataaacaaaaaaatggcattaaagcgaattttaaaattcgccacaacgtggaggatggcagcgtgcagctggctgatcactaccagcaaaacactccaatcggtgatggtcctgttctgctgccagacaatcactatctgagcacgcaaagcgttctgtctaaagatccgaacgagaaacgcgatcatatggttctgctggagttcgtaaccgcagcgggcatcacgcatggtatggatgaactgtacaaatgaccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttataag"
					}

gBlockMaxSize = 3000
PCAMaxSize = 800
gBlockMinSize = 125
oligoAssemblySize = 100
annealingLength = 20
oligoTM = 48

"""
Public Variables:

plasmidName //optional
plasmidSeq
GGtype //1-6, mV, Circular, Custom

inputSeq
inputSeqType //DNA or Protein
partName
partSeq

rsToRemove //list of enzyme names
maxPrimerLength
skipPartPlasmid //True if going straight to BsaI step, else False
method //PCR, gBlocks, Oligo Assembly, or None (will pick best method)

warnings //list of warning messages
errors //list of error messages

maxOverhangConstraints //the constraints used to find the GG overhangs (1-16)
vectorName
vectorSeq
vectorDigest

enzyme //enzyme used in the digests and GG reactions
fragments //list of AssemblyInstructions objects
"""


class GGpart():
	def __init__(self, partName, GGtype, seq, plasmidName="", inputSeqType="DNA", rsToRemove=[], removeHomology=False,
				skipPartPlasmid=False, fiveprime="", threeprime="", method="Any", maxPrimerLength=60, possibleTemplates={}, databaseTemplates=[], addOHs=[]):

		self.partName = partName
		self.GGtype = GGtype
		self.plasmidName = plasmidName

		self.inputSeq = seq
		self.inputSeqType = inputSeqType
		self.rsToRemove = rsToRemove
		self.removeHomology = removeHomology

		self.maxPrimerLength = maxPrimerLength
		self.skipPartPlasmid = skipPartPlasmid
		self.method = method

		self.possibleTemplates = possibleTemplates
		for each in self.possibleTemplates:
			self.possibleTemplates[each]["seq"] = self.possibleTemplates[each]["seq"].replace("\r","").replace("\n","").replace(" ","")

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

		if self.checkInput_beforeOptimization():
			#Generate partSeq by removing RS and codon optimizing
			self.partSeq = ""

			self.partSeq = self.partSeq.lower()
			self.removeInternalRS()

			if self.checkInput_afterOptimization():

				#Initialize GG fragment by making it a single piece
				self._initGGfrag(fiveprime, threeprime)

				#Find the best overhangs
				self.pickGGfrags()

				#Combine adjacent Oligo Assemblies if possible
				if self.method in ["None", "Oligo Assembly"]:
					self.mergeOligoAssemblies()
				try:
					self.vectorDigest = entryVectorDigests[self.vectorName]
					self.vectorSeq = entryVectors[self.vectorName]
				except:
					self.vectorDigest = ""
					self.vectorSeq = ""
				self.plasmidSeq = self.getPlasmidSeq()
				#List of AssemblyInstruction objects
				self.fragments = self.getAssemblyInstructions()

	def __str__(self):
		output = ""
		output += "Name: " + self.partName + "\n"
		output += "GGtype: " + self.GGtype + "\n"
		output += "Enzyme: " + self.enzyme + "\n"
		if self.subclone:
			output += "Method: subclone with " + self.vectorName + "\n"
		else:
			output += "Method: one-pot golden gate with " + self.vectorName + "\n"
		output += "Sequence: " + self.partSeq + "\n"
		output += "\nGolden Gate Fragments:\n\n"
		for each in self.GGfrags:
			output += str(each) + "\n\n"
		return output

	def codonOptimize(self):
		if self.removeHomology == True:
			self.partSeq = codonOptimize(self.inputSeq, windowLen=12, alignmentThreshold=10)
		else:
			self.partSeq = codonOptimize(self.inputSeq, windowLen=9, alignmentThreshold=10)

	def removeInternalRS(self):
		self.partSeq = removeRS(self.inputSeq, self.rsToRemove)

	def mergeOligoAssemblies(self):
		fragsMerged = True
		while fragsMerged:
			fragsMerged = False
			for i in range(len(self.GGfrags)-1):
				currFrag = self.GGfrags[i]
				nextFrag = self.GGfrags[i+1]
				if (currFrag.getLength() + nextFrag.getLength()) < oligoAssemblySize:
					fiveprimeOH = currFrag.fiveprimeOH
					fiveprimeExt = currFrag.fiveprimeExt + currFrag.seq + currFrag.threeprimeExt + currFrag.threeprimeOH + nextFrag.fiveprimeExt
					seq = nextFrag.seq
					threeprimeExt = nextFrag.threeprimeExt
					threeprimeOH = nextFrag.threeprimeOH
					mergedFrag = GGfrag(fiveprimeOH, fiveprimeExt, seq, threeprimeExt, threeprimeOH)
					self.GGfrags.pop(i)
					self.GGfrags[i] = mergedFrag
					fragsMerged = True
					break

	def getAssemblyInstructions(self):
		assemIns = []
		results = self.getPrimersAndMethods()
		primers = results[0]
		methods = results[1]

# 		##### Add unique cutters surrounding each PCA fragment ###########
# 		if self.method == "PCA":
# 			tails[self.enzyme] = (tails[self.enzyme][0][:4] + "GAATTCA" + tails[self.enzyme][0][4:], tails[self.enzyme][1][:-4] + "AGGATCC" + tails[self.enzyme][1][-4:])
# 		##################################################################

		for i in range(len(self.GGfrags)):

			assemIns.append(AssemblyInstructions(methods[i], primers[i], self.GGfrags[i], self.possibleTemplates))
		return assemIns

	def getPrimersAndMethods(self):
		primers = []
		methods = []
		for each in self.GGfrags:
			if self.method == "gBlocks":
				primers.append([])
				methods.append("gBlocks")
			elif self.method == "Oligo Assembly":
				primers.append(each.getOligoAssemPrimers(self.maxPrimerLength))
				methods.append("Oligo Assembly")
			#elif self.method == "PCA":
			#	primers.append(each.getPCAprimers(self.maxPrimerLength))
			#	methods.append("PCA")
            #Check to make sure PCR templates are present
            elif self.method == "PCR":
				primers.append(each.getPCRprimers(self.maxPrimerLength, annealingLength, oligoTM=oligoTM))
				methods.append("PCR")
				if not findTemplate(each.seq, self.possibleTemplates):
					message = "PCR templates could not be found for this assembly. If you have a template, add it to the templates text box. If you don't have a template, consider using an assembly method other than PCR."
					if message not in self.errors:
						self.errors.append(message)
			elif self.method == "None":
				if not findTemplate(each.seq, self.possibleTemplates):
					if each.getLength() <= gBlockMinSize:
						primers.append(each.getOligoAssemPrimers(self.maxPrimerLength))
						methods.append("Oligo Assembly")
					elif each.getLength() + len(tails[self.enzyme][0]) + len(tails[self.enzyme][1]) >  gBlockMaxSize:
						primers.append(each.getPCAprimers(self.maxPrimerLength))
						methods.append("PCA")
					else:
						primers.append([])
						methods.append("gBlocks")
				else:
					if each.getLength() <= oligoAssemblySize:
						primers.append(each.getOligoAssemPrimers(self.maxPrimerLength))
						methods.append("Oligo Assembly")
					else:
						primers.append(each.getPCRprimers(self.maxPrimerLength, annealingLength, oligoTM=oligoTM))
						methods.append("PCR")
		return (primers, methods)


	def getPlasmidSeq(self):
		s = ""
		for index, each in enumerate(self.GGfrags):
			s += each.fiveprimeOH + each.fiveprimeExt + each.seq + each.threeprimeExt
			if index == len(self.GGfrags) - 1 and self.GGtype != "Circular":
				s += each.threeprimeOH
		if self.vectorName == "None":
			#temporary solution until we support skipPartPlasmid better
			if self.skipPartPlasmid:
				s = tails[self.enzyme][0] + s + tails[self.enzyme][1]
		elif self.vectorName in entryVectorDigests.keys():
			s += entryVectorDigests[self.vectorName]
		else:
			vector = ""
			if self.GGtype in ["1","234","2","2a","2b","3","3a","3ab","3ac","3ad","3ae","3b","3bc","3bd","3be","3c","3cd","3ce","3d","3de","3e","4","4a","4b","5","6","7"]:
				vector = "OPL4772"

			s = entry_ohs[self.GGtype][0] + entry_extensions[self.GGtype][0] + s
			s += entry_extensions[self.GGtype][1] + entry_ohs[self.GGtype][1]
			s += entryVectorDigests[vector]
			self.vectorDigest = self.GGfrags[-1].threeprimeOH + entry_extensions[self.GGtype][1] + entry_ohs[self.GGtype][1] + entryVectorDigests[vector] + entry_ohs[self.GGtype][0] + entry_extensions[self.GGtype][0] + self.GGfrags[0].fiveprimeOH
			self.vectorSeq = "NNNNNNNNNN" + self.vectorDigest #mike changed to self.vectorDigest
		return s

	def pickGGfrags(self):
		allowableSize = 0
		forced_method = ""
		if self.method == "gBlocks":
			allowableSize = gBlockMaxSize - 24
		elif self.method == "Oligo Assembly":
			allowableSize = gBlockMinSize - 3
		elif self.method == "PCA":
			allowableSize = PCAMaxSize - 24

		#Divide the initialized GGfrag into chunks based on assembly method
		if self.method in ["gBlocks", "Oligo Assembly", "PCA"]:
			self.GGfrags[0].forced_method = self.method
			self.GGfrags = self.GGfrags[0].divideBySize(allowableSize)
		else:
		##### Run findTemplates to find PCR templates in database####################
        '''Skip this functionality for now - in the future we may want to add the ability to use the API to search Benchling
			dbpath = makeBLASTdb(self.possibleTemplates, self.databaseTemplates)

			results = findTemplates(self.GGfrags[0].seq, dbpath)
			self.GGfrags[0].seq = results[0]
			con = MySQLdb.connect(DATABASE_HOST, DATABASE_USER, DATABASE_PASSWD, DATABASE_NAME);
			cur = con.cursor(MySQLdb.cursors.DictCursor)
			for templateID in results[1].keys():
				try:
					cur.execute("SELECT * FROM Plasmids WHERE (ID=%s)" % templateID)
					plasmid = cur.fetchone()
					seq = plasmid["Sequence"]
					if plasmid["Linear"] == 1:
						self.possibleTemplates[plasmid["Name"]] = {"seq":plasmid["Sequence"], "Linear":True}
					else:
						self.possibleTemplates[plasmid["Name"]] = {"seq":plasmid["Sequence"], "Linear":False}
				except:
					if templateID in ["S288C_Genomic_DNA", "DH10B_Genomic_DNA", "GS115_Genomic_DNA"]:
						self.possibleTemplates[templateID] = {"seq":results[1][templateID], "Linear":True}
        '''
		################### End Find Templates from Database ######################
			self.GGfrags = self.GGfrags[0].divideByCaseChange()
			self._mergeFragments()
			### Split sequence without template into gBlock sized chunks###
			if self.method == "None":
				tempFrags = []
				for index, each in enumerate(self.GGfrags):
					#if not findTemplate(each.seq, self.possibleTemplates):
					if each.getLength() > 3100:
						each.forced_method = "PCA"
						tempFrags += each.divideBySize(PCAMaxSize - 24)
					elif each.getLength() < gBlockMinSize:
						each.forced_method = "Oligo Assembly"
						tempFrags += [each]
					#If length is just over gBlock length, make the assembly include one oligo assembly to save on a gBlock
					elif 0 < each.getLength() % (gBlockMaxSize - 24) < gBlockMinSize - 4 * each.getLength()/(gBlockMaxSize-24):
						oligoAssemSize = max(20, each.getLength() % (gBlockMaxSize - 24) + 4 * each.getLength()/(gBlockMaxSize-24))
						if index == 0:
							oligoFrag = GGfrag(each.getFragSeq()[:4], "", each.getFragSeq()[4:oligoAssemSize], "", "", forced_method="Oligo Assembly")
						else:
							oligoFrag = GGfrag("", "", each.getFragSeq()[:oligoAssemSize], "", "", forced_method="Oligo Assembly")
						if index == len(self.GGfrags)-1:
							gBlockFrag = GGfrag("", "", each.getFragSeq()[oligoAssemSize:-4], "", each.getFragSeq()[-4:], forced_method="gBlocks")
						else:
							gBlockFrag = GGfrag("", "", each.getFragSeq()[oligoAssemSize:], "", "", forced_method="gBlocks")
						tempFrags += [oligoFrag]
						tempFrags += gBlockFrag.divideBySize(gBlockMaxSize - 24)
					else:
						each.forced_method = "gBlocks"
						tempFrags += each.divideBySize(gBlockMaxSize - 24)
					#else:
						#tempFrags.append(each)
				self.GGfrags = tempFrags
		self.pickGGfragOHs()
		if self.maxOverhangConstraints > 12:
			self.warnings.append("Due to its complexity, this assembly uses overhangs that are somewhat similar and may be prone to misannealing.")
		for each in self.GGfrags:
			each.assignTails(tails[self.enzyme], self.enzyme)

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

	#combines items from self.GGfrags where possible based on maxPrimer length
	#works by converting the GGfrags into segments of form [oh+ext, seq, ext+oh]
	#calls mergeSegments to combine them and then turns them back into GGfrags.
	#If I have time, this should be written since it is stupid to convert into segments
	def _mergeFragments(self):
		#Maximum length of primer
		maxPrimerLength = self.maxPrimerLength - 11 #11 = length of type IIs tail
		if len(self.GGfrags) <= 1:
			return
		#a segment is just a ist of 3 sequences, the middle piece is the template and the flanking pieces are extra sequence
		segments = []
		for i in range(len(self.GGfrags)):
			newSeg = [self.GGfrags[i].fiveprimeOH + self.GGfrags[i].fiveprimeExt, self.GGfrags[i].seq, self.GGfrags[i].threeprimeExt + self.GGfrags[i].threeprimeOH]
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
			newFrag = GGfrag(fiveprimeOH, fiveprimeExt, seq, threeprimeExt, threeprimeOH)
			newFrags.append(newFrag)
		self.GGfrags = newFrags

	def initCircularOHs(self):
		dbpath = makeBLASTdb(self.possibleTemplates, self.databaseTemplates)
		self.partSeq = reindexCircular(self.partSeq, dbpath)
		frag1 = GGfrag("","",self.partSeq[-500:],"","")
		frag2 = GGfrag("","",self.partSeq[:500],"","")
		possOHs = findPossOH_byPrimerLength([frag1, frag2], self.maxPrimerLength - 11, annealingLength, gBlockMaxSize, self.enzyme)

		existingOHs = []
		additionalOHs = list(self.additionalOHs)
		constraints = 1
		bestOHs = False
		while not bestOHs:
			if constraints > 16:
				self.errors.append("No valid overhangs could be found")
				return False
			self.maxOverhangConstraints = constraints
			bestOHs = findBestOverhangs(existingOHs, possOHs, additionalOHs, constraints)
			constraints += 1
		ohIndex = bestOHs[0][1]
		if ohIndex >= 0:
			self.partSeq = self.partSeq[ohIndex:] + self.partSeq[:ohIndex]
		else:
			self.partSeq = self.partSeq[ohIndex:] + self.partSeq[:ohIndex]

		frag_seq = self.partSeq[4:]
		fiveprimeOH = self.partSeq[:4].upper()
		threeprimeOH = self.partSeq[:4].upper()
		fiveprimeExt = ""
		threeprimeExt = ""
		newFrag = GGfrag(fiveprimeOH, fiveprimeExt, frag_seq, threeprimeExt, threeprimeOH)
		self.GGfrags = [newFrag]

	#Generate a single GGfrag for the whole part before dividing it up into subfrags
	def _initGGfrag(self, fiveprime, threeprime):
		#only set to true if you have to gel purify the fragments
		self.subclone = False

		frag_seq = self.partSeq
		fiveprimeOH = ""
		threeprimeOH = ""
		fiveprimeExt = ""
		threeprimeExt = ""
		#This are tuples of the various 5' and 3' sequence add-ons (5', 3')
		part_oh = ()
		part_extention = ()
		entry_oh = ()
		entry_extension = ()

        #if custom is selected, assign the proper overhangs
		if self.GGtype == "Custom":
			part_oh = (fiveprime[:4].upper(), threeprime[-4:].upper())
			part_extension = (fiveprime[4:], threeprime[:-4])
			#Custom parts will go into the part entry vector unless skipPartPlasmid is selected
			entry_oh = entry_ohs["3"]
			entry_extension = entry_extensions["3"]
			#could possibly change this default if I want to allow custom mV parts to be made
			self.enzyme = "BsmBI"
		else:
			part_oh = part_ohs[self.GGtype]
			part_extension = part_extensions[self.GGtype]
			entry_oh = entry_ohs[self.GGtype]
			entry_extension = entry_extensions[self.GGtype]
			if self.GGtype in ["1","234","2","2a","2b","3","3a","3ab","3ac","3ad","3ae","3b","3bc","3bd","3be","3c","3cd","3ce","3d","3de","3e","4","4a","4b","5","6","7"]:
				self.enzyme = "BsmBI"
			elif self.GGtype in []:
				self.enzyme = "BbsI"

			fiveprimeOH = entry_oh[0].upper()
			threeprimeOH = entry_oh[1].upper()
			fiveprimeExt = entry_extension[0] + part_oh[0] + part_extension[0]
			threeprimeExt = part_extension[1] + part_oh[1] + entry_extension[1]

			#Set vector
			self.subclone = False
			self.vectorName = "OPL4772"

            newFrag = GGfrag(fiveprimeOH, fiveprimeExt, frag_seq, threeprimeExt, threeprimeOH)
			self.GGfrags = [newFrag]

	def checkInput_beforeOptimization(self):
		seq = Seq(self.inputSeq.upper(), IUPAC.ambiguous_dna)
		if not (_verify_alphabet(seq)):
			self.errors.append("The DNA sequence you entered in invalid.")
		if self.GGtype in ["1", "5"] and "BsmBI" in self.rsToRemove:
			self.warnings.append("You chose to remove BsmBI sites from a type connector part. Doing so will break multigene assembly.")
		if len(self.errors) != 0:
			return False
		else:
			return True

	def checkInput_afterOptimization(self):
		f_seq = FormattedSeq(Seq(self.partSeq), True)

		#Check to make sure BbsI/BsmBI sites aren't present
		if BbsI.search(f_seq) and BsmBI.search(f_seq):
			#comment the line below to let through parts with BbsI and BsmBI
			#be careful, though, these assemblies may be problematic
			#self.errors.append("Your part contains both BbsI and BsmBI sites and cannot be assembled using golden gate.")
			pass
		elif self.GGtype in ["1","234","2","2a","2b","3","3a","3ab","3ac","3ad","3ae","3b","3bc","3bd","3be","3c","3cd","3ce","3d","3de","3e","4","4a","4b","5","6","7"] and BbsI.search(f_seq):
			self.errors.append("Your part contains a BbsI site which must be removed prior to assembly.")
		elif self.GGtype in ["234","2","2a","2b","3","3a","3ab","3ac","3ad","3ae","3b","3bc","3bd","3be","3c","3cd","3ce","3d","3de","3e","4","4a","4b","6","7"] and BsmBI.search(f_seq):
			self.errors.append("Your part contains a BsmBI site which must be removed prior to assembly.")

		#Check to make sure connector parts have BsmBI sites
		bsmBIFor = self.partSeq.upper().count(BsmBI.site)
		bsmBIRev = self.partSeq.upper().count(str(Seq(BsmBI.site).reverse_complement()))
		if str(self.GGtype) == "1":
			if bsmBIFor + bsmBIRev != 1:
				self.warnings.append("Your type 1 part should have exactly 1 BsmBI site. Consider modifying for multigene assembly.")
		if str(self.GGtype) == "5":
			if bsmBIFor + bsmBIRev != 1:
				self.warnings.append("Your type 5 part should have exactly 1 BsmBI site. Consider modifying for multigene assembly.")

		#Warn if ORF has a start codon
		if self.GGtype in ["3","3a","3ab","3ac","3ad","3ae"]:
			if str(self.partSeq.upper())[:3] != "ATG":
				self.warnings.append("Your part is missing a start codon.")
		#Warn if ORF has a stop codon
		if str(self.GGtype) in ["3","3a","3ab","3ac","3ad","3ae","3b","3bc","3bd","3be","3c","3cd","3ce","3d","3de","3e","4a"]:
			if len(self.partSeq) % 3 != 0:
				self.warnings.append("Your part appears to be out of frame (length is not a multiple of 3). If this is a coding sequence, check to make sure it is correct.")
			if Seq(self.partSeq).translate().find("*") > -1:
				self.warnings.append("Your part has a stop codon. If it is not removed, the part cannot be used for making N-terminal fusions.")

		if len(self.errors) != 0:
			return False
		else:
			return True

def _verify_alphabet(sequence):
	letters = sequence.alphabet.letters
	if not letters:
		raise ValueError("Alphabet does not define letters.")
	for letter in sequence:
		if letter not in letters:
			return False
	return True

def main():
	seq = "atggtgagcgagctgattaaggagaacatgcacatgaagctgtacatggagggcaccgtgaacaaccaccacttcaagtgcacatccgagggcgaaggcaagccctacgagggcacccagaccatgagaatcaaggcggtcgagggcggccctctccccttcgccttcgacatcctggctaccagcttcatgtacggcagcaaaaccttcatcaaccacacccagggcatccccgacttctttaagcagtccttccccgagggcttcacatgggagagagtcaccacatacgaagacgggggcgtgctgaccgctacccaggacaccagcctccaggacggctgcctcatctacaacgtcaagatcagaggggtgaacttcccatccaacggccctgtgatgcagaagaaaacactcggctgggaggcctccaccgagaccctgtaccccgctgacggcggcctggaaggcagagccgacatggccctgaagctcgtgggcgggggccacctgatctgcaacttgaagaccacatacagatccaagaaacccgctaagaacctcaagatgcccggcgtctactatgtggacagaagactggaaagaatcaaggaggccgacaaagagacctacgtcgagcagcacgaggtggctgtggccagatactgcgacctccctagcaaactggggcacaga"
	possibleTemplates = {"pGG001":"CGTCTCTGACCAGaccaataaaaaacgcccggcggcaaccgagcgttctgaacaaatccagatggagttctgaggtcattactggatctatcaacaggagtccaagcgagctcgatatcaaattacgccccgccctgccactcatcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatggtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgaaacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagcgatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccgtctttcattgccatacgaaattccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccattttagcttccttagctcctgaaaatctcgataactcaaaaaatacgcccggtagtgatcttatttcattatggtgaaagttggaacctcttacgtgcccgatcaatcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtAGTCGGCGAGACGtccctatcagtgatagagattgacatccctatcagtgatagagatactgagcacggatctgaaagaggagaaaggatctatggcgagtagcgaagacgttatcaaagagttcatgcgtttcaaagttcgtatggaaggttccgttaacggtcacgagttcgaaatcgaaggtgaaggtgaaggtcgtccgtacgaaggtacccagaccgctaaactgaaagttaccaaaggtggtccgctgccgttcgcttgggacatcctgtccccgcagttccagtacggttccaaagcttacgttaaacacccggctgacatcccggactacctgaaactgtccttcccggaaggtttcaaatgggaacgtgttatgaacttcgaagacggtggtgttgttaccgttacccaggactcctccctgcaagacggtgagttcatctacaaagttaaactgcgtggtaccaacttcccgtccgacggtccggttatgcagaaaaaaaccatgggttgggaagcttccaccgaacgtatgtacccggaagacggtgctctgaaaggtgaaatcaaaatgcgtctgaaactgaaagacggtggtcactacgacgctgaagttaaaaccacctacatggctaaaaaaccggttcagctgccgggtgcttacaaaaccgacatcaaactggacatcacctcccacaacgaagactacaccatcgttgaacagtacgaacgtgctgaaggtcgtcactccaccggtgcttaataaggatctccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata",
						"pGG002":"CGTCTCTCTAAtccctatcagtgatagagattgacatccctatcagtgatagagatactgagcacggatctgaaagaggagaaaggatctatggcgagtagcgaagacgttatcaaagagttcatgcgtttcaaagttcgtatggaaggttccgttaacggtcacgagttcgaaatcgaaggtgaaggtgaaggtcgtccgtacgaaggtacccagaccgctaaactgaaagttaccaaaggtggtccgctgccgttcgcttgggacatcctgtccccgcagttccagtacggttccaaagcttacgttaaacacccggctgacatcccggactacctgaaactgtccttcccggaaggtttcaaatgggaacgtgttatgaacttcgaagacggtggtgttgttaccgttacccaggactcctccctgcaagacggtgagttcatctacaaagttaaactgcgtggtaccaacttcccgtccgacggtccggttatgcagaaaaaaaccatgggttgggaagcttccaccgaacgtatgtacccggaagacggtgctctgaaaggtgaaatcaaaatgcgtctgaaactgaaagacggtggtcactacgacgctgaagttaaaaccacctacatggctaaaaaaccggttcagctgccgggtgcttacaaaaccgacatcaaactggacatcacctcccacaacgaagactacaccatcgttgaacagtacgaacgtgctgaaggtcgtcactccaccggtgcttaataaggatctccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttataGCATCGAGACGCTacggttatccacagaatcaggggataacgcaggaaagaacatgtgagcaaaaggccagcaaaaggccaggaaccgtaaaaaggccgcgttgctggcgtttttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaagaacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagaagatcctttgatcttttctacggggtctgacgctcagtggaacgaaaactcacgttaagggattttggtcatgattgatcgggcacgtaagaggttccaactttcaccataatgaaataagatcactaccgggcgtattttttgagttatcgagattttcaggagctaaggaagctaaaatggagaaaaaaatcactggatataccaccgttgatatatcccaatggcatcgtaaagaacattttgaggcatttcagtcagttgctcaatgtacctataaccagaccgttcagctggatattacggcctttttaaagaccgtaaagaaaaataagcacaagttttatccggcctttattcacattcttgcccgcctgatgaatgctcatccggaatttcgtatggcaatgaaagacggtgagctggtgatatgggatagtgttcacccttgttacaccgttttccatgagcaaactgaaacgttttcatcgctctggagtgaataccacgacgatttccggcagtttctacacatatattcgcaagatgtggcgtgttacggtgaaaacctggcctatttccctaaagggtttattgagaatatgtttttcgtttcagccaatccctgggtgagtttcaccagttttgatttaaacgtggccaatatggacaacttcttcgcccccgttttcaccatgggcaaatattatacgcaaggcgacaaggtgctgatgccgctggcgattcaggttcatcatgccgtttgtgatggcttccatgtcggcagaatgcttaatgaattacaacagtactgcgatgagtggcagggcggggcgtaatttgatatcgagctcgcttggactcctgttgatagatccagtaatgacctcagaactccatctggatttgttcagaacgctcggttgccgccgggcgttttttattggtAG"}

	part = GGpart("blah", "3", seq, inputSeqType="DNA", rsToRemove=[], plasmidName="pWCD0510",
	 				skipPartPlasmid=False, method="PCR", fiveprime="", threeprime="",
					possibleTemplates=possibleTemplates)

	print part
	print part.getPrimers()

# 	if len(part.errors) > 0:
# 		print "Errors:"
# 		for error in part.errors: print error
# 	if len(part.warnings) > 0:
# 		print "Warnings:"
# 		for warning in part.warnings: print warning
# 	print "\nPlasmid Sequence:\n" + part.plasmidSeq

	#part.writeXMLfile("test.xml")
	#print getXML(part)


if __name__ == '__main__':
	main()
