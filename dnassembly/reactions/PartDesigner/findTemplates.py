#!/usr/bin/env python
# encoding: utf-8

from Bio.Blast.Applications import NcbiblastnCommandline
import sys, os, warnings
from Bio import SearchIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# BASE_PATH is the absolute path of ../../.. relative to this script location
#BASE_PATH = reduce (lambda l,r: l + os.path.sep + r, os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep )[:-3] )
#sys.path.append( os.path.join( BASE_PATH, "main" ) )
#sys.path.append('../../../main')
#from siteinfo import *

#os.environ['PATH'] += os.pathsep + BLAST_PATH


def runBLAST(infile, outfile, blastdb):
	#Run BLAST and Save Output File
	blastx_cline = NcbiblastnCommandline(query=infile, db=blastdb, outfmt=5, penalty=-9000000, ungapped=True, out=outfile)
	blastx_cline()


def reindexCircular(plasmidSeq, blastdb):
	matches = findTemplates(plasmidSeq+plasmidSeq, blastdb, reindexCircular=True)
	if len(matches) == 0:
		return plasmidSeq
	for match in matches:
		if match.query_range[0] < len(plasmidSeq) < match.query_range[1]:
			reindexed = (plasmidSeq+plasmidSeq)[match.query_range[0]:match.query_range[0]+len(plasmidSeq)]
			return reindexed
	reindexed = (plasmidSeq+plasmidSeq)[matches[0].query_range[0]:matches[0].query_range[0]+len(plasmidSeq)]
	return reindexed

def findTemplates(plasmidSeq, blastdb, reindexCircular=False):
	matches = []
	if blastdb != "":
		redoBLAST = True
		currSeq = plasmidSeq
		while redoBLAST:
			redoBLAST = False
			#Run Blast and save output
			infilename = BASE_PATH + "/plasmids/assets/GGDesigner/blast_output/infile.fa"
			output_handle = open(infilename, "w")
			SeqIO.write(SeqRecord(Seq(currSeq)), output_handle, 'fasta')
			output_handle.close()
			outfilename = BASE_PATH + "/plasmids/assets/GGDesigner/blast_output/out.xml"
			runBLAST(infilename, outfilename, blastdb)

			#Parse Blast output file
			results = SearchIO.read(outfilename, 'blast-xml')
			sort_key = lambda hit: [int(hit.id.split("|")[0]), -int(hit.id.split("|")[1])]
			results = results.sort(key=sort_key, reverse=False, in_place=False)

			#Search through results for all unique sequence matches that aren't completely overlapping
			temp_matches = []
			for hit in results:
				for newMatch in hit:
					newMatch.query_range
					alreadyCovered = False

					#Set alreadyCovered to True if the match is inside another existing match
					for oldMatch in temp_matches:
						if (oldMatch.query_range[0] <= newMatch.query_range[0] <= newMatch.query_range[1] <=  oldMatch.query_range[1]):
							alreadyCovered = True
						if alreadyCovered:
							break
					#If the match covers new sequence, add it and remove any sequences that are completely covered by it
					if not alreadyCovered:
						temp_matches = [oldMatch for oldMatch in temp_matches if not (newMatch.query_range[0] <= oldMatch.query_range[0] <= oldMatch.query_range[1] <= newMatch.query_range[1])]
						temp_matches.append(newMatch)
			#Sort by the starting index of the match
			sort_key = lambda match: int(match.query_range[0])
			temp_matches = sorted(temp_matches, key=sort_key)

			#Remove matches that are covered by the previous and next fragments
			allowable_gap_length = 15 #prev and next match can be up to this far apart from each other and curr will still be removed
			modified = False
			while not modified:
				modified = True
				for i in range(1, len(temp_matches)-1):
					curr = temp_matches[i]
					prev = temp_matches[i-1]
					next = temp_matches[i+1]
					if prev.query_range[1] >= next.query_range[0] - allowable_gap_length:
						temp_matches = temp_matches[:i] + temp_matches[i+1:]
						modified = True
						break
			matches += temp_matches
			sort_key = lambda match: int(match.query_range[0])
			matches = sorted(matches, key=sort_key)


			if len(temp_matches) > 0:
				redoBLAST = True
				for match in temp_matches:
					newSeq = currSeq[:match.query_range[0]]
					for i in range(match.query_range[1] - match.query_range[0]):
						newSeq += "N"
					newSeq += currSeq[match.query_range[1]:]
					currSeq = newSeq

	else:
		pass

	# return matches to reindexCircular()
	if reindexCircular == True:
		return matches


	#Pick the fragment templates to use and adjust capitalization of plasmidSeq to
	#reflect the fragments that were chosen.
	#Note: This algorithm need to be improved to allow for small gaps between parts
	#If you don't allow this it could increase the number of fragments in the assembly
	templates = {}
	index = 0
	fragIndex = 0
	while index < len(plasmidSeq):
		if fragIndex > len(matches)-1:
			rindex = index + 1
		elif index >= matches[fragIndex].query_range[0]:
			rindex = matches[fragIndex].query_range[1]
			hitID = matches[fragIndex].hit.id.split("|")[2]
			hitSeq = str(matches[fragIndex].hit.seq)
			if hitID in ["S288C", "DH10B", "GS115"]:
				hitID = hitID + "_Genomic_DNA"
			if hitID not in templates.keys():
				templates[hitID] = hitSeq
			elif templates[hitID].find(hitSeq) == -1:
				templates[hitID] += "..." + hitSeq
			fragIndex += 1
		else:
			rindex = index + 1
		lindex = index
		leftBase = str(plasmidSeq[index-1:index])
		if leftBase.isupper():
			plasmidSeq = plasmidSeq[:lindex] + plasmidSeq[lindex:rindex].lower() + plasmidSeq[rindex:]
		else:
			plasmidSeq = plasmidSeq[:lindex] + plasmidSeq[lindex:rindex].upper() + plasmidSeq[rindex:]
		index = rindex
	#plasmidSeq has been appropriately capitalized and templates is a list of relevant database ids
	return (plasmidSeq, templates)

def main(args):
	plasmidSeq = "atgtctaaaggtgaagaattattcactggtgttgtcccaattttggttgaattagatggtgatgttaatggtcacaaattttctgtctccggtgaaggtgaaggtgatgctacttacggtaaattgaccttaaaattgatttgtactactggtaaattgccagttccatggccaaccttagtcactactttaggttatggtttgcaatgttttgctagatacccagatcatatgaaacaacatgactttttcaagtctgccatgccagaaggttatgttcaagaaagaactatttttttcaaagatgacggtaactacaagaccagagctgaagtcaagtttgaaggtgataccttagttaatagaatcgaattaaaaggtattgattttaaagaaggtggtaacattttaggtcacaaattggaatacaactataactctcacaatgtttacatcactgctgacaaacaaaagaatggtatcaaagctaacttcaaaattagacacaacattgaagatggtggtgttcaattagctgaccattatcaacaaaatactccaattggtgatggtccagtcttgttaccagacaaccattacttatcctatcaatctgccttatccaaagatccaaacgaaaagagagatcacatggtcttgttagaatttgttactgctgctggtattacccatggtatggatgaattgtacaaatctaaaggtgaagaattattcactggtgttgtcccaattttggttgaattagatggtgatgttaatggtcacaaattttctgtctccggtgaaggtgaaggtgatgctacttacggtaaattgaccttaaaattTATTTGTACTACTGGTAAATTGCCAGTTCCATGGCCAACCTTAGTCACTACTTTAacttggggtgttcaatgtttttcaagatacccagatcatatgaaacaacatgactttttcaagtctgccatgccagaaggttatgttcaagaaagaactatttttttcaaagatgacggtaactacaagaccagagctgaagtcaagtttgaaggtgataccttagttaatagaatcgaattaaaaggtattgattttaaagaagatggtaacattttaggtcacaaattggaatacatttataactctcacaatgtttacatcactgctgacaaacaaaagaatggtatcaaagctaacttcaaaattagacacaacattgaagatggttctgttcaattagctgaccattatcaacaaaatactccaattggtgatggtccagtcttgttaccagacaaccattacttatccactcaatctaggttatccaaagatccaaacgaaaagagggaccacatggtcttgttagaatttgttactgctgctggtattacccatggtatggatgaattgtacaaaatggcttcctccgaagacgttatcaaagagttcatgcgtttcaaagttcgtatggaaggttccgttaacggtcacgagttcgaaatcgaaggtgaaggtgaaggtcgtccgtacgaaggtacccagaccgctaaactgaaagttaccaaaggtggtccgctgccgttcgcttgggacatcctgtccccgcagttccagtacggttccaaagcttacgttaaacacccggctgacatcccggactacctgaaactgtccttcccggaaggtttcaaatgggaacgtgttatgaacttcgaagacggtggtgttgttaccgttacccaggactcctccctgcaagacggtgagttcatctacaaagttaaactgcgtggtaccaacttcccgtccgacggtccggttatgcagaaaaaaaccatgggttgggaagcttccaccgaacgtatgtacccggaagacggtgctctgaaaggtgaaatcaaaatgcgtctgaaactgaaagacggtggtcactacgacgctgaagttaaaaccacctacatggctaaaaaaccggttcagctgccgggtgcttacaaaaccgacatcaaactggacatcacctcccacaacgaagactacaccatcgttgaacagtacgaacgtgctgaaggtcgtcactccaccggtgct"
	blastdb = "blast_dbs/all_users.fa"
	output = findTemplates(plasmidSeq, blastdb)
	print(output[0])
	print(output[1])


if __name__ == "__main__":
	main(sys.argv[1:])