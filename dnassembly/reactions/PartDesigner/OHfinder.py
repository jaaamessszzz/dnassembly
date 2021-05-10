#!/usr/bin/env python
# encoding: utf-8
"""
OHfinder.py

Created by Will DeLoache on 2012-09-04.
"""

from .GG_OHcheck import *
from .GGfrag import *
from Bio.Restriction import BsmBI, BsaI, Restriction, FormattedSeq
from Bio.Seq import Seq

from ...dna.part import Part

def reindexGGfrags(GGfrags, bestOHs):
	indices = []
	for each in bestOHs:
		indices.append(each[1])

	fragments = []
	for each in GGfrags:
		newFrag = [each.fiveprimeOH, each.fiveprimeExt, each.seq, each.threeprimeExt, each.threeprimeOH]
		fragments.append(newFrag)

	#shift fragments based on indicies
	shiftedFrags = []
	for i in range(len(fragments)):
		Lindex = indices[i]
		Rindex = indices[i+1]
		newFrag = ["","","","",""]

		#Set newFrag[0] and newFrag[1]
		if i == 0:
			newFrag[0] = fragments[i][0]
			newFrag[1] = fragments[i][1]
		else:
			if Lindex < 0:
				s = fragments[i-1][2] + fragments[i-1][3]
				s = s[Lindex:]
				newFrag[0] = (s + fragments[i][1] + fragments[i][2])[:4]
				#Think about this
				newFrag[1] = (s + fragments[i][1])[4:]
			else:
				newFrag[0] = (fragments[i][1] + fragments[i][2])[Lindex:Lindex+4]
				#Think about this
				newFrag[1] = fragments[i][1][Lindex+4:]

		#Set newFrag[2]
		newFrag[2] = fragments[i][2]
		if i != 0:
			lSeqindex = 4 + Lindex - len(fragments[i][1])
			if lSeqindex > 0:
				newFrag[2] = newFrag[2][lSeqindex:]
		if i != len(fragments) - 1:
			rSeqindex = Rindex + len(fragments[i][3])
			if rSeqindex < 0:
				newFrag[2] = newFrag[2][:rSeqindex]

		#Set newFrag[3] and newFrag[4]
		if i == len(fragments) - 1:
			newFrag[3] = fragments[i][3]
			newFrag[4] = fragments[i][4]
		else:
			if Rindex + 4 >= 0:
				s = fragments[i+1][1] + fragments[i+1][2]
				s = s[:Rindex + 4]
				newFrag[4] = (fragments[i][2] + fragments[i][3] + s)[-4:]
				#Think about this
				newFrag[3] = (fragments[i][3] + s)[:-4]
			else:
				newFrag[4] = (fragments[i][2] + fragments[i][3])[Rindex:Rindex+4]
				newFrag[3] = fragments[i][3][:Rindex]
		newGGfrag = Part.GGfrag(newFrag[0], newFrag[1], newFrag[2], newFrag[3], newFrag[4])
		shiftedFrags.append(newGGfrag)
	return shiftedFrags

def findBestOverhangs(existingOHs, possOHs, additionalOHs, constraints):
	if len(possOHs) == 0:
		#print "Found overhangs using constraint stringency =", constraints
		return existingOHs
	if len(possOHs[0]) == 0:
		return False
	tempOHs = []
	for oh in existingOHs:
		tempOHs.append(oh[0])
	for oh in additionalOHs:
		tempOHs.append(oh)
	for newOH in possOHs[0]:
		if overhangOK(newOH[0], tempOHs, constraints):
			existingOHs_copy = list(existingOHs)
			existingOHs_copy.append(newOH)
			#if this was the last overhang to check
			if len(possOHs) == 1:
				#print "Found overhangs using constraint stringency =", constraints
				return existingOHs_copy
			return findBestOverhangs(existingOHs_copy, list(possOHs[1:]), list(additionalOHs), constraints)
	return False

#Finds all possible overhangs that can be used for each junction
#where the length of all primers <= maxPrimerlength
#Returns a list of sorted lists [[junction 1 ohs], [junction 2 ohs]...]
#where each ohs is represented as a tuple (oh sequence, index)
#the index denotes the location of the 1st bp in the oh relative to the junction of the fragments
#tuples are sorted by the absolute value of the index (closest to the junction first)
def findPossOH_byPrimerLength(GGfrags, maxPrimerLength, annealingLength, gBlockMaxSize, enzyme):
	segments = []
	forced_methods = []
	for each in GGfrags:
		seg = [each.fiveprimeOH + each.fiveprimeExt, each.seq, each.threeprimeExt + each.threeprimeOH]
		segments.append(seg)
		forced_methods.append(each.forced_method)
	wiggleRoom = []
	for i in range(len(segments)):
		leftWiggle = 0
		rightWiggle = 0
		#If junction is with vector, skip it because it is fixed
		if i == 0:
			pass
		else:
			leftSeg = segments[i-1]
			rightSeg = segments[i]
			#constraints on leftWiggle from primer length
			leftWiggle_primer = maxPrimerLength - len(rightSeg[0]) - annealingLength
			#constraints on leftWiggle from piece length
			leftWiggle_pieceLen = len(leftSeg[1] + leftSeg[2]) - wiggleRoom[i-1][1] - annealingLength
			if forced_methods[i] == "gBlocks":
				leftWiggle_gBlock = gBlockMaxSize - len(rightSeg[0] + rightSeg[1] + rightSeg[2]) - 22
				leftWiggle_pieceLen	= min(leftWiggle_pieceLen, leftWiggle_gBlock)
			elif forced_methods[i] == "Oligo Assembly":
				leftWiggle_oligo = 200 - len(rightSeg[0] + rightSeg[1] + rightSeg[2])
				leftWiggle_pieceLen	= min(leftWiggle_pieceLen, leftWiggle_oligo)
			#assign the minimum constraints to leftWiggle
			leftWiggle = min(leftWiggle_primer, leftWiggle_pieceLen)
			#constraints on rightWiggle from primer length
			rightWiggle_primer = maxPrimerLength - len(leftSeg[2]) - annealingLength - 4
			#constraints on rightWiggle from piece length
			#Don't have to substract 4 here, but it will make primer design easier
			rightWiggle_pieceLen = len(rightSeg[0] + rightSeg[1]) - annealingLength - 4
			#assign the minimum constraints to rightWiggle
			if forced_methods[i-1] == "gBlocks":
				rightWiggle_gBlock = gBlockMaxSize - len(leftSeg[0] + leftSeg[1] + leftSeg[2]) - 22 - 4 + wiggleRoom[i-1][0]
				rightWiggle_pieceLen = min(rightWiggle_pieceLen, rightWiggle_gBlock)
			elif forced_methods[i-1] == "Oligo Assembly":
				rightWiggle_oligo = 200 - len(leftSeg[0] + leftSeg[1] + leftSeg[2]) - 4 + wiggleRoom[i-1][0]
				rightWiggle_pieceLen = min(rightWiggle_pieceLen, rightWiggle_oligo)
			rightWiggle = min(rightWiggle_primer, rightWiggle_pieceLen)
		wiggleRoom.append((-leftWiggle, rightWiggle))



	poss_ohs = []
	for i in range(1,len(segments)):
		leftSeg = segments[i-1]
		rightSeg = segments[i]
		leftPiece = leftSeg[0] + leftSeg[1] + leftSeg[2]
		rightPiece = rightSeg[0] + rightSeg[1] + rightSeg[2]
		combined = leftPiece.upper() + rightPiece.upper()
		oh_possibilities = []
		poss_string = combined[len(leftPiece)+wiggleRoom[i][0]:len(leftPiece)+wiggleRoom[i][1]+1]
		enz = getattr(Restriction, enzyme)
		if len(enz.search(FormattedSeq(Seq(poss_string)))) > 0:
			oh_index = None
			oh_seq = None
			if poss_string.find(enz.site) > -1:
				site_location = poss_string.find(enz.site)
				oh_index = site_location + wiggleRoom[i][0] + 7
				oh_seq = poss_string[site_location + 7 : site_location + 11]
			else:
				site_location = poss_string.find(str(Seq(enz.site).reverse_complement()))
				oh_index = site_location + wiggleRoom[i][0] - 5
				oh_seq = poss_string[site_location - 5 : site_location - 1]
			poss_ohs.append([(oh_seq, oh_index)])
		else:
			for j in range(len(leftPiece)+wiggleRoom[i][0], len(leftPiece)+wiggleRoom[i][1]+1):
				oh_possibilities.append((combined[j:j+4],j-len(leftPiece)))
			#If no overhang options were found for a junction, return false
			if len(oh_possibilities) == 0:
				return False
			else:
				oh_sorted = sorted(oh_possibilities, key=lambda overhang: abs(overhang[1]))
				poss_ohs.append(oh_sorted)
	return poss_ohs


#Finds all possible overhangs that can be used for each junction
#where the length of all GGfrags <= maxFraglength
#Returns a list of sorted lists [[junction 1 ohs], [junction 2 ohs]...]
#where each ohs is represented as a tuple (oh sequence, index)
#the index denotes the location of the 1st bp in the oh relative to the junction of the fragments
#tuples are sorted by the absolute value of the index (closest to the junction first)
def findPossOH_byFragLength(GGfrags, maxFragLength, annealingLength):
	segments = []
	for each in GGfrags:
		seg = [each.fiveprimeOH + each.fiveprimeExt, each.seq, each.threeprimeExt + each.threeprimeOH]
		segments.append(seg)
	wiggleRoom = []
	for i in range(len(segments)):
		leftWiggle = 0
		rightWiggle = 0
		#If junction is with vector, skip it because it is fixed
		if i == 0:
			pass
		else:
			leftSeg = segments[i-1]
			rightSeg = segments[i]
			#constraints on leftWiggle from rightSeg length (<500bps)
			leftWiggle_rightSeg = maxFragLength - len(rightSeg[0] + rightSeg[1] + rightSeg[2])
			#constraints on leftWiggle from leftSeg length (>18bps)
			leftWiggle_leftSeg = len(leftSeg[0] + leftSeg[1] + leftSeg[2]) - wiggleRoom[i-1][1] - annealingLength
			#assign the minimum constraints to leftWiggle
			leftWiggle = min(leftWiggle_leftSeg, leftWiggle_rightSeg)
			#constraints on rightWiggle from leftSeg length (<500bps)
			rightWiggle_leftSeg = maxFragLength - len(leftSeg[0] + leftSeg[1] + leftSeg[2]) - 4
			#constraints on rightWiggle from rightSeg length (>18bps)
			rightWiggle_rightSeg = len(rightSeg[0] + rightSeg[1] + rightSeg[2]) - annealingLength - 4
			#assign the minimum constraints to rightWiggle
			rightWiggle = min(rightWiggle_leftSeg, rightWiggle_rightSeg)
		wiggleRoom.append((-leftWiggle, rightWiggle))
	poss_ohs = []
	for i in range(1,len(segments)):
		leftSeg = segments[i-1]
		rightSeg = segments[i]
		leftPiece = leftSeg[0] + leftSeg[1] + leftSeg[2]
		rightPiece = rightSeg[0] + rightSeg[1] + rightSeg[2]
		combined = leftPiece.upper() + rightPiece.upper()
		oh_possibilities = []
		for j in range(len(leftPiece)+wiggleRoom[i][0], len(leftPiece)+wiggleRoom[i][1]+1):
			oh_possibilities.append((combined[j:j+4],j-len(leftPiece)))
		#If no overhang options were found for a junction, return false
		if len(oh_possibilities) == 0:
			return False
		else:
			poss_ohs.append(sorted(oh_possibilities, key=lambda overhang: abs(overhang[1])))
	return poss_ohs



def main():
	frag1 = GGfrag("ATCG", "gggagacg", "ATCGACTAGCATCAGCTACGACTACGATCAGCATCAG", "ac","")
	frag2 = GGfrag("", "acttacga", "CGGATCAGCTAGGCTACTTAGCTAGATCGATCGACTG", "atggtcgatg", "GGAT")
	GGfrags = [frag1, frag2]
	bestOHs = [("ATCG", 0), ("xxxx", 5), ("GGAT", 0)]
	reindexGGfrags(GGfrags, bestOHs)


if __name__ == '__main__':
	main()
