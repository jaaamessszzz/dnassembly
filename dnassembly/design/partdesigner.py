import os
import sys
import subprocess

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import math
from Bio.Restriction import Restriction
from itertools import chain

from ..dna.part import *
from ..utils import pairwise
"""
GGfrag >> Part class instantiation/function/attribute mapping

GGfrag                  >>      Part

self                    >>      input_part
self.getLength()        >>      len(sequence)
GGfrag()                >>      Part.GGfrag()

"""

def divideBySize(input_part, size):
    """Return a list of Parts divided by size"""
    length = len(input_part)
    numFrags = 0
    if length % size == 0:
        numFrags = int(length/size)
    else:
        numFrags = int(length/size) + 1
    if numFrags <= 1:
        return [input_part]
    else:
        subseqs = _split_seq_byLength(input_part.sequence, numFrags)
        newFrags = []
        for i in range(len(subseqs)):
            newFrag = ""  # ?????
            if i == 0:
                lenFiveSeq = len(input_part.fiveprimeOH) + len(input_part.fiveprimeExt)
                newFrag = Part.GGfrag(fiveprimeOH=input_part.fiveprimeOH, fiveprimeExt=input_part.fiveprimeExt, seq=subseqs[i][lenFiveSeq:], forced_method=input_part.forced_method)
            elif i == len(subseqs) - 1:
                lenThreeSeq = len(input_part.threeprimeExt) + len(input_part.threeprimeOH)
                if lenThreeSeq == 0:
                    newFrag = Part.GGfrag(seq=subseqs[i], forced_method=input_part.forced_method)
                else:
                    newFrag = Part.GGfrag(seq=subseqs[i][:-lenThreeSeq], threeprimeExt=input_part.threeprimeExt, threeprimeOH=input_part.threeprimeOH, forced_method=input_part.forced_method)
            else:
                newFrag = Part.GGfrag(seq=subseqs[i], forced_method=input_part.forced_method)
            newFrags.append(newFrag)
        return newFrags


def divideByIndexTuples(input_part, index_tuples):
    """Split sequence by a list of index tuple, where each tuple represents a slice where a sequence modification has
    been made by removeRS()
    """
    subseqs = _split_seq_byIndicies(input_part.seq, index_tuples)
    subseqs = _merge_single_bases(subseqs)
    if len(subseqs) <= 1:
        input_part.seq = input_part.seq.lower()
        return [input_part]
    else:
        newFrags = []
        for i in range(len(subseqs)):
            newFrag = ""  # ?????
            if i == 0:
                newFrag = Part.GGfrag(fiveprimeOH=input_part.fiveprimeOH, fiveprimeExt=input_part.fiveprimeExt, seq=subseqs[i])
            elif i == len(subseqs) - 1:
                newFrag = Part.GGfrag(seq=subseqs[i], threeprimeExt=input_part.threeprimeExt, threeprimeOH=input_part.threeprimeOH)
            else:
                newFrag = Part.GGfrag(seq=subseqs[i])
            newFrags.append(newFrag)
        return newFrags


def getPCRprimers(input_part, maxPrimerLength, annealingLength, oligoTM):
    # todo: figure out how to deal with tails... these should already be part of a Part sequence
    tails = input_part.tails
    maxAnnealLengthF = maxPrimerLength - len(tails[0] + input_part.fiveprimeOH + input_part.fiveprimeExt)
    maxAnnealLengthR = maxPrimerLength - len(input_part.threeprimeExt + input_part.threeprimeOH + tails[1])
    i = 1
    while get_tm(input_part.seq[:i]) < oligoTM and i < maxAnnealLengthF:
        i += 1
    #fPrimer = tails[0] + input_part.fiveprimeOH + input_part.fiveprimeExt + input_part.seq[:annealingLength]
    fPrimer = tails[0] + input_part.fiveprimeOH + input_part.fiveprimeExt + input_part.seq[:i]

    i = 1
    while get_tm(input_part.seq[-i:]) < oligoTM and i < maxAnnealLengthR:
        i += 1
    #endSeq = input_part.seq[-annealingLength:] + input_part.threeprimeExt + input_part.threeprimeOH + tails[1]
    endSeq = input_part.seq[-i:] + input_part.threeprimeExt + input_part.threeprimeOH + tails[1]
    rPrimer = str(Seq(endSeq).reverse_complement())
    return (fPrimer, rPrimer)


def getOligoAssemPrimers(input_part, maxPrimerLength):

    topStrand = input_part.fiveprimeOH + input_part.fiveprimeExt + input_part.seq + input_part.threeprimeExt
    bottomStrand = str(Seq(input_part.fiveprimeExt + input_part.seq + input_part.threeprimeExt + input_part.threeprimeOH).reverse_complement())
    oligos = []
    #If only 1 oligo is required
    if len(topStrand) <= maxPrimerLength and len(bottomStrand) <= maxPrimerLength:
        oligos.append(topStrand)
        oligos.append(bottomStrand)
    else:
        numOligos = 2
        found = False
        f_divided = []
        r_divided = []
        while not found:
            f_divided = _split_seq_byLength(topStrand[4:], numOligos * 2 - 1)
            r_divided = _split_seq_byLength(bottomStrand[4:], numOligos * 2 - 1)
            if len(f_divided[1]) + len(f_divided[2]) <= maxPrimerLength:
                if len(r_divided[1]) + len(r_divided[2]) <= maxPrimerLength:
                    found = True
            numOligos += 1

        f_oligos = []
        r_oligos = []
        for i in range(0, len(f_divided), 2):
            f_oligo = ""
            r_oligo = ""
            if i == 0:
                f_oligo = topStrand[:4] + f_divided[i]
                r_oligo = bottomStrand[:4] + r_divided[i]
            else:
                f_oligo = f_divided[i-1] + f_divided[i]
                r_oligo = r_divided[i-1] + r_divided[i]
            f_oligos.append(f_oligo)
            r_oligos.append(r_oligo)
            oligos = tuple(f_oligos + r_oligos)
    return oligos

def get_tm(seq):
    salt = .05
    conc = .0000001
    numG = 0.0 + seq.count("G") + seq.count("g")
    numC = 0.0 + seq.count("C") + seq.count("c")
    base = 81.5 + 16.6*math.log10(salt) + 41*((numG + numC)/len(seq))
    return base - (600/len(seq))

#Merges fragments that are 1 base in length
def _merge_single_bases(subseqs):
    new = []
    curr = ""
    for each in subseqs:
        if len(each) == 1:
            curr += each
        else:
            if curr != "":
                new.append(curr.lower())
            new.append(each)
            curr = ""
    if curr != "":
        new.append(curr.lower())
    return new

#Splits a string into an array of equally divided substrings
def _split_seq_byLength(seq, p):
    newseq = []
    n = len(seq) / p    # min items per subsequence
    r = len(seq) % p    # remaindered items
    b,e = 0, n + min(1, r)  # first split
    for i in range(p):
        newseq.append(seq[b:e])
        r = max(0, r-1)  # use up remainders
        b,e = e, e + n + min(1, r)  # min(1,r) is always 0 or 1
    return newseq

#Splits a string into an array of substrings that are divided by case changes
def _split_seq_byCase(seq):
    pieces = []
    leftIndex = 0
    rightIndex = len(seq)
    for index in range (len(seq)-1):
        if seq[index].isupper() != seq[index+1].isupper():
            rightIndex = index + 1
            pieces.append(seq[leftIndex:rightIndex])
            leftIndex = index + 1
    rightIndex = len(seq)
    pieces.append(seq[leftIndex:rightIndex])
    return pieces


def _split_seq_byIndicies(seq, sequence_tuples):
    """Divides a sequence using provided indicies in tuples"""
    index_list = sorted(list(set([idx for idx in chain(*sequence_tuples)] + [0, len(seq)])))
    new_fragments = list()

    for left_cut, right_cut in pairwise(index_list):
        print(left_cut, right_cut)
        new_fragments.append(seq[left_cut:right_cut])
    return new_fragments
