#!/usr/bin/env python
# encoding: utf-8
"""
GGfrag.py

Modified from code written by Will DeLoache.
"""
import os, sys
import subprocess
from functools import reduce

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import math
from Bio.Restriction import Restriction

BASE_PATH = reduce (lambda l,r: l + os.path.sep + r, os.path.dirname( os.path.realpath( __file__ ) ).split( os.path.sep )) + "/"

class GGfrag():
    def __init__(self, fiveprimeOH = "", fiveprimeExt = "", seq = "", threeprimeExt = "",  threeprimeOH = "", forced_method=False):
        self.seq = seq
        self.fiveprimeOH = fiveprimeOH
        self.threeprimeOH = threeprimeOH
        self.fiveprimeExt = fiveprimeExt
        self.threeprimeExt = threeprimeExt
        self.forced_method = forced_method #used to specify assembly method if it is fixed already (e.g. gBlocks)
        self.tails = []

    def __str__(self):
        seq = self.seq
        if seq == "": seq = "_"
        fiveprimeOH = self.fiveprimeOH
        if fiveprimeOH == "": fiveprimeOH = "_"
        threeprimeOH = self.threeprimeOH
        if threeprimeOH == "": threeprimeOH = "_"
        fiveprimeExt = self.fiveprimeExt
        if fiveprimeExt == "": fiveprimeExt = "_"
        threeprimeExt = self.threeprimeExt
        if threeprimeExt == "": threeprimeExt = "_"
        output = fiveprimeOH + " " + fiveprimeExt + " " + seq + " " + threeprimeExt + " " + threeprimeOH
        return output

    def getFragSeq(self):
        output = self.fiveprimeOH + self.fiveprimeExt + self.seq + self.threeprimeExt + self.threeprimeOH
        return output

    def getLength(self):
        length = len(self.fiveprimeOH) + len(self.fiveprimeExt) + len(self.seq) + len(self.threeprimeExt) + len(self.threeprimeOH)
        return length

    def assignTails(self, tails, enzyme):
        enz = getattr(Restriction, enzyme)
        for_site = enz.site
        rev_site = str(Seq(enz.site).reverse_complement())

        #Delete external restriction site on primer if it's already internal
        self.tails = [tails[0], tails[1]]
        if self.getFragSeq()[5:11].upper() == rev_site:
            self.tails[0] = tails[0][0:4]
        if self.getFragSeq()[-11:-5].upper() == for_site:
            self.tails[1] = tails[1][-4:]

    #Divide GGfrag into a list of smaller fragments
    #The length of each fragment is determined by the size variable
    #This is useful when using gBlocks for part synthesis
    def divideBySize(self, size):
        length = self.getLength()
        numFrags = 0
        if length % size == 0:
            numFrags = int(length/size)
        else:
            numFrags = int(length/size) + 1
        if numFrags <= 1:
            return [self]
        else:
            subseqs = _split_seq_byLength(self.getFragSeq(), numFrags)
            newFrags = []
            for i in range(len(subseqs)):
                newFrag = ""
                if i == 0:
                    lenFiveSeq = len(self.fiveprimeOH) + len(self.fiveprimeExt)
                    newFrag = GGfrag(self.fiveprimeOH, self.fiveprimeExt, subseqs[i][lenFiveSeq:], "", "", forced_method=self.forced_method)
                elif i == len(subseqs) - 1:
                    lenThreeSeq = len(self.threeprimeExt) + len(self.threeprimeOH)
                    if lenThreeSeq == 0:
                        newFrag = GGfrag("", "", subseqs[i], "", "", forced_method=self.forced_method)
                    else:
                        newFrag = GGfrag("", "", subseqs[i][:-lenThreeSeq], self.threeprimeExt, self.threeprimeOH, forced_method=self.forced_method)
                else:
                    newFrag = GGfrag(seq=subseqs[i], forced_method=self.forced_method)
                newFrags.append(newFrag)
            return newFrags

    def divideByCaseChange(self):
        subseqs = _split_seq_byCase(self.seq)
        subseqs = _merge_single_bases(subseqs)
        if len(subseqs) <= 1:
            self.seq = self.seq.lower()
            return [self]
        else:
            newFrags = []
            for i in range(len(subseqs)):
                newFrag = ""
                if i == 0:
                    newFrag = GGfrag(self.fiveprimeOH, self.fiveprimeExt, subseqs[i])
                elif i == len(subseqs) - 1:
                    newFrag = GGfrag("", "", subseqs[i], self.threeprimeExt, self.threeprimeOH)
                else:
                    newFrag = GGfrag(seq=subseqs[i])
                newFrags.append(newFrag)
            return newFrags

    def getPCRprimers(self, maxPrimerLength, annealingLength, oligoTM):
        tails = self.tails
        maxAnnealLengthF = maxPrimerLength - len(tails[0] + self.fiveprimeOH + self.fiveprimeExt)
        maxAnnealLengthR = maxPrimerLength - len(self.threeprimeExt + self.threeprimeOH + tails[1])
        i = 1
        while get_tm(self.seq[:i]) < oligoTM and i < maxAnnealLengthF:
            i += 1
        #fPrimer = tails[0] + self.fiveprimeOH + self.fiveprimeExt + self.seq[:annealingLength]
        fPrimer = tails[0] + self.fiveprimeOH + self.fiveprimeExt + self.seq[:i]

        i = 1
        while get_tm(self.seq[-i:]) < oligoTM and i < maxAnnealLengthR:
            i += 1
        #endSeq = self.seq[-annealingLength:] + self.threeprimeExt + self.threeprimeOH + tails[1]
        endSeq = self.seq[-i:] + self.threeprimeExt + self.threeprimeOH + tails[1]
        rPrimer = str(Seq(endSeq).reverse_complement())
        return (fPrimer, rPrimer)


    def getPCAprimers(self, maxPrimerLength):
        tails = self.tails
        path = BASE_PATH + "GeneDesign/"
        product = tails[0] + self.fiveprimeOH + self.fiveprimeExt + self.seq + self.threeprimeExt + self.threeprimeOH + tails[1]
        SeqIO.write(SeqRecord(Seq(product), id="1", description=""), path + "seq.fa", "fasta")
        pipe = subprocess.call(["./Design_Building_Blocks.pl", "--input=seq.fa", "-a=1", "-m=" + str(maxPrimerLength), "-o=" + str(maxPrimerLength-10)], cwd=path, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if pipe != 0:
            raise Exception("Design_Building_Blocks.pl returned with error code " + str(pipe))
        oligos = []
        for seq_record in SeqIO.parse(path + "seq_gdBB_1/1_gdBB_oligos.FASTA", "fasta"):
            oligos.append(str(seq_record.seq))
        return oligos


    def getOligoAssemPrimers(self, maxPrimerLength):

        topStrand = self.fiveprimeOH + self.fiveprimeExt + self.seq + self.threeprimeExt
        bottomStrand = str(Seq(self.fiveprimeExt + self.seq + self.threeprimeExt + self.threeprimeOH).reverse_complement())
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

def main():
    seq = "aaaaaaaaaatgtctaaaggtgaagaattattcactggtgttgtcccaattttggttgaattagatggtgatgttaatggtcacaaattttctgtctccggtgaaggtgaaggtgatgctacttacggtaaattgaccttaaaattgatttgtactactggtaaattgccagttccatggccaaccttagtcactactttaggttatggtttgcaatgttttgctagatacccagatcatatgaaacaacatgactttttcaagtctgccatgccagaaggttatgttcaagaaagaactatttttttcaaagatgacggtaactacaagaccagagctgaagtcaagtttgaaggtgataccttagttaatagaatcgaattaaaaggtattgattttaaagaagGtggtaacattttaggtcacaaattggaatacaactataactctcacaatgtttacatcactgctgacaaacaaaagaatggtatcaaagctaacttcaaaattagacacaacattgaagatggtggtgttcaattagctgaccattatcaacaaaatactccaattggtgatggtccagtcttgttaccagacaaccattacttatcctatcaatctgccttatccaaagatccaaacgaaaagagagaccacatggtcttgttagaatttgttactgctgctggtattacccatggtatggatgaattgtacaaa"
    fiveOH = "GTAA"
    fiveExt = "gggg"
    threeExt = "aaaa"
    threeOH = "AGTT"
    newFrag = GGfrag(fiveOH, fiveExt, seq, threeExt, threeOH)

    oligos = newFrag.getPCAprimers(60, ["GGGGGGG","AAAAAAG"])
    for each in oligos[0]:
        print(each)
        print("\n")
    print("\n\n\n")
    for each in oligos[1]:
        print(each)
        print("\n")
    print("\n")

# 	a = newFrag.divideByCaseChange()


if __name__ == '__main__':
    main()
