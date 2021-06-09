#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
# encoding: utf-8
"""
GG_OHcheck.py

Supplies the overhangOK method:

overhangOK(oh1, existingOHs, constraints)

Given a new overhang, a list of existing overhangs, and a constraints variable,
overhangOK() will return True if the overhang is compatible with the existing
overhangs and False if it is not.

Constraints range from 1-16, where 1 is the strictest and 16 is the most lenient.

Full explanation of constraint levels:

1-4 ensures longest common subsequence of 2
5-8 prevents all overhangs with 3 matching bases
9-12 prevents all overhangs with 3 consecutive matching bases
13-16 prevents exact matches

Within each of the above categories:
	- at least 1 G/C, 1 A/T, no 3-base repeats (1, 5, 9, 12)
	- at least 1 G/C, 1 A/T
	- no 4-base repeats
	- none

Created by Will DeLoache on 2012-09-04.
"""

def _LCS(X, Y):
    m = len(X)
    n = len(Y)
    # An (m+1) times (n+1) matrix
    C = [[0] * (n+1) for i in range(m+1)]
    for i in range(1, m+1):
        for j in range(1, n+1):
            if X[i-1] == Y[j-1]:
                C[i][j] = C[i-1][j-1] + 1
            else:
                C[i][j] = max(C[i][j-1], C[i-1][j])
    return C

def _backTrack(C, X, Y, i, j):
    if i == 0 or j == 0:
        return ""
    elif X[i-1] == Y[j-1]:
        return _backTrack(C, X, Y, i-1, j-1) + X[i-1]
    else:
        if C[i][j-1] > C[i-1][j]:
            return _backTrack(C, X, Y, i, j-1)
        else:
            return _backTrack(C, X, Y, i-1, j)

def _reversecomplement(sequence):
	"""Return the reverse complement of the dna string."""
	#complement = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
	complement = {"A":"T", "C":"G", "G":"C", "T":"A", "R":"Y", "Y":"R",
	 				"S":"S", "W":"W", "K":"M", "M":"K", "B":"V", "D":"H",
					"H":"D", "V":"B", "N":"N"}
	reverse_complement_sequence = ""
	sequence_list = list(sequence)
	sequence_list.reverse()
	for letter in sequence_list:
		reverse_complement_sequence += complement[letter.upper()]
	return reverse_complement_sequence

def _countBaseMatches(oh1, oh2):
	matches = 0
	for index, char in enumerate(oh1):
		if oh1[index]  == oh2[index]:
			matches += 1
	return matches

def _consec3BaseMatch(oh1, oh2):
	if oh1[:3] == oh2[:3]:
		return True
	elif oh1[1:] == oh2[1:]:
		return True
	else:
		return False

def overhangOK(oh1, existingOHs, constraints):
	#Remove overhangs with degenerate bases
	for char in oh1:
		if char not in ["A", "T", "C", "G"]:
			return False

	if constraints in [1, 5, 9, 13]:
		#Remove overhangs with the same base repeated 3 times
		if oh1.count("A") > 2:
			return False
		if oh1.count("T") > 2 :
			return False
		if oh1.count("C") > 2:
			return False
		if oh1.count("G") > 2:
			return False

	if constraints in [1, 2, 5, 6, 9, 10, 13, 14]:
		#Remove overhangs with all A/Ts or all G/Cs
		if oh1.count("G") + oh1.count("C") == 0:
			return False
		if oh1.count("G") + oh1.count("C") == 4:
			return False

	if constraints in [3, 7, 11, 15]:
		#Remove overhangs with the same base repeated 4 times
		if oh1.count("A") > 3:
			return False
		if oh1.count("T") > 3 :
			return False
		if oh1.count("C") > 3:
			return False
		if oh1.count("G") > 3:
			return False

	#Make a list of overhangs to check against
	ohsToCheck = []
	ohsToCheck.append(_reversecomplement(oh1))
	ohsToCheck += existingOHs
	for oh2 in existingOHs:
		ohsToCheck.append(_reversecomplement(oh2))

	#Check overhang against ohsToCheck list

	#Remove overhangs whose LCS > 2 (e.g. AcTG and ATGa)
	if constraints <= 4:
		for oh2 in ohsToCheck:
			longComSubSeq = _backTrack(_LCS(oh1, oh2), oh1, oh2, len(oh1), len(oh2))
			if len(longComSubSeq) > 2:
				return False

	#Remove overhangs who share 3 of 4 bases (e.g. ATaG and ATcG)
	if constraints in range (5, 9):
		for oh2 in ohsToCheck:
			baseMatches = _countBaseMatches(oh1, oh2)
			if baseMatches > 2:
				return False

	#Remove overhangs who share 3 of 4 consecutive bases (e.g. ATGa and ATGt)
	if constraints in range (9, 13):
		for oh2 in ohsToCheck:
			if _consec3BaseMatch(oh1, oh2):
				return False

	#Remove exact matches
	if constraints >= 13:
		for oh2 in ohsToCheck:
			if oh1 == oh2:
				return False
	return True

def main():
	newOH = "TGAA"
	existingOHs = []
	result = overhangOK(newOH, existingOHs, 1)
	print(result)


if __name__ == '__main__':
	main()
