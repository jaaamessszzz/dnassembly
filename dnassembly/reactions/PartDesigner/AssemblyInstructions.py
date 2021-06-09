#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
# encoding: utf-8
"""
AssemIns.py

Created by Will DeLoache on 2012-09-12.
"""
from Bio.Seq import Seq
import pdb

class AssemblyInstructions():
	# possible_templates[name] = {"seq":sequence, "Linear":True}
	def __init__(self, method, primers, GGfrag, possible_templates):
		tails = GGfrag.tails
		self.method = method
		self.primers = primers
		self.product = ""
		self.digest = ""
		if self.method not in ["PCR", "gBlocks", "PCA"]:
			self.product = GGfrag.fiveprimeOH + GGfrag.fiveprimeExt + GGfrag.seq + GGfrag.threeprimeExt + GGfrag.threeprimeOH
		else:
			self.product = tails[0] + GGfrag.fiveprimeOH + GGfrag.fiveprimeExt + GGfrag.seq + GGfrag.threeprimeExt + GGfrag.threeprimeOH + tails[1]
			self.digest = GGfrag.fiveprimeOH + GGfrag.fiveprimeExt + GGfrag.seq + GGfrag.threeprimeExt + GGfrag.threeprimeOH

		self.template = ""
		self.templateSeq = ""
		if self.method == "PCR":
			possTemp = findTemplate(GGfrag.seq, possible_templates)
			if possTemp:
				self.template = possTemp
				self.templateSeq = possible_templates[possTemp]


	def __str__(self):
		return str(self.method) + "\n" + str(self.primers) + "\n" + str(self.template) + "\n" + str(self.product)


# possible_templates[name] = {"seq":sequence, "Linear":True}
def findTemplate(seq, possible_templates):
	s = seq.upper()
	keys = possible_templates.keys()
	#pdb.set_trace()
	#keys.sort()
	for possTemp in keys:
		temp_seq = possible_templates[possTemp].upper()
		if possible_templates[possTemp]:
			temp_seq += possible_templates[possTemp].upper()
		if temp_seq.upper().find(s) > -1:
			return possTemp
		elif Seq(temp_seq.upper()).reverse_complement().find(s) > -1:
			return possTemp
	return False

def main():
	pass


if __name__ == '__main__':
	main()
