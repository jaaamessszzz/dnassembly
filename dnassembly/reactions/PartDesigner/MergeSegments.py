#!/usr/bin/env python
# encoding: utf-8
"""
MergeSegments.py

Created by Will DeLoache on 2012-09-04.
"""

#Takes a set of segments and tries to combine them where possible
#while keeping the primer length below maxPrimerLength
#segments = [5' extension, template, 3' extension]
#primer length = annealingLength + extension region
def mergeSegments(segments, maxPrimerLength, annealingLength):
	if len(segments) <= 1:
		return segments

	#First try to add the smaller part to the bigger part
	#Then try to add the bigger part to the smaller part
	for i in range(len(segments) - 1):
		seg1 = segments[i]
		seg2 = segments[i+1]
		if len(seg1[1]) >= len(seg2[1]):
			if len(seg1[2] + seg2[0] + seg2[1] + seg2[2]) + annealingLength < maxPrimerLength:
				newSeg = [seg1[0], seg1[1], seg1[2] + seg2[0] + seg2[1] + seg2[2]]
				segments[i:i+2] = [newSeg]
				return mergeSegments(segments, maxPrimerLength, annealingLength)
			elif len(seg1[0] + seg1[1] + seg1[2] + seg2[0]) + annealingLength < maxPrimerLength:
				newSeg = [seg1[0] + seg1[1] + seg1[2] + seg2[0], seg2[1], seg2[2]]
				segments[i:i+2] = [newSeg]
				return mergeSegments(segments, maxPrimerLength, annealingLength)
		else:
			if len(seg1[0] + seg1[1] + seg1[2] + seg2[0]) + annealingLength < maxPrimerLength:
				newSeg = [seg1[0] + seg1[1] + seg1[2] + seg2[0], seg2[1], seg2[2]]
				segments[i:i+2] = [newSeg]
				return mergeSegments(segments, maxPrimerLength, annealingLength)
			elif len(seg1[2] + seg2[0] + seg2[1] + seg2[2]) + annealingLength < maxPrimerLength:
				newSeg = [seg1[0], seg1[1], seg1[2] + seg2[0] + seg2[1] + seg2[2]]
				segments[i:i+2] = [newSeg]
				return mergeSegments(segments, maxPrimerLength, annealingLength)

	#Try to split any remaining short segments between the parts to the left and right
	if len(segments) >= 3:
		i = 1
		while i < len(segments) - 1:
			l_seg = segments[i-1]
			m_seg = segments[i]
			r_seg = segments[i+1]
			#substract an extra 5 bases from the roomToLeft and roomToRight in order to ensure
			#at least 2 overhang options will be available.
			roomToLeft = maxPrimerLength - annealingLength - len(l_seg[2]) - 2
			roomToRight = maxPrimerLength - annealingLength - len(r_seg[0]) - 3
			m_seg_combined = m_seg[0] + m_seg[1] + m_seg[2]
			if len(m_seg_combined) < roomToLeft + roomToRight:
				split_index = int(float(roomToLeft)/(roomToRight+roomToLeft)*len(m_seg_combined))
				new_l_seg = [l_seg[0], l_seg[1], l_seg[2] + m_seg_combined[:split_index]]
				new_r_seg = [m_seg_combined[split_index:] + r_seg[0], r_seg[1], r_seg[2]]
				segments[i-1:i+2] = [new_l_seg, new_r_seg]
			else:
				i += 1

	return segments

def main():
	pass


if __name__ == '__main__':
	main()

