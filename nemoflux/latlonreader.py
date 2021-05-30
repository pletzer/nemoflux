import re
import sys
import numpy

PAT = re.compile(r'^\s*\d+\s+\d+\.\d+\s+(\-?\d+\.\d+)\s+(\-?\d+\.\d+)')
class LatLonReader(object):

	def __init__(self, filename):

		self.lonLatTargets = []

		with open(filename) as f:
			for line in f.readlines():
				m = re.match(PAT, line)
				if m:
					lat, lon = float(m.group(1)), float(m.group(2))
					self.lonLatTargets.append((lon, lat, 0.0))

	def getLonLats(self):
		return numpy.array(self.lonLatTargets)

###############################################################################
if __name__ == '__main__':
	llreader = LatLonReader(sys.argv[1])
	latLonTargets = llreader.getLonLats()
	print(latLonTargets)