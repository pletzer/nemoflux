import numpy

EARTH_RADIUS = 6371000. # in metres

def getArcLength(pointA, pointB, radius=EARTH_RADIUS):
    return abs( radius * numpy.arccos( pointA.dot(pointB)/radius**2 ) )