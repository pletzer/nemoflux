import numpy

def getArcLength(pointA, pointB, radius=1):
    return abs( radius * numpy.arccos( pointA.dot(pointB)/radius**2 ) )