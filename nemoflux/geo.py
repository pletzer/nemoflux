import numpy

EARTH_RADIUS = 6371000. # in metres

def getArcLengthArray(xyzA, xyzB, radius=EARTH_RADIUS):
    radiusSquare = radius*radius
    angle = numpy.arccos(  numpy.tensordot(xyzA, xyzB, (1, 1) ) / radiusSquare )
    return abs(radius * angle)

def getArcLength(xyzA, xyzB, radius=EARTH_RADIUS):
    return abs( radius * numpy.arccos( xyzA.dot(xyzB)/radius**2 ) )