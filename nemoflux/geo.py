import numpy

EARTH_RADIUS = 6371000. # in metres
DEG2RAD = numpy.pi/180.

def lonLat2XYZ(p, radius):
    lam = p[0]*DEG2RAD
    the = p[1]*DEG2RAD
    rho = radius*numpy.cos(the)
    return  numpy.array([rho*numpy.cos(lam), 
                         rho*numpy.sin(lam), 
                         radius*numpy.sin(the)])

def lonLat2XYZArray(p, radius):
    lam = p[...,0]*DEG2RAD
    the = p[...,1]*DEG2RAD
    rho = radius*numpy.cos(the)
    xyz = numpy.zeros(p.shape, numpy.float64)
    xyz[..., 0] = rho*numpy.cos(lam) # x
    xyz[..., 1] = rho*numpy.sin(lam) # y
    xyz[..., 2] = radius*numpy.sin(the) # z
    return xyz

def getArcLengthArray(xyzA, xyzB, radius=EARTH_RADIUS):
    radiusSquare = radius*radius
    angle = numpy.arccos(  numpy.sum(xyzA*xyzB, axis=2) / radiusSquare )
    return abs(radius * angle)

def getArcLength(xyzA, xyzB, radius=EARTH_RADIUS):
    return abs( radius * numpy.arccos( xyzA.dot(xyzB)/radius**2 ) )