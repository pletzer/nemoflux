import numpy

EARTH_RADIUS = 1.0 # 6371000. # in metres
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
    angle = numpy.arccos(  numpy.sum(xyzA*xyzB, axis=-1) / radiusSquare )
    return numpy.fabs(radius * angle)

def getArcLength(xyzA, xyzB, radius=EARTH_RADIUS):
    return abs( radius * numpy.arccos( xyzA.dot(xyzB)/radius**2 ) )


class PoleRotator(object):

    def __init__(self, deltaDeg=(0., 0.)):

        # transformation matrix
        alpha = numpy.pi * deltaDeg[1] / 180.
        beta = numpy.pi * deltaDeg[0] / 180.
        cos_alp = numpy.cos(alpha)
        sin_alp = numpy.sin(alpha)
        cos_bet = numpy.cos(beta)
        sin_bet = numpy.sin(beta)

        # http://gis.stackexchange.com/questions/10808/lon-lat-transformation
        rot_alp = numpy.array([[ cos_alp, 0., sin_alp],
                               [ 0.,      1., 0.     ],
                               [-sin_alp, 0., cos_alp]])
        rot_bet = numpy.array([[ cos_bet, sin_bet, 0.],
                               [-sin_bet, cos_bet, 0.],
                               [ 0.     , 0.,      1.]])
        self.transfMatrix = numpy.dot(rot_bet, rot_alp)


    def apply(self, lons, lats):

        lam = lons * DEG2RAD
        the = lats * DEG2RAD

        cos_the = numpy.cos(the)
        sin_the = numpy.sin(the)
        cos_lam = numpy.cos(lam)
        sin_lam = numpy.sin(lam)

        # to Cartesian
        rho = cos_the
        xx = rho * cos_lam
        yy = rho * sin_lam
        zz = sin_the

        # apply rotation
        xyz = numpy.zeros([3] + list(xx.shape), numpy.float64)
        xyz[0, ...] = xx
        xyz[1, ...] = yy
        xyz[2, ...] = zz
        xyzTransformed = numpy.tensordot(self.transfMatrix, xyz, axes=(0, 0))

        # back to lons and lats
        lonsNew = numpy.arctan2(xyzNew[1,...], xyzNew[0,...]) / DEG2RAD
        latsNew = numpy.arcsin(xyzNew[2,...]) / DEG2RAD

        # NOTE may need to apply dateline fix in some cases

        return lonsNew, latsNew
