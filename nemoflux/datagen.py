import netCDF4
import numpy
from math import pi, cos, sin
import math
import defopt
import geo

# precision with which data will be saved in the netCDF file
REAL = 'float64'

class LatLonDataGen(object):

    def __init__(self, prefix=''):
        self.prefix = prefix

    def setBoundingBox(self, xmin, xmax, ymin, ymax, zmin, zmax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax

    def setSizes(self, nx, ny, nz):
        # number of cells in x, y and z
        self.nx = nx
        self.ny = ny
        self.nz = nz

    def build(self, potentialFunction):

        self.uExtensive = numpy.zeros((self.nz, self.ny, self.nx), REAL)
        self.vExtensive = numpy.zeros((self.nz, self.ny, self.nx), REAL)
        self.uArea = numpy.zeros((self.nz, self.ny, self.nx), REAL)
        self.vArea = numpy.zeros((self.nz, self.ny, self.nx), REAL)

        ny1, nx1 = self.ny + 1, self.nx + 1
        dy, dx = (self.ymax - self.ymin)/float(self.ny), (self.xmax - self.xmin)/float(self.nx)

        # cell centres
        self.nav_lat = numpy.zeros((self.ny, self.nx), REAL)
        self.nav_lon = numpy.zeros((self.ny, self.nx), REAL)

        # cell bounds
        self.bounds_lat = numpy.zeros((self.ny, self.nx, 4), REAL)
        self.bounds_lon = numpy.zeros((self.ny, self.nx, 4), REAL)

        #  non-uniform depth
        self.zhalf = numpy.array([self.zmin + (self.zmax - self.zmin)*((k + 0.5)/float(self.nz))**2 for k in range(self.nz)])
        self.ztop  = numpy.array([self.zmin + (self.zmax - self.zmin)*((k + 0  )/float(self.nz))**2 for k in range(self.nz)])
        self.zbot  = numpy.array([self.zmin + (self.zmax - self.zmin)*((k + 1  )/float(self.nz))**2 for k in range(self.nz)])

        A = geo.EARTH_RADIUS # for evaluation

        # iterate over cells
        for k in range(self.nz):
            dz = self.ztop[k] - self.zbot[k]
            z = self.zhalf[k] # z is at half level
            for j in range(self.ny):
                y0 = self.ymin + j*dy
                y1 = y0 + dy
                ym = y0 + 0.5*dy
                for i in range(self.nx):
                    x0 = self.xmin + i*dx
                    x1 = x0 + dx
                    xm = x0 + 0.5*dx

                    self.bounds_lat[j, i, 0] = y0
                    self.bounds_lat[j, i, 1] = y0
                    self.bounds_lat[j, i, 2] = y1
                    self.bounds_lat[j, i, 3] = y1

                    self.bounds_lon[j, i, 0] = x0
                    self.bounds_lon[j, i, 1] = x1
                    self.bounds_lon[j, i, 2] = x1
                    self.bounds_lon[j, i, 3] = x0

                    self.nav_lat[j, i] = ym
                    self.nav_lon[j, i] = xm

                    # east side of the cell
                    x, y = x1, y1,
                    p1 = numpy.array((x, y, 0.))
                    phi1 = eval(potentialFunction)
                    x, y = x1, y0
                    p0 = numpy.array((x, y, 0.))
                    phi0 = eval(potentialFunction)
                    self.uExtensive[k, j, i] = (phi1 - phi0) * dz
                    self.uArea[k, j, i] = geo.getArcLength(p0, p1,  radius=geo.EARTH_RADIUS) * dz

                    # north side of the cell
                    x, y = x1, y1
                    p1 = numpy.array((x, y, 0.))
                    phi1 = eval(potentialFunction)
                    x, y = x0, y1
                    p0 = numpy.array((x, y, 0.))
                    phi0 = eval(potentialFunction)
                    self.vExtensive[k, j, i] = (phi1 - phi0) * dz
                    self.vArea[k, j, i] = geo.getArcLength(p0, p1, radius=geo.EARTH_RADIUS) * dz

    def rotatePole(self, deltaLonDeg=0., deltaLatDeg=0.):

        lats = numpy.zeros((self.ny, self.nx,), numpy.float64)
        lons = numpy.zeros((self.ny, self.nx,), numpy.float64)

        alpha = numpy.pi * deltaLatDeg / 180.
        beta = numpy.pi * deltaLonDeg / 180.
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
        transfMatrix = numpy.dot(rot_bet, rot_alp)

        # original position
        xyz0 = numpy.zeros((3,), numpy.float64)
        # transformed position
        xyz1 = numpy.zeros((3,), numpy.float64)

        for j in range(self.ny):
            for i in range(self.nx):

                the = numpy.pi * self.nav_lat[j, i] / 180.
                lam = numpy.pi * self.nav_lon[j, i] / 180.
                cos_the = numpy.cos(the)
                sin_the = numpy.sin(the)
                rho = cos_the
                cos_lam = numpy.cos(lam)
                sin_lam = numpy.sin(lam)

                xyz0 = rho * cos_lam, rho * sin_lam, sin_the
                xyz1 = numpy.dot(transfMatrix, xyz0)

                self.nav_lat[j, i] = 180. * math.asin(xyz1[2]) / pi
                self.nav_lon[j, i] = 180. * math.atan2(xyz1[1], xyz1[0]) / pi
                if self.nav_lon[j, i] < 0.:
                    self.nav_lon[j, i] += 360.

                for vertex in range(4):

                    the = numpy.pi * self.bounds_lat[j, i, vertex] / 180.
                    lam = numpy.pi * self.bounds_lon[j, i, vertex] / 180.
                    cos_the = numpy.cos(the)
                    sin_the = numpy.sin(the)
                    rho = cos_the
                    cos_lam = numpy.cos(lam)
                    sin_lam = numpy.sin(lam)

                    xyz0 = rho * cos_lam, rho * sin_lam, sin_the
                    xyz1 = numpy.dot(transfMatrix, xyz0)

                    self.bounds_lat[j, i, vertex] = 180. * math.asin(xyz1[2]) / numpy.pi
                    self.bounds_lon[j, i, vertex] = 180. * math.atan2(xyz1[1], xyz1[0]) / numpy.pi
                    if self.bounds_lon[j, i, vertex] < 0:
                        self.bounds_lon[j, i, vertex] += 360.

                    # add/subtract 360 deg to make the cell well behaved
                    if vertex > 0 and self.bounds_lon[j, i, vertex] - self.bounds_lon[j, i, 0] < -270.:
                        self.bounds_lon[j, i, vertex] += 360.
                    if vertex > 0 and self.bounds_lon[j, i, vertex] - self.bounds_lon[j, i, 0] > +270.:
                        self.bounds_lon[j, i, vertex] -= 360.

        # compute the face areas
        for k in range(self.nz):
            dz = self.ztop[k] - self.zbot[k]
            for j in range(self.ny):
                for i in range(self.nx):
                    # east side of the cell
                    p1 = numpy.array((self.bounds_lon[j, i, 2], self.bounds_lon[j, i, 2], 0.))
                    p0 = numpy.array((self.bounds_lon[j, i, 1], self.bounds_lon[j, i, 1], 0.))
                    self.uArea[k, j, i] = geo.getArcLength(p0, p1, radius=geo.EARTH_RADIUS) * dz
                    # north side of the cell
                    p1 = numpy.array((self.bounds_lon[j, i, 3], self.bounds_lon[j, i, 3], 0.))
                    p0 = numpy.array((self.bounds_lon[j, i, 2], self.bounds_lon[j, i, 2], 0.))
                    self.vArea[k, j, i] = geo.getArcLength(p0, p1, radius=geo.EARTH_RADIUS) * dz

    def save(self):

        ncT = netCDF4.Dataset(self.prefix + '_T.nc', 'w')
        ncT.createDimension('y', self.ny)
        ncT.createDimension('x', self.nx)
        ncT.createDimension('nvertex', 4)

        nav_lat = ncT.createVariable('nav_lat', REAL, ('y', 'x'))
        nav_lat.standard_name = 'latitude'
        nav_lat.bounds = 'bounds_lat'
        nav_lat[:] = self.nav_lat

        nav_lon = ncT.createVariable('nv_lon', REAL, ('y', 'x'))
        nav_lon.standard_name = 'longitude'
        nav_lon.bounds = 'bounds_lon'
        nav_lon[:] = self.nav_lon

        bounds_lat = ncT.createVariable('bounds_lat', REAL, ('y', 'x', 'nvertex'))
        bounds_lat[:] = self.bounds_lat

        bounds_lon = ncT.createVariable('bounds_lon', REAL, ('y', 'x', 'nvertex'))
        bounds_lon[:] = self.bounds_lon
        ncT.close()

        ncU = netCDF4.Dataset(self.prefix + '_U.nc', 'w')
        ncU.createDimension('z', self.nz)
        ncU.createDimension('y', self.ny)
        ncU.createDimension('x', self.nx)
        ncU.createDimension('axis_nbounds', 2)
        depthu = ncU.createVariable('depthu', REAL, ('z',))
        depthu.standard_name = 'depth'
        depthu.units = 'm'
        depthu.bounds = 'depthu_bounds'
        depthu[:] = self.zhalf
        depthu_bounds = ncU.createVariable('depthu_bounds', REAL, ('z', 'axis_nbounds'))
        depthu_bounds[:, 0] = self.ztop
        depthu_bounds[:, 1] = self.zbot
        uArea = ncU.createVariable('uArea', REAL, ('z', 'y', 'x'))
        uArea.long_name = 'u face area'
        uArea.units = 'm2'
        uArea[:] = self.uArea
        uo = ncU.createVariable('uo', REAL, ('z', 'y', 'x'), fill_value=1.e20)
        uo.standard_name = 'sea_water_x_velocity'
        uo.units = 'm/s'
        uo[:] = self.uExtensive / self.uArea
        ncU.earthRadius = f'earth radius = {geo.EARTH_RADIUS} in metres'
        ncU.close()

        ncV = netCDF4.Dataset(self.prefix + '_V.nc', 'w')
        ncV.createDimension('z', self.nz)
        ncV.createDimension('y', self.ny)
        ncV.createDimension('x', self.nx)
        ncV.createDimension('axis_nbounds', 2)
        depthv = ncU.createVariable('depthv', REAL, ('z',))
        depthv.standard_name = 'depth'
        depthv.units = 'm'
        depthv.bounds = 'depthv_bounds'
        depthv[:] = self.zhalf
        vArea = ncV.createVariable('vArea', REAL, ('z', 'y', 'x'))
        vArea.long_name = 'v face area'
        vArea.units = 'm2'
        vArea[:] = self.vArea
        depthv_bounds = ncV.createVariable('depthv_bounds', REAL, ('z', 'axis_nbounds'))
        depthv_bounds[:, 0] = self.ztop
        depthv_bounds[:, 1] = self.zbot
        vo = ncV.createVariable('vo', REAL, ('z', 'y', 'x'), fill_value=1.e20)
        vo.standard_name = 'sea_water_y_velocity'
        vo.units = 'm/s'
        vo[:] = self.vExtensive / self.vArea
        ncV.close()


def main(*, potentialFunction: str="(1. - z/self.zmax)*x*2*pi*A/360.", prefix: str, 
            xmin: float=0.0, xmax: float=360., ymin: float=-90., ymax: float=90., 
            nx: int=10, ny: int=4, nz: int=10, deltaLonDeg: float=0., deltaLatDeg: float=0.):
    """Generate data
    :param prefix: file prefix
    :param xmin: min longitude
    :param xmax: max longitude
    :param ymin: min latitude
    :param ymax: max latitude
    :param nx: number of cells in longitude
    :param ny: number of cells in latitude
    :param potentialFunction: potential expression of x and y
    """
    lldg = LatLonDataGen(prefix)
    lldg.setSizes(nx, ny, nz)
    lldg.setBoundingBox(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=0., zmax=1000.)
    lldg.build(potentialFunction)
    lldg.rotatePole(deltaLonDeg=deltaLonDeg, deltaLatDeg=deltaLatDeg)
    lldg.save()

if __name__ == '__main__':
    defopt.run(main)




        

