import netCDF4
import numpy
from math import pi, cos, sin
import math
import defopt
import geo

# precision with which data will be saved in the netCDF file
REAL = 'float64'

class DataGen(object):

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

    def build(self):
        self.buildScaledVertical()
        self.buildUniformHorizontal()

    def buildScaledVertical(self):
        #  non-uniform depth
        self.zhalf = numpy.array([self.zmin + (self.zmax - self.zmin)*((k + 0.5)/float(self.nz))**2 for k in range(self.nz)])
        self.ztop  = numpy.array([self.zmin + (self.zmax - self.zmin)*((k + 0  )/float(self.nz))**2 for k in range(self.nz)])
        self.zbot  = numpy.array([self.zmin + (self.zmax - self.zmin)*((k + 1  )/float(self.nz))**2 for k in range(self.nz)])

    def buildUniformHorizontal(self):

        ny1, nx1 = self.ny + 1, self.nx + 1
        dy, dx = (self.ymax - self.ymin)/float(self.ny), (self.xmax - self.xmin)/float(self.nx)

        # cell centres
        self.nav_lat = numpy.zeros((self.ny, self.nx), REAL)
        self.nav_lon = numpy.zeros((self.ny, self.nx), REAL)

        # cell bounds
        self.bounds_lat = numpy.zeros((self.ny, self.nx, 4), REAL)
        self.bounds_lon = numpy.zeros((self.ny, self.nx, 4), REAL)

        # iterate over the horizontal cells
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

                # cell mid point
                self.nav_lat[j, i] = ym
                self.nav_lon[j, i] = xm

    def applyPotential(self, potentialFunction):
        zmin, zmax = self.zmin, self.zmax
        self.potential = numpy.zeros((self.nz, self.ny, self.nx, 4), numpy.float64)
        A = geo.EARTH_RADIUS
        for k in range(self.nz):
            z = self.zhalf[k]
            for j in range(self.ny):
                for i in range(self.nx):
                    for vertex in range(4):
                        y = self.bounds_lat[j, i, vertex]
                        x = self.bounds_lon[j, i, vertex]
                        self.potential[k, j, i, vertex] = eval(potentialFunction)

    def computeUVFromPotential(self):
        self.u = numpy.zeros((self.nz, self.ny, self.nx), numpy.float64)
        self.v = numpy.zeros((self.nz, self.ny, self.nx), numpy.float64)
        for k in range(self.nz):
            for j in range(self.ny):
                for i in range(self.nx):
                    x0, x1, x2, x3 = self.bounds_lon[j, i, :]
                    y0, y1, y2, y3 = self.bounds_lat[j, i, :]
                    p0 = numpy.array((x0, y0, 0.))
                    p1 = numpy.array((x1, y1, 0.))
                    p2 = numpy.array((x2, y2, 0.))
                    p3 = numpy.array((x3, y3, 0.))
                    dPhi21 = self.potential[k, j, i, 2] - self.potential[k, j, i, 1]
                    dPhi23 = self.potential[k, j, i, 2] - self.potential[k, j, i, 3]
                    ds21 = geo.getArcLength(p2, p1, radius=geo.EARTH_RADIUS)
                    ds23 = geo.getArcLength(p2, p3, radius=geo.EARTH_RADIUS)
                    # east, - d phi/ dy
                    self.u[k, j, i] = - dPhi21 / ds21
                    # north, + d phi/ dx
                    self.v[k, j, i] = + dPhi23 / ds23

    def rotatePole(self, deltaDeg=(0., 0.)):

        lats = numpy.zeros((self.ny, self.nx,), numpy.float64)
        lons = numpy.zeros((self.ny, self.nx,), numpy.float64)

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
                # use the convention 0 <= lon < 360
                self.nav_lon[j, i] %= 360.

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
                    # use the convention 0 <= lon < 360
                    # self.bounds_lon[j, i, vertex] %= 360.

                    # date line fix
                    if self.bounds_lon[j, i, vertex] - self.nav_lon[j, i] > 270.:
                        self.bounds_lon[j, i, vertex] -= 360.
                    if self.bounds_lon[j, i, vertex] - self.nav_lon[j, i] < -270.:
                        self.bounds_lon[j, i, vertex] += 360.

    def save(self):

        ncT = netCDF4.Dataset(self.prefix + '_T.nc', 'w')
        ncT.createDimension('z', self.nz)
        ncT.createDimension('y', self.ny)
        ncT.createDimension('x', self.nx)
        ncT.createDimension('nvertex', 4)
        ncT.createDimension('axis_nbounds', 2)
        deptht = ncT.createVariable('deptht', REAL, ('z',))
        deptht.standard_name = 'depth'
        deptht.units = 'm'
        deptht.bounds = 'depthu_bounds'
        deptht[:] = self.zhalf
        deptht_bounds = ncT.createVariable('deptht_bounds', REAL, ('z', 'axis_nbounds'))
        deptht_bounds[:, 0] = self.ztop
        deptht_bounds[:, 1] = self.zbot

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

        potential = ncT.createVariable('potential', REAL, ('z', 'y', 'x', 'nvertex'))
        potential.long_name = 'potential'
        potential[:] = self.potential
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
        uo = ncU.createVariable('uo', REAL, ('z', 'y', 'x'), fill_value=1.e20)
        uo.standard_name = 'sea_water_x_velocity'
        uo.units = 'm/s'
        uo[:] = self.u
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
        depthv_bounds = ncV.createVariable('depthv_bounds', REAL, ('z', 'axis_nbounds'))
        depthv_bounds[:, 0] = self.ztop
        depthv_bounds[:, 1] = self.zbot
        vo = ncV.createVariable('vo', REAL, ('z', 'y', 'x'), fill_value=1.e20)
        vo.standard_name = 'sea_water_y_velocity'
        vo.units = 'm/s'
        vo[:] = self.v
        ncV.close()


def main(*, potentialFunction: str="y", prefix: str, 
            xmin: float=0.0, xmax: float=360., 
            ymin: float=-90., ymax: float=90.,
            zmin: float=0., zmax: float=1.0,
            nx: int=10, ny: int=4, nz: int=5, deltaDeg: str="(0.,0.)"):
    """Generate data
    :param potentialFunction: potential expression of x and y
    :param prefix: file prefix
    :param xmin: min longitude
    :param xmax: max longitude
    :param ymin: min latitude
    :param ymax: max latitude
    :param zmin: min depth
    :param zmax: max depth
    :param nx: number of cells in longitude
    :param ny: number of cells in latitude
    :param nz: number of vertical cells
    :param deltaDeg: longitude, latitude pole displacement 
    """
    lldg = DataGen(prefix)
    lldg.setSizes(nx, ny, nz)
    lldg.setBoundingBox(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=0., zmax=1.)
    lldg.build()
    lldg.rotatePole(deltaDeg=eval(deltaDeg))
    lldg.applyPotential(potentialFunction)
    lldg.computeUVFromPotential()
    lldg.save()

if __name__ == '__main__':
    defopt.run(main)




        

