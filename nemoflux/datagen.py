import netCDF4
import numpy
import defopt

class LatLonDataGen(object):

    def __init__(self, prefix=''):
        self.prefix = prefix
        self.nav_lat = []
        self.nav_lon = []
        self.bounds_lat = []
        self.bounds_lon = []

    def setBoundingBox(self, xmin, xmax, ymin, ymax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    def setSizes(self, nx, ny):
        self.nx = nx
        self.ny = ny

    def build(self, potentialFunction):

        self.u = numpy.zeros((self.ny, self.nx), numpy.float32)
        self.v = numpy.zeros((self.ny, self.nx), numpy.float32)

        ny1, nx1 = self.ny + 1, self.nx + 1
        dy, dx = (self.ymax - self.ymin)/float(self.ny), (self.xmax - self.xmin)/float(self.nx)
        self.nav_lat = numpy.zeros((self.ny, self.nx), numpy.float32)
        self.nav_lon = numpy.zeros((self.ny, self.nx), numpy.float32)
        self.bounds_lat = numpy.zeros((self.ny, self.nx, 4), numpy.float32)
        self.bounds_lon = numpy.zeros((self.ny, self.nx, 4), numpy.float32)
        # iterate over cells
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
                self.nav_lat[j, i] = xm

                # east side of the cell
                x, y = x1, y1
                phi1 = eval(potentialFunction)
                x, y = x1, y0
                phi0 = eval(potentialFunction)
                self.u[j, i] = phi1 - phi0

                # north side of the cell
                x, y = x1, y1
                phi1 = eval(potentialFunction)
                x, y = x0, y1
                phi0 = eval(potentialFunction)
                self.v[j, i] = phi1 - phi0


    def save(self):

        ncT = netCDF4.Dataset(self.prefix + '_T.nc', 'w')
        ncT.createDimension('y', self.ny)
        ncT.createDimension('x', self.nx)
        ncT.createDimension('nvertex', 4)

        nav_lat = ncT.createVariable('nav_lat', 'float32', ('y', 'x'))
        nav_lat.standard_name = 'latitude'
        nav_lat.bounds = 'bounds_lat'
        nav_lat[:] = self.nav_lat

        nav_lon = ncT.createVariable('nv_lon', 'float32', ('y', 'x'))
        nav_lon.standard_name = 'longitude'
        nav_lon.bounds = 'bounds_lon'
        nav_lon[:] = self.nav_lon

        bounds_lat = ncT.createVariable('bounds_lat', 'float32', ('y', 'x', 'nvertex'))
        bounds_lat[:] = self.bounds_lat

        bounds_lon = ncT.createVariable('bounds_lon', 'float32', ('y', 'x', 'nvertex'))
        bounds_lon[:] = self.bounds_lon

        ncT.close()

        ncU = netCDF4.Dataset(self.prefix + '_U.nc', 'w')
        ncU.createDimension('y', self.ny)
        ncU.createDimension('x', self.nx)
        uo = ncU.createVariable('uo', 'float32', ('y', 'x'))
        uo.standard_name = 'sea_water_x_velocity'
        uo.units = 'm/s'
        uo[:] = self.u
        ncU.close()

        ncV = netCDF4.Dataset(self.prefix + '_V.nc', 'w')
        ncV.createDimension('y', self.ny)
        ncV.createDimension('x', self.nx)
        vo = ncV.createVariable('vo', 'float32', ('y', 'x'))
        vo.standard_name = 'sea_water_y_velocity'
        vo.units = 'm/s'
        vo[:] = self.v
        ncV.close()


def main(*, potentialFunction: str="x", prefix: str, 
            xmin: float=0.0, xmax: float=360., ymin: float=-90., ymax: float=90., 
            nx: int=10, ny: int=4):
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
    lldg.setSizes(nx, ny)
    lldg.setBoundingBox(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
    lldg.build(potentialFunction)
    lldg.save()

if __name__ == '__main__':
    defopt.run(main)




        

