import netCDF4
import geo
import defopt
import numpy
from horizgrid import HorizGrid
import mint
# to read the target points from file
from latlonreader import LatLonReader


EARTH_RADIUS = 6371000.0


class Field(object):

    def __init__(self, tFile, uFile, vFile, lonLatPoints, sverdrup=False):

        self.sverdrup = sverdrup

        # read the cell bounds
        with netCDF4.Dataset(tFile) as nc:
            bounds_lat = nc.variables['bounds_lat'][:]
            bounds_lon = nc.variables['bounds_lon'][:]
            self.bounds_depth = nc.variables['deptht_bounds'][:]

        self.lonmin = bounds_lon.min()
        self.lonmax = bounds_lon.max()
        self.latmin = bounds_lat.min()
        self.latmax = bounds_lat.max()
        print(f'lon-lat box: {self.lonmin}, {self.latmin} -> {self.lonmax}, {self.latmax}')

        self.timeIndex = 0
        self.ncU = netCDF4.Dataset(uFile)
        self.ncV = netCDF4.Dataset(vFile)

        self.nt, self.nz, self.ny, self.nx = self.getSizes()

        self.gr = HorizGrid(tFile)
        self.pli = mint.PolylineIntegral()
        self.pli.build(self.gr.getMintGrid(), lonLatPoints,
                       counterclock=False, periodX=360.0)


        self.thickness = self.bounds_depth[:, 1] - self.bounds_depth[:, 0]

        numCells = self.ny * self.nx
        self.dx = min((self.lonmax - self.lonmin)/float(self.nx), (self.latmax - self.latmin)/float(self.ny))
        self.arcLengths = numpy.zeros((numCells, 4), numpy.float64)
        self.computeArcLengths()

        # compute the edge fluxes from the vector fields
        self.edgeFluxesUArray = numpy.zeros((numCells,), numpy.float64)
        self.edgeFluxesVArray = numpy.zeros((numCells,), numpy.float64)
        self.integratedVelocity = numpy.zeros((numCells, 4), numpy.float64)
        self.maxAbsFlux = 0.

        # read/get the staggered field integrated over the depth
        uVerticallyIntegrated, vVerticallyIntegrated = self.getUV()
        self.computeIntegratedFlux(uVerticallyIntegrated, vVerticallyIntegrated)
        print(f'max vertically integrated edge |flux|: {self.maxAbsFlux}')

        self.buildEdgeUVGrids(bounds_lon, bounds_lat)

        # create a polyline grid along the target points and place vectors on them
        vectorPoints = []
        self.uVectors = []
        for i in range(len(lonLatPoints) - 1):
            begPoint = lonLatPoints[i]
            endPoint = lonLatPoints[i + 1]
            u = endPoint - begPoint
            distance = numpy.sqrt(u.dot(u))
            # normalize
            u /= distance
            nvpts = max(2, int(distance / self.dx))
            vdx = distance / float(nvpts - 1)
            for j in range(nvpts):
                vectorPoints.append(begPoint + u*j*vdx)
                self.uVectors.append(u)
        self.vectorPoints = numpy.array(vectorPoints)

        # compute the vector at the target line
        self.vinterp = mint.VectorInterp()
        self.vinterp.setGrid(self.gr.getMintGrid())
        self.vinterp.buildLocator(numCellsPerBucket=128, periodX=360.)
        self.vinterp.findPoints(self.vectorPoints, tol2=1.e-12)
        self.vectorValues = self.vinterp.getFaceVectors(self.integratedVelocity)

    def update(self):

        # update the data
        uVerticallyIntegrated, vVerticallyIntegrated = self.getUV()
        # this will update self.integratedVelocity
        self.computeIntegratedFlux(uVerticallyIntegrated, vVerticallyIntegrated)

        self.vectorValues[:] = self.vinterp.getFaceVectors(self.integratedVelocity)

    def getSizes(self):
        # get the sizes
        nt, nz, ny, nx = 1, 1, 0, 0
        shapeU = self.ncU.variables['uo'].shape
        try:
            nt, nz, ny, nx = shapeU
        except:
            try:
                nz, ny, nx = shapeU
            except:
                try:
                    ny, nx = shapeU
                except:
                    raise RuntimeError("ERROR: uo's shape does not match (t, z, y, x), (z, y, x) or (y, x)")
        return nt, nz, ny, nx

    def buildEdgeUVGrids(self, bounds_lon, bounds_lat):

        # points, 4 points per cell, 3D
        self.lonlat = numpy.zeros((self.ny, self.nx, 4, 3), numpy.float64)
        self.lonlat[..., 0] = bounds_lon
        self.lonlat[..., 1] = bounds_lat

    def readField(self, nc, fieldName):

        # read u
        try:
            field = nc.variables[fieldName][self.timeIndex, :, :, :]
        except:
            try:
                field = nc.variables[fieldName][...]
            except:
                raise RuntimeError(f'ERROR: could not read {fieldName} field')

        # set to zero where missing
        field = numpy.ma.filled(field, 0.0)

        # integrate vertically, multiplying by the thickness of the layers,
        # this is the sum of multiplying axis 0 of thickness with axis 0 of field...
        fieldVerticallyIntegrated = numpy.tensordot(self.thickness, field, axes=(0, 0))

        return fieldVerticallyIntegrated

    def getUV(self):
        uVerticallyIntegrated = self.readField(self.ncU, 'uo')
        vVerticallyIntegrated = self.readField(self.ncV, 'vo')
        return uVerticallyIntegrated, vVerticallyIntegrated

    def computeArcLengths(self):

        # lon, lat, elev coords
        points = self.gr.getPoints() # array of size (numCells, 4, 3)

        # convert to Cartesian
        xyz = geo.lonLat2XYZArray(points, radius=geo.EARTH_RADIUS) # array of size (numCells, 4, 3)

        # compute the arc lengths
        for i0 in range(4):
            i1 = (i0 + 1) % 4
            self.arcLengths[:, i0] = geo.getArcLengthArray(xyz[:, i0, :], xyz[:, i1, :], radius=geo.EARTH_RADIUS)

    def computeIntegratedFlux(self, uVerticallyIntegrated, vVerticallyIntegrated):

        numCells = self.ny * self.nx
        #        ^
        #        |
        #  3-----V-----2
        #  |           |
        #  |     T     U-->
        #  |           |
        #  0-----------1

        # cross product means that the flux in V points to the negative y
        self.edgeFluxesUArray[:] = + uVerticallyIntegrated.reshape((numCells,)) * self.arcLengths[:, 1]
        self.edgeFluxesVArray[:] = - vVerticallyIntegrated.reshape((numCells,)) * self.arcLengths[:, 2]

        #        2
        #        ^
        #        |
        #  3-----V-----2
        #  |           |
        #  |->3  0     U->1
        #  |     ^     |
        #  |     |     |
        #  0-----------1

        # east
        self.integratedVelocity[:, 1] = self.edgeFluxesUArray
        # north
        self.integratedVelocity[:, 2] = self.edgeFluxesVArray

        # these should provide a different view to the array, NOT A COPY
        eU = self.edgeFluxesUArray.reshape((self.ny, self.nx))
        eV = self.edgeFluxesVArray.reshape((self.ny, self.nx))
        iV = self.integratedVelocity.reshape((self.ny, self.nx, 4))

        # south, (ny, nx, 4)
        iV[1:, :, 0] = eV[:-1, :] # will set self.integratedVelocity on the south side
        # west
        iV[:, 1:, 3] = eU[:, :-1]
        # periodic BCs
        iV[:, 0, 3] = eU[:, -1]

        if self.sverdrup:
            eU *= EARTH_RADIUS / 1.e6
            eV *= EARTH_RADIUS / 1.e6
            iV *= EARTH_RADIUS / 1.e6

        # from now on, edge fluxes are abs values!
        self.edgeFluxesUArray[:] = numpy.fabs(self.edgeFluxesUArray)
        self.edgeFluxesVArray[:] = numpy.fabs(self.edgeFluxesVArray)

        self.maxAbsFlux = max(self.maxAbsFlux, self.edgeFluxesUArray.max(), self.edgeFluxesVArray.max())

def __del__(self):
    self.ncU.close()
    self.ncV.close()


def main(*, tFile: str, uFile: str, vFile: str, lonLatPoints: str='', iFile: str='', sverdrup: bool=False):
    """Visualize fluxes
    :param tFile: netcdf file holding the T-grid
    :param uFile: netcdf file holding u data
    :param vFile: netcdf file holding v data
    :param lonLatPoints: target points "(lon0, lat0), (lon1, lat1),..."
    :param iFile: alternatively read target points from text file
    """
    if lonLatPoints:
        xyVals = numpy.array(eval(lonLatPoints))
    elif iFile:
        llreader = LatLonReader(iFile)
        xyVals = llreader.getLonLats()
    else:
        raise RuntimeError('ERROR must provide either iFile or lonLatPoints!')
    print(f'target points:\n {xyVals}')
    numTargetPoints = xyVals.shape[0]
    lonLatZPoints = numpy.zeros((numTargetPoints, 3), numpy.float64)
    lonLatZPoints[:, :2] = xyVals
    field = Field(tFile, uFile, vFile, lonLatZPoints, sverdrup)


if __name__ == '__main__':
    defopt.run(main)
    