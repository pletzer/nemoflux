import mint
import numpy
import netCDF4
from horizgrid import HorizGrid
import geo
import defopt


class FluxCalc(object):

    def __init__(self, tFile, uFile, vFile, targetPoints):

        self.gr = HorizGrid(tFile)
        self.pli = mint.PolylineIntegral()
        self.pli.build(self.gr.getMintGrid(), targetPoints,
                       counterclock=False, periodX=360.0)

        with netCDF4.Dataset(tFile) as nc:
            depth_half = nc.variables['deptht'][:]
            bounds_depth = nc.variables['deptht_bounds'][:]

        with netCDF4.Dataset(uFile) as nc:
            uo = nc.variables['uo'][:]
        # set the velocity to zero where missing
        uo = numpy.ma.filled(uo, 0.0)

        with netCDF4.Dataset(vFile) as nc:
            vo = nc.variables['vo'][:]
        # set the velocity to zero where missing
        vo = numpy.ma.filled(vo, 0.0)

        nz, ny, nx = uo.shape

        numCells = self.gr.getNumCells()

        # integrate vertically, multiplying by the thickness of the layers
        dz = -(bounds_depth[:, 1] - bounds_depth[:, 0]) # DEPTH HAS OPPOSITE SIGN TO Z
        uo = numpy.tensordot(dz, uo, axes=(0, 0)) # sum of multiplying axis 0 of dz with axis 0 of uo
        vo = numpy.tensordot(dz, vo, axes=(0, 0))

        points = self.gr.getPoints().reshape((ny, nx, 4, 3))

        xyz0 = geo.lonLat2XYZArray(points[:, :, 0, :], radius=geo.EARTH_RADIUS)
        xyz1 = geo.lonLat2XYZArray(points[:, :, 1, :], radius=geo.EARTH_RADIUS)
        xyz2 = geo.lonLat2XYZArray(points[:, :, 2, :], radius=geo.EARTH_RADIUS)
        xyz3 = geo.lonLat2XYZArray(points[:, :, 3, :], radius=geo.EARTH_RADIUS)

        ds01 = geo.getArcLengthArray(xyz0, xyz1, radius=geo.EARTH_RADIUS)
        ds12 = geo.getArcLengthArray(xyz1, xyz2, radius=geo.EARTH_RADIUS)
        ds32 = geo.getArcLengthArray(xyz3, xyz2, radius=geo.EARTH_RADIUS)
        ds03 = geo.getArcLengthArray(xyz0, xyz3, radius=geo.EARTH_RADIUS)

        self.integratedVelocity = numpy.zeros((ny, nx, 4), numpy.float64)
        #        ^
        #        |
        #  3-----V-----2
        #  |           |
        #  |->   T     U->
        #  |     ^     |
        #  |     |     |
        #  0-----------1

        # south
        self.integratedVelocity[1:, :, 0] = vo[:-1, :] * ds01[:-1, :]
        # else:
        #     # folding, assuming even nx
        #     self.integratedVelocity[cellId, 0] = vo[:, j, (i+nx//2)%nx] * ds01

        # east
        self.integratedVelocity[..., 1] = uo * ds12

        # north
        self.integratedVelocity[..., 2] = vo * ds32

        # periodic west
        self.integratedVelocity[:, 1:, 3] = uo[:, :-1] * ds03[:, :-1]
        self.integratedVelocity[:, 0, 3] = uo[:, -1] * ds03[:, -1]

        # array should have shape (numCells, 4)
        self.integratedVelocity = self.integratedVelocity.reshape((numCells, 4))


    def getFlux(self):
        return self.pli.getIntegral(self.integratedVelocity)


def main(*, tFile: str, uFile: str, vFile: str, xyStr: str):
    """Compute flux
    :param tFile: netCDF file containing t grid data
    :param uFile: netCDF file containing U component data
    :param vFile: netCDF file containing V component data
    :param xyStr: array of target points
    """
    xyVals = numpy.array(eval(xyStr))
    numTargetPoints = xyVals.shape[0]
    targetPoints = numpy.zeros((numTargetPoints, 3), numpy.float64)
    targetPoints[:, :2] = xyVals
    fc = FluxCalc(tFile, uFile, vFile, targetPoints)
    print(f'flux: {fc.getFlux()}')

if __name__ == '__main__':
    defopt.run(main)

