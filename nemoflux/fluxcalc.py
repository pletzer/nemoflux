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
        self.integratedVelocity = numpy.zeros((numCells, 4), numpy.float64)

        # integrate vertically, multiplying by the thickness of the layers
        dz = -(bounds_depth[:, 1] - bounds_depth[:, 0]) # DEPTH HAS OPPOSITE SIGN TO Z
        uo = numpy.tensordot(dz, uo, axes=(0, 0)) # sum of multiplying axis 0 of dz with axis 0 of uo
        vo = numpy.tensordot(dz, vo, axes=(0, 0))

        cellId = 0
        for j in range(ny):
            for i in range(nx):

                p0 = self.gr.getPoint(cellId, 0)
                p1 = self.gr.getPoint(cellId, 1)
                p2 = self.gr.getPoint(cellId, 2)
                p3 = self.gr.getPoint(cellId, 3)

                ds01 = geo.getArcLength(p0, p1, radius=geo.EARTH_RADIUS)
                ds12 = geo.getArcLength(p1, p2, radius=geo.EARTH_RADIUS)
                ds32 = geo.getArcLength(p3, p2, radius=geo.EARTH_RADIUS)
                ds03 = geo.getArcLength(p0, p3, radius=geo.EARTH_RADIUS)

                # integrate vertically

                #        ^
                #        |
                #  3-----V-----2
                #  |           |
                #  |->   T     U->
                #  |     ^     |
                #  |     |     |
                #  0-----------1

                # south
                if j > 0:
                    self.integratedVelocity[cellId, 0] = vo[j-1, i+0] * ds01
                # else:
                #     # folding, assuming even nx
                #     self.integratedVelocity[cellId, 0] = vo[:, j, (i+nx//2)%nx] * ds01

                # east
                self.integratedVelocity[cellId, 1] = uo[j+0, i+0] * ds12

                # north
                self.integratedVelocity[cellId, 2] = vo[j+0, i+0] * ds32

                # periodic west
                self.integratedVelocity[cellId, 3] = uo[j+0, (i-1)%nx] * ds03

                cellId += 1


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

