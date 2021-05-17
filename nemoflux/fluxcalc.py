import mint
import numpy
import netCDF4
from grid import Grid
import defopt

class FluxCalc(object):

    def __init__(self, tFile, uFile, vFile, targetPoints):

        self.gr = Grid(tFile)
        self.pli = mint.PolylineIntegral()
        self.pli.build(self.gr.getMintGrid(), targetPoints,
                       counterclock=False, periodX=360.0)

        ncU = netCDF4.Dataset(uFile)
        uo = ncU.variables['uo'][:]
        ncU.close()
        # set the velocity to zero where missing
        uo = numpy.ma.filled(uo, 0.0)

        ncV = netCDF4.Dataset(vFile)
        vo = ncV.variables['vo'][:]
        ncV.close()
        # set the velocity to zero where missing
        vo = numpy.ma.filled(vo, 0.0)
       
        nz, ny, nx = uo.shape

        numCells = self.gr.getNumCells()
        self.integratedVelocity = numpy.zeros((nz, numCells, 4), numpy.float64)

        for k in range(nz):
            cellId = 0
            for j in range(ny):
                for i in range(nx):

                    x0, y0, _ = self.gr.getPoint(cellId, 0)
                    x1, y1, _ = self.gr.getPoint(cellId, 1)
                    x2, y2, _ = self.gr.getPoint(cellId, 2)
                    x3, y3, _ = self.gr.getPoint(cellId, 3)

                    # NEED TO USE ds IN METRES (TO CHANGE)

                    # south
                    ds = numpy.sqrt((x1 - x0)**2 + (y1 - y0)**2)
                    if j >= 1:
                        ds = self.gr.getEdgeArcLength(cellId, 0)
                        self.integratedVelocity[k, cellId, 0] = vo[k, j - 1, i] * ds

                    # east
                    ds = self.gr.getEdgeArcLength(cellId, 1)
                    self.integratedVelocity[k, cellId, 1] = uo[k, j, i] * ds

                    # north
                    ds = self.gr.getEdgeArcLength(cellId, 2)
                    self.integratedVelocity[k, cellId, 2] = vo[k, j, i] * ds

                    # west, periodic boundary
                    ds = self.gr.getEdgeArcLength(cellId, 3)
                    self.integratedVelocity[k, cellId , 3] = uo[k, j, i - 1] * ds

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
    numPoints = xyVals.shape[0]
    targetPoints = numpy.zeros((numPoints, 3), numpy.float64)
    targetPoints[:, :2] = xyVals
    fc = FluxCalc(tFile, uFile, vFile, targetPoints)
    print(f'flux: {fc.getFlux()}')

if __name__ == '__main__':
    defopt.run(main)

