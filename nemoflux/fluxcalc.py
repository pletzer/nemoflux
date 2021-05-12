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

        numCells = self.gr.getNumCells()
        self.integratedVelocity = numpy.zeros((numCells, 4), numpy.float64)

        ncU = netCDF4.Dataset(uFile)
        uo = ncU.variables['uo'][:]
        ncU.close()

        ncV = netCDF4.Dataset(vFile)
        vo = ncV.variables['vo'][:]
        ncV.close()
       
        ny, nx = uo.shape

        k = 0
        for j in range(ny):
            for i in range(nx):

                x0, y0, _ = self.gr.getPoint(k, 0)
                x1, y1, _ = self.gr.getPoint(k, 1)
                x2, y2, _ = self.gr.getPoint(k, 2)
                x3, y3, _ = self.gr.getPoint(k, 3)

                # east
                ds = numpy.sqrt((x2 - x1)**2 + (y2 - y1)**2)
                self.integratedVelocity[k, 1] = uo[j, i] * ds

                # north
                ds = numpy.sqrt((x3 - x2)**2 + (y3 - y2)**2)
                self.integratedVelocity[k, 2] = vo[j, i] * ds

                k += 1

    def getFlux(self):
        return self.pli.getIntegral(self.integratedVelocity)

def main(*, tFile: str, uFile: str, vFile: str, xyStr: str):
    """Compute flux
    :param tFile: netCDF file containing t grid data
    :param uFile: netCDF file containing U component data
    :param vFile: netCDF file containing V component data
    """
    xyVals = numpy.array(eval(xyStr))
    numPoints = xyVals.shape[0]
    targetPoints = numpy.zeros((numPoints, 3), numpy.float64)
    targetPoints[:, :2] = xyVals
    fc = FluxCalc(tFile, uFile, vFile, targetPoints)
    print(f'flux: {fc.getFlux()}')

if __name__ == '__main__':
    defopt.run(main)

