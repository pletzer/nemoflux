import mint
import netCDF4
import numpy
import defopt
import re
import vtk

class HorizGrid(object):

    def __init__(self, tFile):
    
        nc = netCDF4.Dataset(tFile)
        # read the cell bounds
        bounds_lat = nc.variables['bounds_lat'][:]
        bounds_lon = nc.variables['bounds_lon'][:]
        # read the potential values
        potential = nc.variables['potential'][:]
        nc.close()

        ny, nx, nvertex = bounds_lat.shape
        numCells = ny * nx
        self.points = numpy.zeros((ny, nx, nvertex, 3), numpy.float64)
        self.points[..., 0] = bounds_lon
        self.points[..., 1] = bounds_lat
        self.points = self.points.reshape((numCells, nvertex, 3))
        self.grid = mint.Grid()
        self.grid.setPoints(self.points)

    def getMintGrid(self):
        return self.grid

    def getNumCells(self):
        return self.grid.getNumberOfCells()

    def getPoints(self):
        return self.points

    def getPoint(self, cellId, vertex):
        return self.points[cellId, vertex, :]

    def dump(self, fileName):
        """Dump the grid data to a VTK file

        :param fileName: file name
        """
        self.grid.dump(fileName)

def main(*, tFile: str):
    """Create grid
    :param tFile: netCDF file containing t grid data
   """
    gr = HorizGrid(tFile)
    vtkFile = re.sub('.nc', '.vtk', tFile)
    gr.dump(vtkFile)

if __name__ == '__main__':
    defopt.run(main)


