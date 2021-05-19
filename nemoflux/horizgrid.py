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
        nc.close()

        ny, nx, nvertex = bounds_lat.shape
        numCells = ny * nx
        self.points = numpy.zeros((numCells, 4, 3), numpy.float64)

        k = 0
        for j in range(ny):
            for i in range(nx):
                for vertex in range(4):
                    x, y = bounds_lon[j, i, vertex], bounds_lat[j, i, vertex]
                    self.points[k, vertex, :] = x, y, 0.0
                k += 1

        self.grid = mint.Grid()
        self.grid.setPoints(self.points)
        self.grid.computeEdgeArcLengths()

        # compute the areas
        dx10 = self.points[:, 1, 0] - self.points[:, 0, 0]
        dy10 = self.points[:, 1, 1] - self.points[:, 0, 1]
        dx30 = self.points[:, 3, 0] - self.points[:, 0, 0]
        dy30 = self.points[:, 3, 1] - self.points[:, 0, 1]
        self.areas0 = dx10*dy30 - dx30*dy10

        dx32 = self.points[:, 3, 0] - self.points[:, 2, 0]
        dy32 = self.points[:, 3, 1] - self.points[:, 2, 1]
        dx12 = self.points[:, 1, 0] - self.points[:, 2, 0]
        dy12 = self.points[:, 1, 1] - self.points[:, 2, 1]
        self.areas2 = dx32*dy12 - dx12*dy32

        self.grid.attach('areas0', self.areas0)
        self.grid.attach('areas2', self.areas2)
        
    def getEdgeArcLength(self, cellId, edgeIndex):
    	return self.grid.getEdgeArcLength(cellId, edgeIndex)

    def getMintGrid(self):
        return self.grid

    def getNumCells(self):
        return self.points.shape[0]

    def getPoint(self, k, vertex):
        return self.points[k, vertex, :]


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


