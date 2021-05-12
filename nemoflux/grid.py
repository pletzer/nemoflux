import mint
import netCDF4
import numpy
import defopt
import re

class Grid(object):

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

    def dump(self, fileName):
        """Dump the grid data to a VTK file

        :param fileName: file name
        """
        self.grid.dump(fileName)

def main(*, tFile: str):
    """Create grid
    :param tFile: netCDF file containing t grid data
   """
    gr = Grid(tFile)
    vtkFile = re.sub('.nc', '.vtk', tFile)
    gr.dump(vtkFile)

if __name__ == '__main__':
    defopt.run(main)


