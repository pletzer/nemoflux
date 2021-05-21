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
        self.points = numpy.zeros((numCells, 4, 3), numpy.float64)
        self.cellAreas = numpy.zeros((numCells,), numpy.float64)

        cellId = 0
        for j in range(ny):
            for i in range(nx):
                for vertex0 in range(4):
                    vertex1 = (vertex0 + 1) % 4
                    vertex2 = (vertex0 - 1) % 4
                    x0, y0 = bounds_lon[j, i, vertex0], bounds_lat[j, i, vertex0]
                    x1, y1 = bounds_lon[j, i, vertex1], bounds_lat[j, i, vertex1]
                    x2, y2 = bounds_lon[j, i, vertex2], bounds_lat[j, i, vertex2]
                    self.points[cellId, vertex0, :] = x0, y0, 0.0
                    dp10 = numpy.array((x1 - x0, y1 - y0, 0.))
                    dp20 = numpy.array((x2 - x0, y2 - y0, 0.))
                    area = numpy.cross(dp10, dp20)
                    self.cellAreas[cellId] += 0.25*area[2]
                # increment the cell counter
                cellId += 1

        self.grid = mint.Grid()
        self.grid.setPoints(self.points)
        self.grid.computeEdgeArcLengths()
        self.grid.attach('cellAreas', self.cellAreas)

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


