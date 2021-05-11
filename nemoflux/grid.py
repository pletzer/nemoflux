import mint
import netCDF4
import numpy

class Grid(object):

    def __init__(self, cellDataFile, varName)
    """Constructor

    :param cellDataFile: netCDF file containing cell centred data
    :param varName: variable name of the cell centred data
    """
    
    nc = netCDF4.Dataset(cellDataFile)

    # get the variable
    var = nv.variables[varName]
    ny, nx = var.shape
    coordNames = var.coordinates.split()

    # gather the latitues and longitudes
    lonName = ''
    latName = ''
    for cn in coordNames:
        v = nc.variables[cn]
        if v.standard_name == 'longitude':
            lonName = cn
        elif v.standard_name == 'latitude':
            latName = cn

    vLon = nc.variables(lonName)
    vLat = nc.variables(latName)

    # read the cell bounds
    boundsLon = nc.variables[vLon.bounds][:]
    boundsLat = nc.variables[vLon.bounds][:]

    self.points = numpy.zeros((numCells, 4, 3), numpy.float64)
    nc.close()

    self.grid = mint.Grid()
    self.grid.setPoints(self.points)

    def dump(self, fileName):
        """Dump the grid data to a VTK file

        :param fileName: file name
        """
        self.grid.dump(fileName)



