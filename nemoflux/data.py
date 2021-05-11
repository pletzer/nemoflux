import numpy
import netCDF4

class Data(object):

    def __init__(self, uDataFile, vDataFile, uVarName, vVarName)
    """Constructor

    :param uDataFile: netCDF file containing U component data on Arakawa C grid
    :param vDataFile: netCDF file containing V component data on Arakawa C grid
    :param uVarName: name of the U component data in uDataFile
    :param vVarName: name of the V component data in vDataFile
    """
    
    ncU = netCDF4.Dataset(uDataFile)
    ncV = netCDF4.Dataset(vDataFile)

    vU = ncU.variables[uVarName]
    vV = ncV.variables[vVarName]

    # read
    self.u = vU[:]
    self.v = vV[:]

    def getFaceIntegratedData(self, fileName):
        """Dump the grid data to a VTK file

        :param fileName: file name
        """
        self.grid.dump(fileName)



