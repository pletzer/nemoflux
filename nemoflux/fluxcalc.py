import mint
import numpy

class FluxCalc(object):

    def __init__(self, grid, targetPoints, integratedVelocity):
        """Constructor
    
        :param grid: instance of mintGrid
        :param targetPoints: array of target points
        """

        self.pli = mint.PolylineIntegral()
        self.pli.build(grid, targetPoints counterclock=False, periodX=360.0)

    def getFlux(self, integratedVelocity):
        """Compute the flux

        :param integratedVelocity: velocity data integrated over the cell faces
        :returns: flux value
        """        
        return self.pli.getIntegral(integratedVelocity)