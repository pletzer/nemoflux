import vtk
import netCDF4
import geo
import defopt
import numpy
from horizgrid import HorizGrid
import mint

class CallBack(object):

    def __init__(self, fluxviz):
        self.fluxviz = fluxviz

    def execute(self, obj, event):
        key = obj.GetKeySym()
        self.fluxviz.update(key)
        self.fluxviz.renWin.Render()


class FluxViz(object):

    def __init__(self, tFile, uFile, vFile, lonLatPoints):

        self.timeIndex = 0
        self.ncU = netCDF4.Dataset(uFile)
        self.ncV = netCDF4.Dataset(vFile)

        self.nt, self.nz, self.ny, self.nx = self.getSizes()
        numCells = self.ny * self.nx

        self.gr = HorizGrid(tFile)
        self.pli = mint.PolylineIntegral()
        self.pli.build(self.gr.getMintGrid(), lonLatPoints,
                       counterclock=False, periodX=360.0)

        # read the cell bounds
        with netCDF4.Dataset(tFile) as nc:
            bounds_lat = nc.variables['bounds_lat'][:]
            bounds_lon = nc.variables['bounds_lon'][:]
            self.bounds_depth = nc.variables['deptht_bounds'][:]

        self.thickness = -(self.bounds_depth[:, 1] - self.bounds_depth[:, 0]) # DEPTH HAS OPPOSITE SIGN TO Z

        uVerticallyIntegrated, vVerticallyIntegrated = self.getUV()

        self.edgeFluxesUArray = numpy.zeros((numCells,), numpy.float64)
        self.edgeFluxesVArray = numpy.zeros((numCells,), numpy.float64)

        numCells = self.ny * self.nx
        self.integratedVelocity = numpy.zeros((numCells, 4), numpy.float64)
        self.minFlux, self.maxFlux = +float('inf'), -float('inf')

        self.computeIntegratedFlux(uVerticallyIntegrated, vVerticallyIntegrated)

        self.buildTargetLineGrid(lonLatPoints)
        self.buildEdgeUVGrids(bounds_lon, bounds_lat)
        print(f'min/max vertically integrated edge flux: {self.minFlux}/{self.maxFlux}')


    def update(self, key):
        if key == 't':
            # forward in time
            self.timeIndex = (self.timeIndex + 1) % self.nt
        elif key == 'b':
            # backward in time
            self.timeIndex = (self.timeIndex - 1) % self.nt
        # update the data
        uo, vo = self.getUV()
        self.computeIntegratedFlux(uo, vo)
        self.edgeFluxesUArray[:] = self.integratedVelocity[:, 1]
        self.edgeFluxesVArray[:] = self.integratedVelocity[:, 2]
        self.lut.SetTableRange(self.minFlux, self.maxFlux)
        self.lut.Modified()
        self.cbar.Modified()
        totalFlux = self.pli.getIntegral(self.integratedVelocity)
        self.title.SetInput(f'total flux: {totalFlux:10.3f} @ time {self.timeIndex}')
        self.title.Modified()
        self.edgeFluxesU.Modified()
        self.edgeFluxesV.Modified()
        print(f'time index is now {self.timeIndex} min/max flux: {self.minFlux:10.3f}/{self.maxFlux:10.3f} nt = {self.nt}')


    def getSizes(self):
        # get the sizes
        nt, nz, ny, nx = 1, 1, 0, 0
        shapeU = self.ncU.variables['uo'].shape
        try:
            nt, nz, ny, nx = shapeU
        except:
            try:
                nz, ny, nx = shapeU
            except:
                try:
                    ny, nx = shapeU
                except:
                    raise RuntimeError("ERROR: uo's shape does not match (t, z, y, x), (z, y, x) or (y, x)")
        return nt, nz, ny, nx


    def buildTargetLineGrid(self, lonLatPoints):

        ptIds = vtk.vtkIdList()
        ptIds.SetNumberOfIds(2)

        # build the target point mesh
        self.lonLatPoints = lonLatPoints
        numTargetPoints = lonLatPoints.shape[0]
        self.targetPointData = vtk.vtkDoubleArray()
        self.targetPointData.SetNumberOfComponents(3)
        self.targetPointData.SetNumberOfTuples(numTargetPoints)
        self.targetPointData.SetVoidArray(self.lonLatPoints, numTargetPoints*3, 1)
        self.targetPoints = vtk.vtkPoints()
        self.targetPoints.SetData(self.targetPointData)
        self.gridTargetLine = vtk.vtkPolyData()
        self.gridTargetLine.SetPoints(self.targetPoints)
        self.gridTargetLine.Allocate()
        assert(numTargetPoints > 1)
        for i in range(numTargetPoints - 1):
            ptIds.SetId(0, i)
            ptIds.SetId(1, i + 1)
            self.gridTargetLine.InsertNextCell(vtk.VTK_LINE, ptIds)


    def buildEdgeUVGrids(self, bounds_lon, bounds_lat):

        numCells = self.ny * self.nx

        # points, 4 points per cell, 3D
        self.lonlat = numpy.zeros((self.ny, self.nx, 4, 3), numpy.float64)
        self.lonlat[..., 0] = bounds_lon
        self.lonlat[..., 1] = bounds_lat

        self.pointData = vtk.vtkDoubleArray()
        self.pointData.SetNumberOfComponents(3)
        self.pointData.SetNumberOfTuples(4 * numCells)
        self.pointData.SetVoidArray(self.lonlat, 4 * numCells * 3, 1)

        self.points = vtk.vtkPoints()
        self.points.SetData(self.pointData)


        # 1 grid for the U fluxes, 1 grid for the V fluxes
        self.gridU = vtk.vtkPolyData()
        self.gridU.SetPoints(self.points)
        self.edgeFluxesU = vtk.vtkDoubleArray()
        self.edgeFluxesU.SetName('U')
        self.edgeFluxesU.SetNumberOfComponents(1)
        self.edgeFluxesU.SetNumberOfTuples(numCells)
        self.edgeFluxesU.SetVoidArray(self.edgeFluxesUArray, numCells, 1)

        self.gridV = vtk.vtkPolyData()
        self.gridV.SetPoints(self.points)
        self.edgeFluxesV = vtk.vtkDoubleArray()
        self.edgeFluxesV.SetName('V')
        self.edgeFluxesV.SetNumberOfComponents(1)
        self.edgeFluxesV.SetNumberOfTuples(numCells)
        self.edgeFluxesV.SetVoidArray(self.edgeFluxesVArray, numCells, 1)

        self.gridU.Allocate()
        self.gridV.Allocate()

        self.gridU.GetCellData().SetScalars(self.edgeFluxesU)
        self.gridU.GetCellData().SetActiveScalars('U')
        self.gridV.GetCellData().SetScalars(self.edgeFluxesV)
        self.gridV.GetCellData().SetActiveScalars('V')

        ptIds = vtk.vtkIdList()
        ptIds.SetNumberOfIds(2)

        # assemble the cells
        cellId = 0
        for j in range(self.ny):
            for i in range(self.nx):

                ptIds.SetId(0, 4*cellId + 1)
                ptIds.SetId(1, 4*cellId + 2)
                self.gridU.InsertNextCell(vtk.VTK_LINE, ptIds)

                ptIds.SetId(0, 4*cellId + 3)
                ptIds.SetId(1, 4*cellId + 2)
                self.gridV.InsertNextCell(vtk.VTK_LINE, ptIds)

                # increment the cell counter
                cellId += 1


    def getUV(self):

        # read u
        try:
            uo = self.ncU.variables['uo'][self.timeIndex, :, :, :]
        except:
            try:
                uo = self.ncU.variables['uo'][...]
            except:
                raise RuntimeError('ERROR: could not read u field')

        # set the velocity to zero where missing
        uo = numpy.ma.filled(uo, 0.0)

        # read v
        try:
            vo = self.ncV.variables['vo'][self.timeIndex, :, :, :]
        except:
            try:
                vo = self.ncV.variables['vo'][...]
            except:
                raise RuntimeError('ERROR: could not read v field')

        # set the velocity to zero where missing
        vo = numpy.ma.filled(vo, 0.0)

        # integrate vertically, multiplying by the thickness of the layers
        uVerticallyIntegrated = numpy.tensordot(self.thickness, uo, axes=(0, 0)) # sum of multiplying axis 0 of dz with axis 0 of uo
        vVerticallyIntegrated = numpy.tensordot(self.thickness, vo, axes=(0, 0))

        return uVerticallyIntegrated, vVerticallyIntegrated


    def computeIntegratedFlux(self, uVerticallyIntegrated, vVerticallyIntegrated):

        numCells = self.ny * self.nx

        integratedVelocity = self.integratedVelocity.reshape((self.ny, self.nx, 4))

        # now compute the integrated flux, cell by cell
        points = self.gr.getPoints().reshape((self.ny, self.nx, 4, 3))
        xyz0 = geo.lonLat2XYZArray(points[:, :, 0, :], radius=geo.EARTH_RADIUS)
        xyz1 = geo.lonLat2XYZArray(points[:, :, 1, :], radius=geo.EARTH_RADIUS)
        xyz2 = geo.lonLat2XYZArray(points[:, :, 2, :], radius=geo.EARTH_RADIUS)
        xyz3 = geo.lonLat2XYZArray(points[:, :, 3, :], radius=geo.EARTH_RADIUS)

        ds01 = geo.getArcLengthArray(xyz0, xyz1, radius=geo.EARTH_RADIUS)
        ds12 = geo.getArcLengthArray(xyz1, xyz2, radius=geo.EARTH_RADIUS)
        ds32 = geo.getArcLengthArray(xyz3, xyz2, radius=geo.EARTH_RADIUS)
        ds03 = geo.getArcLengthArray(xyz0, xyz3, radius=geo.EARTH_RADIUS)

        #        ^
        #        |
        #  3-----V-----2
        #  |           |
        #  |->   T     U->
        #  |     ^     |
        #  |     |     |
        #  0-----------1

        # south
        integratedVelocity[1:, :, 0] = vVerticallyIntegrated[:-1, :] * ds01[:-1, :]
        # else:
        #     # folding, assuming even nx
        #     self.integratedVelocity[cellId, 0] = vo[:, j, (i+nx//2)%nx] * ds01

        # east
        integratedVelocity[..., 1] = uVerticallyIntegrated * ds12

        # north
        integratedVelocity[..., 2] = vVerticallyIntegrated * ds32

        # periodic west
        integratedVelocity[:, 1:, 3] = uVerticallyIntegrated[:, :-1] * ds03[:, :-1]
        integratedVelocity[:, 0, 3] = uVerticallyIntegrated[:, -1] * ds03[:, -1]

        # array should have shape (numCells, 4)
        self.integratedVelocity = integratedVelocity.reshape((numCells, 4))

        #        ^
        #        |
        #  3-----V-----2
        #  |           |
        #  |     T     U-->
        #  |           |
        #  0-----------1
        self.edgeFluxesUArray[:] = self.integratedVelocity[:, 1]
        self.edgeFluxesVArray[:] = self.integratedVelocity[:, 2]

        self.minFlux = min(self.minFlux, self.edgeFluxesUArray.min(), self.edgeFluxesVArray.min())
        self.maxFlux = max(self.maxFlux, self.edgeFluxesUArray.max(), self.edgeFluxesVArray.max())


    def show(self, npx=1260, npy=960):

        print(f'integrated velocity on edge 0: {self.integratedVelocity[...,0]}')
        print(f'integrated velocity on edge 2: {self.integratedVelocity[...,2]}')
        totalFlux = self.pli.getIntegral(self.integratedVelocity)
        self.title = vtk.vtkTextActor()
        self.title.SetTextScaleMode(0)
        self.title.GetTextProperty().SetFontSize(50)
        self.title.SetInput(f"total flux = {totalFlux:15.8g} @ time {self.timeIndex}")
        self.title.SetPosition((0.6, 0.9))

        self.lut = vtk.vtkLookupTable()
        self.lut.SetHueRange(0.6, 0.07)
        self.lut.SetTableRange(self.minFlux, self.maxFlux)
        self.lut.Build()

        self.cbar = vtk.vtkScalarBarActor()
        self.cbar.SetLookupTable(self.lut)

        radiusMin = 0.05*min(360/self.nx, 180/self.ny)

        self.mapperU = vtk.vtkPolyDataMapper()
        self.mapperU.SetInputData(self.gridU)
        self.mapperU.SetLookupTable(self.lut)
        self.mapperU.SetUseLookupTableScalarRange(1)
        self.actorU = vtk.vtkActor()
        self.actorU.SetMapper(self.mapperU)

        self.mapperV = vtk.vtkPolyDataMapper()
        self.mapperV.SetInputData(self.gridV)
        self.mapperV.SetLookupTable(self.lut)
        self.mapperV.SetUseLookupTableScalarRange(1)
        self.actorV = vtk.vtkActor()
        self.actorV.SetMapper(self.mapperV)

        self.tubePoints = vtk.vtkTubeFilter()
        self.tubePoints.SetRadius(2.0)
        self.tubePoints.SetInputData(self.gridTargetLine)
        self.mapperPoints = vtk.vtkPolyDataMapper()
        self.mapperPoints.SetInputConnection(self.tubePoints.GetOutputPort())
        self.actorPoints = vtk.vtkActor()
        self.actorPoints.SetMapper(self.mapperPoints)

        # Create the graphics structure. The renderer renders into the render
        # window. The render window interactor captures mouse events and will
        # perform appropriate camera or actor manipulation depending on the
        # nature of the events.
        self.ren = vtk.vtkRenderer()
        self.renWin = vtk.vtkRenderWindow()
        self.renWin.AddRenderer(self.ren)
        self.iren = vtk.vtkRenderWindowInteractor()
        self.iren.SetRenderWindow(self.renWin)

        # Add the actors to the renderer, set the background and size
        self.ren.AddActor(self.actorU)
        self.ren.AddActor(self.actorV)
        self.ren.AddActor(self.actorPoints)
        self.ren.AddActor(self.cbar)
        self.ren.AddActor(self.title)
        self.ren.SetBackground((0.1, 0.1, 0.1))
        self.renWin.SetSize(npx, npy)
        self.renWin.SetWindowName('Vertically integrated edge flux')

        # allow the user to interact with the visualisation
        self.callBack = CallBack(self)
        self.iren.AddObserver('KeyPressEvent', self.callBack.execute)
        print('type "t"/"b" to step forward/backward in time')

        # This allows the interactor to initalize itself. It has to be
        # called before an event loop.
        self.iren.Initialize()

        self.renWin.Render()

        # Start the event loop.
        self.iren.Start()

def __del__(self):
    self.ncU.close()
    self.ncV.close()


def main(*, tFile: str, uFile: str, vFile: str, lonLatPointsStr: str):
    """Visualize fluxes
    :param tFile: netcdf file holding the T-grid
    :param uFile: netcdf file holding u data
    :param vFile: netcdf file holding v data
    :param lonLatPointsStr: target points "(lon0, lat0), (lon1, lat1),..."
    """
    xyVals = numpy.array(eval(lonLatPointsStr))
    numTargetPoints = xyVals.shape[0]
    lonLatZPoints = numpy.zeros((numTargetPoints, 3), numpy.float64)
    lonLatZPoints[:, :2] = xyVals
    fv = FluxViz(tFile, uFile, vFile, lonLatZPoints)
    fv.show()


if __name__ == '__main__':
    defopt.run(main)
    