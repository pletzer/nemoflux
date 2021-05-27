import vtk
import netCDF4
import geo
import defopt
import numpy
from horizgrid import HorizGrid
import mint

class FluxViz(object):

    def __init__(self, tFile, uFile, vFile, lonLatPoints):

        self.gr = HorizGrid(tFile)
        self.pli = mint.PolylineIntegral()
        self.pli.build(self.gr.getMintGrid(), lonLatPoints,
                       counterclock=False, periodX=360.0)


        # read the cell bounds
        with netCDF4.Dataset(tFile) as nc:
            bounds_lat = nc.variables['bounds_lat'][:]
            bounds_lon = nc.variables['bounds_lon'][:]
            depth_half = nc.variables['deptht'][:]
            self.bounds_depth = nc.variables['deptht_bounds'][:]

        self.ncU = netCDF4.Dataset(uFile)
        self.ncV = netCDF4.Dataset(vFile)
        uo, vo = self.getUV()

        ny, nx = self.ny, self.nx
        numCells = ny * nx
        self.integratedVelocity = numpy.zeros((numCells, 4), numpy.float64)
        self.computeIntegratedFlux(uo, vo)

        # points, 4 points per cell, 3D
        self.lonlat = numpy.zeros((ny, nx, 4, 3), numpy.float64)
        self.lonlat[..., 0] = bounds_lon
        self.lonlat[..., 1] = bounds_lat

        self.pointData = vtk.vtkDoubleArray()
        self.pointData.SetNumberOfComponents(3)
        self.pointData.SetNumberOfTuples(4 * numCells)
        self.pointData.SetVoidArray(self.lonlat, 4 * numCells * 3, 1)

        self.points = vtk.vtkPoints()
        self.points.SetData(self.pointData)

        # grid for the target points
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

        # 1 grid for the U fluxes, 1 grid for the V fluxes
        self.gridU = vtk.vtkPolyData()
        self.gridU.SetPoints(self.points)
        self.edgeFluxesUArray = numpy.zeros((numCells,), numpy.float64)
        self.edgeFluxesU = vtk.vtkDoubleArray()
        self.edgeFluxesU.SetName('U')
        self.edgeFluxesU.SetNumberOfComponents(1)
        self.edgeFluxesU.SetNumberOfTuples(numCells)
        self.edgeFluxesU.SetVoidArray(self.edgeFluxesUArray, numCells, 1)

        self.gridV = vtk.vtkPolyData()
        self.gridV.SetPoints(self.points)
        self.edgeFluxesVArray = numpy.zeros((numCells,), numpy.float64)
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

        #        ^
        #        |
        #  3-----V-----2
        #  |           |
        #  |     T     U-->
        #  |           |
        #  0-----------1
        self.edgeFluxesUArray[:] = self.integratedVelocity[:, 1]
        self.edgeFluxesVArray[:] = self.integratedVelocity[:, 2]

        self.minFlux = min(self.edgeFluxesUArray.min(), self.edgeFluxesVArray.min())
        self.maxFlux = min(self.edgeFluxesUArray.max(), self.edgeFluxesVArray.max())

        # assemble the cells
        cellId = 0
        for j in range(ny):
            for i in range(nx):

                ptIds.SetId(0, 4*cellId + 1)
                ptIds.SetId(1, 4*cellId + 2)
                self.gridU.InsertNextCell(vtk.VTK_LINE, ptIds)

                ptIds.SetId(0, 4*cellId + 3)
                ptIds.SetId(1, 4*cellId + 2)
                self.gridV.InsertNextCell(vtk.VTK_LINE, ptIds)

                # increment the cell counter
                cellId += 1

        print(f'min/max vertically integrated edge flux: {self.minFlux}/{self.maxFlux}')


    def getUV(self):

        # read u
        uo = self.ncU.variables['uo'][:]
        # set the velocity to zero where missing
        uo = numpy.ma.filled(uo, 0.0)

        # read v
        vo = self.ncV.variables['vo'][:]
        # set the velocity to zero where missing
        vo = numpy.ma.filled(vo, 0.0)

        nz, ny, nx = uo.shape
        numCells = ny * nx
        self.ny, self.nx = ny, nx

        # integrate vertically, multiplying by the thickness of the layers
        dz = -(self.bounds_depth[:, 1] - self.bounds_depth[:, 0]) # DEPTH HAS OPPOSITE SIGN TO Z
        uo = numpy.tensordot(dz, uo, axes=(0, 0)) # sum of multiplying axis 0 of dz with axis 0 of uo
        vo = numpy.tensordot(dz, vo, axes=(0, 0))

        return uo, vo


    def computeIntegratedFlux(self, uo, vo):

        ny, nx = self.ny, self.nx
        numCells = ny*nx

        integratedVelocity = self.integratedVelocity.reshape((ny, nx, 4))

        # now compute the integrated flux, cell by cell
        points = self.gr.getPoints().reshape((ny, nx, 4, 3))
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
        integratedVelocity[1:, :, 0] = vo[:-1, :] * ds01[:-1, :]
        # else:
        #     # folding, assuming even nx
        #     self.integratedVelocity[cellId, 0] = vo[:, j, (i+nx//2)%nx] * ds01

        # east
        integratedVelocity[..., 1] = uo * ds12

        # north
        integratedVelocity[..., 2] = vo * ds32

        # periodic west
        integratedVelocity[:, 1:, 3] = uo[:, :-1] * ds03[:, :-1]
        integratedVelocity[:, 0, 3] = uo[:, -1] * ds03[:, -1]

        # array should have shape (numCells, 4)
        self.integratedVelocity = integratedVelocity.reshape((numCells, 4))


    def show(self, npx=1260, npy=960):

        totalFlux = self.pli.getIntegral(self.integratedVelocity)
        title = vtk.vtkTextActor()
        title.SetTextScaleMode(0)
        title.GetTextProperty().SetFontSize(50)
        title.SetInput(f"total flux = {totalFlux:10.3f}")
        title.SetPosition((0.6, 0.9))

        lut = vtk.vtkLookupTable()
        lut.SetHueRange(0.6, 0.07)
        lut.SetTableRange(self.minFlux, self.maxFlux)
        lut.Build()

        cbar = vtk.vtkScalarBarActor()
        cbar.SetLookupTable(lut)

        radiusMin = 0.1*min(360/self.nx, 180/self.ny)

        tubesU = vtk.vtkTubeFilter()
        tubesU.SetInputData(self.gridU)
        tubesU.SetRadius(radiusMin) # min radius
        tubesU.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
        tubesU.SetRadiusFactor(100)

        tubesV = vtk.vtkTubeFilter()
        tubesV.SetInputData(self.gridV)
        tubesV.SetRadius(radiusMin) # min radius
        tubesV.SetRadiusFactor(100)
        tubesV.SetVaryRadiusToVaryRadiusByAbsoluteScalar()

        mapperU = vtk.vtkPolyDataMapper()
        mapperU.SetInputConnection(tubesU.GetOutputPort())
        mapperU.SetLookupTable(lut)
        mapperU.SetUseLookupTableScalarRange(1)
        actorU = vtk.vtkActor()
        actorU.SetMapper(mapperU)

        mapperV = vtk.vtkPolyDataMapper()
        mapperV.SetInputConnection(tubesV.GetOutputPort())
        mapperV.SetLookupTable(lut)
        mapperV.SetUseLookupTableScalarRange(1)
        actorV = vtk.vtkActor()
        actorV.SetMapper(mapperV)

        tubePoints = vtk.vtkTubeFilter()
        tubePoints.SetRadius(0.2*min(360/self.nx, 180/self.ny))
        tubePoints.SetInputData(self.gridTargetLine)
        mapperPoints = vtk.vtkPolyDataMapper()
        mapperPoints.SetInputConnection(tubePoints.GetOutputPort())
        actorPoints = vtk.vtkActor()
        actorPoints.SetMapper(mapperPoints)

        # Create the graphics structure. The renderer renders into the render
        # window. The render window interactor captures mouse events and will
        # perform appropriate camera or actor manipulation depending on the
        # nature of the events.
        ren = vtk.vtkRenderer()
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)

        # Add the actors to the renderer, set the background and size
        ren.AddActor(actorU)
        ren.AddActor(actorV)
        ren.AddActor(actorPoints)
        ren.AddActor(cbar)
        ren.AddActor(title)
        ren.SetBackground((0.1, 0.1, 0.1))
        renWin.SetSize(npx, npy)
        renWin.SetWindowName('Vertically integrated edge flux')

        # This allows the interactor to initalize itself. It has to be
        # called before an event loop.
        iren.Initialize()

        # We'll zoom in a little by accessing the camera and invoking a "Zoom"
        # method on it.
        # ren.ResetCamera()
        # ren.GetActiveCamera().Zoom(1.5)
        renWin.Render()

        # Start the event loop.
        iren.Start()

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
    