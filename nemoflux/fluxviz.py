import vtk
import netCDF4
import geo
import defopt
import numpy

class FluxViz(object):

    def __init__(self, tFile, uFile, vFile):

        # read the cell bounds
        with netCDF4.Dataset(tFile) as nc:
            bounds_lat = nc.variables['bounds_lat'][:]
            bounds_lon = nc.variables['bounds_lon'][:]
            depth_half = nc.variables['deptht'][:]
            bounds_depth = nc.variables['deptht_bounds'][:]

        # read u
        with netCDF4.Dataset(uFile) as nc:
            uo = nc.variables['uo'][:]
        # set the velocity to zero where missing
        uo = numpy.ma.filled(uo, 0.0)

        # read v
        with netCDF4.Dataset(vFile) as nc:
            vo = nc.variables['vo'][:]
        # set the velocity to zero where missing
        vo = numpy.ma.filled(vo, 0.0)

        nz, ny, nx = uo.shape
        numCells = ny * nx
        self.ny, self.nx = ny, nx

        # integrate vertically, multiplying by the thickness of the layers
        dz = -(bounds_depth[:, 1] - bounds_depth[:, 0]) # DEPTH HAS OPPOSITE SIGN TO Z
        uo = numpy.tensordot(dz, uo, axes=(0, 0)) # sum of multiplying axis 0 of dz with axis 0 of uo
        vo = numpy.tensordot(dz, vo, axes=(0, 0))

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

        # 1 grid for the U fluxes, 1 grid for the V fluxes
        self.gridU = vtk.vtkPolyData()
        self.gridU.SetPoints(self.points)
        self.edgeFluxesU = vtk.vtkDoubleArray()
        self.edgeFluxesU.SetName('U')
        self.edgeFluxesU.SetNumberOfComponents(1)
        self.edgeFluxesU.SetNumberOfTuples(numCells)
        self.edgeFluxesU.Fill(0.)

        self.gridV = vtk.vtkPolyData()
        self.gridV.SetPoints(self.points)
        self.edgeFluxesV = vtk.vtkDoubleArray()
        self.edgeFluxesU.SetName('V')
        self.edgeFluxesV.SetNumberOfComponents(1)
        self.edgeFluxesV.SetNumberOfTuples(numCells)
        self.edgeFluxesV.Fill(0.)

        self.gridU.Allocate()
        self.gridV.Allocate()

        self.gridU.GetCellData().SetScalars(self.edgeFluxesU)
        #self.gridU.GetCellData().SetActiveScalars('U')
        self.gridV.GetCellData().SetScalars(self.edgeFluxesV)
        #self.gridV.GetCellData().SetActiveScalars('verticallyIntegratedVFlux')

        ptIds = vtk.vtkIdList()
        ptIds.SetNumberOfIds(2)

        flxU = numpy.zeros((1,), numpy.float64)
        flxV = numpy.zeros((1,), numpy.float64)
        xyz1 = numpy.zeros((3,), numpy.float64)
        xyz2 = numpy.zeros((3,), numpy.float64)
        xyz3 = numpy.zeros((3,), numpy.float64)

        self.minFlux = float('inf')
        self.maxFlux = -float('inf')

        # compute the vertically integrated lateral fluxes on each horizontal edge
        cellId = 0
        for j in range(ny):
            for i in range(nx):

                #        ^
                #        |
                #  3-----V-----2
                #  |           |
                #  |     T     U-->
                #  |           |
                #  0-----------1

                ptIds.SetId(0, 4*cellId + 1)
                ptIds.SetId(1, 4*cellId + 2)
                self.gridU.InsertNextCell(vtk.VTK_LINE, ptIds)

                ptIds.SetId(0, 4*cellId + 3)
                ptIds.SetId(1, 4*cellId + 2)
                self.gridV.InsertNextCell(vtk.VTK_LINE, ptIds)

                xyz1[:] = geo.lonLat2XYZ(self.lonlat[j, i, 1, :], radius=geo.EARTH_RADIUS)
                xyz2[:] = geo.lonLat2XYZ(self.lonlat[j, i, 2, :], radius=geo.EARTH_RADIUS)
                xyz3[:] = geo.lonLat2XYZ(self.lonlat[j, i, 3, :], radius=geo.EARTH_RADIUS)

                flxU[0] = uo[j, i] * geo.getArcLength(xyz1, xyz2, radius=geo.EARTH_RADIUS)
                self.edgeFluxesU.SetTuple(cellId, flxU)

                flxV[0] = vo[j, i] * geo.getArcLength(xyz3, xyz2, radius=geo.EARTH_RADIUS)
                self.edgeFluxesV.SetTuple(cellId, flxV)

                self.minFlux = min(self.minFlux, flxU[0], flxV[0])
                self.maxFlux = max(self.minFlux, flxU[0], flxV[0])

                # increment the cell counter
                cellId += 1

        print(f'min/max vertically integrated edge flux: {self.minFlux}/{self.maxFlux}')

    def show(self):

        xAxis = vtk.vtkAxisActor()
        xAxis.SetPoint1((0., -90., 0.))
        xAxis.SetPoint2((360., -90., 0.))
        #xAxis.SetTitle('longitude deg.')
        xAxis.SetDeltaRangeMajor(10.)

        yAxis = vtk.vtkAxisActor()
        yAxis.SetPoint1((0., -90., 0.))
        yAxis.SetPoint2((0., 90., 0.))
        #yAxis.SetTitle('latitude deg.')
        yAxis.SetDeltaRangeMajor(10.)

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
        ren.AddActor(cbar)
        ren.AddActor(xAxis)
        ren.AddActor(yAxis)
        ren.SetBackground((0.1, 0.1, 0.1))
        renWin.SetSize(1260, 960)
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

def main(*, tFile: str, uFile: str, vFile: str):
    """Visualize fluxes
    :param tFile: netcdf file holding the T-grid
    :param uFile: netcdf file holding u data
    :param vFile: netcdf file holding v data
    """
    fv = FluxViz(tFile, uFile, vFile)
    fv.show()


if __name__ == '__main__':
    defopt.run(main)
    