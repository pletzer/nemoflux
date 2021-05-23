import vtk
import netCDF4
import geo
import defopt
import numpy

class FluxViz(object):

    def __init__(self, tFile, uFile, vFile):

        # read the cell bounds
        nc = netCDF4.Dataset(tFile)
        bounds_lat = nc.variables['bounds_lat'][:]
        bounds_lon = nc.variables['bounds_lon'][:]
        depth_half = nc.variables['deptht'][:]
        bounds_depth = nc.variables['deptht_bounds'][:]
        nc.close()

        # read u
        nc = netCDF4.Dataset(uFile)
        uo = nc.variables['uo'][:]
        nc.close()
        # set the velocity to zero where missing
        uo = numpy.ma.filled(uo, 0.0)

        # read v
        nc = netCDF4.Dataset(vFile)
        vo = nc.variables['vo'][:]
        nc.close()
        # set the velocity to zero where missing
        vo = numpy.ma.filled(vo, 0.0)

        nz, ny, nx = uo.shape
        numCells = ny * nx

        # points, 4 points per cell, 3D
        self.xyz = numpy.zeros((ny, nx, 4, 3), numpy.float64)
        self.xyz[..., 0] = bounds_lon
        self.xyz[..., 1] = bounds_lat

        self.vPointData = vtk.vtkDoubleArray()
        self.vPointData.SetNumberOfComponents(3)
        self.vPointData.SetNumberOfTuples(4 * numCells)
        self.vPointData.SetVoidArray(self.xyz, 4 * numCells * 3, 1)

        self.vPoints = vtk.vtkPoints()
        self.vPoints.SetData(self.vPointData)

        # 1 grid for the U fluxes, 1 grid for the V fluxes
        self.vGridU = vtk.vtkUnstructuredGrid()
        self.vGridU.SetPoints(self.vPoints)
        self.edgeFluxesU = vtk.vtkDoubleArray()
        self.edgeFluxesU.SetName('verticallyIntegratedUFlux')
        self.edgeFluxesU.SetNumberOfComponents(1)
        self.edgeFluxesU.SetNumberOfTuples(numCells)
        self.edgeFluxesU.Fill(0.)

        self.vGridV = vtk.vtkUnstructuredGrid()
        self.vGridV.SetPoints(self.vPoints)
        self.edgeFluxesV = vtk.vtkDoubleArray()
        self.edgeFluxesU.SetName('verticallyIntegratedVFlux')
        self.edgeFluxesV.SetNumberOfComponents(1)
        self.edgeFluxesV.SetNumberOfTuples(numCells)
        self.edgeFluxesV.Fill(0.)

        self.vGridU.Allocate()
        self.vGridV.Allocate()

        self.vGridU.GetCellData().SetScalars(self.edgeFluxesU)
        self.vGridV.GetCellData().SetScalars(self.edgeFluxesV)

        ptIds = vtk.vtkIdList()
        ptIds.SetNumberOfIds(2)

        flxU = numpy.zeros((1,), numpy.float64)
        flxV = numpy.zeros((1,), numpy.float64)
        p1 = numpy.zeros((3,), numpy.float64)
        p2 = numpy.zeros((3,), numpy.float64)
        p3 = numpy.zeros((3,), numpy.float64)

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
                #  |           U-->
                #  |           |
                #  0-----------1


                ptIds.SetId(0, 4*cellId + 1); ptIds.SetId(1, 4*cellId + 2)
                self.vGridU.InsertNextCell(vtk.VTK_LINE, ptIds)

                ptIds.SetId(0, 4*cellId + 3); ptIds.SetId(1, 4*cellId + 2)
                self.vGridV.InsertNextCell(vtk.VTK_LINE, ptIds)

                p1[:] = self.xyz[j, i, 1, :]
                p2[:] = self.xyz[j, i, 2, :]
                p3[:] = self.xyz[j, i, 3, :]

                # integrate vertically
                for k in range(nz):

                    dz = -(bounds_depth[k, 1] - bounds_depth[k, 0]) # DEPTH HAS OPPOSITE SIGN TO Z

                    self.edgeFluxesU.GetTuple(cellId, flxU)
                    flxU[0] += uo[k, j, i] * geo.getArcLength(p1, p2, radius=geo.EARTH_RADIUS) * dz
                    self.edgeFluxesU.SetTuple(cellId, flxU)

                    self.edgeFluxesV.GetTuple(cellId, flxV)
                    flxV[0] += vo[k, j, i] * geo.getArcLength(p3, p2, radius=geo.EARTH_RADIUS) * dz
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
        lut.SetHueRange(0.666, 0.0)
        lut.SetTableRange(self.minFlux, self.maxFlux)
        lut.Build()

        cbar = vtk.vtkScalarBarActor()
        cbar.SetLookupTable(lut)

        mapperU = vtk.vtkDataSetMapper()
        mapperU.SetInputData(self.vGridU)
        mapperU.SetLookupTable(lut)
        mapperU.SetUseLookupTableScalarRange(1)
        actorU = vtk.vtkActor()
        actorU.SetMapper(mapperU)

        mapperV = vtk.vtkDataSetMapper()
        mapperV.SetInputData(self.vGridV)
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
        ren.SetBackground((0.7, 0.7, 0.7))
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
    