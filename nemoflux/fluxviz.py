import vtk
import netCDF4
import geo
import defopt
import numpy
from horizgrid import HorizGrid
import mint
# to read the target points from file
from latlonreader import LatLonReader
from field import Field


# callback class called when the user interacts with the visualization
class CallBack(object):

    def __init__(self, fluxviz):
        self.fluxviz = fluxviz

    def execute(self, obj, event):
        key = obj.GetKeySym()
        self.fluxviz.update(key)
        self.fluxviz.renWin.Render()


class FluxViz(object):

    def __init__(self, tFile, uFile, vFile, lonLatPoints, sverdrup=False):

        self.field = Field(tFile, uFile, vFile, lonLatPoints, sverdrup)
        self.buildTargetLineGrid(lonLatPoints)
        self.buildEdgeUVGrids()

    def update(self, key):

        camera = self.ren.GetActiveCamera()
        x, y, z = camera.GetPosition()

        if key == 't':
            # forward in time
            self.field.timeIndex = (self.field.timeIndex + 1) % self.field.nt
        elif key == 'T':
            # backward in time
            self.field.timeIndex = (self.field.timeIndex - 1) % self.field.nt
        else:
            if key == 'r':
                # reset 
                lon, lat = 0.5*(self.field.lonmin + self.field.lonmax), \
                           0.5*(self.field.latmin + self.field.latmax)
                camera.SetFocalPoint((lon, lat, 0.))
                camera.SetPosition((lon, lat, 2*self.field.dx*max(self.field.ny, self.field.nx)))
                camera.SetViewUp((0., 1., 0.))
            elif  key == 'x':
                camera.SetPosition(x + self.field.dx, y, z)
                camera.SetFocalPoint(x + self.field.dx, y, 0.)
            elif key == 'X':
                camera.SetPosition(x - self.field.dx, y, z)
                camera.SetFocalPoint(x - self.field.dx, y, 0.)                
            elif key == 'y':
                camera.SetPosition(x, y + self.field.dx, z)
                camera.SetFocalPoint(x, y + self.field.dx, 0.)
            elif key == 'Y':
                camera.SetPosition(x, y - self.field.dx, z)
                camera.SetFocalPoint(x, y - self.field.dx, 0.)
            elif key == 'z':
                camera.SetPosition(x, y, z + self.field.dx)
            elif key == 'Z':
                camera.SetPosition(x, y, z - self.field.dx)
            elif key == 'v':
                self.targetGlyphs.SetScaleFactor(1.3*self.targetGlyphs.GetScaleFactor())
            elif key == 'V':
                self.targetGlyphs.SetScaleFactor(self.targetGlyphs.GetScaleFactor()/1.3)
            elif key == 's':
                # screenshot
                w2f = vtk.vtkWindowToImageFilter()
                w2f.Modified()
                w2f.SetInput(self.renWin)
                w2f.Update()
                wr = vtk.vtkPNMWriter()
                wr.SetInputConnection(w2f.GetOutputPort())
                filename = f'fluxviz_{self.field.timeIndex:04d}.pnm'
                wr.SetFileName(filename)
                wr.Write()
                print(f'saved screen shot for time {self.field.timeIndex} in file {filename}')
            elif key == 'm':
                # save all frames
                for i in range(self.field.nt):
                    self.update(key='t')
                    self.update(key='s')
            camera.Modified()
            return

        # update the data
        self.field.update()

        # update the pipeline
        self.lut.SetTableRange(0., self.field.maxAbsFlux)
        self.lut.Modified()
        self.cbar.Modified()
        totalFlux = self.field.pli.getIntegral(self.field.integratedVelocity)
        if self.field.sverdrup:
            self.title.SetInput(f'flux = {totalFlux:6.3g} (Sv) @ time {self.field.timeIndex}')
        else:
            self.title.SetInput(f'flux = {totalFlux:6.3g} (A m^2/s) @ time {self.field.timeIndex}')
        self.title.Modified()
        self.edgeFluxesU.Modified()
        self.edgeFluxesV.Modified()
        self.targetVectorValues.Modified()
        print(f'time index now {self.field.timeIndex} max |flux|: {self.field.maxAbsFlux:10.3f} nt = {self.field.nt}')

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

    def buildEdgeUVGrids(self):

        numCells = self.field.ny * self.field.nx

        self.pointData = vtk.vtkDoubleArray()
        self.pointData.SetNumberOfComponents(3)
        self.pointData.SetNumberOfTuples(4 * numCells)
        self.pointData.SetVoidArray(self.field.lonlat, 4 * numCells * 3, 1)

        self.points = vtk.vtkPoints()
        self.points.SetData(self.pointData)

        # 1 grid for the U fluxes, 1 grid for the V fluxes
        self.gridU = vtk.vtkPolyData()
        self.gridU.SetPoints(self.points)
        self.edgeFluxesU = vtk.vtkDoubleArray()
        self.edgeFluxesU.SetName('U')
        self.edgeFluxesU.SetNumberOfComponents(1)
        self.edgeFluxesU.SetNumberOfTuples(numCells)
        self.edgeFluxesU.SetVoidArray(self.field.edgeFluxesUArray, numCells, 1)

        self.gridV = vtk.vtkPolyData()
        self.gridV.SetPoints(self.points)
        self.edgeFluxesV = vtk.vtkDoubleArray()
        self.edgeFluxesV.SetName('V')
        self.edgeFluxesV.SetNumberOfComponents(1)
        self.edgeFluxesV.SetNumberOfTuples(numCells)
        self.edgeFluxesV.SetVoidArray(self.field.edgeFluxesVArray, numCells, 1)

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
        for j in range(self.field.ny):
            for i in range(self.field.nx):

                ptIds.SetId(0, 4*cellId + 1)
                ptIds.SetId(1, 4*cellId + 2)
                self.gridU.InsertNextCell(vtk.VTK_LINE, ptIds)

                ptIds.SetId(0, 4*cellId + 3)
                ptIds.SetId(1, 4*cellId + 2)
                self.gridV.InsertNextCell(vtk.VTK_LINE, ptIds)

                # increment the cell counter
                cellId += 1

    def show(self, npx=1260, npy=960):

        totalFlux = self.field.pli.getIntegral(self.field.integratedVelocity)
        self.title = vtk.vtkTextActor()
        self.title.SetTextScaleMode(0)
        self.title.GetTextProperty().SetFontSize(50)
        if self.field.sverdrup:
            self.title.SetInput(f"flux = {totalFlux:10.3g} (Sv) @ time {self.field.timeIndex}")
        else:
            self.title.SetInput(f"flux = {totalFlux:10.3g} (A m^2/s) @ time {self.field}")

        # lookup table
        self.lut = vtk.vtkLookupTable()
        nc1 = 10001
        self.lut.SetNumberOfTableValues(nc1)
        for i in range(nc1):
            x = float(i)/float(nc1-1)
            r = numpy.sqrt(x) # x**2
            g = numpy.sin(numpy.pi*x)
            b = 0.6 + 0.4*numpy.sqrt(x)
            a = 1.0
            if x < 0.0001:
                # land or zero flux
                r, g, b, a = 1.0, 0.9, 0.8, 1.0
            self.lut.SetTableValue(i, r, g, b, a)
        self.lut.SetTableRange(0.0, self.field.maxAbsFlux)
        self.lut.Build()

        # colorbar
        self.cbar = vtk.vtkScalarBarActor()
        self.cbar.SetLookupTable(self.lut)
        if self.field.sverdrup:
            self.cbar.SetTitle('flux (Sv)')
        else:
            self.cbar.SetTitle('flux (A m^2/s)')
        self.cbar.SetBarRatio(0.08)

        # tubes for the u fluxes
        self.tubesU = vtk.vtkTubeFilter()
        self.tubesU.SetRadius(self.field.dx * 0.2)
        self.tubesU.SetNumberOfSides(16)
        self.tubesU.SetInputData(self.gridU)
        self.mapperU = vtk.vtkPolyDataMapper()
        self.mapperU.SetInputConnection(self.tubesU.GetOutputPort())
        self.mapperU.SetLookupTable(self.lut)
        self.mapperU.SetUseLookupTableScalarRange(1)
        self.actorU = vtk.vtkActor()
        self.actorU.SetMapper(self.mapperU)
        self.actorU.GetProperty().SetInterpolationToPhong()

        # tubes for the v fluxes
        self.tubesV = vtk.vtkTubeFilter()
        self.tubesV.SetRadius(self.field.dx * 0.2)
        self.tubesV.SetNumberOfSides(16)
        self.tubesV.SetInputData(self.gridV)
        self.mapperV = vtk.vtkPolyDataMapper()
        self.mapperV.SetInputConnection(self.tubesV.GetOutputPort())
        self.mapperV.SetLookupTable(self.lut)
        self.mapperV.SetUseLookupTableScalarRange(1)
        self.actorV = vtk.vtkActor()
        self.actorV.SetMapper(self.mapperV)
        self.actorV.GetProperty().SetInterpolationToPhong()

        # target line
        self.tubePoints = vtk.vtkTubeFilter()
        self.tubePoints.SetRadius(self.field.dx * 0.3)
        self.tubePoints.SetNumberOfSides(16)
        self.tubePoints.SetInputData(self.gridTargetLine)
        self.mapperPoints = vtk.vtkPolyDataMapper()
        self.mapperPoints.SetInputConnection(self.tubePoints.GetOutputPort())
        self.actorPoints = vtk.vtkActor()
        self.actorPoints.SetMapper(self.mapperPoints)
        self.actorPoints.GetProperty().SetColor(1., 0.7, 0.) # white transect

        # add vector plot to target line
        nvpts = self.field.vectorPoints.shape[0]
        maxVectorLength = max([numpy.sqrt(self.field.vectorValues[i, :].dot(self.field.vectorValues[i, :])) for i in range(nvpts)])

        self.targetGlyphs = vtk.vtkGlyph3D()
        self.targetVectorPointData = vtk.vtkDoubleArray()
        self.targetVectorValues = vtk.vtkDoubleArray()
        self.targetVectorPoints = vtk.vtkPoints()
        self.targetVectorGrid = vtk.vtkStructuredGrid()
        self.targetVectorMapper = vtk.vtkPolyDataMapper()
        self.targetGlyphs = vtk.vtkGlyph3D()
        self.targetArrow = vtk.vtkArrowSource()
        self.targetVectorActor = vtk.vtkActor()

        self.targetVectorPointData.SetNumberOfComponents(3)
        self.targetVectorPointData.SetNumberOfTuples(nvpts)
        # move the point a little up for visualization
        self.field.vectorPoints[:, 2] = 0.5*self.field.dx
        self.targetVectorPointData.SetVoidArray(self.field.vectorPoints, nvpts*3, 1)
        self.targetVectorPoints.SetNumberOfPoints(nvpts)
        self.targetVectorPoints.SetData(self.targetVectorPointData)

        self.targetVectorValues.SetNumberOfComponents(3)
        self.targetVectorValues.SetNumberOfTuples(nvpts)
        self.targetVectorValues.SetVoidArray(self.field.vectorValues, nvpts*3, 1)

        self.targetVectorGrid.SetDimensions(nvpts, 1, 1)
        self.targetVectorGrid.SetPoints(self.targetVectorPoints)
        self.targetVectorGrid.GetPointData().SetVectors(self.targetVectorValues)

        self.targetArrow.SetShaftResolution(8)
        self.targetArrow.SetTipResolution(16)

        self.targetGlyphs.SetVectorModeToUseVector()
        self.targetGlyphs.SetScaleModeToScaleByVector()
        self.targetGlyphs.SetSourceConnection(self.targetArrow.GetOutputPort())
        self.targetGlyphs.SetInputData(self.targetVectorGrid)
        self.targetGlyphs.SetScaleFactor(5*self.field.dx/maxVectorLength)

        self.targetVectorMapper.SetInputConnection(self.targetGlyphs.GetOutputPort())
        self.targetVectorMapper.Update()
        self.targetVectorActor.SetMapper(self.targetVectorMapper)
        self.targetVectorActor.GetProperty().SetColor(1., 0.7, 0.) # yellow arrows

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
        self.ren.AddActor(self.targetVectorActor)
        self.ren.SetBackground((0., 0., 0.))
        camera = self.ren.GetActiveCamera()
        lon, lat = 0.5*(self.field.lonmin + self.field.lonmax), 0.5*(self.field.latmin + self.field.latmax)
        camera.SetFocalPoint((lon, lat, 0.))
        camera.SetPosition((lon, lat, 2*self.field.dx*max(self.field.ny, self.field.nx)))
        camera.SetViewUp((0., 1., 0.))
        self.renWin.SetSize(npx, npy)
        self.renWin.SetWindowName('Vertically integrated edge flux')

        # allow the user to interact with the visualisation
        self.callBack = CallBack(self)
        self.iren.AddObserver('KeyPressEvent', self.callBack.execute)
        print('type "t" to step forward in time')

        # This allows the interactor to initalize itself. It has to be
        # called before an event loop.
        self.iren.Initialize()

        self.renWin.Render()

        # Start the event loop.
        self.iren.Start()

def __del__(self):
    self.ncU.close()
    self.ncV.close()


def main(*, tFile: str, uFile: str, vFile: str, lonLatPoints: str='', iFile: str='', sverdrup: bool=False):
    """Visualize fluxes
    :param tFile: netcdf file holding the T-grid
    :param uFile: netcdf file holding u data
    :param vFile: netcdf file holding v data
    :param lonLatPoints: target points "(lon0, lat0), (lon1, lat1),..."
    :param iFile: alternatively read target points from text file
    :param sverdrup: whether or not to use Sverdrup units (default is A m^2/s)
    """
    if lonLatPoints:
        xyVals = numpy.array(eval(lonLatPoints))
    elif iFile:
        llreader = LatLonReader(iFile)
        xyVals = llreader.getLonLats()
    else:
        raise RuntimeError('ERROR must provide either iFile or lonLatPoints!')
    print(f'target points:\n {xyVals}')
    numTargetPoints = xyVals.shape[0]
    lonLatZPoints = numpy.zeros((numTargetPoints, 3), numpy.float64)
    lonLatZPoints[:, :2] = xyVals
    fv = FluxViz(tFile, uFile, vFile, lonLatZPoints, sverdrup)
    fv.show()


if __name__ == '__main__':
    defopt.run(main)
    