import vtk
import netCDF4
import geo
import defopt
import numpy
from horizgrid import HorizGrid
import mint
from latlonreader import LatLonReader


EARTH_RADIUS = 6371000.0

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

        self.sverdrup = sverdrup

        # read the cell bounds
        with netCDF4.Dataset(tFile) as nc:
            bounds_lat = nc.variables['bounds_lat'][:]
            bounds_lon = nc.variables['bounds_lon'][:]
            self.bounds_depth = nc.variables['deptht_bounds'][:]

        self.lonmin = bounds_lon.min()
        self.lonmax = bounds_lon.max()
        self.latmin = bounds_lat.min()
        self.latmax = bounds_lat.max()
        print(f'lon-lat box: {self.lonmin}, {self.latmin} -> {self.lonmax}, {self.latmax}')

        self.timeIndex = 0
        self.ncU = netCDF4.Dataset(uFile)
        self.ncV = netCDF4.Dataset(vFile)

        self.nt, self.nz, self.ny, self.nx = self.getSizes()

        self.gr = HorizGrid(tFile)
        self.pli = mint.PolylineIntegral()
        self.pli.build(self.gr.getMintGrid(), lonLatPoints,
                       counterclock=False, periodX=360.0)


        self.thickness = self.bounds_depth[:, 1] - self.bounds_depth[:, 0]

        numCells = self.ny * self.nx
        self.dx = min((self.lonmax - self.lonmin)/float(self.nx), (self.latmax - self.latmin)/float(self.ny))
        self.arcLengths = numpy.zeros((numCells, 4), numpy.float64)
        self.computeArcLengths()

        # compute the edge fluxes from the vector fields
        self.edgeFluxesUArray = numpy.zeros((numCells,), numpy.float64)
        self.edgeFluxesVArray = numpy.zeros((numCells,), numpy.float64)
        self.integratedVelocity = numpy.zeros((numCells, 4), numpy.float64)
        self.maxAbsFlux = 0.

        # read/get the staggered field integrated over the depth
        uVerticallyIntegrated, vVerticallyIntegrated = self.getUV()
        self.computeIntegratedFlux(uVerticallyIntegrated, vVerticallyIntegrated)
        print(f'max vertically integrated edge |flux|: {self.maxAbsFlux}')

        self.buildTargetLineGrid(lonLatPoints)
        self.buildEdgeUVGrids(bounds_lon, bounds_lat)

        # create a polyline grid along the target points and place vectors on them
        vectorPoints = []
        self.uVectors = []
        for i in range(len(lonLatPoints) - 1):
            begPoint = lonLatPoints[i]
            endPoint = lonLatPoints[i + 1]
            u = endPoint - begPoint
            distance = numpy.sqrt(u.dot(u))
            # normalize
            u /= distance
            nvpts = max(2, int(distance / self.dx))
            vdx = distance / float(nvpts - 1)
            for j in range(nvpts):
                vectorPoints.append(begPoint + u*j*vdx)
                self.uVectors.append(u)
        self.vectorPoints = numpy.array(vectorPoints)

        # compute the vector at the target line
        self.vinterp = mint.VectorInterp()
        self.vinterp.setGrid(self.gr.getMintGrid())
        self.vinterp.buildLocator(numCellsPerBucket=128, periodX=360.)
        self.vinterp.findPoints(self.vectorPoints, tol2=1.e-12)
        self.vectorValues = self.vinterp.getFaceVectors(self.integratedVelocity)

    def update(self, key):

        camera = self.ren.GetActiveCamera()
        x, y, z = camera.GetPosition()

        if key == 't':
            # forward in time
            self.timeIndex = (self.timeIndex + 1) % self.nt
        elif key == 'T':
            # backward in time
            self.timeIndex = (self.timeIndex - 1) % self.nt
        else:
            if key == 'r':
                # reset 
                lon, lat = 0.5*(self.lonmin + self.lonmax), 0.5*(self.latmin + self.latmax)
                camera.SetFocalPoint((lon, lat, 0.))
                camera.SetPosition((lon, lat, 2*self.dx*max(self.ny, self.nx)))
                camera.SetViewUp((0., 1., 0.))
            elif  key == 'x':
                camera.SetPosition(x + self.dx, y, z)
                camera.SetFocalPoint(x + self.dx, y, 0.)
            elif key == 'X':
                camera.SetPosition(x - self.dx, y, z)
                camera.SetFocalPoint(x - self.dx, y, 0.)                
            elif key == 'y':
                camera.SetPosition(x, y + self.dx, z)
                camera.SetFocalPoint(x, y + self.dx, 0.)
            elif key == 'Y':
                camera.SetPosition(x, y - self.dx, z)
                camera.SetFocalPoint(x, y - self.dx, 0.)
            elif key == 'z':
                camera.SetPosition(x, y, z + self.dx)
            elif key == 'Z':
                camera.SetPosition(x, y, z - self.dx)
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
                filename = f'fluxviz_{self.timeIndex:04d}.pnm'
                wr.SetFileName(filename)
                wr.Write()
                print(f'saved screen shot for time {self.timeIndex} in file {filename}')
            elif key == 'm':
                # save all frames
                for i in range(self.nt):
                    self.update(key='t')
                    self.update(key='s')
            camera.Modified()
            return

        # update the data
        uVerticallyIntegrated, vVerticallyIntegrated = self.getUV()
        # this will update self.integratedVelocity
        self.computeIntegratedFlux(uVerticallyIntegrated, vVerticallyIntegrated)

        self.vectorValues[:] = self.vinterp.getFaceVectors(self.integratedVelocity)

        # update the pipeline
        self.lut.SetTableRange(0., self.maxAbsFlux)
        self.lut.Modified()
        self.cbar.Modified()
        totalFlux = self.pli.getIntegral(self.integratedVelocity)
        if self.sverdrup:
            self.title.SetInput(f'flux = {totalFlux:6.3g} (Sv) @ time {self.timeIndex}')
        else:
            self.title.SetInput(f'flux = {totalFlux:6.3g} (A m^2/s) @ time {self.timeIndex}')
        self.title.Modified()
        self.edgeFluxesU.Modified()
        self.edgeFluxesV.Modified()
        self.targetVectorValues.Modified()
        print(f'time index now {self.timeIndex} max |flux|: {self.maxAbsFlux:10.3f} nt = {self.nt}')

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

    def readField(self, nc, fieldName):

        # read u
        try:
            field = nc.variables[fieldName][self.timeIndex, :, :, :]
        except:
            try:
                field = nc.variables[fieldName][...]
            except:
                raise RuntimeError(f'ERROR: could not read {fieldName} field')

        # set to zero where missing
        field = numpy.ma.filled(field, 0.0)

        # integrate vertically, multiplying by the thickness of the layers
        # sum of multiplying axis 0 of thickness with axis 0 of field...
        fieldVerticallyIntegrated = numpy.tensordot(self.thickness, field, axes=(0, 0))

        return fieldVerticallyIntegrated

    def getUV(self):
        uVerticallyIntegrated = self.readField(self.ncU, 'uo')
        vVerticallyIntegrated = self.readField(self.ncV, 'vo')
        return uVerticallyIntegrated, vVerticallyIntegrated

    def computeArcLengths(self):

        # lon, lat, elev coords
        points = self.gr.getPoints() # array of size (numCells, 4, 3)

        # convert to Cartesian
        xyz = geo.lonLat2XYZArray(points, radius=geo.EARTH_RADIUS) # array of size (numCells, 4, 3)

        # compute the arc lengths
        for i0 in range(4):
            i1 = (i0 + 1) % 4
            self.arcLengths[:, i0] = geo.getArcLengthArray(xyz[:, i0, :], xyz[:, i1, :], radius=geo.EARTH_RADIUS)

    def computeIntegratedFlux(self, uVerticallyIntegrated, vVerticallyIntegrated):

        numCells = self.ny * self.nx

        #        ^
        #        |
        #  3-----V-----2
        #  |           |
        #  |     T     U-->
        #  |           |
        #  0-----------1

        self.edgeFluxesUArray[:] = uVerticallyIntegrated.reshape((numCells,)) * self.arcLengths[:, 1]
        self.edgeFluxesVArray[:] = vVerticallyIntegrated.reshape((numCells,)) * self.arcLengths[:, 2]

        #        2
        #        ^
        #        |
        #  3-----V-----2
        #  |           |
        #  |->3  0     U->1
        #  |     ^     |
        #  |     |     |
        #  0-----------1

        # east
        self.integratedVelocity[:, 1] = self.edgeFluxesUArray
        # north
        self.integratedVelocity[:, 2] = self.edgeFluxesVArray

        # these should provide a different view to the array, NOT A COPY
        eU = self.edgeFluxesUArray.reshape((self.ny, self.nx))
        eV = self.edgeFluxesVArray.reshape((self.ny, self.nx))
        iV = self.integratedVelocity.reshape((self.ny, self.nx, 4))

        # south, (ny, nx, 4)
        iV[1:, :, 0] = eV[:-1, :] # will set self.integratedVelocity on the south side
        # west
        iV[:, 1:, 3] = eU[:, :-1]
        # periodic BCs
        iV[:, 0, 3] = eU[:, -1]

        if self.sverdrup:
            eU *= EARTH_RADIUS / 1.e6
            eV *= EARTH_RADIUS / 1.e6
            iV *= EARTH_RADIUS / 1.e6

        # from now on, edge fluxes are abs values
        self.edgeFluxesUArray[:] = numpy.fabs(self.edgeFluxesUArray)
        self.edgeFluxesVArray[:] = numpy.fabs(self.edgeFluxesVArray)

        self.maxAbsFlux = max(self.maxAbsFlux, self.edgeFluxesUArray.max(), self.edgeFluxesVArray.max())

    def show(self, npx=1260, npy=960):

        totalFlux = self.pli.getIntegral(self.integratedVelocity)
        self.title = vtk.vtkTextActor()
        self.title.SetTextScaleMode(0)
        self.title.GetTextProperty().SetFontSize(50)
        if self.sverdrup:
            self.title.SetInput(f"flux = {totalFlux:10.3g} (Sv) @ time {self.timeIndex}")
        else:
            self.title.SetInput(f"flux = {totalFlux:10.3g} (A m^2/s) @ time {self.timeIndex}")

        # lookup table
        self.lut = vtk.vtkLookupTable()
        nc1 = 10001
        self.lut.SetNumberOfTableValues(nc1)
        for i in range(nc1):
            x = float(i)/float(nc1-1)
            r = x**2
            g = numpy.sin(numpy.pi*x/2.)**2
            b = 0.4 + 0.6*numpy.sqrt(x)
            a = 1.0
            if x < 0.0001:
                # land or zero flux
                r, g, b, a = 0.8, 0.75, 0.7, 1.0
            self.lut.SetTableValue(i, r, g, b, a)
        self.lut.SetTableRange(0.0, self.maxAbsFlux)
        self.lut.Build()

        # colorbar
        self.cbar = vtk.vtkScalarBarActor()
        self.cbar.SetLookupTable(self.lut)
        if self.sverdrup:
            self.cbar.SetTitle('flux (Sv)')
        else:
            self.cbar.SetTitle('flux (A m^2/s)')
        self.cbar.SetBarRatio(0.08)

        # tubes for the u fluxes
        self.tubesU = vtk.vtkTubeFilter()
        self.tubesU.SetRadius(self.dx * 0.2)
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
        self.tubesV.SetRadius(self.dx * 0.2)
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
        self.tubePoints.SetRadius(self.dx * 0.3)
        self.tubePoints.SetNumberOfSides(16)
        self.tubePoints.SetInputData(self.gridTargetLine)
        self.mapperPoints = vtk.vtkPolyDataMapper()
        self.mapperPoints.SetInputConnection(self.tubePoints.GetOutputPort())
        self.actorPoints = vtk.vtkActor()
        self.actorPoints.SetMapper(self.mapperPoints)
        self.actorPoints.GetProperty().SetColor(1., 0.7, 0.) # white transect

        # add vector plot to target line
        nvpts = self.vectorPoints.shape[0]
        maxVectorLength = max([numpy.sqrt(self.vectorValues[i, :].dot(self.vectorValues[i, :])) for i in range(nvpts)])

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
        self.vectorPoints[:, 2] = 0.5*self.dx
        self.targetVectorPointData.SetVoidArray(self.vectorPoints, nvpts*3, 1)
        self.targetVectorPoints.SetNumberOfPoints(nvpts)
        self.targetVectorPoints.SetData(self.targetVectorPointData)

        self.targetVectorValues.SetNumberOfComponents(3)
        self.targetVectorValues.SetNumberOfTuples(nvpts)
        self.targetVectorValues.SetVoidArray(self.vectorValues, nvpts*3, 1)

        self.targetVectorGrid.SetDimensions(nvpts, 1, 1)
        self.targetVectorGrid.SetPoints(self.targetVectorPoints)
        self.targetVectorGrid.GetPointData().SetVectors(self.targetVectorValues)

        self.targetArrow.SetShaftResolution(8)
        self.targetArrow.SetTipResolution(16)

        self.targetGlyphs.SetVectorModeToUseVector()
        self.targetGlyphs.SetScaleModeToScaleByVector()
        self.targetGlyphs.SetSourceConnection(self.targetArrow.GetOutputPort())
        self.targetGlyphs.SetInputData(self.targetVectorGrid)
        self.targetGlyphs.SetScaleFactor(5*self.dx/maxVectorLength)

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
        lon, lat = 0.5*(self.lonmin + self.lonmax), 0.5*(self.latmin + self.latmax)
        camera.SetFocalPoint((lon, lat, 0.))
        camera.SetPosition((lon, lat, 2*self.dx*max(self.ny, self.nx)))
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
    