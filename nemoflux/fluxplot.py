import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import geo
import defopt
import numpy
from horizgrid import HorizGrid
import mint
# to read the target points from file
from latlonreader import LatLonReader
from field import Field
import re
import glob
import os.path
import matplotlib
matplotlib.rcParams.update({'font.size': 20})

def main(*, tFile: str, uFile: str, vFile: str, lonLatPoints: str='', iFiles: str='', sverdrup: bool=False):
    """Visualize fluxes
    :param tFile: netcdf file holding the T-grid
    :param uFile: netcdf file holding u data
    :param vFile: netcdf file holding v data
    :param lonLatPoints: target points "[(lon0, lat0), (lon1, lat1),...],[...]"
    :param iFiles: alternatively read target points from text files "['file1', 'file2',...]"
    :param sverdrup: whether or not to use Sverdrup units (default is A m^2/s)
    """
    if lonLatPoints:
        lonLatZPoints = [numpy.array([(ll[0], ll[1], 0.0) for ll in llp]) for llp in eval(lonLatPoints)] 
    elif iFiles:
        lonLatZPoints = []
        listOfFiles = []
        try:
=            listOfFiles = eval(iFiles)
        except:
            try:
                # can return in any order
                listOfFiles = glob.glob(iFiles)
            except:
                # single target line file?
                listOfFiles.append(iFiles)
        print(f'list of target surfaces: {listOfFiles}')
        for iFile in listOfFiles:
            llreader = LatLonReader(iFile)
            lonLatZPoints.append([(ll[0], ll[1], 0.) for ll in llreader.getLonLats()])
    else:
        raise RuntimeError('ERROR must provide either iFiles (-i) or lonLatPoints (-l)!')
    print(f'target points:\n {lonLatZPoints}')
    fld = Field(tFile, uFile, vFile, lonLatZPoints, sverdrup)
    results = {}
    timeVals = []
    for itime in range(fld.nt):
        fld.update()
        timeVals.append(fld.timeObj.getTimeAsDate(itime))
        lineIndex = 0
        for pli in fld.plis:
            results[lineIndex] = results.get(lineIndex, []) + [pli.getIntegral(fld.integratedVelocity)]
            lineIndex += 1
        fld.timeIndex += 1
    # plotting
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%Y'))
    plt.gcf().autofmt_xdate()
    lineTypes = ['b-', 'm--', 'c-.', 'r:', 'g-', 'k--']
    lgds = []
    for lineIndex, values in results.items():
        plt.plot(timeVals, values, lineTypes[lineIndex % len(lineTypes)])
        lgds.append(os.path.basename(listOfFiles[lineIndex]))
    if len(lgds) > 1:
        plt.legend(lgds)
    plt.title('Water flow')
    plt.xlabel('month-year')
    if sverdrup:
        plt.ylabel('Sv')
    else:
        plt.ylabel('A m^2/s')
    plt.show()


if __name__ == '__main__':
    defopt.run(main)
    