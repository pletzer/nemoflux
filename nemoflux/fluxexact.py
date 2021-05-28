import numpy
from numpy import pi, cos, sin
import math
import defopt

def main(*, potentialFunction: str="(cos(t*2*pi/nt)+2)*(0.5*(y/180)**2 + sin(2*pi*x/360))",
            zmin: float=0., zmax: float=1.0, nz: int=5, nt: int=1, 
            deltaDeg: str="(0.,0.)", lonLatPointsStr: str,
            ):
    """Generate data
    :param potentialFunction: potential expression of x (logical lon), y (logical lat), z (depth) and t (time index)
    :param zmin: min depth
    :param zmax: max depth
    :param nz: number of vertical cells
    :param nt: number of time steps
    :param deltaDeg: longitude, latitude pole displacement 
    :param lonLatPointsStr: target points "(lon0, lat0), (lon1, lat1),..."
    """
    xyVals = numpy.array(eval(lonLatPointsStr))

    dz = (zmax - zmin)/float(nz)
    print(f'zmin/zmax = {zmin}/{zmax}')
    # linear for the time being
    zhalf = numpy.array([zmin + (k + 0.5)*dz for k in range(nz)])
    ztop  = numpy.array([zmin + (k + 0  )*dz for k in range(nz)])
    zbot  = numpy.array([zmin + (k + 1  )*dz for k in range(nz)])
    thickness = -(ztop - zbot) # DEPTH HAS OPPOSITE SIGN TO Z

    xyBeg = xyVals[0, :]
    xyEnd = xyVals[-1, :]
    print(f'beg/end target points: {xyBeg} {xyEnd}')

    print('time_index                 flux')
    for t in range(nt):
        flux = 0
        for k in range(nz):
            z = zhalf[k]
            x, y = xyBeg[:2]
            phiA = eval(potentialFunction)
            x, y = xyEnd[:2]
            phiB = eval(potentialFunction)
            flux += (phiB - phiA) * thickness[k]
        print(f'{t:10d} {flux:20.10g}')


if __name__ == '__main__':
    defopt.run(main)




        

