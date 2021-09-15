import netCDF4
import defopt
import sys
import datetime

def main(*, tfile: str='', 
	        ufile: str='',
	        vfile: str='', 
	        outputdir: str='./', jmin: int, jmax: int, imin: int, imax: int):
    """
    subset nemo data
    :param tfile: name of the netCDF file containing T cell grid data
    :param ufile: name of the netCDF file containing u data
    :param vfile: name of the netCDF file containing v data
    :param outputdir: output directory, the files will be saved as T.nc, U.nc and V.nc
    :param jmin: min j index
    :param jmax: max j index
    :param imin: min i index
    :param imax: max i index
    """

    # T file
    print(f'T file: {tfile}')
    nci = netCDF4.Dataset(tfile)
    nco = netCDF4.Dataset(f'{outputdir}/T.nc', 'w')
    nco.command = sys.argv
    nco.timestamp = f'generated on {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

    for d in nci.dimensions:
        print(f'creating dimension {d}')
        if d == 'x':
            nco.createDimension(d, imax - imin)
        elif d == 'y':
            nco.createDimension(d, jmax - jmin)
        else:
            nco.createDimension(d, nci.dimensions[d].size)

    for vname in 'bounds_lon', 'bounds_lat', 'deptht', 'deptht_bounds':
        print(f'creating variable {vname}')
        vari = nci.variables[vname]
        varo = nco.createVariable(vname, vari.dtype, vari.dimensions)
        for a in vari.ncattrs():
            val = getattr(vari, a)
            print(f'\tattribute {a} has value {val}')
            setattr(varo, a, val)
        # write the variable
        if 'lon' in vname or 'lat' in vname:
            varo[:] = vari[jmin:jmax, imin:imax]
        else:
            varo[:] = vari[:]
    nco.close()
    nci.close()

    # U, V file
    varnameMap = {'U': 'uo', 'V': 'vo'}
    filenameMap = {'U': ufile, 'V': vfile}

    for field in 'U', 'V':
        nci = netCDF4.Dataset(filenameMap[field])
        nco = netCDF4.Dataset(f'{outputdir}/{field}.nc', 'w')
        nco.command = sys.argv
        nco.timestamp = f'generated on {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'

        for d in nci.dimensions:
            print(f'creating dimension {d}')
            if d == 'x':
                nco.createDimension(d, imax - imin)
            elif d == 'y':
                nco.createDimension(d, jmax - jmin)
            else:
                nco.createDimension(d, nci.dimensions[d].size)

        fname = varnameMap[field]
        for vname in 'time_counter', 'time_centered', 'time_centered_bounds', fname:
            print(f'creating variable {vname}')
            vari = nci.variables[vname]
            if vname in ('uo', 'vo'):
                varo = nco.createVariable(vname, vari.dtype, vari.dimensions, fill_value=vari._FillValue, zlib=True)
            else:
                varo = nco.createVariable(vname, vari.dtype, vari.dimensions)
            for a in vari.ncattrs():
                if a == '_FillValue':
                    continue
                val = getattr(vari, a)
                print(f'\tattribute {a} has value {val}')
                setattr(varo, a, val)
            # write
            if vname == 'uo' or vname == 'vo':
                varo[:] = vari[..., jmin:jmax, imin:imax]
            else:
                varo[:] = vari[:]
        nco.close()
        nci.close()



if __name__ == '__main__':
    defopt.run(main)
