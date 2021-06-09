# nemoflux

Scripts to compute the lateral water flow across longitude-latitude lines from ocean NEMO simulation data

## Rerequisites

You'll need the python mint package to be installed. This is most easily done on Linux/Mac OSX with the following commands:
```
conda create -n nemoflux python=3.8
conda activate nemoflux
conda install -c conda-forge python-mint>=1.14.6
```
**NOTE: python-mint version 1.14.6 should become available soon**

## How to download the software

```
git clone https://github.com/pletzer/nemoflux
cd nemoflux/nemoflux
```

## A quick tour of nemoflux

NEMO data can be quite large. The following will generate mock data that comply to NEMO's NetCDF file. The variable names and attributes match those of a NEMO file. This is a good way to get started. 
```
python datagen.py --deltaDeg="20,30" --nz=10 --nt=4 --nx=36 --ny=18 --prefix=test
```
Will produce three files:
```
ls test_?.nc
```
should display test_T.nc, test_U.nc	and test_V.nc. test_T.nc contains the grid information, test_U.nc and test_V.nc contain the u, v velocity components. The following command will display the total flow for the target longitude, latitude points (-129,-80),(-23,-34),(156,78):
```
python fluxviz.py  -t test_T.nc -u test_U.nc -v test_V.nc -l "(-129,-80),(-23,-34),(156,78)"
```

![alt total flow at time 0](https://github.com/pletzer/nemoflux/blob/main/picture/simple.png?raw=true)

You can step in time by typing "t" in the window. Type "q" to quit/exit. Zoom in/out using the mouse/pad. To reset the view type "r".


## How to subset NEMO data

To subset the NEMO data to a smaller domain
```
python python subsetNEMO.py -t $TFILE -u $UFILE -v $VFILE --outputdir=mytest --imin=1200 --imax=1300 --jmin=500 --jmax=600
```
where `$TFILE`, `$UFILE` and `$VFILE` are the names of the T, U, and V netCDF files. 

## How to compute the total flow across an irregular path

```
python fluxviz.py -t ../data/sa/T.nc -u ../data/sa/U.nc -v ../data/sa/V.nc -s ../data/sa/S3_sa.txt
```
![alt total flow at time 0](https://github.com/pletzer/nemoflux/blob/main/picture/sa.png?raw=true)

Feel free to edit the target points in file ../data/sa/S3_sa.txt. 