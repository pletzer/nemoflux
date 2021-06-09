# nemoflux

Scripts to compute the lateral water flow across longitude-latitude lines from ocean NEMO simulation data

## Rerequisites

You'll need the python mint package to be installed. This is most easily done on Linux/Mac OSX with the following commands:
```
conda create -n nemoflux python=3.8
conda activate nemoflux
conda install -c conda-forge python-mint>=1.14.6
```

## How to download the software

```
git clone https://github.com/pletzer/nemoflux

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
should display test_T.nc, test_U.nc	and test_V.nc. test_T.nc contains the grid information, test_U.nc the u velocity component and test_V.nc the v velocity component. The following command will display the total flow for the target longitude, latitude points (-129,-80),(-23,-34),(156,78):
```
python fluxviz.py  -t test_T.nc -u test_U.nc -v test_V.nc -l "(-129,-80),(-23,-34),(156,78)"
```

![alt total flow at time 0](https://github.com/pletzer/nemoflux/blob/main/picture/simple.png?raw=true)

You can step in time by typing "t" in the window. Type "q" to quit/exit.

## How to subset NEMO data

TO DO


## How to compute the total flow across an irregular path

```
cd nemoflux/nemoflux
python fluxviz.py -t 

```

