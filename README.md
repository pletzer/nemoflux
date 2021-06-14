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
cd nemoflux/nemoflux
```

## A simple example

NEMO data can be quite large. The following will generate mock data that mimick NEMO's NetCDF files. The variable names and some attributes match those of a NEMO file in a minimalistic way. This is a good way to get started. 
```
python datagen.py --streamFunction="-x"
```
will produce three files:
```
ls ?.nc
```
should display T.nc, U.nc	and V.nc. Here T.nc contains the grid information, U.nc and V.nc contain the u, v velocity components. The following command will display the total flow for the target longitude-latitude points (-180,-90),(-160,-10),(20,-50),(60,90),(180,40):
```
python fluxviz.py  -t T.nc -u U.nc -v V.nc --lonLatPoints="(-180,-70),(-160,-10),(-35,40),(20,-50),(60,50),(180,40)"
```

![alt total flow at time 0](https://github.com/pletzer/nemoflux/blob/main/pictures/simple.png?raw=true)
The above shows the grid (rectilinear), the fluxes on each grid edge as colour coded tubes and the target line (orange) over which the flux is computed. In this example, the stream function is "-x" and the velocity is a cross product of zHat times grad (-x), which gives a velocity pointing up in the y direction. (Our zHat points down as is expected for ocean depth.) The orange arrows show the flux, perpendicular to the target line. The computed flux is 360, which exactly matches the difference between the end and start longitude coordinates of the target line, as we would expect for a velocity that derives from a stream function.

## A closed contour flux calculation example

The target line can be a closed contour, for instance
```
python fluxviz.py  -t T.nc -u U.nc -v V.nc --lonLatPoints="(-100,-80),(100,-80),(0,80),(-100,-80)"
```
in which case the total flux must be zero. Note that the last point replicates the starting point to close the contour.
![alt total flow at time 0](https://github.com/pletzer/nemoflux/blob/main/pictures/closed.png?raw=true)

## A more complex vector field

Let's increae the resolution and have the stream function vary in a more interesting way
```
python datagen.py --streamFunction="cos(2*pi*y/360) + sin(2*pi*x/360)" --nx=360 --ny=180
python fluxviz.py  -t T.nc -u U.nc -v V.nc --lonLatPoints="(-100,-80),(100,-80),(0,80),(-100,-80)"
```
Type "v" in the visualisation window to increase the vector arrow size ("V" to decrease). Again the total flux for the closed loop is, to within machine accuracy, zero.
![alt total flow at time 0](https://github.com/pletzer/nemoflux/blob/main/pictures/closed2.png?raw=true)


## A curvilinear grid example

An easy way to create a curvilinear grid is by displacing the poles:
```
python datagen.py --streamFunction="cos(2*pi*y/360) + sin(2*pi*x/360)" --nx=360 --ny=180 --deltaDeg="20,30"
python fluxviz.py  -t T.nc -u U.nc -v V.nc --lonLatPoints="(-100,-80),(100,-80),(0,80),(-100,-80)"
```
You can zoom and rotate with the mouse. Type "r" to reset the view, "y"/"Y" to move the camera up and down. 

![alt total flow at time 0](https://github.com/pletzer/nemoflux/blob/main/pictures/rotatedPole.png?raw=true)

## Adding elevation and depth

The stream function can take the additional arguments "z" and "t" for depth and time dependence. Try:
```
python datagen.py --streamFunction="(1+10*z)*(t+1)*(cos(2*pi*y/360) + sin(2*pi*x/360))" --nx=360 --ny=180 --nz=10 --nt=20 --deltaDeg="20,30"
python fluxviz.py  -t T.nc -u U.nc -v V.nc --lonLatPoints="(-100,-80),(100,-80),(0,80)"
```



This shows a "tartan" plot where the u and v velocity fields are integrated vectically across layers and along the T cell edges. The arrows represent the velocity field that is perpendicular to the lateral surface. You can step in time by typing "t" in the window. Type "q" to quit/exit. Zoom in/out using the mouse/pad. To reset the view type "r".


## Computing fluxes from NEMO data

As mentioned previously, NEMO grids can be large and extracting a subset of the global, one that focuses on the domain of interest, makes the application more snappy. 

To subset the NEMO data to a smaller domain
```
python python subsetNEMO.py -t $TFILE -u $UFILE -v $VFILE --outputdir=mytest --imin=1200 --imax=1300 --jmin=500 --jmax=600
```
where `$TFILE`, `$UFILE` and `$VFILE` are the names of the T, U, and V netCDF files, respectively, and `--imin`, `--imax`, `--jmin` and `--jmax` are the start/end indices in the input files. 

## How to compute the total flow across an irregular path

```
python fluxviz.py -t ../data/sa/T.nc -u ../data/sa/U.nc -v ../data/sa/V.nc -s ../data/sa/S3_sa.txt
```
![alt total flow at time 0](https://github.com/pletzer/nemoflux/blob/main/pictures/sa.png?raw=true)

Note that the flux/flow should be multiplied by the earth's radius (A = 6371000 metres) to get m^3/s. Feel free to edit the target points in file ../data/sa/S3_sa.txt. The flux crossing land is set to zero.
