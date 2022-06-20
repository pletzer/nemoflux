# nemoflux

This repository contains scripts that show how to compute the lateral water flow across arbitrary transects. The data are assumed to be stored in 
the NetCDF format. The variable names have been hardwired to match those of ocean NEMO simulations.

## Prerequisites

You'll need the python `mint` package to be installed. This is most easily done on Linux/Mac OSX with the following commands:
```
conda create -n nemoflux
conda activate nemoflux
conda install -c conda-forge python-mint>=1.24.4 xarray defopt
```

## How to download the software

```
git clone https://github.com/pletzer/nemoflux
cd nemoflux/nemoflux
```

## A simple example

NEMO data can be quite large. The following will generate mock data that mimic NEMO's NetCDF files. The variable names and some attributes match those of a NEMO file but in a minimalist way. This is a good way to get started. 
```
python datagen.py --streamFunction="x"
```
will produce three files: T.nc, U.nc and V.nc. Here T.nc contains the grid information (the lon-lat bounds). Files U.nc and V.nc contain the u, v velocity components, respectively. 

The following command will display the total flow for a set of target longitude-latitude points:
```
python fluxviz.py  -t T.nc -u U.nc -v V.nc --lonLatPoints="(-180,-70),(-160,-10),(-35,40),(20,-50),(60,50),(180,40)"
```

![alt simple flux](https://github.com/pletzer/nemoflux/blob/main/pictures/simple.png?raw=true)

The plot shows the rectilinear grid. The edges of each cell are colour coded by the amount of flux. The velocity field is grad(x) x zHat where zHat points out of the screen. With our choice of stream function, the velocity is uniform and points down in the y direction.  

The orange line is the "target" representing the surface extruded in the z direction for which the flux (velocity times area) is computed. For this simple stream function model the total flux is just the difference between the end and start points of the stream function, in our case 360. (If the grid is on the sphere then the units are earth radius `A` times m^2/s).

The velocity on the target line is shown as a set of arrows.

## A singular example

The code is able to recover the exact flux when the target line starts and ends at grid nodes, regardless of the intermediate target points. This is also true when the field is singular. 

To generate a singular field, set the stream function to be proportional to the angle around the singularity

```
python datagen.py --streamFunction="arctan2(y, x+180)/(2*pi)"
python fluxviz.py  -t T.nc -u U.nc -v V.nc --lonLatPoints="(-180,-80), (-10, -80),(-10,80), (-180, 80)"
```

![alt singular flux](https://github.com/pletzer/nemoflux/blob/main/pictures/singular.png?raw=true)

The difference of the stream function taken between the end and the starting points is 0.5, that is the angle difference (`pi`) divided by the normalization factor `2*pi` used in the stream function. 

This flow integral is independent of the path of the target line. Numerically, you should get the exact value if the starting and end points fall on mesh nodes that are at least one cell away from the singularity. Try it!


## A more complex vector field

Let's increase the resolution and have the stream function vary in a more interesting way
```
python datagen.py --streamFunction="cos(2*pi*y/360) + sin(2*pi*x/360)" --nx=360 --ny=180
python fluxviz.py  -t T.nc -u U.nc -v V.nc --lonLatPoints="(-100,-80),(100,-80),(0,80),(-100,-80)"
```
Type "v" in the visualisation window to increase the vector arrow size (or "V" to decrease). The total flux for the closed loop is, to within machine accuracy, zero. This is expected whenever the velocity field derives from a stream function.

![alt total flow at time 0](https://github.com/pletzer/nemoflux/blob/main/pictures/closed2.png?raw=true)


## A curvilinear grid example

An easy way to create a curvilinear grid is by displacing the poles (--deltaDeg option):
```
python datagen.py --streamFunction="cos(2*pi*y/360) + sin(2*pi*x/360)" \
       --nx=360 --ny=180 --deltaDeg="20,30"
python fluxviz.py  -t T.nc -u U.nc -v V.nc --lonLatPoints="(-100,-80),(100,-80),(0,80),(-100,-80)"
```
You can zoom and rotate with the mouse. Type "r" to reset the view, "y"/"Y" to move the camera up and down or "x"/"X" to move the camera to the left/right. 

![rotated pole](https://github.com/pletzer/nemoflux/blob/main/pictures/rotatedPole.png?raw=true)

## Adding elevation and depth

The stream function can take the additional arguments "z" and "t" for depth and time dependence. Try:
```
python datagen.py --streamFunction="(1+10*z)*(t+1)*(cos(2*pi*y/360) + sin(2*pi*x/360))" \
       --nx=360 --ny=180 --nz=10 --nt=20 --deltaDeg="20,30"
python fluxviz.py  -t T.nc -u U.nc -v V.nc --lonLatPoints="(-100,-80),(100,-80),(0,80)"
```

This will show a "tartan" plot where the u and v velocity fields are integrated vertically across the "z" layers and along the T cell edges. The arrows represent the velocity field that is perpendicular to the lateral surface. You can step in time by typing "t" in the window. Type "q" to quit/exit. Zoom in/out using the mouse/pad. To reset the view type "r".


## Computing fluxes from NEMO data

As mentioned previously, NEMO grids can be large and extracting a subset of the global grid, one that focuses on the domain of interest, makes the application more snappy. 

To subset the NEMO data to a smaller domain, type
```
python python subsetNEMO.py -t $TFILE -u $UFILE -v $VFILE --outputdir=mytest \
                            --imin=1200 --imax=1300 --jmin=500 --jmax=600
```
where `$TFILE`, `$UFILE` and `$VFILE` are the names of the T, U, and V netCDF files, respectively, and `--imin`, `--imax`, `--jmin` and `--jmax` are the start/end indices in the input files. This will generate the NetCDF files `fluxviz.py` needs. 
```
python fluxviz.py -t ../data/sa/T.nc -u ../data/sa/U.nc -v ../data/sa/V.nc -s -i ../data/sa/S3_sa.txt
```
![south africa](https://github.com/pletzer/nemoflux/blob/main/pictures/sa.png?raw=true)

Note option "-s" which computes the mass flow in Sverdrup (1 Sv = 1.e6 m^3/s). To create a movie type "m" in the window - this will create PNM files for each frame. You create an MPEG movie form the frames with the command
```
convert -quality 100 -delay 20 fluxviz*.pnm fluxviz.mpg
```
assuming you have ImageMagick installed.

Feel free to edit the target points in file ../data/sa/S3_sa.txt. Note that the target line can intersect land - this will just set the flux to be zero there.
