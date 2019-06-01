
This repository contains code implementing the algorithm introduced in "Cities, Lights, and Skills in Developing Economies" in the _Journal of Urban Economics_ by Jonathan Dingel, Antonio Miscio, and Don Davis.
The `R` code constructs metropolitan areas by aggregating finer geographic units on the basis of contiguous areas of light in nighttime satellite images.
As an example, we apply the algorithm to townships in China in 2000, as in Figure 1 of our article.

We thank Dylan Clarke, who wrote the majority of the `R` code appearing in this repository.


## Software requirements

The algorithm is implemented in [R](https://cran.r-project.org/).
We ran our code using R 3.5.1.
Our R code leverages spatial and measurement packages with additional system requirements, namely [`gdalUtils`](https://cran.r-project.org/web/packages/gdalUtils/index.html), [`rgdal`](https://cran.r-project.org/web/packages/rgdal/index.html), [`rgeos`](https://cran.r-project.org/web/packages/rgeos/index.html), [`sp`](https://cran.r-project.org/web/packages/units/index.html), [`sf`](https://cran.r-project.org/web/packages/sf/index.html), and [`units`](https://cran.r-project.org/web/packages/units/index.html).
We used [GEOS](https://trac.osgeo.org/geos/) 3.7.0, [GDAL](https://www.gdal.org/usergroup0.html) 2.3.2, [PROJ](https://proj4.org/download.html) 4.9, and [udunits](https://www.unidata.ucar.edu/software/udunits/udunits-current/doc/udunits/udunits2.html) 2.2.
We expect the code to work on other versions too.

We automate the downloading of nighttime satellite images and invocation of the `R` script using Unix's [`make`](http://swcarpentry.github.io/make-novice/) utility.
We strongly recommend a computing environment that supports GNU bash,
but this is not necessary to run the `R` code.

## Running the code

### Download

First, download (or clone) this repository by clicking the green `Clone or download` button above.
Uncompress the ZIP file into a working directory on your cluster or local machine.
You will see three folders: `code`, `input`, and `output`.

### Example: Chinese townships in 2000

The `code/params.yaml` file included in the repository contains parameters to produce metropolitan areas for China in 2000 by aggregating townships on the basis of lights at night above a brightness threshold of 30.
The resulting output was used to produce Figure 1 in Dingel, Miscio, and Davis (2019).

At the Unix/Linux/MacOSX command line, navigate to the `code` directory and type `make`.
```
cd code
make
```
This will download the nighttime satellite image from NOAA's website and then execute `calls.R` using the parameters declared in `params.yaml`.
* The `Makefile` assumes that your machine is connected to the Internet and that `Rscript` is a valid command name.
* If you are in a computing environment that supports the [Slurm workload manager](https://slurm.schedmd.com/) (if the `Makefile` detects that the command `sbatch` is valid), tasks will be submitted as jobs to your computing cluster.
* If `sbatch` is not available, the `Makefile` will execute the `Rscript` command locally.

If your environment does not support `make` (e.g., some variants of Windows), follow the instructions in `input/readme.md` to download the [NOAA TIF file](https://ngdc.noaa.gov/eog/data/web_data/v4composites/F152000.v4.tar) to the `input` folder.
Then run `calls.R`.

### Adapting to your use case

Edit `code/params.yaml` file to declare the parameters for the year, geographic area, and projections you desire.
Identify the shapefile for your use case by editing the `geo_shapefile` path in line 12 of `params.yaml`.
Then run `calls.R` (the `Makefile` parses `params.yaml`, so you should be able to just type `make` after editing `params.yaml`).
You should never need to edit `calls.R` nor `functions.R`.
