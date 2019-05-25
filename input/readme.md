This folder needs to contain the NOAA raster TIF file and the relevant shapefile.

If your machine supports GNU bash and the `raster_download` parameter in `params.yaml` is set to `TRUE`, the `R` script will retrieve the raster file from the [NOAA website](https://ngdc.noaa.gov/eog/data/web_data/v4composites/) automatically.
If you use `make`, the `Makefile` will perform this download prior to running the `R` script.

If not, you must download the raster file manually, unpack the TAR file, unzip the TIF file, and store it in this folder.
This folder is referenced by the `prep_raster` function in `functions.R` to load the raw NOAA raster and transform it.

After the file you need is in this `input` folder, set the `raster_download` parameter in `params.yaml` to `FALSE` to ensure the file is not downloaded repeatedly each time the code is run.
