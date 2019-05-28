## Metropolitan area algorithm
# Dingel, Miscio, Davis (2019)
# Chinese townships example


#######################
##### Load resources

rm(list=ls()) # clean directory

# install packages
data_packages     <- c("bit64","data.table","dplyr","gdalUtils","lwgeom","PBSmapping","maptools","raster","rgdal","rgeos","sp","sf","units","yaml")
graphics_packages <- c("CEoptim","devtools","GISTools","igraph","MapColoring","maps","RColorBrewer","viridis")
for (i in c(data_packages,graphics_packages)) {
  if (require(i,character.only=TRUE)==FALSE){
    install.packages(i,repos='http://cran.us.r-project.org')
  } else{ require(i,character.only=TRUE) }
}
devtools::install_github("hunzikp/MapColoring")
print("package installation complete")
source("functions.R") # load functions



#######################
##### Set parameters

params <- yaml::read_yaml("params.yaml") # read in constants from YAML

# if GDAL<2.4 and brightness%%10!=0, then raster values must be scaled by 10
require(dplyr)
gdal_version <- system("gdalinfo --version",intern=T) %>% substr(.,6,8) %>% as.numeric()
gdal_status  <- (gdal_version>=2.4 | params$brightness%%10==0)
threshlevel  <- ifelse(gdal_status,params$brightness,params$brightness*10)



######################################
##### Download lights at night image

if (params$raster_download) paste0("bash download_raster.sh ",params$year) %>% system()
noaa_raster <- paste0("basename ../input/F*",params$year,"*_web.stable_lights.avg_vis.tif .tif") %>% system(.,intern=T)



########################################
##### Build file names from parameters

raster_file      <- paste(noaa_raster,params$country_name,sep="_")
contour_file     <- paste("nightlight_contour",params$year,params$brightness,sep="_")
mapping_file     <- paste("mapping",params$geo_name,params$year,"NTL",sep="_")
assignments_file <- paste(params$geo_name,params$year,"NTL",sep="_")



#############################
#####  Apply DMD Algorithm

# clip, clean, and project NOAA raster
prep_raster(rastr    = noaa_raster,
            indir    = "../input",
            outdir   = "../output",
            country  = params$country_name,
            boundary = params$country_bbox,
            proj4    = params$equidistant_proj4,
            gdal_new = gdal_status,
            histo    = TRUE)

# contour lines -> polygon shapefile
contour_lines(rastr    = paste0(raster_file,".tif"),
              shapeout = contour_file,
              outdir   = "../output",
              thresh   = threshlevel,
              proj4    = params$equidistant_proj4)

# spatial join w/ overlapping area
sp_join(contourpoly = paste0(contour_file,".shp"),
        adminpoly   = params$geo_shapefile,
        outfile     = paste0(mapping_file,params$brightness),
        outdir      = "../output",
        proj4       = params$equalarea_proj4)

# write final output
data_join(infile   = mapping_file,
          outfile  = assignments_file,
          outdir   = "../output",
          threshes = params$brightness,
          join     = FALSE,
          id       = params$geo_key)



############################################
#### Produce Figure 1 if zoom_bbox provided

##If params.yaml specifies "zoom_bbox", then produce three maps as PNG files

zoom_status    <- (length(params$zoom_bbox)==4) # zoom_status=TRUE if zoom_box provided
if (zoom_status) { 
zoom_scale_loc <- c(params$zoom_bbox[2]-0.95,params$zoom_bbox[3]+0.2)
raster_blues <- colorRampPalette(c(RColorBrewer::brewer.pal(9,"Blues"),"black")) # define color

#Figure 1A: Raster layer overlaid on land area
map_raw_raster(rasterfile = paste0("../output/",raster_file,".tif"),
               shapefile  = params$geo_shapefile,
               proj4      = params$equalarea_proj4,
               image_out  = paste("nightlight_raster",params$geo_name,params$year,"map",sep="_"),
               scheme     = rev(raster_blues(63)),
               breakvals  = seq(0,63,1),
               zoom_box   = params$zoom_bbox,
               zoom       = zoom_status)

#Figure 1B: Map of light polygons overlaid on administrative units with scalebar
map_light_shp(admin_shapefile = params$geo_shapefile,
              light_shapefile = paste0("../output/",contour_file,".shp"),
              image_out       = paste("radar_map",params$geo_name,params$year,params$brightness,sep="_"),
              proj4           = params$equalarea_proj4,
              zoom_box        = params$zoom_bbox,
              scale_loc       = zoom_scale_loc,
              zoom            = zoom_status)

#Figure 1C: Map of township-metro assignment
map_metros(metro_shapefile = paste0("../output/",mapping_file,params$brightness,".shp"),
           image_out       = paste0("assignment_map_",params$geo_name,"_",params$year,"_metros_NTL",params$brightness),
           proj4           = params$equalarea_proj4,
           zoom_box        = params$zoom_bbox,
           zoom            = zoom_status)
}
