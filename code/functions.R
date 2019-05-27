
#' Prepare raw NOAA raster input for contour line: crop, recode, project
#' @param rasta Raster* object
#' @param indir input directory
#' @param outdir output directory
#' @param country name of country for output filename
#' @param boundary vector of lat/lon boundary for region of interest c(lon1,lon2,lat1,lat2)
#' @param proj4 projection string. Planar coordinates suggested; WGS84 provided as the default
#' @param histo create histogram of raster cell values? T/F
#' @return Raster* object exported + histogram .png if histo=T

prep_raster <- function(rastr,
                        indir="../input",
                        outdir="../output",
                        country,
                        boundary,
                        proj4="+init=epsg:4326",
                        gdal_new=T,
                        histo=F) {

  library(dplyr)
  raw       <- paste0(rastr,".tif") %>% file.path(indir,.) %>% raster::raster()
  crop      <- raster::crop(raw,boundary)
  no_clouds <- raster::reclassify(crop,cbind(255,Inf,0),right=TRUE)
  if (gdal_new==T){
    recode <- no_clouds
  } else {
    recode <- no_clouds*10 # if GDAL is older than GDAL 2.4, raster must be scaled for gdal_contour to work
  }
  project   <- raster::projectRaster(from=recode,
                                     crs=proj4,
                                     method="ngb")
  raster::writeRaster(project,
                      filename=paste0(rastr,"_",country,".tif") %>% file.path(outdir,.),
                      overwrite=T)

  if (histo==TRUE){
    paste0("histogram_",rastr,"_",country,".png") %>% file.path(outdir,.) %>% png()
    raster::hist(project,
                 maxpixels=1e5,
                 main="Density of 100k pixels",
                 xlab="Pixel value",
                 ylab="Number of pixels")
    dev.off()
  }
}




#' Draw contour lines from prepared raster and convert to polygon
#' This program will only work at levels that are multiples of 10 unless GDAL version >=2.4
#' @param rasta Raster* object
#' @param shapeout file path and name for intermediate output from this routine
#' @param outdir output directory
#' @param thresh contour line brightness threshold (scalar)
#' @param proj4 projection string from prep raster. Planar coordinates suggested; WGS84 provided as the default
#' @return shapefile of light polygons (SpatialPolygonsDataFrame if returned)

contour_lines <- function(rastr,
                          shapeout,
                          outdir="../output",
                          thresh=30,
                          proj4="+init=epsg:4326") {

  system(paste0("rm -r ",file.path(outdir,shapeout),"*"))

  library(dplyr)

    # create contour lines (SpatialLines) file; read back in
    gdalUtils::gdal_contour(src_filename=file.path(outdir,rastr),
                            dst_filename=paste0(shapeout,"_lines.shp") %>% file.path(outdir,.),
                            fl=as.character(thresh),
                            verbose=TRUE)

    # SpatialLines -> SpatialPolygons
    sl_contour <- paste0(shapeout,"_lines.shp") %>% file.path(outdir,.) %>% rgdal::readOGR()
    ps_contour  <- maptools::SpatialLines2PolySet(sl_contour)
    sp_polygon  <- maptools::PolySet2SpatialPolygons(ps_contour,close_polys=T)

    # set up vector IDs
    df <- length(sp_polygon) %>% seq.int() %>% as.matrix() %>% data.frame() %>% data.table::setnames("ntl_id")

    # attach to SpatialPolygonsDataFrame
    sdf_polygon <- sp::SpatialPolygonsDataFrame(sp_polygon,data=as.data.frame(df))

    # re-assign projection
    if ( is.na(sp::proj4string(sp_polygon)) ) {
      sp::proj4string(sdf_polygon) <- sp::CRS(proj4)
    }

    ### dissolve holes
    # if the polygons are few, this is done by exhaustive search fairly quickly
    if (length(sdf_polygon) < 2e3){

        #' @param shp SpatialPolygons* object with holes
        check_dissolve <- function(shp) {
          sapply(1:(length(shp)),function(i)sum(rgeos::gWithin(shp[i,],shp[-i,],byid=T)))
        }

        sdf_polygon_set <- check_dissolve(sdf_polygon)
        sdf_polygon_dissolved <- sdf_polygon[which(sdf_polygon_set==0), ] # keep

    # if polygons are many, use partitions to spatially index polygons
    } else {

        # set up parameters
        n_parts <- 4                                            # number of closed partitions
        raster_bbox  <- raster::extent(sdf_polygon)             # full bounding box (x1,x2,y1,y2)
        x_width  <- (raster_bbox[2]-raster_bbox[1]) / n_parts   # moving window

        # N boxes: moving window partitions
        for (i in 1:n_parts){
          assign(
            paste("box",i,sep="_"),
            as(raster::extent(max(raster_bbox[2]-(n_parts-i+1)*x_width,raster_bbox[1]),
                              min(raster_bbox[2]-(n_parts-i)*x_width,raster_bbox[2]),
                              raster_bbox[3],
                              raster_bbox[4]),"SpatialLines")
          )
        }

        # assign original projections to boxes for geometry operations
        og_proj <- sp::CRS(sp::proj4string(sdf_polygon))
        sp::proj4string(box_1) <- og_proj
        sp::proj4string(box_2) <- og_proj
        sp::proj4string(box_3) <- og_proj
        sp::proj4string(box_4) <- og_proj

        # validate polygon object
        # this will tighten up any self-intersections if they exist
        if (rgeos::gIsValid(sdf_polygon)!=TRUE) {
          sdf_poly_proj <- sp::spTransform(sdf_polygon,sp::CRS(proj4))  # planar coordinates
          sdf_poly_proj <- rgeos::gBuffer(sdf_poly_proj,byid=T,width=0) # requires planar coordinates
          sdf_polygon   <- sp::spTransform(sdf_poly_proj,og_proj)       # return to original projection
        }

        # split up polygon into N parts based on N boxes
        for (i in 1:n_parts){
          assign(
            paste("part",i,sep="_"),
            raster::crop(sdf_polygon,eval(parse(text=paste("box",i,sep="_"))))
          )
        }

        # create N closed partitions + leftover of boundary points
        for (i in 1:n_parts){
          assign(
            paste("result",i,sep="_"),
            eval(parse(text=paste("part",i,sep="_")))
          )
        }
        leftover <- rbind(part_1[which(rgeos::gIntersects(box_1,part_1,byid=T)), ],
                          part_2[which(rgeos::gIntersects(box_2,part_2,byid=T)), ],
                          part_3[which(rgeos::gIntersects(box_3,part_3,byid=T)), ],
                          part_4[which(rgeos::gIntersects(box_4,part_4,byid=T)), ])

        #' @param shp SpatialPolygons* object
        check_for_holes <- function(shp){
          sapply(1:(length(shp)),function(i)
            sum(rgeos::gWithin(shp[i,],shp[-i,],byid=T))
          )
        }

        # check each partition + leftover
        for (i in 1:n_parts) {
          assign(paste0("sdf_hole_set",i),check_for_holes(eval(parse(text=paste("result",i,sep="_")))))
        }
        sdf_hole_left <- check_for_holes(leftover)

        # select polygons that are not holes
        keep_ids <- rbind(result_1[which(sdf_hole_set1==0),]@data,
                          result_2[which(sdf_hole_set2==0),]@data,
                          result_3[which(sdf_hole_set3==0),]@data,
                          result_4[which(sdf_hole_set4==0),]@data,
                          leftover[which(sdf_hole_left==0),]@data)

        # drop holes
        sdf_polygon_dissolved <- sdf_polygon[ sdf_polygon@data$ntl_id %in% unique(keep_ids$ntl_id) , ]

    }

    print(paste0(
      length(sl_contour)," contour lines drawn. ",
      length(sdf_polygon)," polygons formed. ",
      length(sdf_polygon)-length(sdf_polygon_dissolved)," holes found. ",
      length(sdf_polygon_dissolved)," connected spaces produced."
    ))

    # write output
    rgdal::writeOGR(sdf_polygon_dissolved,
                    dsn=outdir,
                    layer=shapeout,
                    driver='ESRI Shapefile',
                    overwrite=TRUE)

}



#' Spatial join between administrative shapefile and polygonized lights shapefile
#' @param contourpoly object of class sf
#' @param adminpoly object of class sf
#' @param outfile name for csv+shp output from this routine
#' @param outdir output directory
#' @param proj4 projection string. Equal area projection suggested; WGS84 provided as the default
#' @return exports admin unit shapefile and csv with -metro area assignments

sp_join <- function(contourpoly,
                    adminpoly,
                    outdir="../output",
                    outfile,
                    proj4="+init=epsg:4326") {

  library(dplyr)
  # load shapefiles
  sf_admin  <- sf::st_read(adminpoly)
  sf_light  <- file.path(outdir,contourpoly) %>% sf::st_read()

  # set up polygons for intersection-area calculations
  sf_admin_proj <- sf::st_transform(sf_admin,proj4)
  sf_light_proj <- sf::st_transform(sf_light,proj4)

  # validate polygons
  if ( sum(sf::st_is_valid(sf_admin_proj)) != nrow(sf_admin_proj) ) {
    sf_admin_proj <- lwgeom::st_make_valid(sf_admin_proj)
  }
  if ( sum(sf::st_is_valid(sf_light_proj)) != nrow(sf_light_proj) ) {
    sf_light_proj <- lwgeom::st_make_valid(sf_light_proj)
  }

  # spatial join
  df_join <- sf::st_join(sf_admin_proj,sf_light_proj,largest=T)

  # re-project out
  base_proj <- sf::st_crs(sf_admin)$proj4string
  sf_out <- sf::st_transform(df_join,base_proj)

  # export town-MSA csv
  paste0(outfile,".csv") %>% file.path(outdir,.) %>% sf::st_write(sf_out,.,delete_dsn=TRUE)
  # export town-MSA shapefile
  paste0(outfile,".shp") %>% file.path(outdir,.) %>% sf::st_write(sf_out,.,delete_dsn=TRUE)

# end of function
}



#' Write data output from shapefile with light polygon metro assignments
#' @param infile Raster* object
#' @param outdir output directory
#' @param outfile name for output of normalized file
#' @param threshes define thresholds to write (individually or combined)
#' @param join single, combined csv output for all metro thresholds? T/F
#' @param id if true, geographic key to join objects on
#' @return exports .csv

data_join <- function(infile,
                      outdir="../output",
                      outfile,
                      threshes=seq(10,60,10),
                      join=T,
                      id){

  for (i in threshes){
    # read files in df_10 df_20 ...
    assign(
      paste("df",i,sep="_"),
      data.table::fread(input=file.path(outdir,paste0(infile,i,".csv")))
    )
    if (join==T) {
      data.table::setnames(eval(parse(text=paste0("df_",i))),
                           "ntl_id",
                           paste0("ntl_",i)) # rename metro ID
    } else {
      data.table::fwrite(eval(parse(text=paste0("df_",i))),
                         file=file.path(outdir,paste0(outfile,i,".csv")),
                         append=F,
                         quote="auto")
    }
  }

  if (join==T) {
    df <- eval(parse(text=paste0("df_",threshes[1])))
    for (j in seq(threshes[2] , threshes[length(threshes)] , by=(threshes[2]-threshes[1]) )){
      df <- dplyr::left_join(df,dplyr::select(eval(parse(text=paste0("df_",j))),paste0("ntl_",j),id),by=id)
    }
    data.table::fwrite(df,
                       file=file.path(outdir,paste0(outfile,".csv")),
                       append=F,
                       quote="auto")
  }


}

## Figure 1 style output

#' Panel A: map two raw input files together: raster and shapefile
#' @param rasterfile Raster* object (input)
#' @param shapefile Spatial* object (input)
#' @param image_out file path and name for graphical output
#' @param proj4 projection string. WGS84 provided as the default, change to your projection
#' @param scheme color for raster. viridis::inferno(9) provided as default
#' @param lcol line color for shapefile. black provided as default
#' @param lwidth line width for shapefile. 0.8 provided as default
#' @param breakvals define thresholds to show on axis scale
#' @param axis_ticks location of ticks on axis scale measure
#' @param axis_labs what to display at those values
#' @param mask mask raster to land area? T/F
#' @param zoom zoom into particular spatial position? T/F
#' @param zoom_box if true, provide vector of coordinates c(lon1,lon2,lat1,lat2) for raster::crop
#' @return exports .png

map_raw_raster <- function(rasterfile,
                           shapefile,
                           image_out,
                           proj4="+init=epsg:4326",
                           scheme=inferno(9),
                           lcol="black",
                           lwidth=0.8,
                           breakvals=c(0,1,10,20,30,40,50,60,63),
                           axis_ticks=seq(0,60,10),
                           axis_labs=c("0","10","20","30","40","50","60+"),
                           mask=TRUE,
                           zoom=FALSE,
                           zoom_box){

  df_raster <- raster::raster(rasterfile)
  sdf_shape <- rgdal::readOGR(shapefile)

  # project shapefile
  sdf_shape <- sp::spTransform(sdf_shape,sp::CRS(proj4))
  if (rgeos::gIsValid(sdf_shape)!=TRUE) {
    sdf_shape <- rgeos::gBuffer(sdf_shape,byid=T,width=0)
  }
  # project raster
  df_raster <- raster::projectRaster(from=df_raster,
                                     crs=proj4,
                                     method="ngb")

  # project bounding box
  if (zoom==T){
    if (proj4!="+init=epsg:4326") {
      d                  <- data.frame(lon=c(zoom_box[1],zoom_box[2]),
                                       lat=c(zoom_box[3],zoom_box[4]))
      sp::coordinates(d) <- c("lon","lat")
      sp::proj4string(d) <- sp::CRS("+init=epsg:4326")
      zoom_box           <- sp::spTransform(d,sp::CRS(proj4))
    }
    boundary <- zoom_box
    df_raster <- raster::crop(df_raster,boundary)
    sdf_shape <- raster::crop(sdf_shape,boundary)
  }

  # mask raster to land area
  if (mask==TRUE) df_plot <- raster::mask(x=df_raster,mask=sdf_shape)

  # plot
  png(paste0("../output/",image_out,".png"))
  par(mar=c(0,0,0,1),oma=c(0,0,0,0),mgp=c(0,0,0))

  raster::plot(sdf_shape,border="black",lwd=0.8)
  raster::plot(df_plot,legend.width=1,alpha=0.85,
              breaks=breakvals,
              axis.args=list(at=axis_ticks,
                         labels=axis_labs),
              col=scheme,
              axes=F,ylab=NA,xlab=NA,add=TRUE,
              xaxs="i",yaxs="i")
  dev.off()
}


# Panel B: radar map of light polygons on administrative land area
#' @param admin_shapefile Spatial* object (input)
#' @param light_shapefile Spatial* object (nightlight_contour_year_??.shp)
#' @param image_out file path and name for graphical output
#' @param proj4 projection string. WGS84 provided as the default
#' @param bubble_fill color for light polygons. red provided as default
#' @param bubble_border color for outline of light polygons. green provided as default
#' @param zoom zoom into particular spatial position? T/F
#' @param zoom_box if true, provide vector pf coordinates c(lon1,lon2,lat1,lat2)
#' @param scale_loc provide vector of coordinates for scale c(lon,lat)
#' @return exports .png

map_light_shp <- function(admin_shapefile,
                          light_shapefile,
                          image_out,
                          proj4="+init=epsg:4326",
                          bubble_fill="Red",
                          bubble_border="Green",
                          zoom=FALSE,
                          zoom_box,
                          scale_loc){

  sdf_admin <- rgdal::readOGR(admin_shapefile)
  sdf_light <- rgdal::readOGR(light_shapefile)

  # project shapefiles
  sdf_admin <- sp::spTransform(sdf_admin,sp::CRS(proj4))
  sdf_light <- sp::spTransform(sdf_light,sp::CRS(proj4))
  if (rgeos::gIsValid(sdf_admin)!=TRUE) {
    sdf_admin <- rgeos::gBuffer(sdf_admin,byid=T,width=0)
  }
  if (rgeos::gIsValid(sdf_light)!=TRUE) {
    sdf_light <- rgeos::gBuffer(sdf_light,byid=T,width=0)
  }

  # project parameters
  if (zoom==T){
    if (proj4!="+init=epsg:4326") {
      # project bounding box
      d                  <- data.frame(lon=c(zoom_box[1],zoom_box[2]),
                                       lat=c(zoom_box[3],zoom_box[4]))
      sp::coordinates(d) <- c("lon","lat")
      sp::proj4string(d) <- sp::CRS("+init=epsg:4326")
      zoom_box           <- sp::spTransform(d,sp::CRS(proj4))
      # project scale location
      s                  <- data.frame(lon=c(scale_loc[1]),
                                       lat=c(scale_loc[2]))
      sp::coordinates(s) <- c("lon","lat")
      sp::proj4string(s) <- sp::CRS("+init=epsg:4326")
      scale_loc          <- sp::spTransform(s,sp::CRS(proj4))
    }
    boundary <- zoom_box
    sdf_admin <- raster::crop(sdf_admin,boundary)
    sdf_light <- raster::crop(sdf_light,boundary)
  }

  png(paste0("../output/",image_out,".png"))
  par(mar=c(1,1,1,1),oma=c(0.25,0,0.25,0),pty="s")

  raster::plot(sdf_admin,col="bisque",border="black",lwd=0.5,xaxs="i",yaxs="i",bg="lightblue")
  raster::plot(sdf_light,col=adjustcolor(c(bubble_fill),alpha.f=0.4),border=bubble_border,lwd=0.5,add=TRUE)
  if (proj4=="+init=epsg:4326") {
    # lat/long
    maps::map.scale(x=scale_loc[1],
                    y=scale_loc[2],
                    ratio=FALSE,
                    relwidth=0.2)
  } else {
    # projected coordinates
    GISTools::map.scale(xc=scale_loc@coords[1],
                        yc=scale_loc@coords[2],
                        len=10e4,
                        ndivs=2,
                        subdiv=50,
                        units="km")
  }
  dev.off()
}


#' Panel C: Map of metro assignments
#' @param metro_shapefile Spatial* object (mapping_geo_year_NTL??.shp)
#' @param image_out file path and name for graphical output
#' @param proj4 projection string. WGS84 provided as the default, change to your projection
#' @param zoom zoom into particular spatial position? T/F
#' @param zoom_box if true, provide vector of coordinates c(lon1,lon2,lat1,lat2) for raster::crop
#' @return .png

map_metros <- function(metro_shapefile,
                       image_out,
                       proj4="+init=epsg:4326",
                       zoom=FALSE,
                       zoom_box){

  sdf_metro <- rgdal::readOGR(metro_shapefile)

  # project admin shapefile
  sdf_metro <- sp::spTransform(sdf_metro,sp::CRS(proj4))
  # validate
  if (rgeos::gIsValid(sdf_metro)!=TRUE) {
    sdf_metro <- rgeos::gBuffer(sdf_metro,byid=T,width=0)
  }

  # project and define bounding box
  if (zoom==T){
    if (proj4!="+init=epsg:4326") {
      d                  <- data.frame(lon=c(zoom_box[1],zoom_box[2]),
                                       lat=c(zoom_box[3],zoom_box[4]))
      sp::coordinates(d) <- c("lon","lat")
      sp::proj4string(d) <- sp::CRS("+init=epsg:4326")
      zoom_box           <- sp::spTransform(d,sp::CRS(proj4))
    }
    boundary <- zoom_box
    sdf_metro <- raster::crop(sdf_metro,boundary)
  }

  png(paste0("../output/",image_out,".png"))
  par(mar=c(1,1,1,1),oma=c(0.25,0,0.25,0),mgp=c(0,0,0))

  # dissolve administrative units into metropolitan areas
  sdf_metro_dissolved <- rgeos::gUnaryUnion(sdf_metro,sdf_metro@data$ntl_id)
  # build coloring for metropolitan areas
  if ("try-error" %in% class(try(
      color_vector <- MapColoring::getOptimalContrast(x=sdf_metro_dissolved,
                                                      col=c("#FF0000","#F2AD00","#00A08A","#9986A5"))
  ))){color_vector <- "red"}

  raster::plot(sdf_metro_dissolved,
               col=color_vector,
               border="black",
               lwd=0.8,xaxs="i",yaxs="i") # metro polygons
  raster::lines(sdf_metro,col="black",lwd=0.4) # administrative borders

  dev.off()
}
