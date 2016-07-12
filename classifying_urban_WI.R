## [[file:utc.org::*Libraries][Libraries:1]]
## install.packages(c("gdalUtils","ascii","rgeos","mlr","broom","rgdal","raster","plyr","ggplot2","dplyr","tidyr","stringr","foreach","doParallel","glcm","randomForest","kernlab","irace","parallelMap"))
##install.packages("e1071")
##  install.packages("FSelector")

  library(gdalUtils)
  library(ascii)
  library(rgeos)
  library(mlr)
  library(broom)
  library(rgdal)
  library(raster)
  library(plyr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(foreach)
  library(doParallel)
  library(glcm)
  library(randomForest)
  library(kernlab)
  library(irace)
  library(parallelMap)
  library(FSelector)
## Libraries:1 ends here

## [[file:utc.org::*Input%20Directories][Input\ Directories:1]]
image.names <- c("NAIP","PanshpSPOT")
image.dirs <- paste0("../RD_",image.names)
pca.dir <- "../RD_PCA_Regions"
training.dir <- "../RD_Training_Regions"
accuracy.dir <- "../RD_Accuracy"
grids.accuracy.dir <- str_c(accuracy.dir, "/Grids")
fieldplots.accuracy.dir<- str_c(accuracy.dir, "/FieldData")
crop.dir <- "../RD_CroplandDataLayer"
water.dir <- "../RD_WI-waterbody-24k"
urban.dir <- "../RD_US_UrbanAreasShapefile"
urban.and.incorporated.dir <- "../RD_merged_WIurbanAreas_and_incorporatedAreas"
## Input\ Directories:1 ends here

## [[file:utc.org::*Variable%20Names%20and%20Paths][Variable\ Names\ and\ Paths:1]]
locations = c("madison","wausau")

image.paths <- expand.grid(image.names,locations) %>% data.frame %>%
      mutate(img.paths = paste0(image.dirs,"/",Var2,Var1,".tif")) %>%
      .$img.paths

ratio.appendage <- "_ratio"
pca.appendage <- "_pca"
model.appendage = "_model"

feature.df.appendage <- "_featureDF"

ModelBuilding.appendage = "_modelBuildingDF"

tile.id.col.nm.for.grid.and.field.accuracy <- c("unq__ID", "Plot")
## Variable\ Names\ and\ Paths:1 ends here

## [[file:utc.org::*Patterns][Patterns:1]]
grid.pattern = "[a-zA-Z]{3}\\.[0-9]+m\\.[0-9]+" #I removed "_" from end. <2016-07-02 Sat>
texture.pattern = "stat-.*_window-.*_angle[-]+[0-9]+"
segmentation.pattern = "Pixel|N-[0-9]+_C-[0-9]+"
target.pattern = "all|grass|impervious|tree"
image.pattern = "[a-zA-Z]{5}[a-zA-Z]+"
model.pattern = "rf_prob|rf_resp|svm_resp"
## Patterns:1 ends here

## [[file:utc.org::*Texture%20Params][Texture\ Params:1]]
band.for.texture.appendage = "_ratio.nir"
window <- list(c(3,3), c(5,5), c(7,7))
statistics = list("homogeneity", "contrast", "correlation", "entropy")
shift = list(c(0,1),c(1,0),c(1,1),c(-1,1))

## band.for.texture.appendage = "_ratio.nir"
## window <- list(c(3,3))
## statistics = list("homogeneity")
## shift = list(c(0,1))

texture.params <- expand.grid(band.appendage = band.for.texture.appendage,window = window, statistics = statistics, shift = shift, stringsAsFactors = F)
## Texture\ Params:1 ends here

## [[file:utc.org::*Segmentation%20Params][Segmentation\ Params:1]]
segment.size <- c(rep(15,3), rep(20,3),rep(30,3),rep(45,3),rep(60,3),rep(100,3))
  compactness <- round(segment.size * c(.3, .5, .6))

  ## segment.size <- c(rep(30,1), rep(100,1))
  ## compactness <- segment.size * c(.5)

## segment.size <- 20
## compactness <- 12

  segment.params <- data.frame(compactness = compactness, segment.size = segment.size)
## Segmentation\ Params:1 ends here

## [[file:utc.org::*Input%20Shapefile%20DSNs%20and%20Layers][Input\ Shapefile\ DSNs\ and\ Layers:1]]
pca.region.dsn <- "../RD_PCA_Regions/"
pca.region.layer.appendage <- "_PCA_regions"

training.region.dsn <- "../RD_Training_Regions/"
training.region.layer.appendage <- "_TrainingPolygons"

grid.accuracy.region.dsn <- "../RD_Accuracy/Grids"
grid.accuracy.region.layer <- "Grids"

field.accuracy.region.dsn <- "../RD_Accuracy/FieldData"
field.accuracy.region.layer <- "fieldPoints"

accuracy.region.dsn <- c(grid.accuracy.region.dsn, field.accuracy.region.dsn)
accuracy.region.layer <- c(grid.accuracy.region.layer, field.accuracy.region.layer)
## Input\ Shapefile\ DSNs\ and\ Layers:1 ends here

## [[file:utc.org::*Derived%20Directories][Derived\ Directories:1]]
# make derived data directory
  derived.dir <- "../DD"

  dd.training.dirs <- str_c(derived.dir, "/",locations,"_Training")

  dd.pca.dirs <- str_c(derived.dir, "/",locations,pca.appendage)

  dd.accuracy.dirs <- str_c(derived.dir, "/",locations,"_Accuracy")

  dd.models.dirs <- paste0(derived.dir,"/",locations,"_Models")

  dd.accuracy.classified.dirs <- str_c(dd.accuracy.dirs, "/ClassifiedTiles")

derived.dirs <- c(derived.dir, dd.training.dirs, dd.pca.dirs, dd.accuracy.dirs, dd.models.dirs, dd.accuracy.classified.dirs)
## Derived\ Directories:1 ends here

## [[file:utc.org::*Make%20Derived%20Directories][Make\ Derived\ Directories:1]]
sapply(derived.dirs, FUN = function(x) dir.create(x))
## Make\ Derived\ Directories:1 ends here

## [[file:utc.org::*Define%20Derived%20Shapefile%20DSNs%20and%20Layers][Define\ Derived\ Shapefile\ DSNs\ and\ Layers:1]]
training.region.imageCRS.dsn <- str_c(derived.dir,"/reprojected.Training_Regions")

pca.region.imageCRS.dsn <- str_c(derived.dir,"/reprojected.PCA_Regions")

accuracy.region.imageCRS.dsn <- str_c(derived.dir,"/reprojected.Accuracy.Regions")


lapply(training.region.imageCRS.dsn, FUN = function(x) dir.create(x))
lapply(pca.region.imageCRS.dsn, FUN = function(x) dir.create(x))
lapply(accuracy.region.imageCRS.dsn, FUN = function(x) dir.create(x))
## Define\ Derived\ Shapefile\ DSNs\ and\ Layers:1 ends here

## [[file:utc.org::*number%20of%20cores][number\ of\ cores:1]]
cores <- detectCores()
cores <- 44
## number\ of\ cores:1 ends here

## [[file:utc.org::*CRS][CRS:1]]
utm16 <- CRS("+init=epsg:32616")
wtm <- CRS("+init=epsg:3071")
## CRS:1 ends here

## [[file:utc.org::*ASCII][ASCII:1]]
options(asciiType = "org")
## ASCII:1 ends here

## [[file:utc.org::*delete?][delete\?:1]]
#  band.names.wRatios <- c("blue","green","red","nir","b_ratio","g_ratio","r_ratio","n_ratio","ndvi")
#  pixel.feature.df.appendage = "_PixelFeatureDF"
#  segmentFeatureDF.appendage = "_SegmentFeatureDF.rds"
#  pca.model.name.appendage = "_pca.rds"

#     mad.grid.id.pattern = "mad.[0-9]+m.[0-9]+"
## delete\?:1 ends here

## [[file:utc.org::*Extract%20Name%20from%20path][Extract\ Name\ from\ path:1]]
# Use basename instead!
extract.name.from.path <- function(path) {
    str_extract(basename(path), "[A-Za-z0-9_]*.") %>%
        str_sub(.,1,-2)
}
## Extract\ Name\ from\ path:1 ends here

## [[file:utc.org::*Reproject%20Shapefile%20to%20Image%20Coordinate%20Reference%20System][Reproject\ Shapefile\ to\ Image\ Coordinate\ Reference\ System:1]]
Reproject_Shapefile_to_Image_CRS <- function(shapefile.dsn,
                                             shapefile.layer,
                                             image.path,
                                             shapefile.out.dsn) {
    r <- stack(image.path)
    shapefile <- readOGR(shapefile.dsn, shapefile.layer)
    shapefile.WimageCRS <- spTransform(shapefile, crs(r))
    image.name <- extract.name.from.path(image.path)
    shapefile.layer  <- str_c(image.name,"_",shapefile.layer)
    writeOGR(shapefile.WimageCRS, shapefile.out.dsn, shapefile.layer, driver = "ESRI Shapefile", overwrite =T)
}
## Reproject\ Shapefile\ to\ Image\ Coordinate\ Reference\ System:1 ends here

## [[file:utc.org::*get%20polygon%20ids][get\ polygon\ ids:1]]
IDs.SpatialPolygonsDataFrame <- function(x,...) {
    vapply(slot(x, "polygons"), function(x) slot(x, "ID"), "")
}
## get\ polygon\ ids:1 ends here

## [[file:utc.org::*Crop%20image%20to%20each%20Shapefile%20polygon][Crop\ image\ to\ each\ Shapefile\ polygon:1]]
Crop_image_to_each_Shapefile_polygon <- function(shapefile.dsn,
                                                 shapefile.layer,
                                                 image.path,
                                                 cores,
                                                 output.dir)  {
    image.name <- extract.name.from.path(image.path)
    shape <- readOGR(shapefile.dsn, str_c(image.name,"_",shapefile.layer))
    polygons <- as(shape, "SpatialPolygons")

    image <- stack(image.path)

    cl <- makeCluster(cores)
    registerDoParallel(cl)

    foreach (i = seq_along(polygons),
             .packages = c("raster")) %dopar% {
                 r <- image
                 r <- crop(r, polygons[i])
                 writeRaster(r, paste0(output.dir,"/",image.name,".",i,".tif"),
                             overwrite = T)
             }
closeAllConnections()
}
## Crop\ image\ to\ each\ Shapefile\ polygon:1 ends here

## [[file:utc.org::*Crop%20image%20to%20regions%20around%20shapefile%20points][Crop\ image\ to\ regions\ around\ shapefile\ points:1]]
# assign the polygon name to the points.
give_polygons_attributes_of_first_point_within <- function(points,
                                                           polygons){
    if (length(points@data$row) >1) {
        points <- points[points@data$row ==1 & points@data$col ==1 ,]
    }
    po <- gIntersects(points, polygons, byid=TRUE)
    out <- foreach(polygon.number = seq_along(polygons), .combine = "rbind") %do% {
        first.point.data <- points[po[polygon.number,],]@data %>%
                                                       slice(1)
        pd <- as(polygons[polygon.number], "SpatialPolygonsDataFrame")
        pd@data <- first.point.data
        pd
    }
}

Crop_image_to_regions_around_points_nameBygrid<- function(shapefile.dsn,
                                                          shapefile.layer,
                                                          image.path,
                                                          cores,
                                                          output.dir,
                                                          column.name = "unq__ID",
                                                          point.buffer.size = 4,
                                                          polygon.buffer.size = 15)  {
    image.name <- extract.name.from.path(image.path)
    points <- readOGR(shapefile.dsn,str_c(image.name,"_",shapefile.layer))
    box <- gBuffer(points, width = point.buffer.size, byid = F)
    box <- disaggregate(box)

    polygons <- as(box, "SpatialPolygons")

    polygons <- give_polygons_attributes_of_first_point_within(points,polygons)

    image <- stack(image.path)

    image.extent <- as(extent(image), "SpatialPolygons")
    proj4string(image.extent) <- proj4string(image)

    polygons.in.image <- foreach(i = seq_along(polygons),.combine = "c") %do% {
        gIntersects(polygons[i,],image.extent)
    }
closeAllConnections()

    polygons <- polygons[polygons.in.image,]

    cl <- makeCluster(cores)
    registerDoParallel(cl)

    foreach (k = seq_along(polygons),
             .packages = c("raster","rgeos")) %do% {
                 r <- image
                 poly <- gBuffer(polygons[k,],width = polygon.buffer.size, byid = T)
                 r <- crop(r, poly)
                 tile.id <- polygons@data[k,column.name]
                 writeRaster(r, paste0(output.dir,"/",image.name,".",tile.id,".tif"),
                             overwrite = T)
             }
closeAllConnections()
}

                                        #  shapefile.dsn = grid.accuracy.region.imageCRS.dsn
                                        #  shapefile.layer = grid.accuracy.region.layer,
                                        #  output.dir = image.cropped.to.grid.accuracy.dir


Crop_image_to_regions_around_points <- function(shapefile.dsn,
                                                shapefile.layer,
                                                image.path,
                                                cores,
                                                output.dir)  {

    points <- readOGR(shapefile.dsn, shapefile.layer)
    box <- gBuffer(points, width = 8)
    box <- disaggregate(box)

    polygons <- as(box, "SpatialPolygons")

    image <- stack(image.path)

    cl <- makeCluster(cores)
    registerDoParallel(cl)

    foreach (i = seq_along(polygons),
             .packages = c("raster")) %dopar% {
                 r <- image
                 r <- crop(r, polygons[i])
                 writeRaster(r, paste0(output.dir,"/",i,".tif"),
                             overwrite = T)
             }
closeAllConnections()
}
## Crop\ image\ to\ regions\ around\ shapefile\ points:1 ends here

## [[file:utc.org::*Save%20each%20band][Save\ each\ band:1]]
save_each_band <- function(tile.path, band.names) {
    tile <- stack(tile.path)
    names(tile) <- band.names
    tile.name <- str_sub(basename(tile.path),1,-5)
    writeRaster(tile, filename = paste0(dirname(tile.path),"/",tile.name,"_",names(tile), ".tif"), bylayer = T, format = "GTiff", overwrite = T)
}
## Save\ each\ band:1 ends here

## [[file:utc.org::*Add%20Texture][Add\ Texture:1]]
named.glcm <- function(tile.dir, tile.basename, band.appendage, window, statistics, shift, na_opt, na_val,...) {

    tile.path <- paste0(tile.dir, "/", tile.basename,band.appendage,".tif")
    x <- raster(tile.path)

    if (statistics == "correlation") {
        texture <- glcm(x, window = window, statistics = statistics, shift = shift, na_opt = na_opt, na_val = na_val)
        texture[texture == -Inf] <- -1
        texture[texture == Inf] <- 1
        texture[is.na(texture)] <- 1
    } else {
        texture <- glcm(x, window = window, statistics = statistics, shift = shift, na_opt = na_opt, na_val = na_val)
    }
    win.size <- paste0("window.",window[1])
    shift.dir <- paste0("angle.",atan(shift[1]/shift[2])*180/pi) # calc shift angle

    tile.dir <- dirname(tile.path)
    tile.name <- str_sub(basename(tile.path),1,-5)
    fn = paste0(tile.dir,"/", tile.basename,band.appendage, "_stat.", statistics, "_", win.size,"_",shift.dir,".tif")
    writeRaster(texture, fn, overwrite = T)
}

calc.texture <- function(texture.params.df,
                         tile.dir,
                         tile.basename) {

    texture <- mapply(named.glcm,
                      tile.dir = tile.dir,
                      tile.basename = tile.basename,
                      band.appendage = texture.params.df$band.appendage,
                      window = texture.params.df$window,
                      statistics = texture.params.df$statistics,
                      shift = texture.params.df$shift,
                      na_opt = "ignore",
                      na_val = NA)
}
## Add\ Texture:1 ends here

## [[file:utc.org::*Make%20new%20ratio%20bands%20from%20image][Make\ new\ ratio\ bands\ from\ image:1]]
calc_ratios <- function(tile.path, band.names, ratio.bands, scale200 = T) {
    tile <- stack(tile.path)
    names(tile) <- band.names

    ratios <- tile[[ratio.bands,drop = F]] / sum(tile)

    if (scale200 == T) {
        ratios <- ratios * 200
    }

    tile.name <- str_sub(basename(tile.path),1,-5)
    names(ratios) <- paste0(tile.name,"_ratio.",ratio.bands)
    writeRaster(ratios, filename= paste0(dirname(tile.path),"/",names(ratios),".tif"),
                bylayer = T, format= "GTiff", overwrite = T,
                datatype = 'INT1U')
}

calc_ndvi <- function(tile.path, band.names, ndvi_appendage = "_ndvi", scale200 = T) {

    tile <- stack(tile.path)
    names(tile) <- band.names

    ndvi <- (tile[["nir"]] - tile[["red"]]) /  (tile[["nir"]] + tile[["red"]])

    ndvi [ndvi < 0] <- 0

    if (scale200 == T) {
        ndvi <- ndvi * 200
    }

    tile.dir <- dirname(tile.path)
    tile.name <- str_sub(basename(tile.path),1,-5)
    writeRaster(ndvi, filename=paste0(tile.dir,"/",tile.name,ndvi_appendage,".tif"), bylayer=TRUE,format="GTiff", overwrite = T,datatype = 'INT1U')
    return(ndvi)
}
## Make\ new\ ratio\ bands\ from\ image:1 ends here

## [[file:utc.org::*Make%20Pixel%20Feature%20DF][Make\ Pixel\ Feature\ DF:1]]
save.pixel.feature.df <- function(tile.dir,
                                  tile.name,
                                  feature.pattern,
                                  feature.df.append = feature.df.appendage ) {
    s <- stack(list.files(tile.dir, pattern = paste0(tile.name,feature.pattern), full.names = T))
    names(s) <- sub(x = names(s), pattern = paste0("(",tile.name,"_)"), replacement = "")
    s.df <- as.data.frame(s, xy = T)
    saveRDS(s.df, file = paste0(tile.dir, "/", tile.name, "_Pixel",feature.df.append, ".rds"))
}
## Make\ Pixel\ Feature\ DF:1 ends here

## [[file:utc.org::*Image%20PCA][Image\ PCA:1]]
pca.transformation <- function(tile.dir,
                               image.name,
                               tile.name,
                               loc,
                               feature.pattern = "_(blue|green|red|nir|ratio.blue|ratio.green|ratio.red|ratio.nir|ndvi).tif",
                               pca.append = pca.appendage,
                               out.image.appendage = pca.appendage,
                               comps.to.use = c(1,2,3),
                               pca.dir = dd.pca.dir) {

    s <- stack(list.files(tile.dir, pattern = paste0(tile.name,feature.pattern), full.names = T))
    names(s) <- sub(x = names(s), pattern = ".*_", replacement = "")

    pca.model <- readRDS(str_c(pca.dir,"/",loc,image.name,pca.append,".rds"))

    r <- predict(s, pca.model, index = comps.to.use)

    min.r <- getRasterMin(r)
    max.r <- getRasterMax(r)
    rescaled.r <- rescale.0.254(r, min.r, max.r)

    out.path <- str_c(tile.dir, "/", tile.name, out.image.appendage, ".tif")
    writeRaster(rescaled.r, filename = out.path, overwrite=TRUE, datatype = 'INT1U', bylayer = F)
}


getRasterMin <- function(t) {
    return(min(cellStats(t, stat = "min")))
}

getRasterMax <- function(t) {
    return(max(cellStats(t, stat = "max")))
}

rescale.0.254 <- function(raster,
                          min,
                          max) {
                              (raster - min)/(max-min) * 254
}

## image.pca <- function(image.name,
##                       pca.model.name.append = pca.model.name.appendage,
##                       tile.dir,
##                       tile.name,
##                       in.image.appendage = ratio.tile.name.append,
##                       out.image.appendage = pca.tile.name.append,
##                       band.names = c("blue","green","red","nir","b_ratio","g_ratio","r_ratio","n_ratio","ndvi"),
##                       comps.to.use = c(1,2,3),
##                       pca.dir = dd.pca.dir) {


##     out.path <- str_c(tile.dir, "/", tile.name, out.image.appendage, ".tif")

##     s <- stack(str_c(tile.dir, "/", tile.name, in.image.appendage,".tif"))
##     names(s) <- band.names

##     pca.model <- readRDS(str_c(pca.dir,"/",image.name,pca.model.name.append))

##     r <- predict(s, pca.model, index = comps.to.use)

##     min.r <- getRasterMin(r)
##     max.r <- getRasterMax(r)
##     rescaled.r <- rescale.0.255(r, min.r, max.r)
##     writeRaster(rescaled.r, filename = out.path, overwrite=TRUE, datatype = 'INT1U')
## }


make.and.save.pca.transformation <- function(image.dir,
                                             image.name,
                                             location,
                                             pca.append = pca.appendage,
                                             max.sample.size = 10000,
                                             core.num = cores,
                                             feature.pattern = ".*_(blue|green|red|nir|ratio.blue|ratio.green|ratio.red|ratio.nir|ndvi).tif",
                                             ratio.appendage = ratio.tile.name.append) {

    tile.paths <- list.files(image.dir, pattern = paste0(image.name,feature.pattern), full.names = T)

    tile.names <- str_match(tile.paths,"(.*\\.[0-9]+)_.*")[,2] %>%  unique() # get the image names of pca regions

    cl <- makeCluster(cores)
    registerDoParallel(cl)

    sr <- foreach (tile.name = tile.names, .packages = c("stringr","raster"), .combine ="rbind") %dopar% {
        t.names <- str_extract(tile.paths, paste0(".*",tile.name,".*")) %>% na.omit()
        tile <- stack(t.names)
        names(tile) <- sub(x = names(tile), pattern = ".*_", replacement = "")
        samp <- sampleRandom(tile, ifelse(ncell(tile) > max.sample.size ,max.sample.size, ncell(tile)))
        colnames(samp) <- names(tile)
        samp
    }
closeAllConnections()

                                        # Perform PCA on sample
    pca <- prcomp(sr, scale = T)
    saveRDS(pca,paste0(image.dir,"/",location,image.name,pca.append,".rds"))
    return(pca)
}


## make.and.save.pca.transformation <- function(image.dir,
##                                              image.name,
##                                              pca.model.name.append = "_pca.rds",
##                                              max.sample.size = 10000,
##                                              core.num = cores,
##                                              band.names = c("blue","green","red","nir","b_ratio","g_ratio","r_ratio","n_ratio","ndvi"),
##                                              ratio.appendage = ratio.tile.name.append) {
##     tile.paths <- list.files(str_c(image.dir), pattern = paste0("*",ratio.appendage), full.names = T)

##     tile.names <- basename(tile.paths)

##     cl <- makeCluster(core.num)
##     registerDoParallel(cl)

##     sr <- foreach (i = seq_along(tile.names), .packages = c("raster"), .combine ="rbind") %dopar% {
##         tile <- stack(tile.paths[i])
##         s <- sampleRandom(tile, ifelse(ncell(tile) > max.sample.size ,max.sample.size, ncell(tile)))
##     }

##     colnames(sr) <- band.names

##                                         # Perform PCA on sample
##     pca <- prcomp(sr, scale = T)
##     saveRDS(pca,paste0(image.dir,"/",image.name,pca.model.name.append))

##     return(pca)
## }


image.pca.forWholeState <- function(pca.model.name.append = pca.model.name.appendage,
                                    tile.dir,
                                    tile.name,
                                    in.image.appendage = ratio.tile.name.append,
                                    out.image.appendage = pca.tile.name.append,
                                    band.names = c("blue","green","red","nir","b_ratio","g_ratio","r_ratio","n_ratio","ndvi"),
                                    comps.to.use = c(1,2,3),
                                    pca.transform) {


    out.path <- str_c(tile.dir, "/", tile.name, out.image.appendage, ".tif")

    s <- stack(str_c(tile.dir, "/", tile.name, in.image.appendage,".tif"))
    names(s) <- band.names

    r <- predict(s, pca.transform, index = comps.to.use)

    min.r <- getRasterMin(r)
    max.r <- getRasterMax(r)
    rescaled.r <- rescale.0.254(r, min.r, max.r)
    writeRaster(rescaled.r, filename = out.path, overwrite=TRUE, datatype = 'INT1U')
}



## image.dir <- image.cropped.to.training.dir
## image.name <- 9
##                         in.image.appendage = ratio.tile.name.append
##                         out.image.appendage = pca.tile.name.append
##                         band.names = c("blue","green","red","nir","b_ratio","g_ratio","r_ratio","n_ratio","ndvi")
##                         max.sample.size = 10000
##                         comps.to.use = c(1,2,3)

##       out.path <- str_c(image.dir, "/", image.name, out.image.appendage, ".tif")

##       s <- stack(str_c(image.dir, "/", image.name, in.image.appendage,".tif"))
##       names(s) <- band.names

##       sr <- sampleRandom(s, ifelse(ncell(s) > max.sample.size, max.sample.size, ncell(s)))
##       pca <- prcomp(sr, scale = T)

##       r <- predict(s, pca, index = comps.to.use)

##       min.r <- getRasterMin(r)
##       max.r <- getRasterMax(r)
##       rescaled.r <- rescale.0.255(r, min.r, max.r)
##       writeRaster(rescaled.r, filename = out.path, overwrite=TRUE, datatype = 'INT1U')









                                        # Function takes raster stack, samples data, performs pca and returns stack of first n_pcomp bands
                                        ## predict_pca_wSampling_parallel <- function(stack, sampleNumber, n_pcomp, nCores = detectCores()-1) {
                                        ##     sr <- sampleRandom(stack,sampleNumber)
                                        ##     pca <- prcomp(sr, scale=T)
                                        ##     beginCluster()
                                        ##     r <- clusterR(stack, predict, args = list(pca, index = 1:n_pcomp))
                                        ##     endCluster()
                                        ##     return(r)
                                        ## }
## Image\ PCA:1 ends here

## [[file:utc.org::*Segment%20image][Segment\ image:1]]
segment.multiple <- function(tile.dir,
                             tile.name,
                             image.name,
                             segment.params.df,
                             krusty  = T) {
    segments <- mapply(segment,
                       tile.dir = tile.dir,
                       image.name = image.name,
                       tile.name = tile.name,
                       compactness = segment.params.df$compactness,
                       segment.size = segment.params.df$segment.size,
                       krusty = krusty)
      }

segment  <- function(tile.dir,
                     image.name,
                     tile.name,
                     compactness,
                     segment.size,
                     krusty = T) {
    pixel_size <- ifelse(image.name == "NAIP", 1, 1.5)
    compactness <- if(image.name == "NAIP") compactness else round(2/3*compactness)
    if (krusty == T) {
        system(paste("/home/erker/.conda/envs/utc/bin/python","fia_segment_cmdArgs.py",pixel_size,segment.size,compactness,tile.name,tile.dir))
    } else {
        system(paste("python","fia_segment_cmdArgs.py",pixel_size,segment.size,compactness,tile.name,tile.dir))
    }
}
## Segment\ image:1 ends here

## [[file:utc.org::*add.features][add\.features:1]]
add.features <- function(tile.dir,
                         tile.name,
                         band.names,
                         ndvi = T,
                         ratio.bands,
                         texture = T,
                         texture.params.df) {

    til.path <- paste0(tile.dir,"/",tile.name,".tif")
    til <- stack(til.path)
    names(til) <- band.names

    save_each_band(tile.path = til.path,
                   band.names = band.names)

    if (ndvi == T) {
        calc_ndvi(tile.path = til.path,
                  band.names = band.names)
    }

    if (length(ratio.bands > 0)) {
        calc_ratios(tile.path = til.path,
                    band.names = band.names,
                    ratio.bands = ratio.bands)
    }

    if (texture == T) {
        calc.texture(texture.params.df = texture.params.df,
                     tile.dir = tile.dir,
                     tile.basename = tile.name)
    }
}
## add\.features:1 ends here

## [[file:utc.org::*segment%20Feature%20DF][segment\ Feature\ DF:1]]
make.segment.feature.df.foreach.segmentation <- function(tile.dir,
                                                         tile.name,
                                                         feature.pattern,
                                                         segmentation.pattern = "_N-[0-9]+_C-[0-9]+.*") {

    segmentation.files <-  list.files(tile.dir, pattern = paste0(tile.name,segmentation.pattern))
    segmentation.param.appendages <- str_match(segmentation.files,paste0(tile.name,"(_.*).tif"))[,2] %>% na.omit()


    out <- lapply(X = segmentation.param.appendages, FUN = function(segmentation.param.appendage) {
        make.segment.feature.df(tile.dir = tile.dir,
                                tile.name = tile.name,
                                segmentation.param.appendage = segmentation.param.appendage,
                                fea.pattern = feature.pattern)
    })

}


make.segment.feature.df <- function(tile.dir,
                                    tile.name,
                                    segmentation.param.appendage,
                                    fea.pattern,
                                    feature.df.append = feature.df.appendage) {

    fea <- stack(list.files(tile.dir, pattern = paste0(tile.name,fea.pattern), full.names = T))
    names(fea) <- sub(x = names(fea), pattern = "(madisonNAIP|madisonPanshpSPOT).*?_", replacement = "")

    seg.path <- paste0(tile.dir,"/",tile.name,segmentation.param.appendage, ".tif")
    seg <- raster(seg.path)

                                        # Create a data_frame where mean and variances are calculated by zone
    x <- as.data.frame(fea, xy = T)
    s <- as.data.frame(seg)
    colnames(s) <- "segment"
    r <- bind_cols(x,s)
    r2 <- r %>%
        group_by(segment)

    mean.and.sd <- r2 %>%
        summarize_each(funs(mean(.,na.rm = T), sd(., na.rm = T))) %>%
        select(-x_mean, -x_sd, -y_mean, -y_sd)

    tile.name.df = data.frame(tile.name = rep(tile.name, nrow(mean.and.sd)))

    out <- bind_cols(mean.and.sd, tile.name.df)


    names <- colnames(out)
    names <- str_replace(names, "\\(",".")
    names <- str_replace(names, "\\)",".")
    names <- str_replace(names, "\\:",".")
    colnames(out) <- names
    saveRDS(out, file = paste0(tile.dir,"/",tile.name,segmentation.param.appendage,feature.df.append,".rds"))
    out
}



                                        #  make.segment.feature.df(dd.training.dir, "madisonNAIP.1", segmentation.param.appendage = "_N-100_C-10", feature.pattern = feature.pattern)
## segment\ Feature\ DF:1 ends here

## [[file:utc.org::*make.feature.df][make\.feature\.df:1]]
make.feature.df <- function(tile.dir,
                            image.name,
                            tile.name,
                            band.names,
                            ndvi = T,
                            ratio.bands,
                            texture = T,
                            texture.params.df,
                            feature.pattern = "_(blue|green|red|nir|ratio.blue|ratio.green|ratio.red|ratio.nir|ndvi|ratio.nir_stat\\.\\w+_window\\.\\d+_angle\\..?\\d+).tif",
                            pixel.df,
                                        #                              pca.features = c("blue","green","red","nir","ndvi","ratio.blue","ratio.green","ratio.red","ratio.nir"),
                            pca.features = c("blue","green","red","nir"),
                            pca.location,
                            segmentation = T,
                            segment.params.df,
                            segment.feature.df = T,
                            using.krusty = T) {

    add.features(tile.dir,
                 tile.name,
                 band.names,
                 ndvi = T,
                 ratio.bands,
                 texture = T,
                 texture.params.df)

    message ( tile.name,"features added")

    if (pixel.df ==T) {

        save.pixel.feature.df(tile.dir = tile.dir,
                              tile.name = tile.name,
                              feature.pattern)}

    message("pixel feature df saved")

    pca.transformation(tile.dir = tile.dir,
                       tile.name = tile.name,
                       image.name = image.name,
                       loc = pca.location)

    message("pca done")

    if (segmentation == T) {

        segment.multiple(tile.dir = tile.dir,
                         tile.name = tile.name,
                         image.name = image.name,
                         segment.params.df = segment.params.df,
                         krusty = using.krusty)}
    message("segmentation done")
    if (segment.feature.df == T) {

        make.segment.feature.df.foreach.segmentation(tile.dir = tile.dir,
                                                     tile.name = tile.name,
                                                     feature.pattern = feature.pattern)}


}
## make\.feature\.df:1 ends here

## [[file:utc.org::*polygonize%20segment%20raster%20with%20gdal%20and%20add%20Class%20to%20shapefile][polygonize\ segment\ raster\ with\ gdal\ and\ add\ Class\ to\ shapefile:1]]
gdal_polygonizeR <- function(x, outshape=NULL, gdalformat = 'ESRI Shapefile',
                             pypath=NULL, readpoly=TRUE, quiet=TRUE) {
    if (isTRUE(readpoly)) require(rgdal)
    if (is.null(pypath)) {
        pypath <- Sys.which('gdal_polygonize.py')
    }
    if (!file.exists(pypath)) stop("Can't find gdal_polygonize.py on your system.")
    owd <- getwd()
    on.exit(setwd(owd))
    setwd(dirname(pypath))
    if (!is.null(outshape)) {
        outshape <- sub('\\.shp$', '', outshape)
        f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
        if (any(f.exists))
            stop(sprintf('File already exists: %s',
                         toString(paste(outshape, c('shp', 'shx', 'dbf'),
                                        sep='.')[f.exists])), call.=FALSE)
    } else outshape <- tempfile()
    if (is(x, 'Raster')) {
        require(raster)
        writeRaster(x, {f <- tempfile(fileext='.asc')})
        rastpath <- normalizePath(f)
    } else if (is.character(x)) {
        rastpath <- normalizePath(x)
    } else stop('x must be a file path (character string), or a Raster object.')
    system2('python', args=(sprintf('"%1$s" "%2$s" -f "%3$s" "%4$s.shp"',
                                    pypath, rastpath, gdalformat, outshape)))
    if (isTRUE(readpoly)) {
        shp <- readOGR(dirname(outshape), layer = basename(outshape), verbose=!quiet)
        return(shp)
    }
    return(NULL)
}


polygonize.and.add.Class <- function(image.dir,
                                     image.name,
                                     segment.appendage = segment.tile.name.append,
                                     no.class = "N") {
    seg <- raster(paste0(image.dir,"/",image.name,segment.appendage,'.tif'))
    segPoly <- gdal_polygonizeR(seg)
    segPoly$Class <- no.class
    writeOGR(obj = segPoly,
             dsn = paste0(image.dir,"/",image.name),
             layer = paste0(image.name,segment.appendage),
             driver = "ESRI Shapefile",
             overwrite = T)
}
## polygonize\ segment\ raster\ with\ gdal\ and\ add\ Class\ to\ shapefile:1 ends here

## [[file:utc.org::*Create%20ModelBuilding%20dataframe][Create\ ModelBuilding\ dataframe:1]]
getSegment.class.and.features.Within.Polygon<-function(SegmentFeatureDF,
                                                       training.sp,
                                                       seg.tiles.dir,
                                                       seg.params){
    seg.files <- list.files(seg.tiles.dir, pattern = str_c(seg.params,".tif$"), full.names = T)
                                        # find number of pixels in each segment
    n.pixels.per.seg <- foreach(seg.file = seg.files, .combine = "rbind") %do% {
        seg <- raster::stack(seg.file)
        s.df <- as.data.frame(seg) %>%
            gather(key = image.name, value = segment.id) %>%
            group_by(segment.id, image.name) %>%
            summarize(n.pixels.per.seg = n())
    }
closeAllConnections()
                                        # find number of pixels in each segment are in a polygon
    n.pixels.per.seg.in.polygon <- foreach(seg.file = seg.files, .combine = "rbind") %do% {

        seg <- raster::stack(seg.file)
        ei <- as(extent(seg), "SpatialPolygons")

        if(gIntersects(ei, as(training.sp,"SpatialPolygons"))) {

            a <- raster::extract(seg, as(training.sp,"SpatialPolygons"), df = T)

            a <- a %>%
                gather(key = image.name, value = segment.id, -ID) %>%
                rename(polygon.id = ID) %>%
                group_by(polygon.id, image.name, segment.id) %>%
                summarize(n.pixels.per.seg.in.polygon = n())
        }
    }
                                        # get pct of segment in a polygon,
                                        # filter segments that have more than 50%,
                                        #join Class information from polygons
    if(!is.null(n.pixels.per.seg.in.polygon)) {

                                        #add 1 because the id from raster's extract is just the order of the polygons
        training.sp@data$polygon.id <- as.numeric(IDs.SpatialPolygonsDataFrame(training.sp))+1

        n.pixels <- left_join(n.pixels.per.seg.in.polygon,n.pixels.per.seg) %>%
            mutate(pct.seg.in.polygon = n.pixels.per.seg.in.polygon/n.pixels.per.seg) %>%
            filter(pct.seg.in.polygon >= .5) %>%
            left_join(.,training.sp@data) %>%
            ungroup() %>%
            mutate(segment = segment.id)


        n.pixels$image.name <- str_match(n.pixels$image.name, "(.*?\\.[0-9]+).*")[,2]

        out <- left_join(n.pixels, SegmentFeatureDF) %>%
            distinct() %>%
            dplyr::select(-id,
                          -segment,
                          -segment.id,
                          -image.name,
                          -image.name,
                          -tile.name,
                          -polygon.id,
                          -n.pixels.per.seg,
                          -n.pixels.per.seg.in.polygon,
                          -pct.seg.in.polygon)        %>%
            filter(complete.cases(.))

        out
    }
}

                                        # returns dataframe of values of pixels within polygon
getPixel.Class.and.Coords.Within.Polygon <- function(PixelFeatureDF,
                                                     training.sp) {
    xy <- dplyr::select(PixelFeatureDF,x,y) %>% data.frame
    PixelFeatureDF <- data.frame(PixelFeatureDF)
    coordinates(PixelFeatureDF) <- xy
    proj4string(PixelFeatureDF) <- utm16

    training.sp <- spTransform(training.sp,utm16)

    pts.in.poly <- sp::over(PixelFeatureDF,training.sp)
    PixelFeatureDF@data <- cbind(PixelFeatureDF@data, pts.in.poly)
    PixelFeatureDF <- PixelFeatureDF[which(complete.cases(pts.in.poly)),]
    PixelFeatureDF@data
}
## Create\ ModelBuilding\ dataframe:1 ends here

## [[file:utc.org::*untuned][untuned:1]]
Build.and.Save.models <- function(dir = dd.training.dir,
                                  modelBuildingData = ModelBuildingRDS,
                                  models.dir = dd.models.dir,
                                  image.name,
                                  location,
                                  model.append = model.appendage){

    dat <- readRDS(paste0(dir,"/",modelBuildingData)) %>%
        as.data.frame() %>%
        filter(complete.cases(.))

    seg.p <- str_extract(modelBuildingData, segmentation.pattern)

    names <- colnames(dat)
    names <- str_replace(names, "\\(",".")
    names <- str_replace(names, "\\)",".")
    names <- str_replace(names, "\\:",".")
    colnames(dat) <- names

                                        # Create Tasks
    tsk <- makeClassifTask(id = paste0(location,image.name,"_all"), data = dat, target = "Class")

                                        # Make Learners
    RF_prob <- makeLearner(id = "rf_prob","classif.randomForest", predict.type = "prob", fix.factors.prediction = TRUE)
                                        #      RF_response <- makeLearner(id = "rf_resp", "classif.randomForest", predict.type = "response", fix.factors.prediction = TRUE)
    SVM_response <- makeLearner(id = "svm_resp", "classif.svm", predict.type = "response", fix.factors.prediction = TRUE)

                                        #      learner.list <- list(RF_prob = RF_prob, RF_response = RF_response, SVM_response = SVM_response)
    learner.list <- list(RF_prob = RF_prob, SVM_response = SVM_response)

                                        # Train Learners on Tasks, Make models
                                        #         cl<-makeCluster(cores)
                                        #         registerDoParallel(cl)

    models <- foreach(lnr = learner.list) %do% {
        mod <- train(lnr, tsk)
        mod
        saveRDS(mod, file = paste0(models.dir,"/",location,image.name,"_",seg.p, "_",lnr$id,"_Untuned",model.append, ".rds"))
    }
closeAllConnections()
}
## untuned:1 ends here

## [[file:utc.org::*features%20selected][features\ selected:1]]
print.feature.importance <- function(dir = dd.training.dir,
                                  modelBuildingData = ModelBuildingRDS,
                                  image.name,
                                  location,
                                  feature.importance.methods = c("information.gain","chi.squared")) {


    dat <- readRDS(paste0(dir,"/",modelBuildingData)) %>%
        as.data.frame() %>%
        filter(complete.cases(.))

    seg.p <- str_extract(modelBuildingData, segmentation.pattern)

    names <- colnames(dat)
    names <- str_replace(names, "\\(",".")
    names <- str_replace(names, "\\)",".")
    names <- str_replace(names, "\\:",".")
    colnames(dat) <- names

                                        # Create Tasks
    tsk <- makeClassifTask(id = paste0(location,image.name,"_all"), data = dat, target = "Class")

    fv = generateFilterValuesData(tsk, method = feature.importance.methods)

    fv$data %>% arrange_(feature.importance.methods[1])
}
## features\ selected:1 ends here

## [[file:utc.org::*features%20selected][features\ selected:2]]
Build.and.Save.FeatureSelected.models <- function(dir = dd.training.dir,
                                  modelBuildingData = ModelBuildingRDS,
                                  models.dir = dd.models.dir,
                                  image.name,
                                  location,
                                  model.append = model.appendage){

    dat <- readRDS(paste0(dir,"/",modelBuildingData)) %>%
        as.data.frame() %>%
        filter(complete.cases(.))

    seg.p <- str_extract(modelBuildingData, segmentation.pattern)

    names <- colnames(dat)
    names <- str_replace(names, "\\(",".")
    names <- str_replace(names, "\\)",".")
    names <- str_replace(names, "\\:",".")
    colnames(dat) <- names


                                        # Create Tasks
    tsk <- makeClassifTask(id = paste0(location,image.name,"_all"), data = dat, target = "Class")


                                        # Make Learners
    RF_prob <- makeLearner(id = "rf_prob","classif.randomForest", predict.type = "prob", fix.factors.prediction = TRUE)
    SVM_response <- makeLearner(id = "svm_resp", "classif.svm", predict.type = "response", fix.factors.prediction = TRUE)


                                        # make filter wrappers
    RF_prob_fil <- makeFilterWrapper(RF_prob, fw.method = "chi.squared")
    SVM_response_fil <- makeFilterWrapper(SVM_response, fw.method = "chi.squared")


                                        # make tune wrapper for feature selection
    # inner
    ps = makeParamSet(makeDiscreteParam("fw.abs", values = seq_len(getTaskNFeats(tsk))))
    ctrl = makeTuneControlGrid()
    inner = makeResampleDesc("CV", iter = 2)

    RF_prob_tunfil = makeTuneWrapper(RF_prob_fil, resampling = inner, par.set = ps, control = ctrl, show.info = FALSE)

    SVM_response_tunfil = makeTuneWrapper(SVM_response_fil, resampling = inner, par.set = ps, control = ctrl, show.info = FALSE)

    learner.list <- list(RF_prob_tunfil = RF_prob_tunfil, SVM_response_tunfil = SVM_response_tunfil)

    # outer
    outer = makeResampleDesc("Subsample", iter = 3)
    res = benchmark(tasks = tsk, learners = learner.list, resampling = outer, show.info = FALSE)

res


    models <- foreach(lnr = learner.list) %do% {
        mod <- train(lnr, tsk)
        mod
        saveRDS(mod, file = paste0(models.dir,"/",location,image.name,"_",seg.p, "_",lnr$id,"_FeatureSelected",model.append, ".rds"))
    }
closeAllConnections()
}
## features\ selected:2 ends here

## [[file:utc.org::*tuned][tuned:1]]

## tuned:1 ends here

## [[file:utc.org::*Classify%20Raster][Classify\ Raster:1]]
classify.segmented.raster <- function(segment.feature.df.dir,
                                      segment.dir,
                                      model.dir,
                                      model.name.rds = "models",
                                      segment.feature.appendage = segment.feature.df.name.append,
                                      segmentation.appendage = segment.tile.name.append,
                                      segmentation.prms,
                                      classify.out.dir,
                                      tile.name = i) {
    df <- readRDS(paste0(segment.feature.df.dir,"/",tile.name,segment.feature.appendage))
    mod <-readRDS(paste0(model.dir,"/",model.name.rds))
                                        #    umod <- unlist(models, recursive = F)
    seg.path <- paste0(segment.dir,"/",tile.name,segmentation.appendage)
    seg <- raster(seg.path)
                                        #       dfRowsWithNA <- which(is.na(df[,2]))
    complete.df <- df[complete.cases(df),] # svm can't predict with NAs

    pred <- predict(mod, newdata = complete.df)
    response <- factor(as.character(pred$data$response), levels = c("g","i","t","o"))
    m <- cbind(zone = complete.df$segment, response)
    m <- left_join(as.data.frame(df["segment"]), as.data.frame(m), by = c("segment" = "zone"))
    seg.df <- as.data.frame(seg, xy = T)
    names(seg.df)[3] <- "segment"
    seg.df <- left_join(seg.df, m)
    seg.df$response <- mapvalues(seg.df$response, from = c(1,2,3,4), to = c("g","i","t","o"))
    seg.df$response <- factor(seg.df$response)
    r <- seg
    values(r) <- seg.df$response

                                        #        x <- data.frame(ID = 1:4, LandCover = c("G","I","T","O")) %>%
                                        #            filter(LandCover %in% levels(factor(response)))
                                        #        levels(r) <- x
                                        # Removing Probability layer because can't have attributes with it.  When I do final classifcaiton I should add back in.

    ## if (ncol(pred$data) > 2) {
    ##     prob <- (pred$data[,grep("prob.*", x = colnames(pred$data))]) # get columns that contain probabilities
    ##     ProbOfClass <- apply(prob, MARGIN = 1, FUN = max)
    ##     m <- cbind(segment = df$segment, ProbOfClass)
    ##     m <- left_join(as.data.frame(df["segment"]), as.data.frame(m))
    ##     p <- reclassify(seg, m)
    ##     r <- stack(r,p)
    ## }
    path <- paste0(dd.accuracy.classified.dir,"/",tile.name,"_",segmentation.prms,"_",mod$task.desc$id,"_",mod$learner$id,".tif")
    writeRaster(r, path, overwrite=TRUE)
    print(path)

}




classify.pixel.raster <- function(tile.dir = dd.accuracy.dir,
                                  tile.name,
                                  pixelFeatureDF.appendage = pixel.feature.df.appendage,
                                  model.dir = Models.dir,
                                  model.rds,
                                  seg.prms = "Pixel") {
    ras <- stack(str_c(tile.dir,"/",tile.name,".tif"))
    pix.mod <- readRDS(str_c(model.dir,"/",model.rds))
                                        #      pix.umods <- unlist(pix.mods, recursive = F)

    pix.feature.df <- readRDS(str_c(tile.dir,"/",tile.name,"_",seg.prms,pixelFeatureDF.appendage,".rds"))

    if(!is.null(pix.feature.df$y)) {
        pix.feature.df <- dplyr::select(pix.feature.df, -x, -y)
    }

                                        # I set NA's to 0 here.  Not the best choice.  Not sure why they exist.
                                        # Maybe because pca transform
                                        # imputing to mean would probably be better

    pix.feature.df <- as.matrix(pix.feature.df)

    pix.feature.df[which(is.na(pix.feature.df))] <- 0

    pix.feature.df <- as.data.frame(pix.feature.df)

    pred <- predict(pix.mod, newdata = pix.feature.df)
    a <- ras[[1]]
    values(a) <- pred$data$response
    path <- paste0(dd.accuracy.classified.dir,"/",tile.name,"_",seg.prms,"_",pix.mod$task.desc$id,"_",pix.mod$learner$id,".tif")
    writeRaster(a, path, overwrite = T)
    print(path)
}
## Classify\ Raster:1 ends here

## [[file:utc.org::*Calculate%20Percent%20Cover%20in%20Classified%20Tiles][Calculate\ Percent\ Cover\ in\ Classified\ Tiles:1]]
get.prcnt.class <- function(points,r) {
    r <- crop(r,points)  # should I do a mask instead??
    g <- cellStats(r == 1, stat = sum)
    im <- cellStats(r == 2, stat = sum)
    tr <- cellStats(r == 3, stat = sum)
    o <-  cellStats(r == 4, stat = sum)
    totC <- ncell(r)
    return(c(pct_g_pred = g/totC, pct_i_pred = im/totC, pct_t_pred = tr/totC, pct_o_pred = o/totC))
}


get_area_convexHull <- function(points) {
    ch <- chull(coordinates(points))
    coords <- coordinates(points)[c(ch,ch[1]),]
    poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)),ID = 1)))
    gArea(poly)
}



calculate.percent.cover.in.classified.tile <- function(pts,
                                                       tile.dir = dd.accuracy.classified.dir,
                                                       tile.pth,
                                                       n.rows.and.columns.subset,
                                                       mod = 1,
                                                       mad.grid.id.pattern = "mad-[0-9]+m-[0-9]+",
                                                       grid.pattern = "[a-zA-Z]{3}-[0-9]+m-[0-9]+_",
                                                       image.pattern = "[a-zA-Z]{5}[a-zA-Z]+",
                                                       target.pattern = "all|grass|impervious|tree",
                                                       model.pattern = "rf_prob|rf_resp|svm_resp",
                                                       seg.prms = "N-[0-9]+_C-[0-9]+|Pixel"
                                                       ) {
    tile.nm <- basename(tile.pth)


    pts.sub <- pts@data  %>%
        filter.by.row.and.col(.,n.rows.and.columns.subset, mod = mod)

    coordinates(pts.sub) <- ~ crds_x1 + crds_x2

    proj4string(pts.sub) <- utm16
    tile.unique.name <- str_extract(tile.pth, mad.grid.id.pattern)
    pts.at.grid <- pts.sub[which(pts.sub@data$unq__ID == tile.unique.name),]
    tile <- raster(tile.pth, proj4string = "+init:epsg=32616")

    area.pts <- get_area_convexHull(pts.at.grid)

    if(!is.null(raster::intersect(extent(tile),bbox(pts.at.grid)))) {

        get.prcnt.class(pts.at.grid,tile) %>%
            t() %>%
            as.data.frame() %>%
            mutate(grid.tile.target.model = tile.nm,
                   grid = str_sub(str_extract(grid.tile.target.model, grid.pattern),1,-2),
                   image =  str_extract(grid.tile.target.model, image.pattern),
                   target.cover = str_extract(grid.tile.target.model, target.pattern),
                   model =  str_extract(grid.tile.target.model, model.pattern),
                   n.points = n.rows.and.columns.subset * n.rows.and.columns.subset,
                   area = area.pts,
                   seg.params = str_extract(grid.tile.target.model, seg.prms),
                   target.type = ifelse(target.cover == "all", "multinomial", "binomial"))
    }
}
## Calculate\ Percent\ Cover\ in\ Classified\ Tiles:1 ends here

## [[file:utc.org::*Calculate%20Percent%20Cover%20of%20Grids,%20subsetted][Calculate\ Percent\ Cover\ of\ Grids\,\ subsetted:1]]
filter.by.row.and.col <- function(df,nrow.and.col, mod) {
    nrow <-df %>%
        group_by(unq__ID) %>%
        summarize(nrow = max(row))

    df <- left_join(df,nrow)

    df %>%
        filter(nrow >= nrow.and.col,   # remove grids that have fewer than the number of rows & columns
               row <= nrow.and.col,    # remove rows greater than the number we are interested in
               col <=nrow.and.col,   # same for columns as rows
               row %% mod == 0,
               col %% mod == 0)
}

add.n.pts.per.grid <- function(df){
    n.pts<-df %>%
        group_by(unq__ID) %>%
        summarize(n.points = n())

    left_join(df,n.pts)
}


get.pct.cvr.typ <- function(df) {
    df %>%
        group_by(unq__ID, cvr_typ,n.points, area) %>%
        summarize(number = n()) %>%
        ungroup() %>%
        mutate(google.truth.pct.cover = number/n.points) %>%
        dplyr::select(-number)
}

combine.classes.to.g.i.t.o <- function(df) {

    df %>%
        mutate(cvr_typ = as.character(cvr_typ),
               cvr_typ = ifelse(cvr_typ == "s",
                                "i",
                                cvr_typ),
               cvr_typ = ifelse(cvr_typ != "g" &
                                cvr_typ != "i" &
                                cvr_typ != "t", "o", cvr_typ)) %>%
        group_by(unq__ID, cvr_typ, n.points, area) %>%
        summarize(google.truth.pct.cover = sum(google.truth.pct.cover))

}


calc.binomial.pct.cvrs <- function(df) {

    out <- foreach(target.cvr.type = c("g","i","t")) %do%{
        df %>%
            mutate(cvr_typ = ifelse(cvr_typ == target.cvr.type, cvr_typ, "o")) %>%
            group_by(unq__ID, n.points, cvr_typ) %>%
            summarize(pct.cover = sum(pct.cover)) %>%
            mutate(target.type = "binomial",
                   target.cover = target.cvr.type,
                   target.cover = ifelse(target.cover == "g", "grass",
                                  ifelse(target.cover == "t", "tree",
                                         "impervious"))) %>%
            spread(key = cvr_typ, value = pct.cover)
    }
    out <- bind_rows(out)
    out %>%
        rename(pct.g.googleEarth = g, pct.i.googleEarth = i, pct.t.googleEarth = t, pct.o.googleEarth = o)
}



get.area.convexHull <- function(x_coord, y_coord) {
    m <- matrix(c(x_coord, y_coord), ncol = 2)
    ch <- chull(m)
    coords <- m[c(ch,ch[1]),]
    poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)),ID = 1)))
    gArea(poly)
}



calc.pct.cvr.for.grid.subset <- function(df,
                                         n.rows.and.columns.for.subset=20,
                                         mod,
                                         gridID = "unq__ID") {


    df <- filter.by.row.and.col(df, n.rows.and.columns.for.subset, mod) %>%
        add.n.pts.per.grid() %>%
        group_by_(gridID)


    area.df <- df %>%
        summarize(area = get.area.convexHull(crds_x1, crds_x2))

    df <- left_join(df, area.df)


    df <- df %>%
        get.pct.cvr.typ() %>%
        combine.classes.to.g.i.t.o() %>%
                                        #               ungroup() %>%
                                        #               dplyr::select(-n.points) %>%
        spread(., key = cvr_typ, value = google.truth.pct.cover, fill = 0)

                                        #         df[is.na(df)] <- 0

    df.multnm <- df %>%
        mutate(target.type = "multinomial") %>%
        rename(pct.g.googleEarth = g, pct.i.googleEarth = i, pct.t.googleEarth = t) %>%
        mutate(target.cover = "all")

    if(!is.null(df.multnm$o)) { df.multnm <- rename(df.multnm, pct.o.googleEarth = o)}

    df <- df %>%
        gather(key = cvr_typ, value = pct.cover, -unq__ID, -n.points)

    df.binm <- df %>%
        calc.binomial.pct.cvrs()


    df.out <- bind_rows(df.binm, df.multnm)
    return(df.out)
}
## Calculate\ Percent\ Cover\ of\ Grids\,\ subsetted:1 ends here

## [[file:utc.org::*Point-wise%20error%20functions][Point-wise\ error\ functions:1]]
calcErrorAllMultinomial <-  function(pts, tile, Pixel = F) {
    classification <- raster::extract(classified.tile, pts)
    if(Pixel == T) {
        lvls <- levels(classified.tile)[[1]]
        classification <- mapvalues(classification, from = lvls[,1], to = as.character(lvls[,2]))
    } else {
        m <- tile@data@attributes[[1]]
        classification <- mapvalues(classification, from = m$ID, to = levels(m$category))
    }
    google = pts@data$cvr_typ
    overall.error <- 1 - mean(classification == google, na.rm = T)
    pct.grass.classified.as.other <- 1 - mean(classification[which(google == "g")] == google[which(google == "g")], na.rm = T)
    pct.impervious.classified.as.other <- 1 - mean(classification[which(google == "i")] == google[which(google == "i")], na.rm = T)
    pct.tree.classified.as.other <- 1 - mean(classification[which(google == "t")] == google[which(google == "t")], na.rm = T)
    error <- c(overall.error = overall.error,
               pct.grass.classified.as.other = pct.grass.classified.as.other,
               pct.impervious.classified.as.other = pct.impervious.classified.as.other,
               pct.tree.classified.as.other = pct.tree.classified.as.other)
    return(error)
}

calcErrorBinomial <-  function(pts, tile, target, Pixel = F) {
    classification <- raster::extract(classified.tile, pts)
    if(Pixel == T) {
        lvls <- levels(classified.tile)[[1]]
        classification <- mapvalues(classification, from = lvls[,1], to = as.character(lvls[,2]))
    } else {
        classification <- mapvalues(classification, from = c(1,2,3,4), to = c("g","i","t","o"))
    }
    classification <- ifelse(classification == target, classification, "o")
    google <- pts@data$cvr_typ
    google <- ifelse(google == target, google, "o")
    overall.error <- 1 - mean(classification == google)
    pct.grass.classified.as.other <- 1 - mean(classification[which(google == "g")] == google[which(google == "g")])
    pct.impervious.classified.as.other <- 1 - mean(classification[which(google == "i")] == google[which(google == "i")])
    pct.tree.classified.as.other <- 1 - mean(classification[which(google == "t")] == google[which(google == "t")])
    error <- c(overall.error = overall.error,
               pct.grass.classified.as.other = pct.grass.classified.as.other,
               pct.impervious.classified.as.other = pct.impervious.classified.as.other,
               pct.tree.classified.as.other = pct.tree.classified.as.other)
    return(error)
}




calcConfusionMat <- function(pts, tile) {
    classification <- raster::extract(classified.tile, pts)
    classification <- mapvalues(classification, from = c(1,2,3,4), to = c("g","i","t","o"))
    table(classification, google = pts@data$cvr_typ)
}
## Point-wise\ error\ functions:1 ends here

## [[file:utc.org::*Plot%20points%20on%20classifed%20tile][Plot\ points\ on\ classifed\ tile:1]]
pts.on.classified.tile.plot.ErrorinTitle <- function(error, grd.pts, classified.tile.path, fig.dir, target = NULL) {

      grid.name <- str_match(classified.tile.path, ".*([a-z]{3}\\.[0-9]+m\\.[0-9]+)_.*")[,2]
      pts <- grd.pts[grd.pts@data$unq__ID == grid.name,]
  pts@data <- pts@data %>%
        mutate(x = coordinates(pts)[,1],
               y = coordinates(pts)[,2])

    if(target == "a") {
        pts.plot <- geom_point(data = pts@data, aes(x = x, y = y, color = cvr_typ))
        pts.plot <- geom_point(data = pts@data, aes(x = x, y = y, fill = cvr_typ), shape = 21, color = "black", size =2, stroke = .2)
    } else {
        pts@data <- pts@data %>%
            mutate(cvr_typ = ifelse(cvr_typ == target, cvr_typ, "o"))
        pts.plot <- geom_point(data = pts@data, aes(x = x, y = y, color = cvr_typ))
    }
    r.df <- as.data.frame(raster(classified.tile.path), xy = T)
    names(r.df) <- c("x","y","cvr_typ")
                                        #        r.df <- r.df %>%
                                        #            mutate(cvr_typ = mapvalues(cvr_typ, from = c(1,2,3,4), to = c("g","i","t","o")))
    pxls.plot <- ggplot() + geom_raster(data = r.df, aes(x = x, y = y, fill = cvr_typ))
    title <- ggtitle(label = paste0("err:",round(error,2),"_",names(raster(classified.tile.path))))
    UTC_pal <- c(g = "#ffff99", i = "#f0027f", t = "#7fc97f", o = "#666666")
    plt <- pxls.plot + pts.plot + title + scale_fill_manual(values = UTC_pal)+ scale_color_manual(values = UTC_pal) + coord_equal()

    dir.create(fig.dir)

    png(filename = paste0(fig.dir,"/","Err.",round(error,2),"_",names(raster(classified.tile.path)),".png"))
    print(plt)
    dev.off()
#    plt
}

pts.on.classified.tile.plot <- function(pts, classified.tile, fig.dir, target = NULL) {

    if(target == "a") {
        pts.plot <- geom_point(data = pts@data, aes(x = x, y = y, color = cvr_typ))
        pts.plot <- geom_point(data = pts@data, aes(x = x, y = y, fill = cvr_typ), shape = 21, color = "black", size =2, stroke = .2)
    } else {
        pts@data <- pts@data %>%
            mutate(cvr_typ = ifelse(cvr_typ == target, cvr_typ, "o"))
        pts.plot <- geom_point(data = pts@data, aes(x = x, y = y, color = cvr_typ))
    }
    r.df <- as.data.frame(classified.tile, xy = T)
    names(r.df) <- c("x","y","cvr_typ")
                                        #        r.df <- r.df %>%
                                        #            mutate(cvr_typ = mapvalues(cvr_typ, from = c(1,2,3,4), to = c("g","i","t","o")))
    pxls.plot <- ggplot() + geom_raster(data = r.df, aes(x = x, y = y, fill = cvr_typ))
    title <- ggtitle(label = names(classified.tile))
    UTC_pal <- c(g = "#ffff99", i = "#f0027f", t = "#7fc97f", o = "#666666")
    plt <- pxls.plot + pts.plot + title + scale_fill_manual(values = UTC_pal)+ scale_color_manual(values = UTC_pal) + coord_equal()

    dir.create(fig.dir)

    png(filename = paste0(fig.dir,"/",names(classified.tile),".png"))
    print(plt)
    dev.off()
    plt
}
## Plot\ points\ on\ classifed\ tile:1 ends here

## [[file:utc.org::*other%20Functions][other\ Functions:1]]
image_to_classified_image <- function()





                                        # contained urban, don't intersect water = as is
                                        # contained urban, intersect water = mask water
                                        # intersect urban, don't intersect water = mask urban
                                        # intersect urban, intersect water = mask urban & water
                                        # if none of the above, don't write the raster



    Mask_water_crops_urban <- function(image.full.path, water, crops, urban) {

    }




Water_Urban_mask <- function(tile.path, tile.name, urban, water) {
                                        # load image tile
    tile <- stack(tile.path)
                                        # get extent image and make sp object
    et <- as(extent(tile), "SpatialPolygons")
    proj4string(et) <- "+init=epsg:26916"
                                        # Mask out non-urban areas
    if(gContainsProperly(urban,et) & !gIntersects(water,et)){
        writeRaster(tile, filename = str_c(masked.tiles.directory,"/",tile.name), overwrite = T)
    } else if (gContainsProperly(urban,et) & gIntersects(water,et)) {
        tile <- mask(tile, water, inverse = T)
        writeRaster(tile, filename = str_c(masked.tiles.directory,"/",tile.name), overwrite = T)
    } else if (gIntersects(urban, et) & !gIntersects(water,et)) {
        tile <- mask(tile, urban)
        writeRaster(tile, filename = str_c(masked.tiles.directory,"/",tile.name), overwrite = T)
    } else if (gIntersects(urban, et) & gIntersects(water,et)) {
        tile <- mask(tile, urban)
        tile <- mask(tile, water, inverse = T)
        writeRaster(tile, filename = str_c(masked.tiles.directory,"/",tile.name), overwrite = T)
    }
}

Crop_mask <- function(tile.path, tile.name, CDL_stack, n_years){

    tile <- stack(tile.path)
    crops <- crop(CDL_stack, tile)

                                        # These are the values in the CDL that correspond to non crop cover types and not water
    NonCroppedValues <- c(0,63:65, 81:83, 87:88, 112, 121:124, 131, 141:143, 152, 176, 190, 195)
                                        # open water is 111

    NonCroppedValues <- c(0,63:65, 81:83, 87:88, 112, 121:124, 131, 141:143, 152, 176, 190, 195)
                                        # open water is 111. I don't include it in the above list so that it gets masked

                                        # I'm going to add 37, Other Hay/Non-alfalfa, to the non crop cover types
    NonCroppedValues <- c(NonCroppedValues, 37)
                                        # I'm going to add 36, Alfalfa, to the non crop cover types
    NonCroppedValues <- c(NonCroppedValues, 36)

                                        # find cells that have been assigned crop all three years
    crops[crops %in% NonCroppedValues] <- 0
    crops[!(crops %in% NonCroppedValues)] <- 1
    cropsum <- overlay(crops, fun = sum)

    dis.cropsum <- disaggregate(cropsum, fact = 20)
    dis.cropsum <- resample(dis.cropsum, tile, "ngb")
    masked_tile <- mask(tile, dis.cropsum, maskvalue = n_years)

                                        #               Save Image
    writeRaster(masked_tile, paste0(crop.masked.tiles.directory, "/", tile.name), overwrite = T)
}
## other\ Functions:1 ends here

## [[file:utc.org::*Set%20location%20to%20Madison][Set\ location\ to\ Madison:1]]
location <- "madison"
image.paths <- str_extract(image.paths, paste0(".*",location,".*")) %>% na.omit
dd.pca.dir <-  str_extract(dd.pca.dirs, paste0(".*",location,".*")) %>% na.omit
dd.training.dir <- str_extract(dd.training.dirs, paste0(".*",location,".*")) %>% na.omit
dd.models.dir <- str_extract(dd.models.dirs, paste0(".*",location,".*")) %>% na.omit
dd.accuracy.dir <- str_extract(dd.accuracy.dirs, paste0(".*",location,".*")) %>% na.omit
dd.accuracy.classified.dir <-str_extract(dd.accuracy.classified.dirs, paste0(".*",location,".*")) %>% na.omit
## Set\ location\ to\ Madison:1 ends here

## [[file:utc.org::*read%20in%20pca%20model%20if%20it%20exists.%20If%20I%20run%20this,%20don't%20run%20rest%20of%20pca%20code%20in%20this%20subtre][read\ in\ pca\ model\ if\ it\ exists\.\ \ If\ I\ run\ this\,\ don\'t\ run\ rest\ of\ pca\ code\ in\ this\ subtre:1]]
## pca <- foreach(i = seq_along(image.names)) %do% {
##    readRDS(str_c(dd.pca.dir,"/madisonNAIP_pca.rds"))
## }
## read\ in\ pca\ model\ if\ it\ exists\.\ \ If\ I\ run\ this\,\ don\'t\ run\ rest\ of\ pca\ code\ in\ this\ subtre:1 ends here

## [[file:utc.org::*Reproject%20and%20Crop%20PCA%20Region%20Shapefile%20to%20Image][Reproject\ and\ Crop\ PCA\ Region\ Shapefile\ to\ Image:1]]
foreach(img.pth = image.paths) %do% {

         Reproject_Shapefile_to_Image_CRS(pca.region.dsn,
                                         str_c(location,pca.region.layer.appendage),
                                         img.pth,
                                         pca.region.imageCRS.dsn)

       Crop_image_to_each_Shapefile_polygon(pca.region.imageCRS.dsn,
                                         str_c(location,pca.region.layer.appendage),
                                        img.pth,
                                        cores = cores,
                                        output.dir = dd.pca.dir)
}
  closeAllConnections()
## Reproject\ and\ Crop\ PCA\ Region\ Shapefile\ to\ Image:1 ends here

## [[file:utc.org::*Add%20Features%20(ratios%20and%20ndvi)][Add\ Features\ \(ratios\ and\ ndvi\):1]]
cl <- makeCluster(cores)
   registerDoParallel(cl)

    tile.names <- list.files(dd.pca.dir) %>%
        str_extract(., pattern = ".*[0-9]+.tif") %>%
            str_extract(., pattern = ".*[0-9]+") %>%
                na.omit()

   ratios <- foreach (j = tile.names,
            .packages = c("raster","stringr")) %dopar% {
                add.features(tile.dir = dd.pca.dir,
                             tile.name = j,
                             band.names = c("blue","green","red","nir"),
                             ratio.bands = c("blue","green","red","nir"),
                             texture = F)
            }

closeAllConnections()
## Add\ Features\ \(ratios\ and\ ndvi\):1 ends here

## [[file:utc.org::*Create%20and%20Save%20PCA%20model/rotation][Create\ and\ Save\ PCA\ model/rotation:1]]
pca <- foreach(img.nm = image.names) %do% {
                make.and.save.pca.transformation(image.dir = dd.pca.dir,
                                                 image.name = img.nm,
                                                 location = location)
    }
closeAllConnections()
## Create\ and\ Save\ PCA\ model/rotation:1 ends here

## [[file:utc.org::*Make%20Training%20Tiles][Make\ Training\ Tiles:1]]
foreach(img.pth = image.paths) %do% {

      Reproject_Shapefile_to_Image_CRS(training.region.dsn,
                                       str_c(location,training.region.layer.appendage),
                                       img.pth,
                                       training.region.imageCRS.dsn)

      Crop_image_to_each_Shapefile_polygon(training.region.imageCRS.dsn,
                                       str_c(location,training.region.layer.appendage),
                                           img.pth,
                                           cores = cores,
                                           output.dir = dd.training.dir)
}
closeAllConnections()
## Make\ Training\ Tiles:1 ends here

## [[file:utc.org::*Make%20Feature%20data%20frames,%20for%20Each%20Training%20Tile][Make\ Feature\ data\ frames\,\ for\ Each\ Training\ Tile:1]]
cl <- makeCluster(cores)
   registerDoParallel(cl)

   pixel.added.features.raster.list <- foreach(img.nm = image.names) %do% {

       tile.names <- list.files(dd.training.dir) %>%
           str_extract(., pattern = str_c(location,img.nm,".[0-9]+.tif")) %>%
           str_extract(., pattern = str_c(location,img.nm,".[0-9]+")) %>%
           na.omit()

       foreach (i = tile.names,
                .packages = c("glcm","raster","stringr","dplyr")) %dopar% {


                    feature.dfs <- make.feature.df(tile.dir = dd.training.dir,
                                                   tile.name = i,
                                                   image.name = img.nm,
                                                   band.names = c("blue","green","red","nir"),
                                                   ndvi = T,
                                                   ratio.bands = c("blue","green","red","nir"),
                                                   texture.params.df = texture.params,
                                                   pixel.df = T,
                                                   pca.location = location,
                                                   segmentation = T,
                                                   segment.params.df = segment.params)

                }
   }
closeAllConnections()
## Make\ Feature\ data\ frames\,\ for\ Each\ Training\ Tile:1 ends here

## [[file:utc.org::*Combine%20Feature%20Dataframes][Combine\ Feature\ Dataframes:1]]
cl <- makeCluster(cores)
    registerDoParallel(cl)

  feature.dfs <- list.files(dd.training.dir, full.names = T) %>%
      str_extract(paste0(".*(",feature.df.appendage,").*")) %>%
      na.omit()

  foreach(img.nm = image.names) %do% {
      img.feature.dfs <- str_extract(feature.dfs, str_c(".*",img.nm,".*")) %>%
          na.omit()
      SegParams <- unique(str_extract(img.feature.dfs, segmentation.pattern)) %>%
          na.omit()

      foreach(seg.param.set = SegParams, .packages = c("dplyr","stringr")) %dopar% {
          img.seg.feature.dfs = str_extract(img.feature.dfs, str_c(".*",seg.param.set,".*")) %>%
              na.omit()
          dfs <- lapply(img.seg.feature.dfs, readRDS)
          combined.dfs <- bind_rows(dfs)
          saveRDS(combined.dfs, file = str_c(dd.training.dir, "/", location,img.nm, "_",seg.param.set, feature.df.appendage,".rds"))
      }
  }
closeAllConnections()
## Combine\ Feature\ Dataframes:1 ends here

## [[file:utc.org::*Create%20Model%20Building%20Dataframes,%20assign%20Class%20to%20feature%20dfs][Create\ Model\ Building\ Dataframes\,\ assign\ Class\ to\ feature\ dfs:1]]
cl <- makeCluster(cores)
  registerDoParallel(cl)


  model.building.dfs <-  foreach(img.nm = image.names) %do% {

      featureDF.files <- list.files(dd.training.dir) %>%
          str_extract(., str_c(location,img.nm,"_(",segmentation.pattern,")", feature.df.appendage,".rds$")) %>%
          na.omit()

      training.polygon.layer <- list.files(training.region.dsn) %>%
          str_extract(.,str_c(".*",location,img.nm, ".*")) %>%
          na.omit() %>%
          extract.name.from.path() %>%
          unique()

      training.polygons <- readOGR(dsn = training.region.dsn, layer = training.polygon.layer)

      foreach(feature.df.rds = featureDF.files, .packages = c("mlr","foreach","doParallel", "stringr", "raster","rgeos","dplyr","sp","tidyr")) %dopar% {

          feature.df <- readRDS(file = str_c(dd.training.dir,"/",feature.df.rds))

          if(complete.cases(str_extract(feature.df.rds, "Pixel"))) {
              model.building.df <- getPixel.Class.and.Coords.Within.Polygon(PixelFeatureDF = feature.df,
                                                                            training.sp = training.polygons)

              model.building.df <- model.building.df %>%
                  dplyr::select(-x, -y, -id)

              saveRDS(object = model.building.df, file = paste0(dd.training.dir,"/",location,img.nm,"_Pixel",ModelBuilding.appendage,".rds"))
          } else          {
              segment.parameters <- str_extract(feature.df.rds, segmentation.pattern)
              model.building.df <- getSegment.class.and.features.Within.Polygon(SegmentFeatureDF = feature.df,
                  training.sp = training.polygons,
                  seg.tiles.dir = dd.training.dir,
                  seg.params = segment.parameters)
              saveRDS(model.building.df, file = str_c(dd.training.dir,"/",location,img.nm,"_",segment.parameters,ModelBuilding.appendage,".rds"))
          }
      }
  }
closeAllConnections()
## Create\ Model\ Building\ Dataframes\,\ assign\ Class\ to\ feature\ dfs:1 ends here

## [[file:utc.org::*Plot%20Model%20Building%20Dataframes??%20Visualize%20discriminating%20features][Plot\ Model\ Building\ Dataframes\?\?\ Visualize\ discriminating\ features:1]]
seg.p <- "_N-13_C-8"
img.nm <- "PanshpSPOT"
mod.df <- readRDS(paste0(dd.training.dir, "/",location,img.nm, seg.p, ModelBuilding.appendage, ".rds"))

    ggplot(mod.df, aes(color = factor(Class), y = ndvi_mean, x = red_sd)) + geom_point(alpha = .9)
#    ggplot(out, aes(color = factor(Class), y = ndvi_mean, x = red_sd)) + geom_point(alpha = .9)
#  ggplot(model.building.df, aes(color = factor(Class), y = ndvi_mean, x = red_sd)) + geom_point(alpha = .5)
## Plot\ Model\ Building\ Dataframes\?\?\ Visualize\ discriminating\ features:1 ends here

## [[file:utc.org::*Plot%20Model%20Building%20Dataframes??%20Visualize%20discriminating%20features][Plot\ Model\ Building\ Dataframes\?\?\ Visualize\ discriminating\ features:2]]
mod.df <- readRDS(paste0(dd.training.dir, "/",location,img.nm, "_Pixel", ModelBuilding.appendage,".rds"))
ggplot(mod.df, aes(color = factor(Class), y = ndvi, x = ratio.nir_stat.homogeneity_window.3_angle.0)) + geom_point(alpha = .5)
## Plot\ Model\ Building\ Dataframes\?\?\ Visualize\ discriminating\ features:2 ends here

## [[file:utc.org::*untuned%20models][untuned\ models:1]]
cl <- makeCluster(cores)
  registerDoParallel(cl)

  foreach(img.nm = image.names) %do% {

              ModelBuildingRDSs <- list.files(dd.training.dir) %>%
                  str_extract(., str_c(location,img.nm,".*",ModelBuilding.appendage, ".rds")) %>%
                  na.omit()

              foreach(ModelBuildingRDS = ModelBuildingRDSs,
          .packages = c("parallelMap","randomForest","kernlab","irace","mlr","stringr","dplyr","foreach","doParallel")) %dopar% {

                  Build.and.Save.models(dir = dd.training.dir,
                                        modelBuildingData = ModelBuildingRDS,
                                        models.dir = dd.models.dir,
                                        image.name = img.nm,
                                        location = location)
              }
          }
closeAllConnections()
## untuned\ models:1 ends here

## [[file:utc.org::*Show%20importance%20of%20features][Show\ importance\ of\ features:1]]
cl <- makeCluster(cores)
registerDoParallel(cl)

fv <- foreach(img.nm = image.names, .combine = "rbind") %do% {

    ModelBuildingRDSs <- list.files(dd.training.dir) %>%
        str_extract(., str_c(location,img.nm,".*",ModelBuilding.appendage, ".rds")) %>%
        na.omit()

     foreach(ModelBuildingRDS = ModelBuildingRDSs,
            .packages = c("parallelMap","randomForest","kernlab","irace","mlr","stringr","dplyr","foreach","doParallel"),
            .combine = "rbind") %do% {

                fv <- print.feature.importance(dir = dd.training.dir,
                                         modelBuildingData = ModelBuildingRDS,
                                         image.name = img.nm,
                                         location = location,
                                         feature.importance.methods = c("information.gain","chi.squared"))

              fv$modelBuildingDF <- ModelBuildingRDS
              fv
            }
}
## Show\ importance\ of\ features:1 ends here

## [[file:utc.org::*Show%20importance%20of%20features][Show\ importance\ of\ features:2]]
options(warn = -1)
## fv.u <- unlist(fv, recursive = F)
## fv.u <- unlist(fv.u, recursive = F)
## sapply(fv.u, ascii)
fv %>% ascii
options(warn = 1)
## Show\ importance\ of\ features:2 ends here

## [[file:utc.org::*build%20models][build\ models:1]]
cl <- makeCluster(cores)
registerDoParallel(cl)

foreach(img.nm = image.names) %do% {

    ModelBuildingRDSs <- list.files(dd.training.dir) %>%
        str_extract(., str_c(location,img.nm,".*",ModelBuilding.appendage, ".rds")) %>%
        na.omit()

    foreach(ModelBuildingRDS = ModelBuildingRDSs,
            .packages = c("parallelMap","randomForest","kernlab","irace","mlr","stringr","dplyr","foreach","doParallel")) %dopar% {

                Build.and.Save.FeatureSelected.models(dir = dd.training.dir,
                                                      modelBuildingData = ModelBuildingRDS,
                                                      models.dir = dd.models.dir,
                                                      image.name = img.nm,
                                                      location = location)
            }
}
## build\ models:1 ends here

## [[file:utc.org::*tuned][tuned:1]]
Build.and.Save.Tuned.models <- function( dir = dd.training.dir,
                                   modelBuildingData = ModelBuildingRDS,
                                   models.dir = Models.dir,
                                   image.name){

     dat <- readRDS(paste0(dir,"/",modelBuildingData)) %>%
         as.data.frame()

     image.and.segmentation.stem = str_replace(modelBuildingData, ModelBuilding.appendage,"")

     names <- colnames(dat)
     names <- str_replace(names, "\\(",".")
     names <- str_replace(names, "\\)",".")
     names <- str_replace(names, "\\:",".")
     colnames(dat) <- names

                                         # Create Task
     utc.task <- makeClassifTask(id = image.name, data = dat, target = "Class")

                                         # make parameter set for tuning

     rf.ps <- makeParamSet(makeIntegerParam("ntree", lower = 1L, upper = 500L),
                           makeIntegerParam("mtry", lower = 1L, upper = 50L))

     svm.ps <- makeParamSet(makeNumericParam("C", lower = -12, upper = 12, trafo = function(x) 2^x),
                            makeDiscreteParam("kernel", values = c("vanilladot", "polydot", "rbfdot")),
                            makeNumericParam("sigma", lower = -12, upper = 12, trafo = function(x) 2^x,
                                             requires = quote(kernel == "rbfdot")),
                            makeIntegerParam("degree", lower = 2L, upper = 5L,
                                             requires = quote(kernel == "polydot")))

                                         # tune
                                         # inner

     ctrl = makeTuneControlIrace(maxExperiments = 200L)
     inner = makeResampleDesc("CV", iters = 2L)
     svm.lrn = makeTuneWrapper("classif.ksvm", resampling = inner, par.set = svm.ps, control = ctrl, show.info = T)
     rf.lrn = makeTuneWrapper("classif.randomForest", resampling = inner, par.set = rf.ps, control = ctrl, show.info = T)

                                         #outer
     lrnrs = list(svm.lrn, rf.lrn)
     outer = makeResampleDesc("CV", iters = 3L)

 #    parallelStartMulticore(cores)

     res = benchmark(lrnrs, utc.task, outer, measures = acc, show.info = FALSE)

 #   parallelStop()

     saveRDS(res, file = paste0(models.dir,"/",image.and.segmentation.stem, models.appendage))
 }
## tuned:1 ends here

## [[file:utc.org::*Look%20at%20models][Look\ at\ models:1]]
df <- readRDS(paste0(dd.training.dir, "/madisonNAIP_N-30_C-15.ModelBuilding.rds"))
mod <- readRDS(paste0(Models.dir, "/madisonNAIP_N-100_C-50.models.rds"))



getBMRModels(mod)
getBMRLearners(mod)
getBMRPerformances(mod)
getBMRTuneModults(mod, as.df = T)

getBMRTuneModults(mod, as.df = T) %>%
    group_by(learner.id) %>%
    summarize_each(funs = "mean")


mods<-getBMRModels(mod)
## Look\ at\ models:1 ends here

## [[file:utc.org::*Make%20tiles%20at%20accuracy%20regions][Make\ tiles\ at\ accuracy\ regions:1]]
foreach(i = 1) %do% {

    foreach(img.pth = image.paths) %do% {

        Reproject_Shapefile_to_Image_CRS(accuracy.region.dsn[i],
                                         accuracy.region.layer[i],
                                         img.pth,
                                         accuracy.region.imageCRS.dsn)

        Crop_image_to_regions_around_points_nameBygrid(shapefile.dsn = accuracy.region.imageCRS.dsn,
                                                       shapefile.layer = accuracy.region.layer[i],
                                                       image.path = img.pth,
                                                       cores = cores,
                                                       output.dir = dd.accuracy.dir,
                                                       column.name = tile.id.col.nm.for.grid.and.field.accuracy[i])

    }
}
## Make\ tiles\ at\ accuracy\ regions:1 ends here

## [[file:utc.org::*Make%20Feature%20data%20frames,%20for%20each%20Accuracy%20Region%20tile][Make\ Feature\ data\ frames\,\ for\ each\ Accuracy\ Region\ tile:1]]
cl <- makeCluster(cores)
registerDoParallel(cl)

pixel.added.features.raster.list <- foreach(img.nm = image.names) %do% {

                                        #img.nm <- image.names[1]

    tile.names <- list.files(dd.accuracy.dir) %>%
        str_match(., pattern = str_c("(",location,img.nm,".",grid.pattern,")(.tif)"))

    tile.names <- tile.names[,2] %>% na.omit()

    foreach (i = tile.names,
             .packages = c("glcm","raster","stringr","dplyr")) %dopar% {

                 feature.dfs <- make.feature.df(tile.dir = dd.accuracy.dir,
                                                tile.name = i,
                                                image.name = img.nm,
                                                band.names = c("blue","green","red","nir"),
                                                ndvi = T,
                                                ratio.bands = c("blue","green","red","nir"),
                                                texture.params.df = texture.params,
                                                pixel.df = T,
                                                pca.location = location,
                                                segmentation = T,
                                               segment.params.df = segment.params)

             }
}
## Make\ Feature\ data\ frames\,\ for\ each\ Accuracy\ Region\ tile:1 ends here

## [[file:utc.org::*Classify%20Tiles%20at%20accuracy%20regions][Classify\ Tiles\ at\ accuracy\ regions:1]]
cl <- makeCluster(cores)
 registerDoParallel(cl)


 classified.grid.tiles <-       foreach(img.nm = image.names) %do% {

         models <- list.files(dd.models.dir) %>%
             str_extract(., str_c(".*",location,img.nm,".*")) %>%
             na.omit()

         tile.names <- list.files(dd.accuracy.dir) %>%
             str_match(., pattern = str_c("(",location,img.nm,".*?)_.*\\.tif$"))

         tile.names <- tile.names[,2] %>% na.omit() %>% unique()


         foreach(tile.nm = tile.names,
                 .packages = c("plyr","dplyr","raster","stringr","mlr","foreach","doParallel")) %do% {

             foreach(model = models) %do% {

                 seg.p <- str_extract(model, segmentation.pattern)

                 if(grepl("N-[0-9]+_C-[0-9]+",seg.p)) {
                        segment.tile.name.append <- paste0("_",seg.p,".tif")
                        segment.feature.df.name.append <- paste0("_",seg.p,feature.df.appendage,".rds")

                        classify.segmented.raster(segment.feature.df.dir = dd.accuracy.dir,
                                        model.dir = dd.models.dir,
                                        segment.dir = dd.accuracy.dir,
                                        classify.out.dir = dd.accuracy.dir,
                                        tile.name = tile.nm,
                                        segmentation.appendage = segment.tile.name.append,
                                        model.name.rds = model,
                                        segment.feature.appendage = segment.feature.df.name.append,
                                        segmentation.prms = seg.p)

                 } else {
                     classify.pixel.raster(tile.dir = dd.accuracy.dir,
                                           tile.name = tile.nm,
                                           pixelFeatureDF.appendage = feature.df.appendage,
                                           model.dir = dd.models.dir,
                                           model.rds = model,
                                           seg.prms = seg.p)
                 }
             }
         }
     }


stopCluster(cl)
## Classify\ Tiles\ at\ accuracy\ regions:1 ends here

## [[file:utc.org::*Point-wise%20accuracy.%20regular%20confusion%20matrix%20thing.%20I%20should%20do%20this%20for%20the%20grids%20and%20the%20field%20plot%20data][Point-wise\ accuracy\.\ \ regular\ confusion\ matrix\ thing\.\ \ I\ should\ do\ this\ for\ the\ grids\ and\ the\ field\ plot\ data:1]]
grd <- readOGR(dsn = grid.accuracy.region.dsn, layer = grid.accuracy.region.layer, stringsAsFactors = F)

    xy <- coordinates(grd)
    grd@data$x <- xy[,1]
    grd@data$y <- xy[,2]

classified.tile.paths <- list.files(str_c(dd.accuracy.classified.dir), full.names = T) %>%
    str_extract(., pattern = ".*.tif$") %>%
        str_extract(., pattern = str_c(".*",grid.pattern, ".*")) %>%
        na.omit()


grid.names <- classified.tile.paths %>%
    str_match(., paste0(".*(",grid.pattern,").*"))

grid.names <- grid.names[,2] %>%
    unique() %>%
    na.omit()

## grid.name = str_extract(grid.names, ".*150m-[56].*") %>% na.omit()



    cl <- makeCluster(cores)
    registerDoParallel(cl)


    error.df <- foreach(grid.name = grid.names, .combine = "rbind") %do% {

        pts <- grd[grd@data$unq__ID== grid.name,]

        classified.tile.paths.at.grid <- str_extract(classified.tile.paths, str_c(".*",grid.name,"_.*")) %>%
            na.omit()

        ## classified.tile.paths.at.grid2 = classified.tile.paths.at.grid %>%
        ##      str_extract(., ".*madisonNAIP.*N-105.*svm_.*") %>%
        ##      na.omit()

#         classified.tile.path.at.grid = classified.tile.paths.at.grid[1]



        foreach(classified.tile.path.at.grid = classified.tile.paths.at.grid,
                .combine = "rbind",
                .packages = c("plyr","raster","dplyr", "stringr","ggplot2")) %dopar% {

                    classified.tile.name.at.grid <- basename(classified.tile.path.at.grid)
                    classified.tile <- raster(classified.tile.path.at.grid)

                    tgt <- str_extract(classified.tile.name.at.grid, "tree|grass|impervious|all")
                    tgt <- mapvalues(tgt, c("tree","grass","impervious","all"), c("t","g","i","a"))

                   ##  png(str_c("figs/","ClassifiedVersusGrid","/",names(classified.tile),".png"))
                   ## print(pts.on.classified.tile.plot(pts, classified.tile, target = tgt))
                   ## dev.off()

                    PixBool <- !is.na((str_extract(classified.tile.path.at.grid, "_Pixel_")))

                    if(!is.na(str_extract(classified.tile.path.at.grid, "_all_"))) {
                        error <- calcErrorAllMultinomial(pts, classified.tile, Pixel = PixBool)
                        error <- error %>%
                            t() %>%
                            data.frame() %>%
                            mutate(grid = grid.name,
                                   image =  str_extract(classified.tile.name.at.grid, image.pattern),
                                   target.cover = str_extract(classified.tile.name.at.grid, target.pattern),
                                   model =  str_extract(classified.tile.name.at.grid, model.pattern),
                                   seg.params = str_extract(classified.tile.name.at.grid, segmentation.pattern))
                        error
                    } else {
                        target = str_extract(classified.tile.name.at.grid, "tree|grass|impervious")
                        target <- mapvalues(target, c("tree","grass","impervious"), c("t","g","i"))
                        error <- calcErrorBinomial(pts, classified.tile, target, Pixel = PixBool)
                        error <- error %>%
                            t() %>%
                            data.frame() %>%
                            mutate(grid = grid.name,
                                   image =  str_extract(classified.tile.name.at.grid, image.pattern),
                                   target.cover = str_extract(classified.tile.name.at.grid, target.pattern),
                                   model =  str_extract(classified.tile.name.at.grid, model.pattern),
                                   seg.params = str_extract(classified.tile.name.at.grid, seg.prms))

                        error
                    }
                }
    }



    saveRDS(error.df, str_c(derived.dir, "/point2pixel.error.df.rds"))
## Point-wise\ accuracy\.\ \ regular\ confusion\ matrix\ thing\.\ \ I\ should\ do\ this\ for\ the\ grids\ and\ the\ field\ plot\ data:1 ends here

## [[file:utc.org::*Summarize%20Accuracy%20Assessment%20Results][Summarize\ Accuracy\ Assessment\ Results:1]]
error.df <- readRDS(str_c(derived.dir, "/point2pixel.error.df.rds"))

    error.df %>%
        arrange(overall.error) %>%
        head()

    error.df %>%
        arrange(desc(overall.error)) %>%
        head()

    error.df %>%
        filter(seg.params != "Pixel") %>%
        arrange(desc(overall.error)) %>%
        head()

error.df <- error.df %>%
    mutate(segment.size = as.numeric(ifelse(!is.na(str_match(seg.params, "N-([0-9]+)_C-[0-9]+")[,2]), str_match(seg.params, "N-([0-9]+)_C-[0-9]+")[,2], 1)),
           segment.size = ifelse(image == "panshpSPOT", segment.size * 1.5, segment.size),
           compactness = as.numeric(str_match(seg.params, "N-[0-9]+_C-([0-9]+)")[,2]))
## Summarize\ Accuracy\ Assessment\ Results:1 ends here

## [[file:utc.org::*Table%20showing%20performance%20of%20classifiers][Table\ showing\ performance\ of\ classifiers:1]]

## Table\ showing\ performance\ of\ classifiers:1 ends here

## [[file:utc.org::*load%20grid.points][load\ grid\.points:1]]
grid.points <- readOGR(dsn = accuracy.region.imageCRS.dsn,
                       layer = "madisonNAIP_Grids")
## load\ grid\.points:1 ends here

## [[file:utc.org::*Plots%20of%2020%20best%20classified%20grids%20with%20points%20superimposed][Plots\ of\ 20\ best\ classified\ grids\ with\ points\ superimposed:1]]
best.classified.grids <- error.df %>%
        ungroup() %>%
        group_by(grid) %>%
        top_n(1, desc(overall.error)) %>%
        ungroup() %>%
        arrange(overall.error) %>%
        select(overall.error, grid,image, target.cover, model, seg.params) %>%
        mutate(path = paste0(dd.accuracy.classified.dir,"/",location,image,".",grid,"_",seg.params,"_",image,"_",target.cover,"_",model,".tif")) %>%
        head(n = 20)

options(warn = -1)
  best.classified.grids %>% ascii
options(warn = 1)
## Plots\ of\ 20\ best\ classified\ grids\ with\ points\ superimposed:1 ends here

## [[file:utc.org::*Plots%20of%2020%20best%20classified%20grids%20with%20points%20superimposed][Plots\ of\ 20\ best\ classified\ grids\ with\ points\ superimposed:2]]
lapply(1:nrow(best.classified.grids), function(i){
    pts.on.classified.tile.plot.ErrorinTitle(error = best.classified.grids$overall.error[i],
                                         grd.pts = grid.points,
                                         classified.tile.path = best.classified.grids$path[i],
                                         fig.dir = "figs/bestgrids",
                                         target = "a")
})

## plts <- lapply(best.classified.grids$path, function(path) {
##   grid.name <- str_match(path, ".*([a-z]{3}\\.[0-9]+m\\.[0-9]+)_.*")[,2]
##   points <- grid.points[grid.points@data$unq__ID == grid.name,]
##   points@data <- points@data %>%
##       mutate(x = coordinates(points)[,1],
##              y = coordinates(points)[,2])
##   ras <- raster(path)
##   pts.on.classified.tile.plot(fig.dir = "figs/bestgrids",points, ras, target = "a")
## })
## Plots\ of\ 20\ best\ classified\ grids\ with\ points\ superimposed:2 ends here

## [[file:utc.org::*Plots%20of%2020%20best%20classified%20grids%20with%20points%20superimposed][Plots\ of\ 20\ best\ classified\ grids\ with\ points\ superimposed:3]]
best.grid.paths <- list.files("figs/bestgrids", full.names = T)

a <- sapply(best.grid.paths, function(x) message("[[file:",x,"]]"))
## Plots\ of\ 20\ best\ classified\ grids\ with\ points\ superimposed:3 ends here

## [[file:utc.org::*Plots%20of%2020%20worst%20classified%20grids%20with%20points%20superimposed.%20NONE%20SHOULD%20BE%20>50%25%20wrong!][Plots\ of\ 20\ worst\ classified\ grids\ with\ points\ superimposed\.\ \ NONE\ SHOULD\ BE\ >50%\ wrong!:1]]
worst.classified.grids <- error.df %>%
        ungroup() %>%
        group_by(grid) %>%
        top_n(1, overall.error) %>%
        ungroup() %>%
        arrange(desc(overall.error)) %>%
        select(overall.error, grid,image, target.cover, model, seg.params) %>%
        mutate(path = paste0(dd.accuracy.classified.dir,"/",image,".",grid,"_",seg.params,"_",image,"_",target.cover,"_",model,".tif")) %>%
        head(n = 20)

options(warn = -1)
  worst.classified.grids %>% ascii
options(warn = 1)
## Plots\ of\ 20\ worst\ classified\ grids\ with\ points\ superimposed\.\ \ NONE\ SHOULD\ BE\ >50%\ wrong!:1 ends here

## [[file:utc.org::*Plots%20of%2020%20worst%20classified%20grids%20with%20points%20superimposed.%20NONE%20SHOULD%20BE%20>50%25%20wrong!][Plots\ of\ 20\ worst\ classified\ grids\ with\ points\ superimposed\.\ \ NONE\ SHOULD\ BE\ >50%\ wrong!:2]]
worst.grid.paths <- list.files("figs/worstgrids", full.names = T)

a <- sapply(worst.grid.paths, function(x) message("[[file:",x,"]]"))
## Plots\ of\ 20\ worst\ classified\ grids\ with\ points\ superimposed\.\ \ NONE\ SHOULD\ BE\ >50%\ wrong!:2 ends here

## [[file:utc.org::*Table%20showing%20performace%20of%20classifiers,%20average%20over%20all%20grids,%20increasing%20accuracy][Table\ showing\ performace\ of\ classifiers\,\ average\ over\ all\ grids\,\ increasing\ accuracy:1]]
error.df.avg.class <- error.df %>%
      select(-grid) %>%
      group_by(image, target.cover, model, seg.params, segment.size, compactness) %>%
      summarize_each(funs(mean(.,na.rm = T))) %>%
      ungroup() %>%
      arrange(overall.error)


options(warn = -1)
  error.df.avg.class %>% ascii
options(warn = 1)
## Table\ showing\ performace\ of\ classifiers\,\ average\ over\ all\ grids\,\ increasing\ accuracy:1 ends here

## [[file:utc.org::*Table%20showing%20performace%20of%20classifiers,%20average%20over%20all%20grids,%20decreasing%20accuracy][Table\ showing\ performace\ of\ classifiers\,\ average\ over\ all\ grids\,\ decreasing\ accuracy:1]]
options(warn = -1)
  error.df.avg.class %>%
      arrange(desc(overall.error)) %>%
      ascii
options(warn = 1)
## Table\ showing\ performace\ of\ classifiers\,\ average\ over\ all\ grids\,\ decreasing\ accuracy:1 ends here

## [[file:utc.org::*Plots%20of%2020%20*best*%20classified%20grids%20by%20*best*%20classifier%20with%20points%20superimposed][Plots\ of\ 20\ *best*\ classified\ grids\ by\ *best*\ classifier\ with\ points\ superimposed:1]]
best.classif.overall <- error.df.avg.class %>%
    arrange(overall.error) %>%
    slice(1) %>%
    data.frame()

best.classif.best.grids <- best.classif.overall %>%
    select(image, target.cover, model, seg.params) %>%
    left_join(., error.df) %>%
    arrange(overall.error) %>%
    select(overall.error, grid,image, target.cover, model, seg.params) %>%
    mutate(path = paste0(dd.accuracy.classified.dir,"/",image,".",grid,"_",seg.params,"_",image,"_",target.cover,"_",model,".tif"))

    options(warn = -1)

          best.classif.best.grids %>%
          ascii

    options(warn = 1)
## Plots\ of\ 20\ *best*\ classified\ grids\ by\ *best*\ classifier\ with\ points\ superimposed:1 ends here

## [[file:utc.org::*Plots%20of%2020%20*best*%20classified%20grids%20by%20*best*%20classifier%20with%20points%20superimposed][Plots\ of\ 20\ *best*\ classified\ grids\ by\ *best*\ classifier\ with\ points\ superimposed:2]]
best.classif.best.grids <- best.classif.best.grids %>%  head(n=20)


  lapply(1:nrow(best.classif.best.grids), function(i){
      pts.on.classified.tile.plot.ErrorinTitle(error = best.classif.best.grids$overall.error[i],
                                           grd.pts = grid.points,
                                           classified.tile.path = best.classif.best.grids$path[i],
                                           fig.dir = "figs/bestclassif.bestgrids",
                                           target = "a")
  })



  ## plts <- lapply(best.classif.best.grids$path, function(path) {
  ##     grid.name <- str_match(path, ".*([a-z]{3}\\.[0-9]+m\\.[0-9]+)_.*")[,2]
  ##     points <- grid.points[grid.points@data$unq__ID == grid.name,]
  ##     points@data <- points@data %>%
  ##         mutate(x = coordinates(points)[,1],
  ##                y = coordinates(points)[,2])
  ##     ras <- raster(path)
  ##     pts.on.classified.tile.plot(fig.dir = "figs/bestclassif.bestgrids",points, ras, target = "a")
  ## })
## Plots\ of\ 20\ *best*\ classified\ grids\ by\ *best*\ classifier\ with\ points\ superimposed:2 ends here

## [[file:utc.org::*Plots%20of%2020%20*worst*%20classified%20grids%20by%20*best*%20classifier%20with%20points%20superimposed][Plots\ of\ 20\ *worst*\ classified\ grids\ by\ *best*\ classifier\ with\ points\ superimposed:1]]
best.classif.overall <- error.df.avg.class %>%
      arrange(overall.error) %>%
        slice(1) %>%
      data.frame()

best.classif.worst.grids <- best.classif.overall %>%
  select(image, target.cover, model, seg.params) %>%
    left_join(., error.df) %>%
    arrange(desc(overall.error)) %>%
    select(overall.error, grid,image, target.cover, model, seg.params) %>%
    mutate(path = paste0(dd.accuracy.dir,"/",ClassifiedTilesDirName,"/",image,".",grid,"_",seg.params,"_",image,"_",target.cover,"_",model,".tif")) %>%
    head(n=20)


lapply(1:nrow(best.classif.worst.grids), function(i){
    pts.on.classified.tile.plot.ErrorinTitle(error = best.classif.worst.grids$overall.error[i],
                                         grd.pts = grid.points,
                                         classified.tile.path = best.classif.worst.grids$path[i],
                                         fig.dir = "figs/bestclassif.worstgrids",
                                         target = "a")
})

## plts <- lapply(best.classif.worst.grids$path, function(path) {
##   grid.name <- str_match(path, ".*([a-z]{3}\\.[0-9]+m\\.[0-9]+)_.*")[,2]
##   points <- grid.points[grid.points@data$unq__ID == grid.name,]
##   points@data <- points@data %>%
##       mutate(x = coordinates(points)[,1],
##              y = coordinates(points)[,2])
##   ras <- raster(path)
##   pts.on.classified.tile.plot(fig.dir = "figs/bestclassif.worstgrids",points, ras, target = "a")
## })
## Plots\ of\ 20\ *worst*\ classified\ grids\ by\ *best*\ classifier\ with\ points\ superimposed:1 ends here

## [[file:utc.org::*plot%201][plot\ 1:1]]
ggplot(error.df, aes(y = overall.error, x = segment.size, color = compactness, group = model)) + geom_point() + facet_grid(model~image)
## plot\ 1:1 ends here

## [[file:utc.org::*Table%20showing%20the%20*best*%20classified%20grids,%20averaged%20across%20all%20classifiers][Table\ showing\ the\ *best*\ classified\ grids\,\ averaged\ across\ all\ classifiers:1]]
error.df.avg.grids <- error.df %>%
      select(-image, -model, -target.cover, -seg.params, -segment.size, -compactness) %>%
      group_by(grid) %>%
      summarize_each(funs(mean(.,na.rm = T))) %>%
      ungroup() %>%
      arrange(overall.error)


options(warn = -1)
  error.df.avg.grids %>% ascii
options(warn = 1)
## Table\ showing\ the\ *best*\ classified\ grids\,\ averaged\ across\ all\ classifiers:1 ends here

## [[file:utc.org::*Table%20showing%20the%20*worst*%20classified%20grids,%20averaged%20across%20all%20classifiers][Table\ showing\ the\ *worst*\ classified\ grids\,\ averaged\ across\ all\ classifiers:1]]
error.df.avg.grids <- error.df %>%
      select(-image, -model, -target.cover, -seg.params, -segment.size, -compactness) %>%
      group_by(grid) %>%
      summarize_each(funs(mean(.,na.rm = T))) %>%
      ungroup() %>%
      arrange(overall.error)

options(warn = -1)
  error.df.avg.grids %>% arrange(desc(overall.error)) %>% ascii
options(warn = 1)
## Table\ showing\ the\ *worst*\ classified\ grids\,\ averaged\ across\ all\ classifiers:1 ends here

## [[file:utc.org::*Create%20DF%20of%20%25%20cover%20from%20grids%20cropped%20to%20different%20extents][Create\ DF\ of\ %\ cover\ from\ grids\ cropped\ to\ different\ extents:1]]
grd <- readOGR(dsn = grid.accuracy.region.dsn, layer = grid.accuracy.region.layer)
     grd.df <- grd@data

n.rows.and.columns.for.subset = c(15)

     out <- foreach(n.rows.and.columns.for.sub = n.rows.and.columns.for.subset) %do% {
         calc.pct.cvr.for.grid.subset(grd.df, n.rows.and.columns.for.sub)
     }

     Percent.Cover.Grids.with.diff.targettypes.and.diff.subsets <- bind_rows(out)

Percent.Cover.Grids.with.diff.targettypes.and.diff.subsets <- Percent.Cover.Grids.with.diff.targettypes.and.diff.subsets %>%
    rename(grid = unq__ID)

  saveRDS(Percent.Cover.Grids.with.diff.targettypes.and.diff.subsets, str_c(derived.dir,"/","Wausau.Percent.Cover.Grids.with.diff.targettypes.and.diff.subsets.dataframe",".rds"))
## Create\ DF\ of\ %\ cover\ from\ grids\ cropped\ to\ different\ extents:1 ends here

## [[file:utc.org::*Create%20DF%20of%20%25%20cover%20from%20classified%20rasters%20cropped%20to%20different%20extents][Create\ DF\ of\ %\ cover\ from\ classified\ rasters\ cropped\ to\ different\ extents:1]]
grd <- readOGR(dsn = grid.accuracy.region.dsn, layer = grid.accuracy.region.layer)


  # get path of grid tiles (not interested in fieldplot classified tiles)
      classified.tile.paths <- list.files(str_c(dd.accuracy.classified.dir), full.names = T) %>%
          str_extract(., pattern = ".*.tif$") %>%
          str_extract(., pattern = str_c(".*",wau.grid.id.pattern, ".*")) %>%
            na.omit()


n.rows.and.columns.for.subset = c(15)


cl <- makeCluster(cores)
registerDoParallel(cl)


    out <- foreach(n.rows.and.columns.for.sub = n.rows.and.columns.for.subset) %do% {
           pct.class.cover <- foreach(tile.path = classified.tile.paths, .packages = c("raster","dplyr","stringr")) %dopar% {
             calculate.percent.cover.in.classified.tile(pts = grd,
                                                         tile.pth = tile.path,
                                                         n.rows.and.columns.subset = n.rows.and.columns.for.sub)

          }
              saveRDS(pct.class.cover, str_c(derived.dir,"/","Wausau.Percent.Cover.Classified.Tiles.nPoints",n.rows.and.columns.for.sub, ".rds"))
    }


class.cover.files <- list.files(derived.dir, pattern = "Wausau.Percent.Cover.Classified.Tiles.nPoints*", full.names = T)

class.cover.dfs <- lapply(class.cover.files, readRDS)

out <- unlist(class.cover.dfs,recursive = F)

     Percent.Cover.Classified.Tiles.dataframe <- bind_rows(out)





# delete this line if I run it again.
## Percent.Cover.Classified.Tiles.dataframe <-rename(Percent.Cover.Classified.Tiles.dataframe,
##                                                   image = tile,
##                                                   pct_g_pred = pct_g,
##                                                   pct_i_pred = pct_i,
##                                                   pct_t_pred = pct_t,
##                                                   pct_o_pred = pct_o)

  ## saveRDS(Percent.Cover.Classified.Tiles.dataframe, str_c(derived.dir,"/","Percent.Cover.Classified.Tiles.dataframe",".rds"))
## Create\ DF\ of\ %\ cover\ from\ classified\ rasters\ cropped\ to\ different\ extents:1 ends here

## [[file:utc.org::*Join%20Cover%20from%20Grids%20with%20predicted%20Cover%20from%20images][Join\ Cover\ from\ Grids\ with\ predicted\ Cover\ from\ images:1]]
Percent.Cover.Grids.with.diff.targettypes.and.diff.subsets <- readRDS(str_c(derived.dir,"/","Wausau.Percent.Cover.Grids.with.diff.targettypes.and.diff.subsets.dataframe",".rds"))

  str(Percent.Cover.Grids.with.diff.targettypes.and.diff.subsets)
  str(Percent.Cover.Classified.Tiles.dataframe)

Percent.Cover.Classified.Tiles.dataframe %>%
    filter(seg.params == "Pixel") %>%
    data.frame() %>%
    head()

  Percent.Cover.Grids.with.diff.targettypes.and.diff.subsets %>%
      filter(n.points == 400)


#Percent.Cover.Classified.Tiles.dataframe <- Percent.Cover.Classified.Tiles.dataframe %>%
#    rename(pct_g_pred = pct_g, pct_t_pred = pct_t, pct_i_pred = pct_i, pct_o_pred = pct_o)


  grid.master.df <- left_join(Percent.Cover.Classified.Tiles.dataframe, Percent.Cover.Grids.with.diff.targettypes.and.diff.subsets)

  # Should join by Joining by: c("grid", "target.cover", "n.points", "target.type")

  str(grid.master.df)

  grid.master.df %>%
#      filter(n.points == 400) %>%
      data.frame() %>%
      head(n=40)
## Join\ Cover\ from\ Grids\ with\ predicted\ Cover\ from\ images:1 ends here

## [[file:utc.org::*Make%20RMSE%20plots][Make\ RMSE\ plots:1]]
sub.for.rmse.plot <- grid.master.df %>%
      filter(target.type == "multinomial",
             image == "wausauNAIP",
             n.points == 225)


  ggplot(sub.for.rmse.plot, aes( x = pct.t.googleEarth, y = pct_t_pred, color = model)) +
geom_point() + geom_smooth() + theme_classic() +
geom_line(data = data.frame(pct.t.googleEarth = c(0,1), pct_t_pred = c(0,1), seg.params = "1:1"),
color = "black", size = 1) +
ggtitle("NAIP, n.pts: 225")
## Make\ RMSE\ plots:1 ends here

## [[file:utc.org::*Calc%20Error%20Column][Calc\ Error\ Column:1]]
error_tree <- grid.master.df %>%
    filter(target.cover == "tree" | target.cover == "all") %>%
    select(-target.cover) %>%
    group_by(image, model, n.points, seg.params, target.type) %>%
    mutate(t_error = (pct_t_pred - pct.t.googleEarth))

error_tree %>%
    select(image, model, n.points, seg.params, target.type, grid, t_error) %>%
    filter(n.points == 225) %>%
    ungroup() %>%
    arrange(desc(abs(t_error))) %>%
    data.frame() %>%
    head(n=50)
## Calc\ Error\ Column:1 ends here

## [[file:utc.org::*Calc%20Error%20Column][Calc\ Error\ Column:2]]
RMSE_tree <- grid.master.df %>%
      filter(target.cover == "tree" | target.cover == "all") %>%
      select(-target.cover) %>%
      group_by(image, model, n.points, seg.params, target.type) %>%
      summarize(RMSE_t = sqrt( mean( (pct_t_pred - pct.t.googleEarth)^2, na.rm =T ) ) )

RMSE_tree <- RMSE_tree %>%
    mutate(segment.size = ifelse(!is.na(str_extract(seg.params, ".*105.*")), 105,
                          ifelse(!is.na(str_extract(seg.params, ".*60.*")), 60,
                          ifelse(!is.na(str_extract(seg.params, ".*30.*")), 30,
                          ifelse(!is.na(str_extract(seg.params, ".*70.*")), 105,
                          ifelse(!is.na(str_extract(seg.params, ".*40.*")), 60,
                          ifelse(!is.na(str_extract(seg.params, ".*20.*")), 30,1)))))))
## Calc\ Error\ Column:2 ends here

## [[file:utc.org::*RMSE%20analysis][RMSE\ analysis:1]]
options(asciiType = "org")
options(warn = -1)
  RMSE_tree %>%
      ungroup() %>%
      arrange(RMSE_t) %>%
      head(n = 30) %>%
      ascii()
## RMSE\ analysis:1 ends here

## [[file:utc.org::*RMSE%20analysis][RMSE\ analysis:2]]
ggplot(RMSE_tree, aes(x = n.points, y = RMSE_t, color = model)) + geom_point() +
    facet_grid(segment.size~image)
## RMSE\ analysis:2 ends here

## [[file:utc.org::*RMSE%20analysis][RMSE\ analysis:3]]
RMSE_tree.sub <- RMSE_tree%>%
    filter(segment.size == 60, image == "madisonNAIP", target.type == "binomial", model == "svm_resp") %>%
    mutate(area_meters_squared = ((sqrt(n.points) - 1) * 7)^2)


ggplot(RMSE_tree.sub, aes(x = area_meters_squared, y = RMSE_t), color = "blue") + geom_point() +
    labs(y = "Root Mean Squared Prediction Error \n for Percent Tree Cover") +
    theme_classic() +
    theme(axis.title = element_text(size = 24),
          axis.text =  element_text(size = 22)) +
    xlim(0,45000)
## RMSE\ analysis:3 ends here

## [[file:utc.org::*RMSE%20analysis][RMSE\ analysis:4]]
ggplot(RMSE_tree, aes(x = segment.size, y = RMSE_t, color = n.points, group = interaction(n.points,target.type))) + geom_line() +
    facet_grid(model~image)
## RMSE\ analysis:4 ends here

## [[file:utc.org::*RMSE%20analysis][RMSE\ analysis:5]]
m1 <-lm(RMSE_t*100 ~ image * (model +  target.type + n.points * segment.size), data = RMSE_tree)
  tidy(m1, digits = 2) %>%
ascii()
## RMSE\ analysis:5 ends here

## [[file:utc.org::*Save%20Best%20Model][Save\ Best\ Model:1]]
best.classif.overall <- error.df.avg.class %>%
    arrange(overall.error) %>%
    slice(1) %>%
    left_join(.,error.df) %>%
  mutate(path = paste0(image,"_",seg.params,"_",model,"_model.rds"))


best.model <- readRDS(paste0(dd.models.dir,"/",best.classif.overall$path))

saveRDS(best.model, paste0(dd.models.dir,"/best_mad_model_",best.classif.overall$image,"_",best.classif.overall$seg.params,"_",best.model$learner$id,".rds"))
## Save\ Best\ Model:1 ends here
