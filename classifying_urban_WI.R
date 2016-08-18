## [[file:utc.org::*For%20the%20larger%20urban%20areas][For\ the\ larger\ urban\ areas:1]]
i_areas_high_quant <- which(areas[,1] >= high.quant)
## For\ the\ larger\ urban\ areas:1 ends here

## [[file:utc.org::*For%20the%20larger%20urban%20areas][For\ the\ larger\ urban\ areas:2]]
cl <- makeCluster(cores)
    registerDoParallel(cl)


out <- foreach(i = i_areas_high_quant[23],
              .packages = c("sp","raster","rgdal","rgeos", "stringr","doParallel","gdalUtils","plyr","dplyr","mlr","glcm")) %do% {

    urb.poly <- urb.polys[i]
## For\ the\ larger\ urban\ areas:2 ends here

## [[file:utc.org::*Set%20temp%20dir%20for%20this%20urban%20area][Set\ temp\ dir\ for\ this\ urban\ area:1]]
temp_i <- paste0(R_raster_temp,"/",i)
dir.create(temp_i)
rasterOptions(tmpdir=temp_i)
## Set\ temp\ dir\ for\ this\ urban\ area:1 ends here

## [[file:utc.org::*Make%20output%20dir%20for%20this%20urban%20area][Make\ output\ dir\ for\ this\ urban\ area:1]]
urb.path <- paste0(classified.urban.areas.dir, i)
dir.create(urb.path)
## Make\ output\ dir\ for\ this\ urban\ area:1 ends here

## [[file:utc.org::*Get%20names%20of%20NAIP%20tiles%20that%20intersect%20with%20Urban%20Area][Get\ names\ of\ NAIP\ tiles\ that\ intersect\ with\ Urban\ Area:1]]
tiles.in.urban <-  lapply(naip.extents, function(naip.extent) {
    inter <- raster::intersect(naip.extent, urb.poly)
    ifelse(is.null(inter), F, T)
})

tile.index <- which(unlist(tiles.in.urban))


tiles.inter.urb.poly <- lapply(tile.index, function(i) {
    naip.extent <- as(naip.extents[[i]], "SpatialPolygons")
    proj4string(naip.extent) <- wtm
    inter <- gIntersects(naip.extent, urb.poly)
})

tile.index <- tile.index[which(unlist(tiles.inter.urb.poly))]

tiles.names.at.urb.poly <- naip.tif.names[tile.index]
## Get\ names\ of\ NAIP\ tiles\ that\ intersect\ with\ Urban\ Area:1 ends here

## [[file:utc.org::*For%20each%20NAIP%20tile%20that%20intersects%20with%20Urban%20Area][For\ each\ NAIP\ tile\ that\ intersects\ with\ Urban\ Area:1]]
foreach(t.p = tiles.names.at.urb.poly,
                .packages = c("sp","raster","rgdal","rgeos", "stringr","doParallel","gdalUtils","plyr","dplyr","mlr","glcm")) %dopar% {
## For\ each\ NAIP\ tile\ that\ intersects\ with\ Urban\ Area:1 ends here

## [[file:utc.org::*Make%20output%20dir%20for%20this%20tile][Make\ output\ dir\ for\ this\ tile:1]]
tile.name <- basename(t.p) %>%
    str_sub(start = 1, end = -5)  # remove .tif
  tile.urb.path <- paste0(urb.path,"/",tile.name)
  dir.create(tile.urb.path)
message("make tile output dir", tile.urb.path)
## Make\ output\ dir\ for\ this\ tile:1 ends here

## [[file:utc.org::*Crop%20to%20intersection%20of%20image%20and%20Urban%20Extent][Crop\ to\ intersection\ of\ image\ and\ Urban\ Extent:1]]
# Crop image
eu <- extent(urb.poly)
ei <- extent(raster(t.p))
e <- raster::intersect(ei,eu)

inFile <- t.p
outFile <- str_c(tile.urb.path,"/urbanExtent.tif")

gdal_translate(inFile, outFile,
               projwin = c(xmin(e), ymax(e), xmax(e), ymin(e)))


message("Crop to Urban Extent")
## Crop\ to\ intersection\ of\ image\ and\ Urban\ Extent:1 ends here

## [[file:utc.org::*make%20sure%20the%20intersection%20image%20has%20datavalues][make\ sure\ the\ intersection\ image\ has\ datavalues:1]]
r.test <- raster(paste0(tile.urb.path,"/urbanExtent.tif"))

if (sum(!is.na(values(r.test))) > 500) {
## make\ sure\ the\ intersection\ image\ has\ datavalues:1 ends here

## [[file:utc.org::*Generate%20Feature%20data%20frame][Generate\ Feature\ data\ frame:1]]
feature.dfs <- make.feature.df(tile.dir = tile.urb.path,
                                 tile.name = "urbanExtent",
                                 image.name = "NAIP",
                                 band.names = c("blue","green","red","nir"),
                                 ndvi = T,
                                 ratio.bands = c("blue","green","red","nir"),
                                 texture.params.df = texture.params,
                                 pixel.df = F,
                                 pca.location = location,
                                 segmentation = T,
                                 segment.params.df = segment.params)
## Generate\ Feature\ data\ frame:1 ends here

## [[file:utc.org::*Delete%20intermediate%20files][Delete\ intermediate\ files:1]]
intermediate.work <- list.files(tile.urb.path, full.names = T, recursive = T)
  ratio.intermediate.work <- str_extract(intermediate.work, ".*_ratio.*")
  band.intermediate.work <- str_extract(intermediate.work, ".*_(red|blue|green|nir).*")
unlink(c(ratio.intermediate.work, band.intermediate.work))
## Delete\ intermediate\ files:1 ends here

## [[file:utc.org::*Classify][Classify:1]]
model <- list.files(dd.models.dir) %>%
          str_extract("best_.*") %>% na.omit()

      seg.p <- str_extract(model, segmentation.pattern)

      if(grepl("N-[0-9]+_C-[0-9]+",seg.p)) {
          segment.tile.name.append <- paste0("_",seg.p,".tif")
          segment.feature.df.name.append <- paste0("_",seg.p,feature.df.appendage,".rds")

          classify.segmented.raster(segment.feature.df.dir = tile.urb.path,
                                    model.dir = dd.models.dir,
                                    segment.dir = tile.urb.path,
                                    classify.out.dir = tile.urb.path,
                                    tile.name = "urbanExtent",
                                    segmentation.appendage = segment.tile.name.append,
                                    model.name.rds = model,
                                    segment.feature.appendage = segment.feature.df.name.append,
                                    segmentation.prms = seg.p)

      } else {
          classify.pixel.raster(tile.dir = tile.urb.path,
                                tile.name = "urbanExtent",
                                pixelFeatureDF.appendage = feature.df.appendage,
                                model.dir = dd.models.dir,
                                model.rds = model,
                                classify.out.dir = tile.urb.path,
                                seg.prms = seg.p)
      }
      message("Image Classified")

      message("Done with",tile.urb.path)

}
}
## Classify:1 ends here

## [[file:utc.org::*Merge%20NAIP%20Tiles%20if%20there%20is%20more%20than%20one%20over%20an%20urban%20area%20and%20Save%20Classified%20image%20as%20<UrbanArea>.tif][Merge\ NAIP\ Tiles\ if\ there\ is\ more\ than\ one\ over\ an\ urban\ area\ and\ Save\ Classified\ image\ as\ <UrbanArea>\.tif:1]]
image.name <- "NAIP"

classified.tiles <- list.files(urb.path, recursive = T, full.names = T) %>%
    str_extract(.,     pattern = paste0(".*(",segmentation.pattern,")_",location,image.name,"_(", target.pattern,")_(", tuned.pattern, ")_(",model.pattern,").*.tif$")) %>%
    na.omit()

rlist <- lapply(classified.tiles, stack)

if(length(rlist) > 1) {
    out <- do.call(mosaic, c(rlist,list(fun = sample.with.na.rm, tolerance = 0.5)))
} else {
    out <- rlist[[1]]
}

writeRaster(x = out, filename = paste0(urb.path,"_ClassifiedUrbanArea.tif"), overwrite = T, datatype = 'INT1U')

message("Wrote ","ClassifiedUrbanArea_",i,".tif")


rgbn.tiles <- list.files(urb.path, recursive = T, full.names = T) %>%
    str_extract(.,     pattern = ".*urbanExtent.tif") %>%
    na.omit()

rlist <- lapply(rgbn.tiles, stack)

if(length(rlist) > 1) {
out <- do.call(mosaic, c(rlist,list(fun = mean, tolerance = 0.5)))
} else {
    out <- rlist[[1]]
}



writeRaster(x = out, filename = paste0(urb.path,"_rgbn.tif"), overwrite = T)
message("Wrote ","RGBN_",i,".tif")
## Merge\ NAIP\ Tiles\ if\ there\ is\ more\ than\ one\ over\ an\ urban\ area\ and\ Save\ Classified\ image\ as\ <UrbanArea>\.tif:1 ends here

## [[file:utc.org::*Delete%20intermediate%20steps][Delete\ intermediate\ steps:1]]
intermediate.work <- list.files(urb.path, full.names = T)
intermediate.work <- intermediate.work[!grepl(intermediate.work,pattern = ".*(tif)$", perl = T)]
unlink(intermediate.work, recursive = T)
unlink(urb.path)
unlink(temp_i, recursive = T)
## Delete\ intermediate\ steps:1 ends here

## [[file:utc.org::*End%20if%20statement%20and%20Loop%20for%20all%20urban%20areas][End\ if\ statement\ and\ Loop\ for\ all\ urban\ areas:1]]
}
}
}
## End\ if\ statement\ and\ Loop\ for\ all\ urban\ areas:1 ends here
