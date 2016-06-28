

# This should be run from the virtual environment in the project directory
import sys
from skimage.segmentation import slic
import gdal
import os
import numpy as np


def build_tiff(source_raster_path, target_raster_path, image_array):
    source_raster = gdal.Open(source_raster_path)
    columns = source_raster.RasterXSize
    rows = source_raster.RasterYSize
    driver = gdal.GetDriverByName("GTiff")
    new_image = driver.Create(target_raster_path, columns, rows, 1, gdal.GDT_Float64)
    new_image.SetGeoTransform(source_raster.GetGeoTransform())
    new_image.SetProjection(source_raster.GetProjectionRef())
    new_image.GetRasterBand(1).WriteArray(image_array)
#    print 'Tiff made successfully'
    return


def create_array_forSeg_fromTiff(source_raster_path, bandsToUse):
    source_raster = gdal.Open(source_raster_path)
    columns = source_raster.RasterXSize
    rows = source_raster.RasterYSize
    array = source_raster.ReadAsArray()
    array = array.transpose(1, 2, 0)
    array = array[:, :, bandsToUse]  # slice out bandsToUse
    array[array < 0] = 0 # Set values less than zero to zero
    return (rows, columns, array)


def seg_tiff(RowsColumnsArray, PixPerSeg, compactness):
    rows = RowsColumnsArray[0]
    columns = RowsColumnsArray[1]
    num_segs = rows * columns / PixPerSeg
    array = RowsColumnsArray[2].astype(np.uint8)
#    print 'Making segments'
    segments = slic(array, n_segments=num_segs, compactness=compactness, enforce_connectivity=True)
#    print 'Made dem segments'
    return segments

PixelSize = sys.argv[1]
AreaForSegment = sys.argv[2]
compactnessValue = sys.argv[3]
image_name = sys.argv[4]

PixPerSegValue = float(AreaForSegment) / float(PixelSize)


print 'average number of pixels per segment is ' + str(PixPerSegValue)
print 'compactness parameter is ' + compactnessValue


bandsToUse = slice(0,3) # Use only principal components 1-3

image_directory = os.getcwd()

files = os.listdir(image_directory)
files = [x for x in files if x.endswith('_pca.tif')]
files = [x for x in files if x.startswith(image_name)]

print (files)


slic_params = np.array([[int(round(PixPerSegValue)), int(compactnessValue)]])



for f in files :
    f_out = f[:-8]
    for pc in slic_params :
        source_raster_path = os.path.expanduser(f)
        target_raster_path = os.path.expanduser(f_out + '_N-%s_C-%s.tif') %(pc[0],pc[1])
        RowsColumnsArray = create_array_forSeg_fromTiff(f, bandsToUse)
        segments = seg_tiff(RowsColumnsArray, PixPerSeg = pc[0], compactness = pc[1])
        build_tiff(source_raster_path, target_raster_path, segments)






# # it would probably be better to create a directory for each combination of pix and compactness values
# # then save tiles into them.  Instead I append parameter values to the filename and it's clunky

# for f in files :
#     for pc in slic_params :
#         source_raster_path = os.path.expanduser(ScaledPCADirectory + '/' + f)
#         target_raster_path = os.path.expanduser(SegTileDirectory + '/' + f + '_N-%s_C-%s.tif') %(pc[0],pc[1])
#         RowsColumnsArray = create_array_forSeg_fromTiff(ScaledPCADirectory + "/" + f, bandsToUse)
#         segments = seg_tiff(RowsColumnsArray, PixPerSeg = pc[0], compactness = pc[1])
#         build_tiff(source_raster_path, target_raster_path, segments)



