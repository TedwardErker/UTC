# This should be run from the virtual environment in the project directory

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
    print 'Tiff made successfully'
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
    print 'Making segments'
    segments = slic(array, n_segments=num_segs, compactness=compactness, enforce_connectivity=True)
    print 'Made dem segments'
    return segments

# SegDirectory = '~/mydata/DD_NAIP-imagery/madison-SegTiles'

# os.mkdir(SegDirectory)

bandsToUse = slice(0,3) # Use only principal components 1-3

files = os.listdir("../DD_NAIP-imagery/madison-ScaledPCATiles")

#num_pix = np.array([[15,10], [30, 15], [60, 30], [105,32]])

#num_pix = np.array([[15,10],[105,32]])

#num_pix = np.array([[16,10]])
num_pix = np.array([[30, 15]])
#num_pix = np.array([[60, 30]])
#num_pix = np.array([[105,32]])
for f in files :
    for pc in num_pix :
        source_raster_path = os.path.expanduser('~/mydata/DD_NAIP-imagery/madison-ScaledPCATiles/' + f)
        target_raster_path = os.path.expanduser('~/mydata/DD_NAIP-imagery/madison-SegTiles/madison-' + f + '_N-%s_C-%s.tif') %(pc[0],pc[1])
        RowsColumnsArray = create_array_forSeg_fromTiff("../DD_NAIP-imagery/madison-ScaledPCATiles/" + f, bandsToUse)
        segments = seg_tiff(RowsColumnsArray, PixPerSeg = pc[0], compactness = pc[1])
        build_tiff(source_raster_path, target_raster_path, segments)



