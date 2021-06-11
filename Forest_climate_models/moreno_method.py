#!/usr/bin/env python

"""project_extent: this code creates an extent from the boundary, it includes
    functions for creating buffer and creating the extent with and without
    buffer
"""

__author__ = "Pariya Pourmohammadi"
__date__ = "05/14/2021"
__credits__ = ["Pariya Pourmohammadi"]
__version__ = "01.0"
__maintainer__ = "Pariya Pourmohammadi"
__email__ = "ppourmohammadi@ucmerced.edu"
__status__ = "Draft"
#!/usr/bin/env python

"""project_extent: this code creates an extent from the boundary, it includes
    functions for creating buffer and creating the extent with and without
    buffer
"""

__author__ = "Pariya Pourmohammadi"
__date__ = "05/06/2021"
__credits__ = ["Pariya Pourmohammadi"]
__version__ = "01.0"
__maintainer__ = "Pariya Pourmohammadi"
__email__ = "ppourmohammadi@ucmerced.edu"
__status__ = "Finalized"

import os
import geopandas as gpd
from geopandas import GeoSeries
from shapely.geometry import Polygon
import rioxarray as rxr
from shapely.geometry import mapping
import rasterio as rio
from rasterio.features import rasterize
from rasterio.transform import from_bounds
import affine
from rasterio.warp import reproject, Resampling, aligned_target
import numpy as np


class SetExtent:
    """ SetExtent includes creates extent shapefiles
    Attributes
    ----------
    shp: shapefile
        the boundary shapefile
    """

    def __init__(self, shp=None):
        self.shp = shp
        self.crs = None


    def create_buffer(self, buffer_distance=10000):
        """
        Parameters
        ----------
        buffer_distance: num
            distance of buffer

        Returns
        -------
        buffered: shp
            buffered shapefile
        """
        buffered = self.shp.buffer(buffer_distance)
        return buffered


    def create_extent(self, epsg=3310):
        """

        Parameters
        ----------
        epsg: int optional
            espg of the projection scr, if this is not identified the
            projection of out file will be the same as infile

        Returns
        -------
        out_geo_poly: shp
            rectangular vector shape of the extent
        """
        shape_file = self.shp
        bounds = shape_file.total_bounds
        new_poly = Polygon([(bounds[0], bounds[1]),
                            (bounds[2], bounds[1]),
                            (bounds[2], bounds[3]),
                            (bounds[0], bounds[3])])
        geo_poly = GeoSeries([new_poly])
        if epsg is None:
            epsg = self.crs
        out_geo_poly = geo_poly.set_crs(epsg=epsg)

        return out_geo_poly


    def create_buffered_extent(self, buffer_distance, epsg=None):
        """

        Parameters
        ----------
        buffer_distance: num
            distance of buffer
        epsg: int, optional
            espg of the projection scr, if this is not identified the
            projection of out file will be the same as infile

        Returns
        -------
        extended: shp
            rectangular vector shape of the extent with the designated buffer
        """

        buffered_shp = self.create_buffer(buffer_distance)
        if epsg is None:
            epsg = self.crs
        buffered = SetExtent(buffered_shp)
        extended = buffered.create_extent(epsg)
        return extended


def project_file(area, ref=None, epsg=None):
    """
    Parameters
    ----------
    area
    ref
    epsg

    Returns
    -------
    projected
    """

    global projected
    if epsg is not None:
        projected = area.rio.reproject('EPSG:' + epsg)

    else:
        try:
            projected = area.rio.reproject(ref.crs)
        except ValueError:
            print("Set the reference file or epsg! Try again ...")

    return projected


def clip_to_extent(area_file, extent_file):
    """
    This function reprojects the file to the extent and clips it to the polygon
     that is passed to the function

    Parameters
    ----------
    area_file: rasterio object
        the area that is desired to be clipped
    extent_file: poly object
        the extent based on which the area_file gets clipped


    Returns
    -------
    target_file_clipped: rio object
        the clipped rio object
    """

    target_file = rxr.open_rasterio(area_file, masked=True)
    clip_extent = gpd.read_file(extent_file)
    target_file_proj = target_file.rio.reproject(clip_extent.crs)
    target_file_clipped = \
        target_file_proj.rio.clip(clip_extent.geometry.apply(mapping))

    return target_file_clipped


def write_raster(crs='EPSG:3310',
                 file=None,
                 driver='GTiff',
                 f_name=None,
                 dir_name=None,
                 no_data=-9999):
    """
    Parameters
    ----------
    crs: str
        string of desired a coordinate reference system identifier
        or description , this string should be in format of EPSG:value
    driver: str
        the name of the desired format driver
    file: rio file
        rasterio file with attributes of raster data
    f_name: str
        the desired name to be assigned to the file
    dir_name: str
        the desired directory name for the file to be written to
    no_data: num
        Value assigned to nodata

    Returns
    -------

    """
    dat_type = rio.float64

    if file.dtype == 'float64':
        dat_type = rio.float64
    elif file.dtype == 'float32':
        dat_type = rio.float32
    elif file.dtype == 'int32':
        dat_type = rio.int32
    elif file.dtype == 'int16':
        dat_type = rio.int16
    elif file.dtype == 'uint16':
        dat_type = rio.uint16
    elif file.dtype == 'int8':
        dat_type = rio.int8
    elif file.dtype == 'uint8':
        dat_type = rio.uint8

    out_meta = {'count': file.shape[0],
                'crs': crs,
                'dtype': dat_type,
                'height': file.shape[1],
                'width': file.shape[2],
                'driver': driver,
                'transform': file.rio.transform(),
                'nodata': no_data}

    if f_name not in dir_name:
        with rio.open(dir_name + '/' + f_name, 'w', **out_meta) \
                as dst:
            dst.write(file)
            print(f_name + ' is written to ' + dir_name)


def raster_extent(boundary,
                  extent_out,
                  driver='GTiff',
                  dtype=rio.uint8,
                  width=1000,
                  height=1000,
                  crs='EPSG:3310',
                  shape_from_bounds=True,
                  res=(30.0, 30.0)):
    """

    Parameters
    ----------
    boundary: shp
        boundary vector data     
    extent_out: str
        name of out file
    driver: str
        driver
    dtype: rio dtype
        out data type
    width: int
        width of the temp/adjustable
    height: int
        height of the temp/adjustable
    crs: str
        crs of the out file
    shape_from_bounds: boolean
        if the shepe should be created using the boundary
    res: tuple
        our resolution

    Returns
    -------

    """""

    # Load boundary
    df = gpd.read_file(boundary)

    in_shape = width, height

    # use bounds to designate a shape of 1m resolution
    if shape_from_bounds:
        bounds = df.total_bounds
        in_shape = round((bounds[3] - bounds[1]) / 30), \
                   round((bounds[2] - bounds[0]) / 30)

        height = in_shape[0]
        width = in_shape[1]

    transform = from_bounds(*df.total_bounds, width, height)

    rasterized_dat = rasterize(
        [(shape, 1) for shape in df['geometry']],
        out_shape=in_shape,
        fill=0,
        transform=transform,
        all_touched=False,
        dtype=rio.uint8)

    with rio.open(
            'temp.tif',
            'w',
            driver=driver,
            dtype=dtype,
            count=1,
            width=width,
            height=height,
            transform=transform,
            crs=crs
    ) as dst:
        dst.write(rasterized_dat, indexes=1)

    target_file_clipped = clip_to_extent('temp.tif',
        boundary)
    write_raster(crs='EPSG:3310',
        driver='GTiff',
        file=target_file_clipped,
        f_name='new_temp.tif',
        dir_name='.')

    resample(in_file='new_temp.tif',
        dst_file=extent_out,
        new_res=res,
        crs=crs)

    os.remove('temp.tif')
    os.remove('new_temp.tif')


def resample(in_file,
             dst_file,
             src_file=None,
             new_res=(30, 30),
             crs='EPSG:3310',
             snap=False):
    """

    Parameters
    ----------
    in_file: str
        name of the tif file to be projected
    dst_file: str
        name of the destination file
    src_file
    new_res: tuple
        new resolution
    crs: str
        String in form of EPSG:3310
    snap

    Returns
    -------

    """

    # in_file = 'Climate_limits/prism_tmin_Dec_Feb.tif'
    # dst_file = 'Climate_limits/resampled_py_prism_tmin.tif'
    # src_file = 'Extent/buffered_extent.tif'
    # snap = True
    # read the source raster
    with rio.open(in_file) as src:
        array = src.read()
        old_resolution = src.res
        aff = src.transform

        if snap:
            with rio.open(src_file) as src_dat:
                new_aff = aligned_target(transform=src.transform,
                    width=src_dat.shape[1],
                    height=src_dat.shape[0],
                    resolution=src_dat.res)
                # x_res_ratio = old_resolution[0] / src_dat.res[0]
                # y_res_ratio = old_resolution[1] / src_dat.res[1]
                new_aff = new_aff[0]
                new_array = np.empty(shape=(
                    array.shape[0], *src_dat.shape))

        else:
            new_res = new_res
            # setup the transform to change the resolution
            x_res_ratio = old_resolution[0] / new_res[0]
            y_res_ratio = old_resolution[1] / new_res[1]
            new_aff = affine.Affine(aff.a / x_res_ratio, aff.b, aff.c, aff.d,
                                    aff.e / y_res_ratio, aff.f)

            new_array = np.empty(shape=(
                array.shape[0], int(round(array.shape[1] * x_res_ratio)),
                int(round(array.shape[2] * y_res_ratio))))

        new_array = reproject(array,
            new_array,
            src_transform=aff,
            dst_transform=new_aff,
            src_crs=src.crs,
            dst_crs=crs,
            count=1,
            resample=Resampling.bilinear)

        # write results to file
        with rio.open(dst_file,
                'w',
                driver=src.driver,
                height=new_array[0].shape[1],
                width=new_array[0].shape[2],
                nodata=src.nodata,
                dtype=src.dtypes[0],
                count=1,
                crs=crs,
                transform=new_aff) as dst:
            dst.write(new_array[0])


if __name__ == '__main__':
    shp_file = 'AB2551Watersheds/AB2551Watersheds.shp'
    shape = gpd.read_file(shp_file)
    data = SetExtent(shape)
    extent_w_buffer = data.create_buffered_extent(10000, epsg=3310)
    extent = data.create_extent(epsg=3310)
    if 'Extent' not in os.listdir():
        os.mkdir('Extent')
    extent_w_buffer.to_file('Extent/buffered_extent.shp')
    extent.to_file('Extent/extent.shp')
    raster_extent(boundary='Extent/buffered_extent.shp',
        extent_out='Extent/buffered_extent.tif')
    raster_extent(boundary='Extent/extent.shp',
        extent_out='Extent/extent.tif')
    
