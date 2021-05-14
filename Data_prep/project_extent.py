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


class SetExtent:
    """ SetExtent includes creates extent shapefiles
    Attributes
    ----------
    shp: shapefile
        the boundary shapefile
    """

    def __init__(self, shp):
        self.shp = shp
        self.crs = None


    def create_buffer(self, buffer_distance):
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


    def create_extent(self, epsg=None):
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


if __name__ == 'project_extent':
    shp_file = '../AB2551Watersheds/AB2551Watersheds.shp'
    shape = gpd.read_file(shp_file)
    data = SetExtent(shape)
    extent_w_buffer = data.create_buffered_extent(10000, epsg=3310)
    extent = data.create_extent(epsg=3310)
    if 'Extent' not in os.listdir('..'):
        os.mkdir('../Extent')
    extent_w_buffer.to_file("../Extent/buffered_extent.shp")
    extent.to_file("../Extent/extent.shp")

    # base = shape.plot(color='white', edgecolor='black')
    # shape.plot(ax=base)
