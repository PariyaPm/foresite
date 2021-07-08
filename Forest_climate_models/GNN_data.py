import sys
import os
import gc
import geopandas as gpd
sys.path.append('Data_prep')
import project_extent as pe


__author__ = "Pariya Pourmohammadi"
__date__ = "06/07/2021"
__credits__ = ["Pariya Pourmohammadi"]
__version__ = "01.0"
__maintainer__ = "Pariya Pourmohammadi"
__email__ = "ppourmohammadi@ucmerced.edu"
__status__ = "Finalized"


class GNN:
    """ GNN creates specific GNN layer from a list of layers
    and visualized them
    Attributes
    ----------

    """

    def __init__(self,
                 GNN_files=None,
                 gnn_extent=None,
                 gnn_boundary=None,
                 crs=None,
                 filetype=None,
                 boundary=None,
                 in_dir=None,
                 out_dir=None):
        """

        Parameters
        ----------
        GNN_files
        gnn_extent
        gnn_boundary
        crs
        filetype
        boundary
        in_dir
        out_dir
        """


        self.GNN_files = GNN_files
        self.gnn_extent = gnn_extent
        self.gnn_boundary = gnn_boundary
        self.crs = crs
        self.filetype = filetype
        self.boundary = boundary
        self.in_dir = in_dir
        self.out_dir = out_dir


    def read_gnn(self):
        """

        Returns
        -------

        """
        for file in self.GNN_files:
            if file[-4:] == '.'+self.filetype:
                temp = pe.clip_to_extent(
                    area_file=os.path.join(self.path, file),
                    extent_file=self.gnn_extent)
                pe.write_raster(crs=self.crs,
                    driver=self.driver,
                    file=temp,
                    f_name=file,
                    dir_name=self.out_dir)
            print(file)
            temp = None
            gc.collect()


    def visualize_gnn(self):
        """

        Returns
        -------

        """
        for new_file in self.GNN_files:
            if new_file[-4:] == '.'+self.filetype:
                pe.visualize_raster(raster_lyr=self.in_dir + new_file,
                    bound=self.boundary,
                    fig_name=new_file[:-4],
                    path_to_fig=self.out_dir + new_file[:-4])
            gc.collect()


def rescale_LEMMA(file_name, out_name, scaler):
    """

    Parameters
    ----------
    file_name
    out_name
    scaler

    Returns
    -------

    """
    import rasterio

    file = pe.read_raster_array(file_name)
    array_out = file / scaler
    import rasterio

    with rasterio.open(file_name) as src:
        meta = src.meta

    meta.update(dtype=rasterio.float32, crs='EPSG:3310')

    # Write output file
    with rasterio.open(out_name, 'w', **meta) as dst:
        dst.write(array_out.astype(rasterio.float32), 1)

    print('file %s is rescales to 1/%d and saved to %s' % (
    file_name, scaler,
    out_name))


if __name__ == 'GNN_data':
    if 'GNN_AB2551' not in os.listdir():
        os.mkdir('GNN_AB2551')
    GNN_dat = GNN(GNN_files=os.listdir('GNN_AB2551'),
        gnn_extent='Extent/buffered_extent.shp',
        gnn_boundary=gpd.read_file('AB2551Watersheds/AB2551Watersheds.shp'),
        filetype='tif',
        boundary=gpd.read_file('AB2551Watersheds/AB2551Watersheds.shp'),
        crs='EPSG:3310',
        in_dir='GNN_AB2551/',
        out_dir='GNN_AB2551/')


