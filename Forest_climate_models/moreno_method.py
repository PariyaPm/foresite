#!/usr/bin/env python

"""moreno_method: this code creates analyzes the climate forest structure
relationships
"""

__author__ = "Pariya Pourmohammadi"
__date__ = "05/14/2021"
__credits__ = ["Pariya Pourmohammadi"]
__version__ = "01.0"
__maintainer__ = "Pariya Pourmohammadi"
__email__ = "ppourmohammadi@ucmerced.edu"
__status__ = "Draft"

import sys
import geopandas as gpd
import matplotlib.pyplot as plt
import rioxarray as rxr
import rasterio as rio
import gc
import pandas as pd
import numpy as np
import os
from glob import glob

sys.path.append('Data_prep')

import project_extent as pe

# TODO add the params
# TODO  complete the code documentation

def read_raster_array(raster_file):
    """

    Parameters
    ----------
    raster_file

    Returns
    -------

    """
    with rio.open(raster_file) as src:
        return np.array(src.read(1))


def avg_calc(path_list, out_name):
    """

    Parameters
    ----------
    path_list
    out_name

    Returns
    -------

    """
    import rasterio

    lyr_list = [read_raster_array(x) for x in path_list]
    array_out = np.mean(lyr_list, axis=0)

    with rasterio.open(path_list[0]) as src:
        meta = src.meta

    meta.update(dtype=rasterio.float32)

    # Write output file
    with rasterio.open(out_name, 'w', **meta) as dst:
        dst.write(array_out.astype(rasterio.float32), 1)

    return array_out


def visualize_raster(raster_lyr,
                     bound,
                     fig_name,
                     path_to_fig,
                     axes=False):
    """
    Parameters
    ----------
    raster_lyr
    bound
    fig_name
    path_to_fig
    axes

    Returns
    -------

    """

    # array_out_rio = rxr.open_rasterio('Climate_limits/resampled_py_height.tif')
    array_out_rio = rxr.open_rasterio(raster_lyr)
    #
    f, ax = plt.subplots(figsize=(10, 4))

    array_out_rio.plot(ax=ax)
    bound.boundary.plot(ax=ax, color='black')
    # bound.boundary.plot(ax=ax, color='black')

    ax.set(title=fig_name)
    if not axes:
        ax.set_axis_off()
    plt.show()

    plt.savefig(path_to_fig)
    plt.close()
    print(fig_name + ' is plotted and saved')


def calc_bp(data):
    d = np.array(data)
    dict_res = {'median': np.median(d), 'upper_quartile': np.percentile(d, 75),
                'lower_quartile': np.percentile(d, 25)}

    iqr = np.percentile(d, 75) - np.percentile(d, 25)

    dict_res['iqr'] = iqr
    dict_res['upper_whisker'] = d[d <= np.percentile(d, 75) + 1.5 * iqr].max()
    dict_res['lower_whisker'] = d[d >= np.percentile(d, 25) - 1.5 * iqr].min()
    return dict_res


def plot_lines(x,
               y,
               x_lab,
               y_lab,
               title,
               name, ):
    plt.plot(x, y)
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)

    plt.grid()
    # displaying the title
    plt.title(title)
    plt.savefig(name)
    plt.show()
    plt.close()


def plt_scatter(var1, var2, var3,dat1, dat2, dat3,title,path_to_fig):
    one = pd.DataFrame({var1: dat1.flatten(),
                        var2: dat2.flatten(),
                        var3: dat3.flatten()})
    one.dropna()
    one.round(decimals=1)
    one.drop_duplicates()
    t = one.groupby([var1, var2]).mean().reset_index()
    fig, ax = plt.subplots()
    plt.scatter(x=t[var1], y=t[var2], c=t[var3],
        marker='s', s=0.1)
    plt.xlabel(var1)
    plt.ylabel(var2)
    plt.title(title)
    ax.set_xlabel(var1)
    ax.set_ylabel(var2)
    cbar = plt.colorbar()
    cbar.set_label(var3)
    plt.show()
    plt.savefig(path_to_fig)
    plt.close()

# dirs = os.listdir('../../Documents/ForesiteProject/Data/GNN')
# files= Forest_structure_lst
files = os.listdir('GNN_AB2551')

if 'GNN_AB2551' not in os.listdir():
    os.mkdir('GNN_AB2551')

extent = 'Extent/buffered_extent.shp'

for file in files:
    if file[-4:] == '.tif':
        temp = pe.clip_to_extent(
            area_file=os.path.join(
                '../../Documents/ForesiteProject/Data/GNN/rasters', file),
            extent_file=extent)
        pe.write_raster(crs='EPSG:3310',
            driver='GTiff',
            file=temp,
            f_name=file,
            dir_name='GNN_AB2551')
    print(file)
    temp = None
    gc.collect()


# plot each of the forest structure layers
GNN_files = os.listdir('GNN_AB2551')

boundary = gpd.read_file('AB2551Watersheds/AB2551Watersheds.shp')
for new_file in GNN_files:
    if new_file[-4:] == '.tif':
        visualize_raster(raster_lyr='GNN_AB2551/' + new_file,
            bound=boundary,
            fig_name=new_file[:-4],
            path_to_fig='GNN_AB2551/' + new_file[:-4])
    gc.collect()

# call tmin avg
tmin1 = glob(os.path.join(
    '../../Documents/ForesiteProject/Data/Climate_data/',
    'prism_tmin*12.tif'))

tmin2 = glob(os.path.join(
    '../../Documents/ForesiteProject/Data/Climate_data/',
    'prism_tmin*01.tif'))
tmin3 = glob(os.path.join(
    '../../Documents/ForesiteProject/Data/Climate_data/',
    'prism_tmin*02.tif'))

tmin_list = tmin1 + tmin2 + tmin3
tmin_avg = avg_calc(tmin_list, 'Climate_limits/prism_tmin_Dec_Feb.tif')

pe.resample(in_file='Climate_limits/prism_tmin_Dec_Feb.tif',
    dst_file='Climate_limits/resampled_py_prism_tmin.tif',
    src_file='Extent/buffered_extent.tif',
    snap=True)

visualize_raster(raster_lyr='Climate_limits/resampled_py_prism_tmin.tif',
    bound=boundary,
    fig_name='Average Minimum Temperature(°C)',
    path_to_fig='Graphics/resampled_py_prism_tmin')

# call tmax avg
tmax1 = glob(os.path.join(
    '../../Documents/ForesiteProject/Data/Climate_data/',
    'prism_tmax*06.tif'))

tmax2 = glob(os.path.join(
    '../../Documents/ForesiteProject/Data/Climate_data/',
    'prism_tmax*07.tif'))
tmax3 = glob(os.path.join(
    '../../Documents/ForesiteProject/Data/Climate_data/',
    'prism_tmax*08.tif'))

tmax_list = tmax1 + tmax2 + tmax3
tmax_avg = avg_calc(tmax_list, 'Climate_limits/prism_tmax_jun_aug.tif')
pe.resample(in_file='Climate_limits/prism_tmax_jun_aug.tif',
    dst_file='Climate_limits/resampled_py_prism_tmax.tif',
    src_file='Extent/buffered_extent.tif',
    snap=True)
visualize_raster(raster_lyr='Climate_limits/resampled_py_prism_tmax.tif',
    bound=boundary,
    fig_name='Average Maximum Temperature(°C)',
    path_to_fig='Graphics/resampled_py_prism_tmax')

# call precip avg
pecip_list = glob(os.path.join(
    '../../Documents/ForesiteProject/Data/Climate_data/',
    'prism_ppt*.tif'))
precip_avg = avg_calc(pecip_list, 'Climate_limits/prism_ppt.tif')
pe.resample(in_file='Climate_limits/prism_ppt.tif',
    dst_file='Climate_limits/resampled_py_prism_ppt.tif',
    src_file='Extent/buffered_extent.tif',
    snap=True)
visualize_raster(raster_lyr='Climate_limits/resampled_py_prism_ppt.tif',
    bound=boundary,
    fig_name='Average Precipitation(mm)',
    path_to_fig='Graphics/resampled_py_prism_ppt')

"""
stndhgt_2017.tif: "average" Stand height, computed as average of heights of 
all dominant
and codominant trees (m/ha)
--Sum tree height multiplied by TPH for dominant and codominant trees 
(identified by the CROWN_CLASS field)
HT_SUM = ∑ [TREE_LIVE.HT_M * TPH_FC] (for DBH_CM >= 2.5 and CROWN_CLASS <= 3)
--Sum TPH for dominant and codominant trees
TPH_SUM = ∑ TPH_FC (for DBH_CM >= 2.5 and CROWN_CLASS <= 3)
--Divide sum of tree heights by TPH sum 


ba_ge_3_2017.tif: Basal area of live trees >=2.5 cm dbh
∑ TREE_LIVE.BAPH_FC (for DBH_CM >= 2.5) 

mndbhba_2017.tif: Basal-area weighted mean diameter of all live trees(cm)
--Sum basal area per hectare multiplied by DBH_CM for all live trees
BAPH_DBH_SUM = ∑ [TREE_LIVE.BAPH_FC * DBH_CM] (for DBH_CM >= 2.5)
--Sum basal area per hectare for all live trees
BAPH_SUM = ∑ BAPH_FC (for DBH_CM >= 2.5)
--Calculate basal-area weighted mean diameter
MNDBHBA = BAPH_DBH_SUM / BAPH_SUM 
Ref: LEMMA
"""
Forest_structure_lst = ['stndhgt_2017.tif', 'ba_ge_3_2017.tif',
                        'mndbhba_2017.tif']

pe.resample(in_file=os.path.join('GNN_AB2551',
    Forest_structure_lst[0]),
    dst_file='Climate_limits/resampled_py_height.tif',
    src_file='Extent/buffered_extent.tif',
    snap=True)
visualize_raster(raster_lyr='Climate_limits/resampled_py_height.tif',
    bound=boundary,
    fig_name='Stand Height(m)',
    path_to_fig='Graphics/resampled_py_height')

pe.resample(in_file=os.path.join('GNN_AB2551',
    Forest_structure_lst[1]),
    dst_file='Climate_limits/resampled_py_basal.tif',
    src_file='Extent/buffered_extent.tif',
    snap=True)
visualize_raster(raster_lyr='Climate_limits/resampled_py_basal.tif',
    bound=boundary,
    fig_name='Tree Basal Area(m^2/ha)',
    path_to_fig='Graphics/resampled_py_basal')

pe.resample(in_file=os.path.join('GNN_AB2551',
    Forest_structure_lst[2]),
    dst_file='Climate_limits/resampled_py_diameter.tif',
    src_file='Extent/buffered_extent.tif',
    snap=True)
visualize_raster(raster_lyr='Climate_limits/resampled_py_diameter.tif',
    bound=boundary,
    fig_name='Basal-area weighted mean diameter(cm)',
    path_to_fig='Graphics/resampled_py_diameter')

max_temp = read_raster_array('Climate_limits/resampled_py_prism_tmax.tif')
min_temp = read_raster_array('Climate_limits/resampled_py_prism_tmin.tif')
precip = read_raster_array('Climate_limits/resampled_py_prism_ppt.tif')

bhd = read_raster_array('Climate_limits/resampled_py_diameter.tif')
basal = read_raster_array('Climate_limits/resampled_py_basal.tif')
height = read_raster_array('Climate_limits/resampled_py_height.tif')

dict_climate_struct = {'max_temp': max_temp.flatten(),
                       'min_temp': min_temp.flatten(),
                       'precipitation': precip.flatten(),
                       'tree_diameter': bhd.flatten(),
                       'basal_area': basal.flatten(),
                       'tree_height': height.flatten()}

climate_struct_df = pd.DataFrame(dict_climate_struct)

dat = climate_struct_df.query('max_temp !=0 and min_temp!=0 and '
                              'precipitation !=0 ')
dat = dat.dropna()
dat = dat[dat['tree_diameter'] > 0]
dat = dat[dat['basal_area'] > 0]
dat = dat[dat['tree_height'] > 0]

new_dat = dat.round(decimals=2)
unique_vals = new_dat.T.apply(lambda x: x.nunique(), axis=1)

p_range_val = np.array(np.arange(round(np.min(new_dat.precipitation)),
    round(np.max(new_dat.precipitation) + 1), 10))
precip = {}
per_vals = {}

_t_dbh = {'median': 0.0, 'upper_quartile': 0.0, 'lower_quartile': 0.0,
          'iqr': 0.0, 'upper_whisker': 0.0, 'lower_whisker': 0.0}
_t_diam = {'median': 0.0, 'upper_quartile': 0.0, 'lower_quartile': 0.0,
           'iqr': 0.0, 'upper_whisker': 0.0, 'lower_whisker': 0.0}
_t_h = {'median': 0.0, 'upper_quartile': 0.0, 'lower_quartile': 0.0,
        'iqr': 0.0, 'upper_whisker': 0.0, 'lower_whisker': 0.0}

for i in p_range_val[:-1]:
    tmp = new_dat[new_dat['precipitation'].between(i, i + 1, inclusive=True)]
    # print(tmp)

    precip['dbh_{}'.format(i)] = tmp.tree_diameter
    precip['basal_{}'.format(i)] = tmp.basal_area
    precip['h_{}'.format(i)] = tmp.tree_height

    if len(tmp.tree_diameter) > 0:
        _t_dbh = calc_bp(tmp.tree_diameter)
        _t_diam = calc_bp(tmp.basal_area)
        _t_h = calc_bp(tmp.tree_height)

    else:
        _t_dbh['upper_whisker'] = 0
        _t_diam['upper_whisker'] = 0
        _t_h['upper_whisker'] = 0

    if 'precip' in per_vals:
        per_vals['precip'].append(i)
        per_vals['upper_whisk_dbh'].append(_t_dbh['upper_whisker'])
        per_vals['upper_whisk_basal'].append(_t_diam['upper_whisker'])
        per_vals['upper_whisk_height'].append(_t_h['upper_whisker'])

    else:
        per_vals['precip'] = [i]
        per_vals['upper_whisk_dbh'] = [_t_dbh['upper_whisker']]
        per_vals['upper_whisk_basal'] = [_t_diam['upper_whisker']]
        per_vals['upper_whisk_height'] = [_t_h['upper_whisker']]
    print(i)

plot_lines(x=per_vals['precip'],
    y=per_vals['upper_whisk_dbh'],
    x_lab='Precipitation(mm)',
    y_lab='Upper whisker of basal-area weighted mean diameter(cm)',
    title='Basal-area weighted mean diameter limits vs. precipitation',
    name='Graphics/dbh_prcp')

plot_lines(x=per_vals['precip'],
    y=per_vals['upper_whisk_basal'],
    x_lab='Precipitation(mm)',
    y_lab='Upper whisker of basal areas (m^2/ha)',
    title='Basal area limits vs. precipitation',
    name='Graphics/basal_prcp')

plot_lines(x=per_vals['precip'],
    y=per_vals['upper_whisk_height'],
    x_lab='Precipitation(mm)',
    y_lab='Upper whisker of average stand height in one hectare(m)',
    title='Stand Height limits vs. precipitation',
    name='Graphics/height_prcp')

min_temp_range_val = np.arange(round(np.min(new_dat.min_temp)),
    round(np.max(new_dat.min_temp) + 1))
min_t = {}
min_t_vals = {}
for i in min_temp_range_val[:-1]:
    tmp = new_dat[new_dat['min_temp'].between(i, i + 1, inclusive=True)]
    # print(tmp)

    min_t['dbh_{}'.format(i)] = tmp.tree_diameter
    min_t['basal_{}'.format(i)] = tmp.basal_area
    min_t['h_{}'.format(i)] = tmp.tree_height

    _t_dbh = calc_bp(tmp.tree_diameter)
    _t_diam = calc_bp(tmp.basal_area)
    _t_h = calc_bp(tmp.tree_height)

    if 'temp' in min_t_vals:
        min_t_vals['temp'].append(i)
        min_t_vals['upper_whisk_dbh'].append(_t_dbh['upper_whisker'])
        min_t_vals['upper_whisk_basal'].append(_t_diam['upper_whisker'])
        min_t_vals['upper_whisk_height'].append(_t_h['upper_whisker'])

    else:
        min_t_vals['temp'] = [i]
        min_t_vals['upper_whisk_dbh'] = [_t_dbh['upper_whisker']]
        min_t_vals['upper_whisk_basal'] = [_t_diam['upper_whisker']]
        min_t_vals['upper_whisk_height'] = [_t_h['upper_whisker']]
    print(i)

plot_lines(x=min_t_vals['temp'],
    y=min_t_vals['upper_whisk_dbh'],
    x_lab='Minimum temperature(°C)',
    y_lab='Upper whisker of basal-area weighted mean diameter(cm)',
    title='Basal-area weighted mean diameter limits in max temp',
    name='Graphics/dbh_min_tmp')

plot_lines(x=min_t_vals['temp'],
    y=min_t_vals['upper_whisk_basal'],
    x_lab='Minimum temperature(°C)',
    y_lab='Upper whisker of basal areas (m^2/ha)',
    title='Basal area limits in max temp',
    name='Graphics/basal_min_tmp')

plot_lines(x=min_t_vals['temp'],
    y=min_t_vals['upper_whisk_height'],
    x_lab='Minimum temperature(°C)',
    y_lab='Upper whisker of average stand height in one hectare(m)',
    title='Stand Height limits in max temp',
    name='Graphics/height_min_tmp')

max_temp_range_val = np.arange(round(np.min(new_dat.max_temp)),
    round(np.max(new_dat.max_temp) + 1))
max_t = {}
max_t_vals = {}
for i in max_temp_range_val[:-1]:
    tmp = new_dat[new_dat['max_temp'].between(i, i + 1, inclusive=True)]
    # print(tmp)

    max_t['dbh_{}'.format(i)] = tmp.tree_diameter
    max_t['basal_{}'.format(i)] = tmp.basal_area
    max_t['h_{}'.format(i)] = tmp.tree_height

    _t_dbh = calc_bp(tmp.tree_diameter)
    _t_diam = calc_bp(tmp.basal_area)
    _t_h = calc_bp(tmp.tree_height)

    if 'temp' in max_t_vals:
        max_t_vals['temp'].append(i)
        max_t_vals['upper_whisk_dbh'].append(_t_dbh['upper_whisker'])
        max_t_vals['upper_whisk_basal'].append(_t_diam['upper_whisker'])
        max_t_vals['upper_whisk_height'].append(_t_h['upper_whisker'])

    else:
        max_t_vals['temp'] = [i]
        max_t_vals['upper_whisk_dbh'] = [_t_dbh['upper_whisker']]
        max_t_vals['upper_whisk_basal'] = [_t_diam['upper_whisker']]
        max_t_vals['upper_whisk_height'] = [_t_h['upper_whisker']]
    print(i)

plot_lines(x=max_t_vals['temp'],
    y=max_t_vals['upper_whisk_dbh'],
    x_lab='Maximum temperature(°C)',
    y_lab='Upper whisker of basal-area weighted mean diameter(cm)',
    title='Basal-area weighted mean diameter limits in max temp',
    name='Graphics/dbh_max_tmp')

plot_lines(x=max_t_vals['temp'],
    y=max_t_vals['upper_whisk_basal'],
    x_lab='Maximum temperature(°C)',
    y_lab='Upper whisker of basal areas (m^2/ha)',
    title='Basal area limits in max temp',
    name='Graphics/basal_max_tmp')

plot_lines(x=max_t_vals['temp'],
    y=max_t_vals['upper_whisk_height'],
    x_lab='Maximum temperature(°C)',
    y_lab='Upper whisker of average stand height in one hectare(m/ha)',
    title='Stand Height limits in max temp',
    name='Graphics/height_max_tmp')


plt_scatter(var1='minimum temperature(°C)',
    var2='precipitation(mm)',
    var3='tree diameter(cm)',
    dat1=np.array(min_temp),
    dat2=np.array(precip),
    dat3=np.array(bhd),
    title='precipitation and minimum temperature vs tree diameter',
    path_to_fig='Graphics/p_mint_d')

plt_scatter(var1='max_temp',
    var2='precipitation',
    var3='tree_diameter(cm)',
    dat1=max_temp, dat2=precip, dat3=bhd,
    title='precipitation and maximum temperature vs tree diameter',
    path_to_fig='Graphics/p_maxt_d')

plt_scatter(var1='minimum temperature(°C)',
    var2='precipitation(mm)',
    var3='basal area(m^2/ha)',
    dat1=np.array(min_temp),
    dat2=np.array(precip),
    dat3=basal,
    title='Precipitation and minimum temperature vs tree basal area',
    path_to_fig='Graphics/p_mint_b')

plt_scatter(var1='maximum temperature(°C)',
    var2='precipitation(mm)',
    var3='basal area(m^2/ha)',
    dat1=max_temp, dat2=precip, dat3=basal,
    title='precipitation and maximum temperature vs basal area',
    path_to_fig='Graphics/p_maxt_b')

plt_scatter(var1='minimum temperature(°C)',
    var2='precipitation(mm)',
    var3='tree height(m/ha)',
    dat1=np.array(min_temp), dat2=np.array(precip), dat3=height,
    title='precipitation and minimum temperature vs tree height',
    path_to_fig='Graphics/p_mint_h')

plt_scatter(var1='maximum temperature(°C)',
    var2='precipitation(mm)',
    var3='tree height(m/ha)',
    dat1=max_temp, dat2=precip, dat3=height,
    title='precipitation and maximum temperature vs tree height',
    path_to_fig='Graphics/p_maxt_h')

