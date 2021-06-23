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
import rasterio
import rioxarray as rxr
import rasterio as rio
import gc
import pandas as pd
import numpy as np
import os
from glob import glob

sys.path.append('Data_prep')

import project_extent as pe

boundary = gpd.read_file('AB2551Watersheds/AB2551Watersheds.shp')

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

    meta.update(dtype=rasterio.float32,crs='EPSG:3310')

    # Write output file
    with rasterio.open(out_name, 'w', **meta) as dst:
        dst.write(array_out.astype(rasterio.float32), 1)

    return array_out


def annual_avg(path_list,years,month_list_1,month_list_2,out_name):
    '''

    Parameters
    ----------
    path_list
    years
    month_list_1
    month_list_2
    out_name

    Returns
    -------

    '''

    import rasterio

    lyrs=[]

    for year in years:
        globals()[f'list_{year}'] = [x for x in path_list if
                                    (x[-10:-6] == '%d' % (year-1) and
                                    x[-6:-4] in month_list_1) or
                                    (x[-10:-6] == '%d' % year and
                                    x[-6:-4] in month_list_2)]
        lyr_list = [read_raster_array(x) for x in globals()['list_{year}']]
        globals()[f'array_out_{year}'] = np.sum(lyr_list, axis=0)
        lyrs.append(globals()[f'array_out_{year}'])
        print(year)

    average_precip = np.mean(lyrs, axis=0)
    with rasterio.open(path_list[0]) as src:
        meta = src.meta

    meta.update(dtype=rasterio.float32,crs='EPSG:3310')

    # Write output file
    with rasterio.open(out_name, 'w', **meta) as dst:
        dst.write(average_precip.astype(rasterio.float32), 1)

    return average_precip


def visualize_raster(raster_lyr,
                     bound,
                     fig_name,
                     path_to_fig,
                     axes=False,
                     min_lim=float("NaN"),
                     max_lim=float("NaN"),
                     cmap=None,
                     dpi=300):
    """
    Parameters
    ----------
    raster_lyr
    bound
    fig_name
    path_to_fig
    axes
    min_lim
    max_lim
    cmap
    dpi

    Returns
    -------

    """

    # array_out_rio = rxr.open_rasterio('Climate_limits/resampled_py_height.tif')
    array_out_rio = rxr.open_rasterio(raster_lyr)

    if ~np.isnan(min_lim):
        array_out_rio.data[array_out_rio.data < min_lim] = min_lim
    if ~np.isnan(max_lim):
        array_out_rio.data[array_out_rio.data < max_lim] = max_lim

    f, ax = plt.subplots(figsize=(10, 4))

    if cmap is not None:
        array_out_rio.plot(ax=ax, cmap=cmap)
    else:
        array_out_rio.plot(ax=ax)
    bound.boundary.plot(ax=ax, color='black')

    ax.set(title=fig_name)
    if not axes:
        ax.set_axis_off()
    plt.show()

    plt.savefig(path_to_fig, dpi=dpi)
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
               name,
               dpi=300):
    """
    Parameters
    ----------
    x
    y
    x_lab
    y_lab
    title
    name
    dpi
    Returns
    -------

    """
    plt.plot(x, y)
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)

    plt.grid()
    # displaying the title
    plt.title(title)
    plt.savefig(name,dpi=dpi)
    plt.show()
    plt.close()


def plt_scatter(var1, var2, var3, dat1, dat2, dat3,
                title, path_to_fig, dpi=300):
    """

    Parameters
    ----------
    var1
    var2
    var3
    dat1
    dat2
    dat3
    title
    path_to_fig
    dpi

    Returns
    -------

    """
    one = pd.DataFrame({var1: dat1.flatten(),
                        var2: dat2.flatten(),
                        var3: dat3.flatten()})
    one.dropna()
    one.round(decimals=1)
    one.drop_duplicates()
    t = one.groupby([var1, var2]).mean().reset_index()
    fig, ax = plt.subplots()
    plt.grid(linestyle=':', linewidth='0.5', color='black')
    ax.grid(which='minor', linestyle=':', linewidth='0.5')
    plt.scatter(x=t[var1], y=t[var2], c=t[var3],
        marker='o', s=0.5, edgecolors='none')
    plt.xlabel(var1)
    plt.ylabel(var2)
    plt.title(title)
    ax.set_xlabel(var1)
    ax.set_ylabel(var2)
    cbar = plt.colorbar()
    cbar.set_label(var3)
    plt.show()
    plt.savefig(path_to_fig, dpi=dpi)
    plt.close()

    print(title + ' plotted to '+ path_to_fig)



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

    file = read_raster_array(file_name)
    array_out = file / scaler
    import rasterio

    with rasterio.open(file_name) as src:
        meta = src.meta

    meta.update(dtype=rasterio.float32, crs='EPSG:3310')

    # Write output file
    with rasterio.open(out_name, 'w', **meta) as dst:
        dst.write(array_out.astype(rasterio.float32), 1)

    print('file %s is rescales to 1/%d and saved to %s'%(file_name ,scaler
    ,out_name))

# call tmin avg
tmin = glob(os.path.join(
    '../../Documents/ForesiteProject/Data/Climate_data/',
    'prism_tmin*.tif'))

tmin_list = tmin
tmin_avg = avg_calc(tmin_list, 'Climate_limits/prism_tmin.tif')

pe.resample(in_file='Climate_limits/prism_tmin.tif',
    dst_file='Climate_limits/resampled_py_prism_tmin.tif',
    src_file='Extent/buffered_extent.tif',
    snap=True)

visualize_raster(raster_lyr='Climate_limits/resampled_py_prism_tmin.tif',
    bound=boundary,
    fig_name='Average Minimum Temperature(°C)',
    path_to_fig='Graphics/resampled_py_prism_tmin',
    cmap='coolwarm')

# call tmax avg
tmax = glob(os.path.join(
    '../../Documents/ForesiteProject/Data/Climate_data/',
    'prism_tmax*.tif'))
tmax_list = tmax
tmax_avg = avg_calc(tmax_list, 'Climate_limits/prism_tmax.tif')

pe.resample(in_file='Climate_limits/prism_tmax.tif',
    dst_file='Climate_limits/resampled_py_prism_tmax.tif',
    src_file='Extent/buffered_extent.tif',
    snap=True)

visualize_raster(raster_lyr='Climate_limits/prism_tmax.tif',
    bound=boundary,
    fig_name='Average Maximum Temperature(°C)',
    path_to_fig='Graphics/resampled_py_prism_tmax',
    cmap='YlOrRd')

# call precip avg
percip_list = glob(os.path.join(
    '../../Documents/ForesiteProject/Data/Climate_data/',
    'prism_ppt*.tif'))

month_list_1 = ['10', '11', '12']
month_list_2 = ['01', '02','03', '04', '05', '06', '07', '08', '09']
years = np.arange(1982, 2020, 1)

precip = annual_avg(percip_list,
                    years,
                    month_list_1,
                    month_list_1,
                    'Climate_limits/prism_ppt_annual.tif')

pe.resample(in_file='Climate_limits/prism_ppt_annual.tif',
    dst_file='Climate_limits/resampled_py_prism_ppt_annual.tif',
    src_file='Extent/buffered_extent.tif',
    snap=True)

visualize_raster(raster_lyr='Climate_limits/resampled_py_prism_ppt_annual.tif',
    bound=boundary,
    fig_name='Average Precipitation(mm)',
    path_to_fig='Graphics/resampled_py_prism_ppt_wy',
    cmap='YlGnBu')


"""
  Attribute        Year  Filename                  Scalar  Number of neighbors
  BA_GE_3          2017  ba_ge_3_2017.tif          100.0   7    
  MNDBHBA          2017  mndbhba_2017.tif          10.0    7    
  STNDHGT          2017  stndhgt_2017.tif          100.0   7  
  
stndhgt_2017.tif: "average" Stand height, computed as average of heights of 
all dominant and codominant trees (m)
--Sum tree height multiplied by TPH for dominant and codominant trees 
(identified by the CROWN_CLASS field)
HT_SUM = ∑ [TREE_LIVE.HT_M * TPH_FC] (for DBH_CM >= 2.5 and CROWN_CLASS <= 3)
--Sum TPH for dominant and codominant trees
TPH_SUM = ∑ TPH_FC (for DBH_CM >= 2.5 and CROWN_CLASS <= 3)
--Divide sum of tree heights by TPH sum 
STNDHGT = HT_SUM / TPH_SUM 

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


rescale_LEMMA(file_name='Climate_limits/resampled_py_height.tif',
              out_name='Climate_limits/tree_height.tif',
              scaler=100)

visualize_raster(raster_lyr='Climate_limits/tree_height.tif',
    bound=boundary,
    fig_name='Stand Height(m)',
    path_to_fig='Graphics/tree_height',
    min_lim=0,
    cmap='Greens')


pe.resample(in_file=os.path.join('GNN_AB2551',
    Forest_structure_lst[1]),
    dst_file='Climate_limits/resampled_py_basal.tif',
    src_file='Extent/buffered_extent.tif',
    snap=True)
rescale_LEMMA(file_name='Climate_limits/resampled_py_basal.tif',
              out_name='Climate_limits/basal_area.tif',
              scaler=100)
visualize_raster(raster_lyr='Climate_limits/basal_area.tif',
    bound=boundary,
    fig_name='Tree Basal Area(m^2/ha)',
    path_to_fig='Graphics/basal_area',
    min_lim=0,
    cmap='YlGn')

pe.resample(in_file=os.path.join('GNN_AB2551',
    Forest_structure_lst[2]),
    dst_file='Climate_limits/resampled_py_diameter.tif',
    src_file='Extent/buffered_extent.tif',
    snap=True)
rescale_LEMMA(file_name='Climate_limits/resampled_py_diameter.tif',
              out_name='Climate_limits/tree_diameter.tif',
              scaler=10)
visualize_raster(raster_lyr='Climate_limits/tree_diameter.tif',
    bound=boundary,
    fig_name='Basal-area weighted mean diameter(cm)',
    path_to_fig='Graphics/tree_diameter',
    min_lim=0,
    cmap='BuGn')


max_temp = read_raster_array('Climate_limits/resampled_py_prism_tmax.tif')
min_temp = read_raster_array('Climate_limits/resampled_py_prism_tmin.tif')
precip = read_raster_array('Climate_limits/resampled_py_prism_ppt_annual.tif')

bhd = read_raster_array('Climate_limits/tree_diameter.tif')
basal = read_raster_array('Climate_limits/basal_area.tif')
height = read_raster_array('Climate_limits/tree_height.tif')

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

#t = Averaging moving window
t = 20

# p_range_val = np.array(np.arange(round(np.min(new_dat.precipitation)),
#     round(np.max(new_dat.precipitation))))

p_range_val = np.array(np.arange(round(np.min(new_dat.precipitation)),
    round(np.max(new_dat.precipitation)+1), 10))

precip = {}
per_vals = {}

_t_dbh = {'median': 0.0, 'upper_quartile': 0.0, 'lower_quartile': 0.0,
          'iqr': 0.0, 'upper_whisker': 0.0, 'lower_whisker': 0.0}
_t_diam = {'median': 0.0, 'upper_quartile': 0.0, 'lower_quartile': 0.0,
           'iqr': 0.0, 'upper_whisker': 0.0, 'lower_whisker': 0.0}
_t_h = {'median': 0.0, 'upper_quartile': 0.0, 'lower_quartile': 0.0,
        'iqr': 0.0, 'upper_whisker': 0.0, 'lower_whisker': 0.0}

for i in p_range_val:
    tmp = new_dat[new_dat['precipitation'].between(i - t/2, i + t/2,
        inclusive=True)]
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
    # print(_t_dbh['upper_whisker'])
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
    title='Basal-area weighted mean diameter limits in min temp',
    name='Graphics/dbh_min_tmp')

plot_lines(x=min_t_vals['temp'],
    y=min_t_vals['upper_whisk_basal'],
    x_lab='Minimum temperature(°C)',
    y_lab='Upper whisker of basal areas (m^2/ha)',
    title='Basal area limits in min temp',
    name='Graphics/basal_min_tmp')

plot_lines(x=min_t_vals['temp'],
    y=min_t_vals['upper_whisk_height'],
    x_lab='Minimum temperature(°C)',
    y_lab='Upper whisker of average stand height in one hectare(m)',
    title='Stand Height limits in min temp',
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


max_temp = read_raster_array('Climate_limits/resampled_py_prism_tmax.tif')
min_temp = read_raster_array('Climate_limits/resampled_py_prism_tmin.tif')
precip = read_raster_array('Climate_limits/resampled_py_prism_ppt_annual.tif')

bhd = read_raster_array('Climate_limits/tree_diameter.tif')
basal = read_raster_array('Climate_limits/basal_area.tif')
height = read_raster_array('Climate_limits/tree_height.tif')

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

