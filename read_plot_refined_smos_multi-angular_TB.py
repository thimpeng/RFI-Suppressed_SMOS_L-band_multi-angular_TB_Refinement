'''
@Description: 
@Author: Peng
@Date: 2023-04-11 16:27:00
LastEditTime: 2023-04-11 17:18:22
LastEditors: Peng
'''
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def read_refined_dgg_file(file):
    with nc.Dataset(file, 'r') as ds:
        inc = ds.variables['inc'][:]
        dgg_seconds = ds.variables['dgg_time'][:]
        lon = ds.variables['dgg_lon'][:]
        lat = ds.variables['dgg_lat'][:]
        TBh = ds.variables['TBh'][:]
        TBv = ds.variables['TBv'][:]
        TBh_flag = ds.variables['TBh_flag'][:]
        TBv_flag = ds.variables['TBv_flag'][:]
        TBh.mask = TBh_flag
        TBv.mask = TBv_flag
        # TBh[TBh_flag!=1] = np.nan
        # TBv[TBv_flag!=1] = np.nan
    return inc, dgg_seconds, lon, lat, TBh, TBv

def read_refined_easegrid_file(file):
    with nc.Dataset(file, 'r') as ds:
        inc = ds.variables['inc'][:]
        utc_seconds = ds.variables['utc_seconds'][:]
        lon = ds.variables['longitude'][:]
        lat = ds.variables['latitude'][:]
        TBh = ds.variables['TBh'][:]
        TBv = ds.variables['TBv'][:]
    return inc, utc_seconds, lon, lat, TBh, TBv

def plot_data(lon, lat, outangles, TBh, TBv, plot_angle=40, fontsize=14):
    x, y = m(lon, lat)
    fig, ax = plt.subplots(2,1,figsize=(16,10), dpi=300)
    plt.subplot(2,1,1)
    m.scatter(x,y,c=TBh[:,outangles==plot_angle], s=0.1, cmap='jet', vmin=100, vmax=320)
    m.drawparallels(np.arange(-45, 90, 45), labels=[True,False,False,False], linewidth=0.1, fontsize=fontsize)
    m.drawmeridians(np.arange(-180, 180, 45), labels=[False,False,False,True], linewidth=0.1, fontsize=fontsize)
    plt.colorbar(pad=0.02)
    plt.title('Horizontal polarized TB at {} degree'.format(plot_angle))
    plt.subplot(2,1,2)
    m.scatter(x,y,c=TBv[:,outangles==plot_angle], s=0.1, cmap='jet', vmin=100, vmax=320)
    m.drawparallels(np.arange(-45, 90, 45), labels=[True,False,False,False], linewidth=0.1, fontsize=fontsize)
    m.drawmeridians(np.arange(-180, 180, 45), labels=[False,False,False,True], linewidth=0.1, fontsize=fontsize)
    plt.colorbar(pad=0.02)
    plt.title('Vertical polarized TB at {} degree'.format(plot_angle))
    return fig

if __name__ == '__main__':
    m = Basemap(epsg=6933, resolution='h', llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180)
    
    # read refined data in ISEA4H9 15-KM DGG
    file = r'N:\ISEA4H9_15KM\2016\20160601\SM_REPR_MIR_SCLF1C_A_20160601T002552_20160601T011912_724_200_1.nc'
    ds = nc.Dataset(file, 'r')
    print(ds.variables)
    ds.close()
    outangles, utc_seconds, lon, lat, TBh, TBv = read_refined_dgg_file(file)
    
    # read refined data in EASE-GRID 2.0 Global 25-KM
    file = r'N:\EASEGRID2_25KM\2016\20160601\SM_REPR_MIR_SCLF1C_A_20160601T002552_20160601T011912_724_200_1_Hamming_25KM.nc'
    ds = nc.Dataset(file, 'r')
    print(ds.variables)
    ds.close()
    outangles, utc_seconds, lon, lat, TBh, TBv = read_refined_easegrid_file(file)
    # convert seconds to time
    smos_j2000 = np.datetime64('2000-01-01T00')
    utc_time = smos_j2000+np.timedelta64(1,'s')*utc_seconds
    print(outangles) #[2.5,  7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 40. , 42.5, 47.5, 52.5, 57.5, 62.5, 67.5]
    # plt data at 40 degree
    fig = plot_data(lon, lat, outangles, TBh, TBv, 40)
    plt.savefig(file.replace('.nc','_TB_40.tif'), dpi=300, bbox_inches='tight', pil_kwargs = {"compression": "tiff_lzw"})