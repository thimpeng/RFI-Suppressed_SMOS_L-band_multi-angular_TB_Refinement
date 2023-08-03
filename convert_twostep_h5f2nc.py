'''
Description: 
Author: Peng
Date: 2022-07-22 10:52:33
LastEditTime: 2022-08-24 19:37:07
LastEditors: Peng
'''
import numpy as np
import h5py
import sys
from datetime import datetime, timedelta, timezone
import os
import netCDF4 as nc
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--numparts', default=30, type=int)
parser.add_argument('--part', default=1, type=int)
args = parser.parse_args()

def iterbrowse(path):
    for home, dirs, files in os.walk(path):
        for filename in files:
            yield os.path.join(home, filename)

def FilterFileList(path, keywords, *keywords_not_include):

    fileList = []

    for fullname in iterbrowse(path):
        if keywords_not_include:
            if all([keyword in fullname for keyword in keywords]) & all([keyword_not_include not in fullname for keyword_not_include in keywords_not_include]):
                fileList.append(fullname)
        else:
            if all([keyword in fullname for keyword in keywords]):
                fileList.append(fullname)
    return fileList

def traverse_datasets(hdf_file):

    """Traverse all datasets across all groups in HDF5 file."""

    def h5py_dataset_iterator(g, prefix=''):
        for key in g.keys():
            item = g[key]
            path = '{}/{}'.format(prefix, key)
            if isinstance(item, h5py.Dataset): # test for dataset
                yield (path, item)
            elif isinstance(item, h5py.Group): # test for group (go down)
                yield from h5py_dataset_iterator(item, path)

    with h5py.File(hdf_file, 'r') as f:
        for (path, dset) in h5py_dataset_iterator(f):
            # yield path.split('/')[-1]
            yield path

def compute_scale_and_offset(da, n=16):
    """Calculate offset and scale factor for int conversion

    Based on Krios101's code above.
    """

    vmin = np.nanmin(da).item()
    vmax = np.nanmax(da).item()

    # stretch/compress data to the available packed range
    scale_factor = (vmax - vmin) / (2 ** n - 1)

    # translate the range to be symmetric about zero
    add_offset = vmin + 2 ** (n - 1) * scale_factor

    return scale_factor, add_offset

def convert_h5f2nc(h5file):

    kwargs = {'zlib': True}
    file_time = os.path.basename(h5file).split('_')[5:7]
    scen = os.path.basename(h5file).split('_')[4]
    raw_filename = FilterFileList(os.path.join(r'/data2/public2/Satellite/Observation/SMOS/SCLF1C/v724', file_time[0][:4], file_time[0][4:6], file_time[0][6:8]), file_time, *['.HDR','.hdr'])
    assert len(raw_filename) == 1
    new_filename = os.path.basename(raw_filename).replace('_SCLF1C_','_SCLF1C_{}_'.format(scen)).replace('.DBL', '.nc')
    ncfile = os.path.path.join(os.path.dirname(h5file), new_filename)
    
    with h5py.File(h5file, 'r') as f:
        with nc.Dataset(ncfile, 'w', format='NETCDF4') as ds:

            out_angles = np.insert(np.arange(2.5, 68, 5), 8, 40)
            n_inc = out_angles.size

            # create dimension
            ds.createDimension('inc', n_inc)
            ds.createDimension('n_grid_points', f['dgg_id'].shape[0])

            # create variables
            incidence_angles = ds.createVariable('inc', 'float32', ('inc',))
            incidence_angles[:] = out_angles

            dgg_lon = ds.createVariable('dgg_lon', 'float32', ('n_grid_points',), **kwargs)
            dgg_lon.units = 'degrees_east'
            dgg_lon.long_name = "Longitude of the DGG cell's center identified by dgg_id."
            dgg_lon[:] = f['dgg_lon'][:]

            dgg_lat = ds.createVariable('dgg_lat', 'float32', ('n_grid_points',), **kwargs)
            dgg_lat.units = 'degrees_north'
            dgg_lat.long_name = "Latitude of the DGG cell's center identified by dgg_id."
            dgg_lat[:] = f['dgg_lat'][:]

            dgg_id = ds.createVariable('dgg_id', 'int32', ('n_grid_points',), **kwargs)
            dgg_id.units = 'dimensionless'
            dgg_id.coordinates = 'dgg_lat dgg_lon'
            dgg_id.long_name = "Unique identifier for Earth fixed grid point of ISEA 4H9 equal-area discrete global gird (DGG) system."
            dgg_id[:] = f['dgg_id'][:]

            dgg_mask = ds.createVariable('dgg_mask', 'uint8', ('n_grid_points',), **kwargs)
            dgg_mask.units = 'dimensionless'
            dgg_mask.coordinates = 'dgg_lat dgg_lon'
            dgg_mask.long_name = "Flag indicating land/sea USGS content, coastline distance, and Ice content, from SCLF1C"
            dgg_mask[:] = f['dgg_mask'][:]

            smos_start_date = np.datetime64('2000-01-01T00')
            dgg_time = smos_start_date+np.timedelta64(1,'s')*f['dgg_time'][:]
            days = (dgg_time.astype('datetime64[D]')-smos_start_date.astype('datetime64[D]')).astype(int)
            seconds = (dgg_time - dgg_time.astype('datetime64[D]')).astype(int)

            dgg_time = ds.createVariable('dgg_time', 'uint32', ('n_grid_points',), **kwargs)
            dgg_time.units = 'seconds since 2000-01-01 00:00:00'
            dgg_time.long_name = 'Median of the acquisition time of all the scens were taken.'
            dgg_time.coordinates = 'dgg_lat dgg_lon'
            # dgg_time[:] = nc.num2date(f['dgg_time'][:], dgg_time.units)
            dgg_time[:] = f['dgg_time'][:]

            Tbh = ds.createVariable('TBh', 'float32', ('n_grid_points','inc'), **kwargs)
            Tbh.coordinates = 'dgg_lat dgg_lon'
            Tbh.units = 'Kelvin'
            Tbh.long_name = 'TB of horizontal polarization at incidence angle (inc) identified by dgg_id.'
            Tbh[:,:-1] = f['TBh'][:]

            Tbv = ds.createVariable('TBv', 'float32', ('n_grid_points','inc'), **kwargs)
            Tbv.coordinates = 'dgg_lat dgg_lon'
            Tbv.units = 'Kelvin'
            Tbv.long_name = 'TB at vertical polarization at incidence angle (inc) identified by dgg_id.'
            Tbv[:,:-1] = f['TBv'][:]

            dof = ds.createVariable('dof', 'short', ('n_grid_points',), **kwargs)
            dof.coordinates = 'dgg_lat dgg_lon'
            dof.units = 'dimensionless'
            dof.long_name = 'Degree of freedom of horizontal polarization in the second step regression, identified by dgg_id.'
            dof[:] = f['dof'][:]

            flag = ds.createVariable('flag', 'uint32', ('n_grid_points',), **kwargs)
            flag.coordinates = 'dgg_lat dgg_lon'
            flag.units = 'dimensionless'
            flag.long_name = 'Bits flags that represent quality of the regression result, identified by dgg_id.'
            flag[:] = f['flag'][:]

            A_1st = ds.createVariable('A_first_step', 'float64', ('n_grid_points',), **kwargs)
            A_1st.coordinates = 'dgg_lat dgg_lon'
            A_1st.units = 'dimensionless'
            A_1st.long_name = 'Parameter A of the first step results, identified by dgg_id.'
            A_1st[:] = f['params_1st_step'][:,0]

            C = ds.createVariable('C', 'float64', ('n_grid_points',), **kwargs)
            C.coordinates = 'dgg_lat dgg_lon'
            C.units = 'Kelvin'
            C.long_name = 'Parameter C (the first Stokes parameter) of the first step results, identified by dgg_id.'
            C[:] = f['params_1st_step'][:,1]

            assert np.allclose(f['params_1st_step'][:,1], f['params_h'][:,3]*2, equal_nan=True)
            assert np.allclose(f['params_1st_step'][:,1], f['params_v'][:,3]*2, equal_nan=True)

            redChi2_1st = ds.createVariable('redChi2_first_step', 'float64', ('n_grid_points',))
            redChi2_1st.coordinates = 'dgg_lat dgg_lon'
            redChi2_1st.units = 'dimensionless'
            redChi2_1st.long_name = 'The Reduced Chi-squared of the first step regression, identified by dgg_id.'
            redChi2_1st[:] = f['params_1st_step'][:,2]

            aic_1st = ds.createVariable('aic_first_step', 'float64', ('n_grid_points',), **kwargs)
            aic_1st.coordinates = 'dgg_lat dgg_lon'
            aic_1st.units = 'Kelvin'
            aic_1st.long_name = 'Akaike information criterion of the first step regression, identified by dgg_id.'
            aic_1st[:] = f['params_1st_step'][:,3]

            bic_1st = ds.createVariable('bic_first_step', 'float64', ('n_grid_points',), **kwargs)
            bic_1st.coordinates = 'dgg_lat dgg_lon'
            bic_1st.units = 'Kelvin'
            bic_1st.long_name = 'Bayesian information criterion of the first step regression, identified by dgg_id.'
            bic_1st[:] = f['params_1st_step'][:,4]

            ah_2nd = ds.createVariable('a_H_second_step', 'float64', ('n_grid_points',), **kwargs)
            ah_2nd.coordinates = 'dgg_lat dgg_lon'
            ah_2nd.units = 'dimensionless'
            ah_2nd.long_name = 'Parameter a of the second step results at horizontal polarization, identified by dgg_id.'
            ah_2nd[:] = f['params_h'][:,0]

            bh_2nd = ds.createVariable('b_H_second_step', 'float64', ('n_grid_points',), **kwargs)
            bh_2nd.coordinates = 'dgg_lat dgg_lon'
            bh_2nd.units = 'dimensionless'
            bh_2nd.long_name = 'Parameter b of the second step results at horizontal polarization, identified by dgg_id.'
            bh_2nd[:] = f['params_h'][:,1]

            redChi2h_2nd = ds.createVariable('redChi2_H_second_step', 'float64', ('n_grid_points',), **kwargs)
            redChi2h_2nd.coordinates = 'dgg_lat dgg_lon'
            redChi2h_2nd.units = 'dimensionless'
            redChi2h_2nd.long_name = 'The Reduced Chi-squared of the second step regression at horizontal polarization, identified by dgg_id.'
            redChi2h_2nd[:] = f['params_h'][:,4]

            aich_2nd = ds.createVariable('aic_H_second_step', 'float64', ('n_grid_points',), **kwargs)
            aich_2nd.coordinates = 'dgg_lat dgg_lon'
            aich_2nd.units = 'dimensionless'
            aich_2nd.long_name = 'Akaike information criterion of the second step regression at horizontal polarization, identified by dgg_id.'
            aich_2nd[:] = f['params_h'][:,5]

            bich_2nd = ds.createVariable('bic_H_second_step', 'float64', ('n_grid_points',), **kwargs)
            bich_2nd.coordinates = 'dgg_lat dgg_lon'
            bich_2nd.units = 'dimensionless'
            bich_2nd.long_name = 'Bayesian information criterion of the second step regression at horizontal polarization, identified by dgg_id.'
            bich_2nd[:] = f['params_h'][:,6]

            av_2nd = ds.createVariable('a_V_second_step', 'float64', ('n_grid_points',), **kwargs)
            av_2nd.coordinates = 'dgg_lat dgg_lon'
            av_2nd.units = 'dimensionless'
            av_2nd.long_name = 'Parameter a of the second step results at vertical polarization, identified by dgg_id.'
            av_2nd[:] = f['params_v'][:,0]

            bv_2nd = ds.createVariable('b_V_second_step', 'float64', ('n_grid_points',), **kwargs)
            bv_2nd.coordinates = 'dgg_lat dgg_lon'
            bv_2nd.units = 'dimensionless'
            bv_2nd.long_name = 'Parameter b of the second step results at vertical polarization, identified by dgg_id.'
            bv_2nd[:] = f['params_v'][:,1]

            dv_2nd = ds.createVariable('d_V_second_step', 'float64', ('n_grid_points',), **kwargs)
            dv_2nd.coordinates = 'dgg_lat dgg_lon'
            dv_2nd.units = 'dimensionless'
            dv_2nd.long_name = 'Parameter d of the second step results at vertical polarization, identified by dgg_id.'
            dv_2nd[:] = f['params_v'][:,2]

            redChi2v_2nd = ds.createVariable('redChi2_V_second_step', 'float64', ('n_grid_points',), **kwargs)
            redChi2v_2nd.coordinates = 'dgg_lat dgg_lon'
            redChi2v_2nd.units = 'dimensionless'
            redChi2v_2nd.long_name = 'The Reduced Chi-squared of the second step regression at vertical polarization, identified by dgg_id.'
            redChi2v_2nd[:] = f['params_v'][:,4]

            aicv_2nd = ds.createVariable('aic_V_second_step', 'float64', ('n_grid_points',), **kwargs)
            aicv_2nd.coordinates = 'dgg_lat dgg_lon'
            aicv_2nd.units = 'dimensionless'
            aicv_2nd.long_name = 'Akaike information criterion of the second step regression at vertical polarization, identified by dgg_id.'
            aicv_2nd[:] = f['params_v'][:,5]

            bicv_2nd = ds.createVariable('bic_V_second_step', 'float64', ('n_grid_points',), **kwargs)
            bicv_2nd.coordinates = 'dgg_lat dgg_lon'
            bicv_2nd.units = 'dimensionless'
            bicv_2nd.long_name = 'Bayesian information criterion of the second step regression at vertical polarization, identified by dgg_id.'
            bicv_2nd[:] = f['params_v'][:,6]

            MeanBias_H_Fitted = ds.createVariable('MeanBias_H_Fitted', 'float32', ('n_grid_points',), **kwargs)
            MeanBias_H_Fitted.coordinates = 'dgg_lat dgg_lon'
            MeanBias_H_Fitted.units = 'Kelvin'
            MeanBias_H_Fitted.long_name = 'Averaged bias (horizontal polarization) between filtered brightness temperature used to fit and the corresponding fitted ones at same incidence angles.'
            MeanBias_H_Fitted[:] = f['MeanBias_H_Fitted'][:]

            MeanBias_H_Fitted_Bins = ds.createVariable('MeanBias_H_Fitted_Bins', 'float32', ('n_grid_points','inc'), **kwargs)
            MeanBias_H_Fitted_Bins.coordinates = 'dgg_lat dgg_lon'
            MeanBias_H_Fitted_Bins.units = 'Kelvin'
            MeanBias_H_Fitted_Bins.long_name = 'Angle class averaged bias (horizontal polarization, bin center ±2.5 degree) between filtered brightness temperature used to fit and the corresponding fitted ones at same incidence angles.'
            MeanBias_H_Fitted_Bins[:] = f['MeanBias_H_Fitted_Bins'][:]

            MeanBias_H_L1C_Fit = ds.createVariable('MeanBias_H_L1C_Fit', 'float32', ('n_grid_points',), **kwargs)
            MeanBias_H_L1C_Fit.coordinates = 'dgg_lat dgg_lon'
            MeanBias_H_L1C_Fit.units = 'Kelvin'
            MeanBias_H_L1C_Fit.long_name = 'Averaged bias (horizontal polarization) between SCLF1C brightness temperature (directly interpolated and convert to ground reference) and corresponding fitted ones at same incidence angles.'
            MeanBias_H_L1C_Fit[:] = f['MeanBias_H_L1C_Fit'][:]

            MeanBias_H_L1C_Fit_Bins = ds.createVariable('MeanBias_H_L1C_Fit_Bins', 'float32', ('n_grid_points','inc'), **kwargs)
            MeanBias_H_L1C_Fit_Bins.coordinates = 'dgg_lat dgg_lon'
            MeanBias_H_L1C_Fit_Bins.units = 'Kelvin'
            MeanBias_H_L1C_Fit_Bins.long_name = 'Angle class averaged bias (horizontal polarization, bin center ±2.5 degree) between SCLF1C brightness temperature (directly interpolated and convert to ground reference) and corresponding fitted ones at same incidence angles.'
            MeanBias_H_L1C_Fit_Bins[:] = f['MeanBias_H_L1C_Fit_Bins'][:]

            MeanBias_V_Fitted = ds.createVariable('MeanBias_V_Fitted', 'float32', ('n_grid_points',), **kwargs)
            MeanBias_V_Fitted.coordinates = 'dgg_lat dgg_lon'
            MeanBias_V_Fitted.units = 'Kelvin'
            MeanBias_V_Fitted.long_name = 'Averaged bias (vertical polarization) between filtered brightness temperature used to fit and the corresponding fitted ones at same incidence angles.'
            MeanBias_V_Fitted[:] = f['MeanBias_V_Fitted'][:]

            MeanBias_V_Fitted_Bins = ds.createVariable('MeanBias_V_Fitted_Bins', 'float32', ('n_grid_points','inc'), **kwargs)
            MeanBias_V_Fitted_Bins.coordinates = 'dgg_lat dgg_lon'
            MeanBias_V_Fitted_Bins.units = 'Kelvin'
            MeanBias_V_Fitted_Bins.long_name = 'Angle class averaged bias (vertical polarization, bin center ±2.5 degree) between filtered brightness temperature used to fit and the corresponding fitted ones at same incidence angles.'
            MeanBias_V_Fitted_Bins[:] = f['MeanBias_V_Fitted_Bins'][:]

            MeanBias_V_L1C_Fit = ds.createVariable('MeanBias_V_L1C_Fit', 'float32', ('n_grid_points',), **kwargs)
            MeanBias_V_L1C_Fit.coordinates = 'dgg_lat dgg_lon'
            MeanBias_V_L1C_Fit.units = 'Kelvin'
            MeanBias_V_L1C_Fit.long_name = 'Averaged bias (vertical polarization) between SCLF1C brightness temperature (directly interpolated and convert to ground reference) and corresponding fitted ones at same incidence angles.'
            MeanBias_V_L1C_Fit[:] = f['MeanBias_V_L1C_Fit'][:]

            MeanBias_V_L1C_Fit_Bins = ds.createVariable('MeanBias_V_L1C_Fit_Bins', 'float32', ('n_grid_points','inc'), **kwargs)
            MeanBias_V_L1C_Fit_Bins.coordinates = 'dgg_lat dgg_lon'
            MeanBias_V_L1C_Fit_Bins.units = 'Kelvin'
            MeanBias_V_L1C_Fit_Bins.long_name = 'Angle class averaged bias (vertical polarization, bin center ±2.5 degree) between SCLF1C brightness temperature (directly interpolated and convert to ground reference) and corresponding fitted ones at same incidence angles.'
            MeanBias_V_L1C_Fit_Bins[:] = f['MeanBias_V_L1C_Fit_Bins'][:]

            RMSD_H_Fitted = ds.createVariable('RMSD_H_Fitted', 'float32', ('n_grid_points',), **kwargs)
            RMSD_H_Fitted.coordinates = 'dgg_lat dgg_lon'
            RMSD_H_Fitted.units = 'Kelvin'
            RMSD_H_Fitted.long_name = 'Root mean square deviation (horizontal polarization) between filtered brightness temperature used to fit and the corresponding fitted ones at same incidence angles.'
            RMSD_H_Fitted[:] = f['RMSD_H_Fitted'][:]

            RMSD_H_Fitted_Bins = ds.createVariable('RMSD_H_Fitted_Bins', 'float32', ('n_grid_points','inc'), **kwargs)
            RMSD_H_Fitted_Bins.coordinates = 'dgg_lat dgg_lon'
            RMSD_H_Fitted_Bins.units = 'Kelvin'
            RMSD_H_Fitted_Bins.long_name = 'Angle class root mean square deviation (horizontal polarization, bin center ±2.5 degree) between filtered brightness temperature used to fit and the corresponding fitted ones at same incidence angles.'
            RMSD_H_Fitted_Bins[:] = f['RMSD_H_Fitted_Bins'][:]

            RMSD_H_L1C_Fit = ds.createVariable('RMSD_H_L1C_Fit', 'float32', ('n_grid_points',), **kwargs)
            RMSD_H_L1C_Fit.coordinates = 'dgg_lat dgg_lon'
            RMSD_H_L1C_Fit.units = 'Kelvin'
            RMSD_H_L1C_Fit.long_name = 'Root mean square deviation (horizontal polarization) between SCLF1C brightness temperature (directly interpolated and convert to ground reference) and corresponding fitted ones at same incidence angles.'
            RMSD_H_L1C_Fit[:] = f['RMSD_H_L1C_Fit'][:]

            RMSD_H_L1C_Fit_Bins = ds.createVariable('RMSD_H_L1C_Fit_Bins', 'float32', ('n_grid_points','inc'), **kwargs)
            RMSD_H_L1C_Fit_Bins.coordinates = 'dgg_lat dgg_lon'
            RMSD_H_L1C_Fit_Bins.units = 'Kelvin'
            RMSD_H_L1C_Fit_Bins.long_name = 'Angle class root mean square deviation (horizontal polarization, bin center ±2.5 degree) between SCLF1C brightness temperature (directly interpolated and convert to ground reference) and corresponding fitted ones at same incidence angles.'
            RMSD_H_L1C_Fit_Bins[:] = f['RMSD_H_L1C_Fit_Bins'][:]

            RMSD_V_Fitted = ds.createVariable('RMSD_V_Fitted', 'float32', ('n_grid_points',), **kwargs)
            RMSD_V_Fitted.coordinates = 'dgg_lat dgg_lon'
            RMSD_V_Fitted.units = 'Kelvin'
            RMSD_V_Fitted.long_name = 'Root mean square deviation (vertical polarization) between filtered brightness temperature used to fit and the corresponding fitted ones at same incidence angles.'
            RMSD_V_Fitted[:] = f['RMSD_V_Fitted'][:]

            RMSD_V_Fitted_Bins = ds.createVariable('RMSD_V_Fitted_Bins', 'float32', ('n_grid_points','inc'), **kwargs)
            RMSD_V_Fitted_Bins.coordinates = 'dgg_lat dgg_lon'
            RMSD_V_Fitted_Bins.units = 'Kelvin'
            RMSD_V_Fitted_Bins.long_name = 'Angle class root mean square deviation (vertical polarization, bin center ±2.5 degree) between filtered brightness temperature used to fit and the corresponding fitted ones at same incidence angles.'
            RMSD_V_Fitted_Bins[:] = f['RMSD_V_Fitted_Bins'][:]

            RMSD_V_L1C_Fit = ds.createVariable('RMSD_V_L1C_Fit', 'float32', ('n_grid_points',), **kwargs)
            RMSD_V_L1C_Fit.coordinates = 'dgg_lat dgg_lon'
            RMSD_V_L1C_Fit.units = 'Kelvin'
            RMSD_V_L1C_Fit.long_name = 'Root mean square deviation (vertical polarization) between SCLF1C brightness temperature (directly interpolated and convert to ground reference) and corresponding fitted ones at same incidence angles.'
            RMSD_V_L1C_Fit[:] = f['RMSD_V_L1C_Fit'][:]

            RMSD_V_L1C_Fit_Bins = ds.createVariable('RMSD_V_L1C_Fit_Bins', 'float32', ('n_grid_points','inc'), **kwargs)
            RMSD_V_L1C_Fit_Bins.coordinates = 'dgg_lat dgg_lon'
            RMSD_V_L1C_Fit_Bins.units = 'Kelvin'
            RMSD_V_L1C_Fit_Bins.long_name = 'Angle class root mean square deviation (vertical polarization, bin center ±2.5 degree) between SCLF1C brightness temperature (directly interpolated and convert to ground reference) and corresponding fitted ones at same incidence angles.'
            RMSD_V_L1C_Fit_Bins[:] = f['RMSD_V_L1C_Fit_Bins'][:]

            STD_of_X_Processed = ds.createVariable('STD_of_X_Processed', 'float32', ('n_grid_points',), **kwargs)
            STD_of_X_Processed.coordinates = 'dgg_lat dgg_lon'
            STD_of_X_Processed.units = 'Kelvin'
            STD_of_X_Processed.long_name = 'Standard deviation of X-pol brightness temperature at antenna reference frame after filtering and interpolated (before transformed to ground reference frame).'
            STD_of_X_Processed[:] = f['STD_of_X_Processed'][:]

            STD_of_X_Processed_Bins = ds.createVariable('STD_of_X_Processed_Bins', 'float32', ('n_grid_points','inc'), **kwargs)
            STD_of_X_Processed_Bins.coordinates = 'dgg_lat dgg_lon'
            STD_of_X_Processed_Bins.units = 'Kelvin'
            STD_of_X_Processed_Bins.long_name = 'Angle class standard deviation (bin center ±2.5 degree) of X-pol brightness temperature at antenna reference frame after filtering and interpolated (before transformed to ground reference frame).'
            STD_of_X_Processed_Bins[:] = f['STD_of_X_Processed_Bins'][:]

            STD_of_X_L1C = ds.createVariable('STD_of_X_L1C', 'float32', ('n_grid_points',), **kwargs)
            STD_of_X_L1C.coordinates = 'dgg_lat dgg_lon'
            STD_of_X_L1C.units = 'Kelvin'
            STD_of_X_L1C.long_name = 'Standard deviation of SCLF1C X-pol brightness temperature at antenna reference frame.'
            STD_of_X_L1C[:] = f['STD_of_X_L1C'][:]

            STD_of_X_L1C_Bins = ds.createVariable('STD_of_X_L1C_Bins', 'float32', ('n_grid_points','inc'), **kwargs)
            STD_of_X_L1C_Bins.coordinates = 'dgg_lat dgg_lon'
            STD_of_X_L1C_Bins.units = 'Kelvin'
            STD_of_X_L1C_Bins.long_name = 'Angle class standard deviation (bin center ±2.5 degree) of SCLF1C X-pol brightness temperature at antenna reference frame.'
            STD_of_X_L1C_Bins[:] = f['STD_of_X_L1C_Bins'][:]

            STD_of_Y_Processed = ds.createVariable('STD_of_Y_Processed', 'float32', ('n_grid_points',), **kwargs)
            STD_of_Y_Processed.coordinates = 'dgg_lat dgg_lon'
            STD_of_Y_Processed.units = 'Kelvin'
            STD_of_Y_Processed.long_name = 'Standard deviation of Y-pol brightness temperature at antenna reference frame after filtering and interpolated (before transformed to ground reference frame).'
            STD_of_Y_Processed[:] = f['STD_of_Y_Processed'][:]

            STD_of_Y_Processed_Bins = ds.createVariable('STD_of_Y_Processed_Bins', 'float32', ('n_grid_points','inc'), **kwargs)
            STD_of_Y_Processed_Bins.coordinates = 'dgg_lat dgg_lon'
            STD_of_Y_Processed_Bins.units = 'Kelvin'
            STD_of_Y_Processed_Bins.long_name = 'Angle class standard deviation (bin center ±2.5 degree) of Y-pol brightness temperature at antenna reference frame after filtering and interpolated (before transformed to ground reference frame).'
            STD_of_Y_Processed_Bins[:] = f['STD_of_Y_Processed_Bins'][:]

            STD_of_Y_L1C = ds.createVariable('STD_of_Y_L1C', 'float32', ('n_grid_points',), **kwargs)
            STD_of_Y_L1C.coordinates = 'dgg_lat dgg_lon'
            STD_of_Y_L1C.units = 'Kelvin'
            STD_of_Y_L1C.long_name = 'Standard deviation of SCLF1C Y-pol brightness temperature at antenna reference frame.'
            STD_of_Y_L1C[:] = f['STD_of_Y_L1C'][:]

            STD_of_Y_L1C_Bins = ds.createVariable('STD_of_Y_L1C_Bins', 'float32', ('n_grid_points','inc'), **kwargs)
            STD_of_Y_L1C_Bins.coordinates = 'dgg_lat dgg_lon'
            STD_of_Y_L1C_Bins.units = 'Kelvin'
            STD_of_Y_L1C_Bins.long_name = 'Angle class standard deviation (bin center ±2.5 degree) of SCLF1C Y-pol brightness temperature at antenna reference frame.'
            STD_of_Y_L1C_Bins[:] = f['STD_of_Y_L1C_Bins'][:]

            Mean_accuracy_H_Obs = ds.createVariable('Mean_accuracy_H_Obs', 'float32', ('n_grid_points',), **kwargs)
            Mean_accuracy_H_Obs.coordinates = 'dgg_lat dgg_lon'
            Mean_accuracy_H_Obs.units = 'Kelvin'
            Mean_accuracy_H_Obs.long_name = 'Averaged accuracy measurement at horizontal polarization in SCLF1C.'
            Mean_accuracy_H_Obs[:] = f['Mean_accuracy_H_Obs'][:]

            Mean_accuracy_V_Obs = ds.createVariable('Mean_accuracy_V_Obs', 'float32', ('n_grid_points',), **kwargs)
            Mean_accuracy_V_Obs.coordinates = 'dgg_lat dgg_lon'
            Mean_accuracy_V_Obs.units = 'Kelvin'
            Mean_accuracy_V_Obs.long_name = 'Averaged accuracy measurement at vertical polarization in SCLF1C.'
            Mean_accuracy_V_Obs[:] = f['Mean_accuracy_V_Obs'][:]

            N_Alias_Bins = ds.createVariable('N_Alias_Bins', 'uint8', ('n_grid_points','inc'), **kwargs)
            N_Alias_Bins.coordinates = 'dgg_lat dgg_lon'
            N_Alias_Bins.units = 'dimensionless'
            N_Alias_Bins.long_name = 'Number of flagged observations (bin center ±2.5 degree) not inside the exclusive zone of Alias free in SCLF1C.'
            N_Alias_Bins[:] = f['N_Alias_Bins'][:]

            N_L1C_Bins = ds.createVariable('N_L1C_Bins', 'uint8', ('n_grid_points','inc'), **kwargs)
            N_L1C_Bins.coordinates = 'dgg_lat dgg_lon'
            N_L1C_Bins.units = 'dimensionless'
            N_L1C_Bins.long_name = 'Number of SCL1FC observations in bin center ±2.5 degree.'
            N_L1C_Bins[:] = f['N_L1C_Bins'][:]

            N_RFI_Point_Source_Bins = ds.createVariable('N_RFI_Point_Source_Bins', 'uint8', ('n_grid_points','inc'), **kwargs)
            N_RFI_Point_Source_Bins.coordinates = 'dgg_lat dgg_lon'
            N_RFI_Point_Source_Bins.units = 'dimensionless'
            N_RFI_Point_Source_Bins.long_name = 'Number of flagged observations (bin center ±2.5 degree) affected by point source radio frequency interference in SCLF1C.'
            N_RFI_Point_Source_Bins[:] = f['N_RFI_Point_Source_Bins'][:]

            N_RFI_Point_Source_Tail_Bins = ds.createVariable('N_RFI_Point_Source_Tail_Bins', 'uint8', ('n_grid_points','inc'), **kwargs)
            N_RFI_Point_Source_Tail_Bins.coordinates = 'dgg_lat dgg_lon'
            N_RFI_Point_Source_Tail_Bins.units = 'dimensionless'
            N_RFI_Point_Source_Tail_Bins.long_name = 'Number of flagged observations (bin center ±2.5 degree) affected by tails from point source radio frequency interference in SCLF1C.'
            N_RFI_Point_Source_Tail_Bins[:] = f['N_RFI_Point_Source_Tail_Bins'][:]

            N_RFI_Point_Source_or_Tail_Bins = ds.createVariable('N_RFI_Point_Source_or_Tail_Bins', 'uint8', ('n_grid_points','inc'), **kwargs)
            N_RFI_Point_Source_or_Tail_Bins.coordinates = 'dgg_lat dgg_lon'
            N_RFI_Point_Source_or_Tail_Bins.units = 'dimensionless'
            N_RFI_Point_Source_or_Tail_Bins.long_name = 'Number of flagged observations (bin center ±2.5 degree) affected by BOTH of point source and tails from point source radio frequency interference in SCLF1C.'
            N_RFI_Point_Source_or_Tail_Bins[:] = f['N_RFI_Point_Source_or_Tail_Bins'][:]

            N_RFI_Strong_Fitted_Bins = ds.createVariable('N_RFI_Strong_Fitted_Bins', 'uint8', ('n_grid_points','inc'), **kwargs)
            N_RFI_Strong_Fitted_Bins.coordinates = 'dgg_lat dgg_lon'
            N_RFI_Strong_Fitted_Bins.units = 'dimensionless'
            N_RFI_Strong_Fitted_Bins.long_name = 'Number of deleted brightness temperature (bin center ±2.5 degree) in the fitted data due to strong RFI (outside of acceptable geophysical ranges (TBh or Tbv: 50 to 340 K; TB3 or TB4: -50 to 50 K)), these deleted brightness temperatures were filled by linear b-spline interpolation.'
            N_RFI_Strong_Fitted_Bins[:] = f['N_RFI_Strong_Fitted_Bins'][:]

            N_RFI_other_Fitted_Bins = ds.createVariable('N_RFI_other_Fitted_Bins', 'uint8', ('n_grid_points','inc'), **kwargs)
            N_RFI_other_Fitted_Bins.coordinates = 'dgg_lat dgg_lon'
            N_RFI_other_Fitted_Bins.units = 'dimensionless'
            N_RFI_other_Fitted_Bins.long_name = 'Number of deleted brightness temperature (bin center ±2.5 degree) in the fitted data due to other suspected RFI or anomaly, these deleted brightness temperatures were filled by linear b-spline interpolation.'
            N_RFI_other_Fitted_Bins[:] = f['N_RFI_other_Fitted_Bins'][:]

            N_RFI_Strong_L1C_Bins = ds.createVariable('N_RFI_Strong_L1C_Bins', 'uint8', ('n_grid_points','inc'), **kwargs)
            N_RFI_Strong_L1C_Bins.coordinates = 'dgg_lat dgg_lon'
            N_RFI_Strong_L1C_Bins.units = 'dimensionless'
            N_RFI_Strong_L1C_Bins.long_name = 'Number of brightness temperature (bin center ±2.5 degree) due to strong RFI (outside of acceptable geophysical ranges (TBx or Tby: 50 to 340 K; real or image part of TBxy: -50 to 50 K)) at ground reference frame before filter and interpolations.'
            N_RFI_Strong_L1C_Bins[:] = f['N_RFI_Strong_L1C_Bins'][:]

            N_RFI_other_L1C_Bins = ds.createVariable('N_RFI_other_L1C_Bins', 'uint8', ('n_grid_points','inc'), **kwargs)
            N_RFI_other_L1C_Bins.coordinates = 'dgg_lat dgg_lon'
            N_RFI_other_L1C_Bins.units = 'dimensionless'
            N_RFI_other_L1C_Bins.long_name = 'Number of brightness temperature (bin center ±2.5 degree) due to other suspected RFI or anomaly at ground reference frame before filter and interpolations.'
            N_RFI_other_L1C_Bins[:] = f['N_RFI_other_L1C_Bins'][:]

            N_RFI_X_Bins = ds.createVariable('N_RFI_X_Bins', 'uint8', ('n_grid_points','inc'), **kwargs)
            N_RFI_X_Bins.coordinates = 'dgg_lat dgg_lon'
            N_RFI_X_Bins.units = 'dimensionless'
            N_RFI_X_Bins.long_name = 'Number of flagged observations (bin center ±2.5 degree) affected by radio frequency interference at X polarization in SCLF1C.'
            N_RFI_X_Bins[:] = f['N_RFI_X_Bins'][:]

            N_RFI_Y_Bins = ds.createVariable('N_RFI_Y_Bins', 'uint8', ('n_grid_points','inc'), **kwargs)
            N_RFI_Y_Bins.coordinates = 'dgg_lat dgg_lon'
            N_RFI_Y_Bins.units = 'dimensionless'
            N_RFI_Y_Bins.long_name = 'Number of flagged observations (bin center ±2.5 degree) affected by radio frequency interference at Y polarization in SCLF1C.'
            N_RFI_Y_Bins[:] = f['N_RFI_Y_Bins'][:]

            N_Unfiltered_Fitted_Bins = ds.createVariable('N_Unfiltered_Fitted_Bins', 'uint8', ('n_grid_points','inc'), **kwargs)
            N_Unfiltered_Fitted_Bins.coordinates = 'dgg_lat dgg_lon'
            N_Unfiltered_Fitted_Bins.units = 'dimensionless'
            N_Unfiltered_Fitted_Bins.long_name = 'Number of Unfiltered observations (bin center ±2.5 degree) in the fitted data.'
            N_Unfiltered_Fitted_Bins[:] = f['N_Unfiltered_Fitted_Bins'][:]

if __name__=='__main__':
    files_folder = r'/data2/public2/Satellite/Product/SMOS/SCLF1C_Refinement'
    files = FilterFileList(files_folder, ['SM_','MIR','.h5'])
    files = np.array_split(files, args.numparts)[args.part]
    for file in files:
        print(file)
        convert_h5f2nc(file)
        os.remove(file)