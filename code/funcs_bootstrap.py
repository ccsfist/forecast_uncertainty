# Functions for block-bootstrapping

from numba import jit as njit
import numpy as np


@njit
def get_boot_threshs_1d(ts, start_years, samp_size, q):
    ''' numbaized block bootstrap quantile calculation code ''' 
    n_run, n_block = start_years.shape

    ts = ts.astype(np.float32)
    # Create empty 
    assert ts.ndim==1
    out = np.empty(n_run,dtype=ts.dtype)

    # Process by run
    if not np.all(np.isnan(ts)):
        for run_idx in range(n_run):
            draws = np.empty((n_block, samp_size), dtype=ts.dtype)
            for b in range(n_block):
                sy = start_years[run_idx, b]
                draws[b, :] = ts[sy:sy + samp_size]
            # Ensure final time series is same length as original, without overhangs
            draws = draws.flatten()[0:ts.shape[0]]
            out[run_idx] = np.quantile(draws, q)
    else:
        out = (out*np.nan).astype(np.float32)
        
    out = out.astype(np.float32)
    
    return out

@njit
def get_boot_threshs_2d(ts, start_years, samp_size, q):
    ''' numbaized block bootstrap quantile calculation code ''' 
    n_run, n_block = start_years.shape

    ts = ts.astype(np.float32)
    # Create empty 
    assert ts.ndim==2
    out = np.empty((ts.shape[1],n_run),dtype=ts.dtype)

    # Process by run
    if not np.all(np.isnan(ts)):
        for run_idx in range(n_run):
            for idv_idx in range(ts.shape[1]):
                draws = np.empty((n_block, samp_size), dtype=ts.dtype)
                for b in range(n_block):
                    sy = start_years[run_idx, b]
                    draws[b, :] = ts[sy:sy + samp_size,idv_idx]
                # Ensure final time series is same length as original, without overhangs
                draws = draws.flatten()[0:ts.shape[0]]
                out[idv_idx,run_idx] = np.quantile(draws, q)
    else:
        out = (out*np.nan).astype(np.float32)
    out = out.astype(np.float32)
    
    return out



def boot_threshs_block(ds,samp_size,aa_trigger,
                       var='pr',boot_dim = 'year',
                       idx_var = 'start_years',
                       ):
    ''' quantile calculation of block bootstrap samples, optimized for xr.map_blocks() 

    Performs block boostrapping along of variable `var` along dimension `boot_dim`

    Parameters
    ---------------
    ds : xr.Dataset
        Dataset containing both `var` (below) _and_ 
        `idx_var` (below), a `draw` x `block` datarray
        giving the indices of the start years of 
        each bootstrap block. Note that the `idx_var` 
        must be integers of block starts along the 
        dimension `boot_dim`. 

    var : str, by default 'pr'
        The variable in `ds` to process

    idx_var : str, by default 'start_years'
        The `draw` x `block` variable giving block starts,
        with `draw` representing separate time series draws,
        and `block` representing each individual block of the
        time series. 

    boot_dim : str, by default 'year'
        The dimension of `var` across which to bootstrap

    aa_trigger : float
        The quantile to calculate (e.g., 0.2)

    samp_size : int
        The bootstrap block size (e.g., 5)

    Returns
    ---------------
    threshs : xr.DataArray
        The quantiles 
    
    
    '''
    # Transpose main variable ensure correct orders for numba / 
    # raw numpy below
    main_dim_order = ['lat','lon',boot_dim]
    other_dims = [dim for dim in ds[var].dims if dim not in main_dim_order]
    ds[var] = ds[var].transpose(*main_dim_order,*other_dims)

    # Stack dimensions that aren't lat, lon, and the 
    # dimension to bootstrap over, to make np.quantile
    # work in numba (the axis argument isn't supported)
    if len(other_dims)>1:
        ds = ds.stack(idv = other_dims)
        stacked = True
    else:
        stacked = False
        
    # Transpose start years to ensure correct orders
    # for numba / raw numpy below
    ds['start_years'] = ds['start_years'].transpose('draw','block')

    # Create empty array for output
    if len(other_dims)==0:
        threshs = np.empty((ds.sizes['lat'], ds.sizes['lon'], ds.sizes['draw']), dtype=ds[var].dtype)
    else:
        # (flexible, since it's only "idv" in the case of more than one "other_dim",
        # otherwise, it's the original name of that dimension)
        odim_name = [dim for dim in ds[var].dims if dim not in main_dim_order][0]
        threshs = np.empty((ds.sizes['lat'],ds.sizes['lon'],ds.sizes[odim_name],ds.sizes['draw']),
                       dtype=ds[var].dtype)

    # Run block bootstrap threshold calculation
    for i in range(ds.sizes['lat']):
        for j in range(ds.sizes['lon']):
            ts = ds[var].isel(lat=i,lon=j).values
            if ts.ndim == 1:
                threshs[i,j,:] = get_boot_threshs_1d(ts,ds['start_years'].values,samp_size,q=aa_trigger)
            elif ts.ndim == 2:
                threshs[i,j,:] = get_boot_threshs_2d(ts,ds['start_years'].values,samp_size,q=aa_trigger)
            else:
                raise Error

    # Assign to ds to allow unstacking, if necessary
    if len(other_dims) == 0:
        ds['threshs'] = (('lat','lon','draw'),threshs)
    else:
        ds['threshs'] = (('lat','lon',odim_name,'draw'),threshs)

    # Unstack, if necessary
    if stacked:
        out_da = ds['threshs'].unstack()
    else:
        out_da = ds['threshs']

    out_da = out_da.transpose('lat','lon',*other_dims,'draw')

    
    return out_da