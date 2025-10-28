import xarray as xr
import xagg as xa
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import re
import os
import glob
import warnings
import datetime

class NotUniqueFile(Exception):
    """ Exception for when one file needs to be loaded, but the search returned multiple files """
    pass

# Function to convert integer to Roman values
def printRoman(number):
    # from https://www.geeksforgeeks.org/python-program-to-convert-integer-to-roman/
    num = [1, 4, 5, 9, 10, 40, 50, 90,
        100, 400, 500, 900, 1000]
    sym = ["I", "IV", "V", "IX", "X", "XL",
        "L", "XC", "C", "CD", "D", "CM", "M"]
    i = 12

    output = ''
    while number:
        div = number // num[i]
        number %= num[i]
 
        while div:
            output = output + sym[i]
            div -= 1
        i -= 1

    return output

def get_params():
    ''' Get parameters 
    
    Outputs necessary general parameters. 
    
    Parameters:
    ----------------------
    (none)
    
    
    Returns:
    ----------------------
    dir_list : dict()
        a dictionary of directory names for file system 
        managing purposes: 
            - 'raw':   where raw climate files are stored, in 
                        subdirectories by model/product name
            - 'proc':  where processed climate files are stored,
                        in subdirectories by model/product name
            - 'aux':   where aux files (e.g. those that transcend
                        a single data product/model) are stored
    '''

    # Dir_list
    dir_list = pd.read_csv('dir_list.csv')
    dir_list = {d:dir_list.set_index('dir_name').loc[d,'dir_path'] for d in dir_list['dir_name']}


    # Return
    return dir_list

dir_list = get_params()

def get_filepaths(source_dir = 'raw',
                  mod = None,
                  dir_list = dir_list,
                  col_namer = {'(hadley$)|(CMIP[0-9]$)':'forcing_dataset',
                                 'PDO$':'pdo_state',
                                  'AMO$':'amo_state',}):
    ''' Get filepaths of climate data, split up by CMIP filename component
    
    
    Uses modified CMIP5/6 filename standards used by Kevin Schwarzwald's 
    filesystem - in other words, with the additional optional "suffix" 
    between the daterange and the filetype extension. 
    
    Returns
    ------------
    df : pd.DataFrame
        A dataframe containing information for all files in 
        `dir_list[source_dir]/mod/*.nc/.zarr`, with the full filepath in the
        column `path`, and filename components `varname`, `freq`, 
        `model`, `exp`, and optionally `run`, `grid`, `time`, 
        'gwl', 'proj_method'/'proj_target' for  method',`suffix`, in their own
        columns. `grid` may be Nones if files use CMIP5 conventions, 
        `suffix` may be Nones if no suffixes are found. Last time modified
        (from :py:meth:`os.path.getmtime()`) is listed as 'mtime'. 
        
        If `exp` has a match for the regex r"-", then additionally
        extra columns for each experiment name component will be 
        created, if possible, using the `col_namer` input.
    
    
    '''
    
    def id_fncomps(comps,col_namer=col_namer):
        # Make sure there are enough components 
        if len(comps)<6:
            # For now - but there has to be a better way to 
            # flag this
            slots = {'varname':None}
        else:
            try:
                # Figure out filetype
                filetype = re.split(r'\.',comps[-1])[-1]                    
                
                # Prepopulate set components
                slots = {s:comps[n] for n,s in zip(np.arange(0,4),['varname','freq','model','exp'])}
                proc_comps = list(np.arange(0,4))
                
                # Find slot for run, which will be of the form r##i## or "reanalysis" or "ALLRUNS"
                run_match = [re.search(r'(r[0-9]{1,3}([ipf][0-9]{1,2}){1,3})|(ALLRUNS)|(allruns)|(reanalysis)|(obs)|(stations)|(run[0-9]{1,3})|(RUNS[0-9]{1,3})|(ens[0-9]{1,3})',comp) 
                                   for comp in comps]
                if np.any(run_match):
                    slots['run'] = [k.group() for k in run_match
                                     if k is not None][0]
                    proc_comps.append(np.where(run_match)[0][0])
                
                # Get which slot is the timeframe (fx have "na" as timeframe)
                time_match = [re.search('([0-9]{4,8}'+r'-'+'[0-9]{4,8})|(^na($|'+r'.'+'))|(ALLPERIODS)',comp) 
                                   for comp in comps]
                if np.any(time_match):
                    slots['time'] = [k.group() for k in time_match
                                     if k is not None][0]
                    proc_comps.append(np.where(time_match)[0][0])
                
                # Determine whether there's a grid slot 
                # (assuming of the form 'gX(X)')
                grid_match = [re.search('^g[a-z0-9]{1,2}$',comp) for comp in comps]
                if np.any(grid_match): 
                    slots['grid'] = [k.group() for k in grid_match
                                     if k is not None][0]
                    proc_comps.append(np.where(grid_match)[0][0])

                # GWL slot
                gwl_match = [re.search('(^GWL)|(ALLGWLs)|(ALLGWLS)',comp) 
                           for comp in comps]
                if np.any(gwl_match):
                    slots['gwl'] = np.where(gwl_match)[0][0]
                    proc_comps.append(slots['gwl'])

                    if (comps[slots['gwl']] == 'ALLGWLs') or (comps[slots['gwl']] == 'ALLGWLS'):
                        slots['gwl'] = 'ALLGWLs'
                    else:
                        slots['gwl'] = re.sub(r'\-', '.', comps[slots['gwl']][3:None])

                # Projection slot
                proj_match = [re.search(r'^proj[a-zA-Z0-9]*\-base[a-zA-Z0-9\-]*',comp) 
                                   for comp in comps]
                if np.any(proj_match):
                    proj_info = [k.group() for k in proj_match
                                 if k is not None][0]
                    slots['proj_method'] = re.split(r'\-',re.split('proj',proj_info)[1])[0]
                    slots['proj_base'] = re.split('base',proj_info)[1]
                    proc_comps.append(np.where(proj_match)[0][0])

                # Downscaling slot
                dwscl_match = [re.search(r'^dwnscl[a-zA-Z0-9]*\-target[a-zA-z0-9\-]*',comp) 
                                   for comp in comps]
                if np.any(dwscl_match):
                    dwnscl_info = [k.group() for k in dwscl_match
                                 if k is not None][0]
                    slots['dwnscl_method'] = re.split(r'\-',re.split('dwnscl',dwnscl_info)[1])[0]
                    slots['dwnscl_target'] = re.split('target',dwnscl_info)[1]
                    proc_comps.append(np.where(dwscl_match)[0][0])

                # Seasstats slot
                seasstats_match = [re.search(r'seasstats',comp)
                                   for comp in comps]
                if np.any(seasstats_match):
                    seasstats_match = np.where(seasstats_match)[0][0]
                    slots['seasstats'] = comps[seasstats_match+1]
                    proc_comps.append(seasstats_match)
                    proc_comps.append(seasstats_match+1)


                # Get remaining, otherwise unlabeled components
                unlabeled_idxs = np.array([idx for idx in range(len(comps)) if idx not in np.array(proc_comps)])
                if len(unlabeled_idxs) > 0:
                    comps = np.array(comps)[unlabeled_idxs]
                else:
                    comps = []

                # If a remaining, unlabeled slot has the file ending in it, it's the "suffix"
                if np.any([re.search(r'[a-zA-z0-9]*\.'+filetype+r'$',comp) for comp in comps]):
                    slots['suffix'] = [k.group() for k in [re.search(r'[a-zA-z0-9]*\.'+filetype+r'$',comp) for comp in comps]
                                       if k is not None][0]
                    slots['suffix'] = re.split(r'\.',slots['suffix'])[0]
                    # Remove it from the list
                    comps = [comp for comp in comps if not re.search(r'[a-zA-z0-9]*\.'+filetype+r'$',comp)]

                # Add remaining unlabeled slots to just a generic catch-all 
                if len(comps)>0:
                    slots['unlabeled'] = str(comps)

                # If the experiment slot has multiple sub-experiments,
                # save them seperately using the column namer dict
                exp_comps = re.split(r'-',slots['exp'])
                if len(exp_comps)>1:
                    for exp_comp in exp_comps:
                        if np.any([re.search(k,exp_comp) for k in col_namer]):
                            match_type = [v for k,v in col_namer.items() if re.search(k,exp_comp)]
                            if len(match_type) > 1:
                                warnings.warn('More than one column match found for '+exp_comp+
                                              '. Check col_namer, no exp has been split.')
                            else:
                                slots[match_type[0]] = exp_comp

                # Add filetype as a slot
                slots['filetype'] = filetype
            except:
                # Assuming that if there's an error it's because the
                # file in question isn't in a standard respected form
                # For now - but there has to be a better way to 
                # flag this
                slots = {'varname':None}

        return slots

    #---------- Get list of files ----------
    if mod is None:
        # Get all mods
        mods = [re.split('/',mod)[-1] for mod in glob.glob(dir_list[source_dir]+'*')]
    else:
        mods = [mod]
        
    fns_all = [None]*len(mods)
    for mod,mod_idx in zip(mods,np.arange(0,len(mods))):
        # Get list of subdirectories (nc and zarr)
        fns = [*glob.glob(dir_list[source_dir]+mod+'/*.nc'),
               *glob.glob(dir_list[source_dir]+mod+'/*.zarr')]

        # Split up filename by components
        fn_comps = [re.split(r'_',re.split(r'/',fn)[-1]) for fn in fns]
        # Identify components, concatenate with path
        fns_all[mod_idx] = pd.DataFrame([id_fncomps(comps) for comps in fn_comps])
        fns_all[mod_idx] = pd.concat([fns_all[mod_idx],pd.DataFrame([{'path':fn} for fn in fns])],axis=1)

    # Concatentate into single df
    df = pd.concat(fns_all)

    # Add last time modified
    df.loc[:,'mtime'] = [datetime.datetime.fromtimestamp(os.path.getmtime(fn),
                            datetime.UTC)
                         for fn in df.path.values]

    #---------- Return ----------
    return df

# The next two are from https://towardsdatascience.com/the-correct-way-to-average-the-globe-92ceecd172b7

def earth_radius(lat):
    '''
    calculate radius of Earth assuming oblate spheroid
    defined by WGS84
    
    Input
    ---------
    lat: vector or latitudes in degrees  
    
    Output
    ----------
    r: vector of radius in meters
    
    Notes
    -----------
    WGS84: https://earth-info.nga.mil/GandG/publications/tr8350.2/tr8350.2-a/Chapter%203.pdf
    '''
    from numpy import deg2rad, sin, cos

    # define oblate spheroid from WGS84
    a = 6378137
    b = 6356752.3142
    e2 = 1 - (b**2/a**2)

    # convert from geodecic to geocentric
    # see equation 3-110 in WGS84
    lat = deg2rad(lat)
    lat_gc = np.arctan( (1-e2)*np.tan(lat) )

    # radius equation
    # see equation 3-107 in WGS84
    r = (
        (a * (1 - e2)**0.5)
         / (1 - (e2 * np.cos(lat_gc)**2))**0.5
        )

    return r

def area_grid(lat, lon):
    """
    Calculate the area of each grid cell
    Area is in square meters
    
    Input
    -----------
    lat: vector of latitude in degrees
    lon: vector of longitude in degrees
    
    Output
    -----------
    area: grid-cell area in square-meters with dimensions, [lat,lon]
    
    Notes
    -----------
    Based on the function in
    https://github.com/chadagreene/CDT/blob/master/cdt/cdtarea.m
    """
    from numpy import meshgrid, deg2rad, gradient, cos
    from xarray import DataArray

    xlon, ylat = meshgrid(lon, lat)
    R = earth_radius(ylat)

    dlat = deg2rad(gradient(ylat, axis=0))
    dlon = deg2rad(gradient(xlon, axis=1))

    dy = dlat * R
    dx = dlon * R * cos(deg2rad(ylat))

    area = dy * dx

    xda = DataArray(
        area,
        dims=["lat", "lon"],
        coords={"lat": lat, "lon": lon},
        attrs={
            "long_name": "area_per_pixel",
            "description": "area per pixel",
            "units": "m^2",
        },
    )
    return xda

def area_mean(ds,assume_rectangular=True):
    """ Calculate area-weighted mean of all variables in a  dataset
    
    Mean over lat / lon, weighted by the relative size of each
    pixel, dependent on latitude. Only weights by latitude, does
    not take into account lat/lon bounds, if present. 
    
    Parameters
    ------------------
    ds : xr.Dataset
    
    Returns
    ------------------
    dsm : xr.Dataset
        The input dataset, `ds`, averaged.
    
    """
    
    if (ds.sizes['lat'] == 1) and (ds.sizes['lon'] == 1):
        # If just one pixel, return that one pixel
        ds = ds.isel(lat=0,lon=0).drop(['lat','lon'])
        
    elif (ds.sizes['lat'] == 1) and (assume_rectangular):
        # If only one lat row, but multiple long rows, 
        # just get the cartesian mean, if assuming rectangular
        # grids. 
        ds = ds.mean(('lat','lon'))
        
    else:
        # Calculate area in each pixel
        weights = area_grid(ds.lat,ds.lon)

        # Remove nans, to make weight sum have the right magnitude
        weights = weights.where(~np.isnan(ds))

        # Calculate mean
        ds = ((ds*weights).sum(('lat','lon'))/weights.sum(('lat','lon')))

    # Return 
    return ds

def utility_save(ds,output_fn,dir_list=None,raw_overwrite_flag=False,create_dir=True,
                 keep_chunk_encoding = True, save_kwargs = {}):
    ''' Save xarray dataset as netcdf or zarr file, with safeguards
    By default overwrites `output_fn`, *unless* `output_fn` is in the 
    raw data directory as defined by `dir_list`. Creates the implied
    directory if it does not already exist.

    Parameters
    ---------------
    ds : :py:class:`xr.Dataset`

    output_fn : :py:class:`str`

    dir_list : :py:class:`dict` or `None`, default `None`
        If None, then directories are grabbed using `get_params()`.
        Otherwise, put in a manual `dir_list` - which only requires
        `['raw']` as a field (to test whether anything in the 
        `raw` directory is being touched

    raw_overwrite_flag : :py:class:`bool`, default False
        If False, then if `output_fn` already exists in the `dir_list['raw']` 
        directory, an error is raised instead of overwriting the file

    create_dir : :py:class`bool`, default True
        If True, then creates the implied directory (using 
        `os.path.dirname(output_fn)`) if it does not yet exist. 
        
    '''
    
    if dir_list is None:
        dir_list = get_params()

    if not os.path.exists(os.path.dirname(output_fn)):
        os.mkdir(os.path.dirname(output_fn))
        print(os.path.dirname(output_fn)+' created!')

    if os.path.exists(output_fn):
        if not raw_overwrite_flag:
            if re.search(r'^'+dir_list['raw'],output_fn):
                raise FileExistsError('Trying to overwrite a file in the "raw" data directory '+dir_list['raw']+'. '+
                                      'If this is on purpose, set `raw_overwrite_flag=True`.\n'+
                                      'Attempted output filename: '+output_fn)

        os.remove(output_fn)
        print(output_fn+' removed to allow overwrite!')

    if not keep_chunk_encoding:
        from funcs_aux import _remove_chunk_encoding
        ds = _remove_chunk_encoding(ds)

    if re.search(r'\.zarr$',os.path.basename(output_fn)):
        # If saving zarr
        ds.to_zarr(output_fn,**save_kwargs)
    elif re.search(r'\.nc$',os.path.basename(output_fn)):
        # If saving netcdf
        ds.to_netcdf(output_fn,**save_kwargs)
    print(output_fn+' saved!')
    

def utility_print(output_fn,formats=['pdf','png']):
    if 'pdf' in formats:
        plt.savefig(output_fn+'.pdf')
        print(output_fn+'.pdf saved!')

    if 'png' in formats:
        plt.savefig(output_fn+'.png',dpi=300)
        print(output_fn+'.png saved!')

    if 'svg' in formats:
        plt.savefig(output_fn+'.svg')
        print(output_fn+'.svg saved!')


def get_varlist(source_dir=None,var=None,varsub='all',
                experiment=None,freq=None,
                empty_warnings=False):
    ''' Get a list of which models have which variables
    
    Searches the filesystem for all models (directory names) and 
    all variables (first part of filenames, before the first 
    underscore), and returns either that information for all 
    models and variables, or an array of models that have 
    files for specified variables. 
    
    NB: if no experiment or frequency is specified, and the
    full dataframe is returned (`var=None`), then the fields
    have True whenever any file with that variable in the filename
    for that model is present (and potentially more than one). 
    In general, the code does not differentiate between multiple
    files for a single model/variable combination. 
    
    Parameters
    ---------------
    source_dir : str; default dir_list['raw']
        a path to the directory with climate data (all 
        subdirectories are assumed to be models, all files in
        these directories are assumed to be climate data files
        in rough CMIP format).
        
    var : str, list; default `None`
        one variable name or a list of variables for which to 
        subset the model list of. If not `None`, then only a list
        of models for which this variable(s) is present is returned
        (instead of the full Dataframe).
        
    varsub : str; default 'all'
        - if 'all', then if `var` has multiple variables, 
          only models that have files for all of the variables 
          are returned
        - if 'any', then if `var` has multiple variables, 
          models that have files for any of the variables are 
          returned
          
    experiment : str; default `None`
        if not None, then only returns models / True if files
        for the given 'experiment' (in CMIP6 parlance, the 
        fourth filename component) are found. If not None, the
        variable is piped into re.search(), allowing for re
        searches for the experiment. 
        
    freq : str; default `None`
        if not None, then only returns models / True if files
        for the given 'frequency' (in CMIP6 parlance, the 
        second filename component) are found. If not None, the
        variable is piped into re.search(), allowing for re
        searches for the frequency. 
        
    empty_warnings : bool; default `False`
        if True, a warning is thrown if no files at all (before 
        subsetting) are found for a model. 
    
    
    Returns
    ---------------
    varindex : pd.DataFrame()
        if `var` is None, then a models x variables pandas
        DataFrame is returned, with `True` if that model has 
        a file with that variable, and `False` otherwise.
        
    mods : list
        if `var` is not None, then a list of model names 
        that have the variables, subject to the subsetting above
    
    
    '''
    if source_dir is None:
        dir_list = get_params()
        source_dir = dir_list['raw']
    
    
    ##### Housekeeping
    # Ensure the var input is a list of strings, and not a string
    if type(var) == str:
        var = [var]
    
    ##### Identify models
    # Figure out in which position of the filename path the model name
    # directory is located (based on how many directory levels there 
    # are in the parent directory)
    modname_idx = len(re.split('/',source_dir)) - 1
    # Get list of all the models (the directory names in source_dir)
    all_mods = [re.split('/',x)[modname_idx] for x in [x[0] for x in os.walk(source_dir)] if re.split('/',x)[modname_idx]!='']
    all_mods = [mod for mod in list(np.unique(all_mods)) if 'ipynb' not in mod]
    
    ##### Identify variables
    # Get list of all variables used and downloaded
    # Make this a pandas dataarray - mod x var
    varlist = []
    for mod in all_mods[:]:
        varlist.append([re.split(r'_',fn)[0] for fn in [x for x in os.walk(source_dir+mod+'/')][0][2]])
    varlist = [item for sublist in varlist for item in sublist]

    varlist = list(np.unique(varlist))

    # Remove "README" and ".nc" files 
    varlist = [var for var in [var for var in varlist if 'READ' not in var] if '.nc' not in var]
    
    ##### Populate dataframe
    # Create empty dataframe to populate with file existence
    varindex = pd.DataFrame(columns=['model',*varlist])

    # Populate the model column
    varindex['model'] = all_mods

    # Actually, just set the models as the index
    varindex = varindex.set_index('model')
    
    # Now populate the dataframe with Trues if that model has that variable as a file
    for mod in all_mods:
        # Get variable name of each file 
        file_varlist = [re.split(r'_',fn)[0] for fn in [x for x in os.walk(source_dir+mod+'/')][0][2]]

        if len(file_varlist) == 0:
            if empty_warnings:
                warnings.warn('No relevant files found for model '+mod)
            varindex.loc[mod] = False
        else:
            # Subset by frequency, or experiment, if desired
            if freq is not None:
                try:
                    freq_bools = [(re.search(freq,re.split(r'_',fn)[1]) != None) for fn in [x for x in os.walk(source_dir+mod+'/')][0][2]]
                except IndexError:
                    freq_bools = [False]*len(file_varlist)
                    if empty_warnings:
                        warnings.warn('Model '+mod+' has files not in CMIP format.')
                    continue
            else:
                freq_bools = [True]*len(file_varlist)

            if experiment is not None:
                try:
                    exp_bools = [(re.search(experiment,re.split(r'_',fn)[3]) != None) for fn in [x for x in os.walk(source_dir+mod+'/')][0][2]]
                except IndexError:
                    exp_bools = [False]*len(file_varlist)
                    if empty_warnings:
                        warnings.warn('Model '+mod+' has files not in CMIP format.')
                    continue
            else:
                exp_bools = [True]*len(file_varlist)

            # Remove from list if it doesn't fit the frequency/experiment subset
            file_varlist = list(np.asarray(file_varlist)[np.asarray(freq_bools) & np.asarray(exp_bools)])

            # Add to dataframe
            varindex.loc[mod] = [var in file_varlist for var in varlist]

    # Fill NaNs with False
    varindex = varindex.fillna(False)

    ##### Return
    if var is None: 
        return varindex
    else:
        if type(var) == str:
            var = [var]
        if varsub == 'all':
            # (1) is to ensure the `all` is across variables/columns, not rows/models
            return list(varindex.index[varindex[var].all(1)].values)
        elif varsub == 'any':
            return list(varindex.index[varindex[var].any(1)].values)
        else:
            raise KeyError(str(varsub) + ' is not a supported variable subsetting method, choose "all" or "any".')

def id_timeframe(r,cond = 'longest',out = 'timestr'):
    ''' Choose between filepaths based on a temporal condition

    Subset output from :py:meth:`get_filepaths()` based on a condition
    on the timelength of a file. 

    Parameters
    -----------------
    r : :py:class:`pd.Dataframe`, output from :py:meth:`get_filepaths`
        Crucially, needs the "time" column

    cond : :py:class:`str`
        One of:
            - 'longest': choose the file with the longest timeframe
            - 'shortest': choose the file with the shortest timeframe
            - 'earliest': choose the file that starts the earliest
            - 'latest': choose the file that ends the latest

    out : :py:class:`str`
        Determines the return, one of: 
            - 'timestr': the timestring of the relevant file
            - 'df': the full dataframe with just that file's row

    Returns
    -----------------
    depends on `out` above

    TODO: allow arbitrary nested conditions (so, if two have the same
    length, choose one using a different condition)
    
    '''
    # Allowable conds
    conds = ['longest','shortest','earliest','latest']
    
    # Get timeframes of each file
    ts = r.time.values

    # Get timeframe in times
    def convert_to_dt(t):
        t_out = []
        for dt in re.split(r'-',t):
            try:
                t_out.append(pd.to_datetime(dt,format='%Y%m%d'))
            except pd.errors.OutOfBoundsDatetime as e1:
                # Dates past 2262 can't fit into pandas, gotta go to 
                # base datetime
                t_out.append(datetime.datetime(int(dt[0:4]),int(dt[4:6]),int(dt[6:8])))
            except ValueError as e2:
                # Some files got saved with incorrect timeframe slots, 
                # i.e., 31 in any month, even in months with fewer days.
                # If the days slot is larger than the # days that should 
                # be in that month, convert to the last day of that month 
                # before turning to datetime format. 
                if int(dt[-2:]) > pd.to_datetime(dt[0:6]+'01',format='%Y%m%d').daysinmonth:
                    t_out.append(pd.to_datetime(dt[0:6]+str(pd.to_datetime(dt[0:6]+'01',format='%Y%m%d').daysinmonth),
                                            format='%Y%m%d'))
                else:
                    raise e2
        return t_out

    tdts = [convert_to_dt(t) for t in ts]

    # Get length of each timeframe
    tlengths = [np.diff(tdt) for tdt in tdts]

    # Get timestring based on condition above
    if cond == 'longest':
        idx = np.argmax(tlengths)
    elif cond == 'shortest':
        idx = np.argmin(tlengths)
    elif cond == 'earliest':
        idx = np.argmin([t[0] for t in tdts])
    elif cond == 'latest':
        idx = np.argmax([t[1] for t in tdts])
    else:
        raise KeyError('cond must be one of '+', '.join(conds)+'.')

    # Return
    if out == 'timestr':
        return ts[idx]
    elif out == 'df':
        return r.iloc[idx,:]
    else:
        raise KeyError('out must be one of "timestr", "df"')


# This whole  business is because np.nanargmax (and its 
# derivatives, including ds.argmax in xr) can't deal if 
# the whole column of months is NaNs. So in that case, we 
# just keep the value NaN, and run nanargmax on the columns
# that do have values (land pixels in GPCC, for example)

def nan_argmax_xr(x,val=0,dim='month'):
    """ Get the index of each maximum month in the 'month'
    dimension of an arbitrary dataarray with dimensions
    'month' and others. This spits out NaN for any 
    row of months that's entirely NaNs, and therefore 
    provides a workaround for np.argmax() and 
    np.nanargmax(), which both fail in this situation.
    Furthermore, it automatically stacks/unstacks for 
    the calculation, so the input can have an arbitrary
    number of dimensions. 
    """
    
    # Stack to have [__ x month]
    input_dims = list(x.dims)
    
    if dim not in input_dims:
        raise LookupError("no '"+dim+"' dimension found.")
    
    input_dims.remove(dim)
    
    if len(input_dims)>1:
        x = x.stack(alld=(tuple(input_dims)))
        unstack = True
    else: 
        unstack = False
    
    if x.ndim>1:
        # Pre-build np.nan
        out_vals = np.zeros((np.shape(x)[np.argmax([key!=dim for key in x.dims])]))*np.nan
    
        #out_vals[~np.isnan(x[0,:])] = x[:,~np.isnan(x.values[0,:])].argsort(0).isel({dim:-1-val})
        #out_vals = xr.DataArray(out_vals,dims=x.dims[1],coords={x.dims[1]:x[x.dims[1]]})
        nan_idxs = np.isnan(x).all(dim).compute()
        
        if not np.all(nan_idxs): #else keep it nan
            out_vals[~nan_idxs] = x[:,~nan_idxs].argmax(dim)
        out_vals = xr.DataArray(out_vals,dims=x.dims[1],coords={x.dims[1]:x[x.dims[1]]})
    else:
        if np.all(np.isnan(x)):
            out_vals = np.nan
        else:
            out_vals = x.argmax(dim)
        # Pretty sure this just needs to be one value... to be consistent
        #out_vals[~np.isnan(x)] = x[~np.isnan(x.values)].argsort(0).isel({dim:-1-val})
        #out_vals = xr.DataArray(out_vals,dims=x.dims[0],coords={x.dims[0]:x[x.dims[0]]})
        out_vals = xr.DataArray(out_vals)
    
    if unstack:
        out_vals = out_vals.unstack()
    
    return out_vals

def nan_argmin_xr(x,val=0,dim='month'):
    """ Get the index of each maximum month in the 'month'
    dimension of an arbitrary dataarray with dimensions
    'month' and others. This spits out NaN for any 
    row of months that's entirely NaNs, and therefore 
    provides a workaround for np.argmax() and 
    np.nanargmax(), which both fail in this situation.
    Furthermore, it automatically stacks/unstacks for 
    the calculation, so the input can have an arbitrary
    number of dimensions. 
    """
    
    # Stack to have [__ x month]
    input_dims = list(x.dims)
    
    if dim not in input_dims:
        raise LookupError("no '"+dim+"' dimension found.")
    
    input_dims.remove(dim)
    
    if len(input_dims)>1:
        x = x.stack(alld=(tuple(input_dims)))
        unstack = True
    else: 
        unstack = False
    
    if x.ndim>1:
        # Pre-build np.nan
        out_vals = np.zeros((np.shape(x)[np.argmax([key!=dim for key in x.dims])]))*np.nan
        
        #out_vals[~np.isnan(x[0,:])] = x[:,~np.isnan(x.values[0,:])].argsort(0).isel({dim:-1-val})
        #out_vals = xr.DataArray(out_vals,dims=x.dims[1],coords={x.dims[1]:x[x.dims[1]]})
        nan_idxs = np.isnan(x).all(dim).compute()
        
        if not np.all(nan_idxs): #else keep it nan
            out_vals[~nan_idxs] = x[:,~nan_idxs].argmin(dim)
        out_vals = xr.DataArray(out_vals,dims=x.dims[1],coords={x.dims[1]:x[x.dims[1]]})
    else:
        if np.all(np.isnan(x)):
            out_vals = np.nan
        else:
            out_vals = x.argmin(dim)
        # Pretty sure this just needs to be one value... to be consistent
        #out_vals[~np.isnan(x)] = x[~np.isnan(x.values)].argsort(0).isel({dim:-1-val})
        #out_vals = xr.DataArray(out_vals,dims=x.dims[0],coords={x.dims[0]:x[x.dims[0]]})
        out_vals = xr.DataArray(out_vals)
    
    if unstack:
        out_vals = out_vals.unstack()
    
    return out_vals