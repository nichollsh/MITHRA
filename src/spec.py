import requests, os, shutil
import numpy as np
import xarray as xr

import src.utils as utils

def download_grid(use_cache=True):
    flag:int = 0

    # Path to full data
    url = "https://zenodo.org/records/8015969/files/BT_SETTL_full.nc?download=1"
    out = os.path.join(utils.dirs["data"] , "btsettl_full.nc")

    # Check if path exists
    if os.path.exists(out) and (not use_cache):
        os.remove(out)

    # Download if required
    if not os.path.exists(out):
        print("Downloading BT_SETTL grid...")
        with requests.get(url, stream=True) as r:
            with open(out, 'wb') as f:
                shutil.copyfileobj(r.raw, f)

    return flag

def read_file():
    src = os.path.join(utils.dirs["data"],"btsettl_full.nc")
    if not os.path.exists(src):
        print("WARNING: Cannot find BT_SETTL file!")
        return None
    return xr.open_dataset(src)

def get_axes(ds):
    arr_teff = ds.coords['par1']
    arr_logg = ds.coords['par2']
    print(arr_teff)
    arr_wave = ds.coords["wavelength"]

    return np.array(arr_wave), np.array(arr_teff), np.array(arr_logg)

def get_spec(ds, teff, logg):
    
    # get axes
    arr_wave, arr_teff, arr_logg = get_axes(ds)

    # get indices
    i = np.argmin(np.abs(arr_teff - teff))
    j = np.argmin(np.abs(arr_logg - logg))

    print(arr_teff[i])
    print(arr_logg[j])

    # get data 
    flux_grid = ds.data_vars["grid"]

    # return spectrum
    return arr_wave, np.array(flux_grid[:,i,j])



