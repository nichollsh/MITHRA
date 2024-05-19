import requests, os
import numpy as np
import xarray as xr

import src.utils as utils

def download_spec(use_cache=True):
    flag:int = 0

    # Path to full data
    url = "https://zenodo.org/records/8015969/files/BT_SETTL_full.nc?download=1"
    out = os.path.join(utils.dirs["data"] , "btsettl_full.nc")

    # Check if path exists
    if os.path.exists(out) and (not use_cache):
        os.remove(out)

    # Download if required
    if not os.path.exists(out):

        # Check server
        response = requests.get(url)
        
        # Download file if request is ok
        if response.status_code == 200:
            utils.rmsafe(out)
            with open(out, 'wb') as file:
                file.write(response.content)
        else:
            return 1
    
    return flag


def read_spec():
    src = os.path.join(utils.dirs["utils"],"btsettl_full.nc")
    ds_disk = xr.open_dataset(src)

