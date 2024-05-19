import requests, os, shutil
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
        print("Downloading BT_SETTL grid...")
        with requests.get(url, stream=True) as r:
            with open(out, 'wb') as f:
                shutil.copyfileobj(r.raw, f)

    return flag


def read_spec():
    src = os.path.join(utils.dirs["utils"],"btsettl_full.nc")
    ds_disk = xr.open_dataset(src)

