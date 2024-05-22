import requests, os
import numpy as np
import scipy.interpolate as interp

import mithra.utils as utils

def download_bhac(use_cache=True)->int:
    flag:int = 0

    # Path to full data
    url = "https://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/BHAC15_tracks+structure"
    out = os.path.join(utils.dirs["data"] , "bhac15_full.dat")

    # Check if path exists
    if not use_cache:
        utils.rmsafe(out)

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


def read_bhac()->dict:

    # Source file 
    src:str = os.path.join(utils.dirs["data"], "bhac15_full.dat")

    # Read file
    with open(src,"r") as hdl:
        raw = hdl.readlines()

    # skip header 
    lines = raw[43:]
    nlines = len(lines)

    # mass tracks 
    tracks = []

    # loop through file
    i = 0
    arr_keys = ["age",          
                "Teff",      
                "L/Ls",       
                "logg",      
                "R/Rs",     
                "logli",
                "tcentral", 
                "rhocentral", 
                "mcore", 
                "rcore", 
                "k2c", 
                "k2r"
                ]
    track = {}
    while i < nlines:

        l = lines[i][1:-1].strip()
        if len(l) == 0:
            i += 1
            continue 

        if "----------" in l:
            # skip header
            i += 3
            tracks.append(track)
            track = {"mass":0.0}
            for k in arr_keys:
                track[k] = []
            continue

        # read data
        lsplit = l.split(" ")
        values = []
        for l in lsplit:
            l = l.replace(' ','')
            if len(l) > 0:
                values.append(float(l))
        track["mass"] = values[0]
        for j,key in enumerate(arr_keys):
            track[key].append(values[j+1])

        i += 1

    # Drop dummy value
    tracks = tracks[1:]

    # Rescale units
    for t in tracks:
        for k in arr_keys:
            t[k] = np.array(t[k], dtype=float)

        # Convert logt to years
        t["age"] = 10.0 ** t["age"]

    return tracks
    
def get_params_bhac(tracks:dict, mass:float, age:float, params:list)->list:

    # Find closest track
    itrack:int   = -1
    dtrack:float = 1e99
    for i,t in enumerate(tracks):
        d = abs(t["mass"] - mass)
        if d < dtrack:
            itrack = i
            dtrack = d
    track = tracks[itrack]

    # Interpolate over time, for each param
    age = max(age, np.amin(track["age"]))
    age = min(age, np.amax(track["age"]))

    out = []
    for p in params:
        out.append(interp.pchip_interpolate(track["age"],track[p],age))
    return out

