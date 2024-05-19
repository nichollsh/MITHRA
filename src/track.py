import requests, os
import numpy as np

import src.utils as utils

def download_track(use_cache=True):
    flag:int = 0

    # Path to full data
    url = "https://perso.ens-lyon.fr/isabelle.baraffe/BHAC15dir/BHAC15_tracks+structure"
    out = os.path.join(utils.dirs["data"] , "bhac15_full.dat")

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


def read_track():

    # Source file 
    src = os.path.join(utils.dirs["data"], "bhac15_full.dat")

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
    arr_keys = ["arr_age",          
                "arr_Teff",      
                "arr_L/Ls",       
                "arr_logg",      
                "arr_R/Rs",     
                "arr_logli",
                "arr_tcentral", 
                "arr_rhocentral", 
                "arr_mcore", 
                "arr_rcore", 
                "arr_k2c", 
                "arr_k2r"
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
        t["arr_age"] = 10.0 ** t["arr_age"]

    return tracks
    
