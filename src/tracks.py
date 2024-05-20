import requests, os
import numpy as np

import src.utils as utils

def download_bhac(use_cache=True)->int:
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
    
def get_params_bhac(tracks:dict, mass:float, age:float)->dict:

    # Find closest track
    itrack:int   = -1
    dtrack:float = 1e99
    for i,t in enumerate(tracks):
        d = abs(t["mass"] - mass)
        if d < dtrack:
            itrack = i
            dtrack = d
    track = tracks[itrack]
    print("Best track = %d (%.2e)"%(itrack,track["mass"]))

    # Find closest age 
    iage:int   = -1
    dage:float = 1e99
    for i,t in enumerate(track["arr_age"]):
        d = abs(age-t)
        if d < dage:
            iage = i
            dage = d 

    print("Best age   = %d (%.2e yr)"%(iage,track["arr_age"][iage]))
    
    # Get parameters at this mass,age
    out = {}
    for key in track.keys():
        if "arr" in key:
            key_short = key.replace("arr_","")
            out[key_short] = track[key][iage]
        else:
            out[key] = track[key]

    # Return parameters
    return out 

