import glob, os
import numpy as np

import src.utils as utils

def get_params_from_name(fpath:str):
    name = fpath.split("/")[-1].split(".")[0]
    splt = name.split("_")
    teff = float(splt[-2])
    logg = float(splt[-1])/100.0
    return teff,logg

def get_spec(teff, logg, extend_planck=False, xmax=1.0e9):

    # get all files
    files = glob.glob(os.path.join(utils.dirs["data"], "btsettl-cifist*.npy"))
    names = [f.split("/")[-1].split(".")[0] for f in files]

    # get all teff,logg
    arr_teff = []
    arr_logg = []
    for f in names:
        t,l = get_params_from_name(f)
        arr_teff.append(t)
        arr_logg.append(l)
    arr_teff = np.unique(arr_teff)
    arr_logg = np.unique(arr_logg)
    
    # get best indices
    grid_i = np.argmin(np.abs(arr_teff - teff))
    grid_j = np.argmin(np.abs(arr_logg - logg))

    # get data 
    best = os.path.join(utils.dirs["data"], "btsettl-cifist_%04d_%04d.npy"%(arr_teff[grid_i],100*arr_logg[grid_j]))
    data = np.load(best)
    wl = data[0]
    fl = data[1]

    # truncate?
    if wl[-1] > xmax:
        imax = np.argmin(np.abs(wl - xmax))
        wl = wl[:imax]
        fl = fl[:imax]

    print(len(wl))

    # extend with planck function?
    if extend_planck and (wl[-1] < xmax):
        dx = wl[-1]*0.01
        print(dx)
        xp = np.arange(wl[-1]+dx, xmax, dx)
        yp = np.zeros(np.shape(xp))
        for i,x in enumerate(xp):
            lam = x * 1.0e-9  # nm -> m

            # Calculate planck function value [W m-2 sr-1 m-1]
            # http://spiff.rit.edu/classes/phys317/lectures/planck.html
            yp[i] = 2.0 * utils.h_pl * utils.c_vac**2.0 / lam**5.0   *   1.0 / ( np.exp(utils.h_pl * utils.c_vac / (lam * utils.k_B * arr_teff[grid_i])) - 1.0)

            # Integrate solid angle (hemisphere), convert units
            yp[i] = yp[i] * np.pi * 1.0e-9 # [W m-2 nm-1]
            yp[i] = yp[i] * 1000.0 # [erg s-1 cm-2 nm-1]

        wl = np.concatenate((wl, xp))
        fl = np.concatenate((fl, yp))

    # return spectrum
    return wl,fl


