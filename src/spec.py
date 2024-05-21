import glob, os
import numpy as np
import scipy.interpolate as interp

import src.utils as utils

def get_params_from_name(fpath:str):
    name = fpath.split("/")[-1].split(".")[0]
    splt = name.split("_")
    teff = float(splt[-2])
    logg = float(splt[-1])/100.0
    return teff,logg

def list_files():
    return glob.glob(os.path.join(utils.dirs["data"], "btsettl-cifist*.npy"))

def get_axes():
    # get all files
    files = list_files()
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

    return arr_teff, arr_logg

def create_interp(num_wl=100):

    # get wavelength data from first file 
    data = np.load(list_files()[0])
    min_wave, max_wave = np.amin(data[0]), np.amax(data[1])

    # limit wl range
    max_wave = min(max_wave, 1e5)
    min_wave = max(min_wave, 1.0)

    target_wave = np.logspace(np.log10(min_wave+0.1), np.log10(max_wave-0.1), num_wl)
    print(len(target_wave))

    # flattened data from files
    flat_teff = []
    flat_logg = []
    flat_wave = []
    flat_flux = []
    for i,f in enumerate(list_files()):
        # print("Read file %d"%i)
        # get data
        data = np.load(f)
        t,l = get_params_from_name(f)
        w,f = data[0],data[1]
        # print("    %d"%len(data[0]))

        # interpolate flux
        interp_flux = np.interp(target_wave, w, f)

        # store
        flat_wave.extend(list(target_wave))
        flat_flux.extend(list(interp_flux))
        flat_teff.extend(list(np.ones(len(target_wave))*t))
        flat_logg.extend(list(np.ones(len(target_wave))*l))

    print(len(flat_flux))
    print(len(flat_teff))

    # unique
    uniq_teff = np.unique(flat_teff)
    uniq_logg = np.unique(flat_logg)

    # target grid to interpolate to
    xi = np.linspace(np.amin(uniq_teff), np.amax(uniq_teff), len(uniq_teff))
    yi = np.linspace(np.amin(uniq_logg), np.amax(uniq_logg), len(uniq_logg))
    zi = target_wave


    print(len(xi))
    print(len(yi))
    print(len(zi))
    print("meshgrid...",flush=True)
    xi,yi,zi = np.meshgrid(xi,yi,zi)

    print(np.shape(xi))
    print(np.shape(yi))
    print(np.shape(zi))
    print("interpolate...",flush=True)

    # interpolate
    vi = interp.griddata((flat_teff,flat_logg,flat_wave),flat_flux,(xi,yi,zi),method='linear')

    print("done")
    return vi, xi, yi, zi

def get_spec(teff, logg, extend_planck=False, xmax=1.0e9):

    # get axes 
    arr_teff, arr_logg = get_axes()
    
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


