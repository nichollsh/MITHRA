import glob, os, gc, subprocess
import xarray as xr
import numpy as np
import scipy.interpolate as interp

import src.utils as utils

def download_npy(url:str):
    '''
    Scrape a stellar spectrum from the Centro de Astrobiología website.
    
    This spectrum is saved to an npy file in the data folder of this package. 
    The npy file contains a 2D array. The first element is the wavelength array,
    and the second is the spectral flux array [erg s-1 cm-2 nm-1].

    Parameters
        url (str): url from which to download the file

    Returns
        None
    '''

    # download file 
    tmp = "/tmp/btsettl"
    if os.path.exists(tmp):
        os.remove(tmp)
    subprocess.run(["wget","-q","-O",tmp,url])

    # get name 
    with open(tmp,'r') as hdl:
        lines = hdl.readlines()
    teff = float(lines[1].split(" ")[3])
    logg = float(lines[2].split(" ")[3]) * 100.0
    out = os.path.join(utils.dirs["data"],"btsettl-cifist_%04d_%04d"%(teff,logg))
    if os.path.exists(out):
        os.remove(out)

    wl = []
    fl = []
    for i,l in enumerate(lines[9:]):
        s = l.split()
        if len(s) == 2:
            wl.append(float(s[0]))
            fl.append(float(s[1]))
        else:
            break

    # skip lines, to reduce file size
    pitch = 2
        
    # save data
    wl = np.array(wl[::pitch])*0.1   # convert to nm
    fl = np.array(fl[::pitch])/0.1   # convert to erg s-1 cm-2 nm-1
    np.save(out, np.array([wl,fl]))

    del wl
    del fl
    del lines
    return 

def download_all():
    '''
    Iteratively download all stellar spectra from the BT-SETTL CIFIST grid from
    the Centro de Astrobiología website. Saved to the data folder.

    Parameters
        None

    Returns
        None
    '''

    for i in range(1,447,1):
        print("Downloading %03d..."%i)
        url = "http://svo2.cab.inta-csic.es/theory/newov2/ssap.php?model=bt-settl-cifist&fid=%d&format=ascii"%i
        download_npy(url)
    return 

def get_params_from_name(fpath:str):
    '''
    Get model gridpoint parameters (Teff, logg) from its file name.

    Parameters
        fpath (str): path to npy file

    Returns
        teff (float): effective temperature [K]
        logg (float): log surface gravity [log cm/s^2]
    '''

    name = fpath.split("/")[-1].split(".")[0]
    splt = name.split("_")
    teff = float(splt[-2])
    logg = float(splt[-1])/100.0
    return teff,logg

def list_files():
    '''
    List all BT-SETTL npy files in the data directory

    Parameters
        None

    Returns
        files (list): list of absolute file paths
    '''
    files = glob.glob(os.path.join(utils.dirs["data"], "btsettl-cifist*.npy")) 
    return list(files)

def create_interp(num_wl=40, num_teff=0, num_logg=0, teff_lims=(1.0, 2e5), logg_lims=(1.0, 50.0)):
    '''
    Create interpolated grid of teff, logg, and wavelength from npy files.

    Reads all npy files in the data folder, creating a single data structure of 
    flux versus wavelength, teff, and logg. This grid may have missing values.
    The data are interpolated within the same (or reduced) teff and logg grid, 
    so no extrapolation is performed.

    Parameters
        num_wl (int): number of interpolatedwavelength points
        num_teff (int): number of interpolated teff points (0 for same as input)
        num_logg (int): number of interpolated logg points (0 for same as input)
        teff_lims (tuple): minimum and maximum teff values to use
        logg_lims (tuple): minimum and maximum logg values to use

    Returns
        vi (np.ndarray): array of fluxes in a 3D array (teff, logg, wavelength)
        xi (np.ndarray): array of teff [K] on same axes as above
        yi (np.ndarray): array of logg [log cm/s^2] on same axes as above
        zi (np.ndarray): array of wavelength [nm] on same axes as above
    '''

    # get wavelength data from first file 
    data = np.load(list_files()[0])
    min_wave, max_wave = np.amin(data[0]), np.amax(data[1])

    # limit wl range
    max_wave = min(max_wave, 1e5)
    min_wave = max(min_wave, 1.0)
    target_wave = np.linspace(np.log10(min_wave+0.1), np.log10(max_wave-0.1), num_wl)

    # flattened data from files
    flat_teff = []
    flat_logg = []
    flat_wave = []
    flat_flux = []
    for i,f in enumerate(list_files()):
        # get data
        data = np.load(f)
        t,l = get_params_from_name(f)
        w,f = np.log10(data[0]),data[1]

        # skip this point?
        if not(teff_lims[0] <= t <= teff_lims[1]):
            continue 
        if not(logg_lims[0] <= l <= logg_lims[1]):
            continue

        # drop duplicate values and ensure sorted
        _,mask = np.unique(w, return_index=True)
        w = w[mask]
        f = f[mask]

        # interpolate (downsample) flux and wavelength arrays...
        # this is to ensure that all spectra cover the same wavelength range
        # and also to significantly improve the performance of the unstructured
        # interpolation step by reducing the number of points
        w_ds = np.linspace(w[0],w[-1], num_wl*2)
        f_ds = interp.pchip_interpolate(w,f,w_ds)

        # store
        flat_wave.extend(list(w_ds))
        flat_flux.extend(list(f_ds))
        flat_teff.extend(list(np.ones(len(w_ds))*t))
        flat_logg.extend(list(np.ones(len(w_ds))*l))

        del w_ds 
        del f_ds 
        del w 
        del f 
        del data 
        gc.collect()

    # unique
    uniq_teff = np.unique(flat_teff)
    uniq_logg = np.unique(flat_logg)

    # output sizes
    if not(1 <= num_teff <= len(uniq_teff)):
        num_teff = len(uniq_teff)
    if not(1 <= num_logg <= len(uniq_logg)):
        num_logg = len(uniq_logg)

    # target grid to interpolate to
    xi = np.linspace(np.amin(uniq_teff), np.amax(uniq_teff), num_teff)
    yi = np.linspace(np.amin(uniq_logg), np.amax(uniq_logg), num_logg)
    zi = target_wave

    print(len(xi))
    print(len(yi))
    print(len(zi))
    print("meshgrid...",flush=True)
    xi,yi,zi = np.meshgrid(xi,yi,zi, indexing='ij')

    print(np.shape(xi))
    print("interpolate...",flush=True)

    # interpolate
    vi = interp.griddata((flat_teff,flat_logg,flat_wave),flat_flux,(xi,yi,zi),method='linear')

    print("done")
    return (vi, xi, yi, 10.0**zi)


def create_dataset(itp:tuple):
    '''
    Write interpolated grid to NetCDF file using Xarray.

    Parameters
       itp (tuple): tuple of four 3D output arrays from create_interp()

    Returns
        ds (xr.Dataset): dataset of interpolated data
        
    '''

    # get coords
    v,x,y,z = itp[0],itp[1],itp[2],itp[3]

    arr_teff = np.unique(x)
    arr_logg = np.unique(y)
    arr_wave = z[0,0,:]

    coords={
            "teff": arr_teff,
            "logg": arr_logg,
            "wave": arr_wave,
            }
    
    
    # get data
    data = {"flux":xr.DataArray(data=v, coords=coords, dims=["teff","logg","wave"])}
    
    # create dataset
    ds = xr.Dataset(data_vars=data, coords=coords)

    return ds

def write_dataset(itp:tuple):
    '''
    Write interpolated grid to NetCDF file using Xarray.

    Parameters
       ds (xr.Dataset): dataset of interpolated data

    Returns
        None
        
    '''

    # create dataset
    ds = create_dataset(itp)    

    # save
    fpath = os.path.join(utils.dirs["data"], "btsettl_interp.nc")
    if os.path.exists(fpath):
        os.remove(fpath)
    ds.to_netcdf(fpath)
    return 

def get_spec_from_interp(itp:tuple, teff:float, logg:float):
    '''
    Get stellar spectrum from interpolated grid, searching for best point.

    Parameters
       itp (tuple): tuple of four 3D output arrays from create_interp()
       teff (float): target teff
       logg (float): target logg

    Returns
        wl (np.ndarray): wavelengths [nm]
        fl (np.ndarray): spectral fluxes [erg s-1 cm-2 nm-1]
    '''

    v,x,y,z = itp[0], itp[1], itp[2], itp[3]
    sh = np.shape(v)

    iclose = -1 
    jclose = -1
    dclose = 1e99
    for i in range(sh[0]):
        for j in range(sh[1]):
            t = x[i,j,0]
            l = y[i,j,0]
            d = np.sqrt(((teff-t)/teff)**2 + ((logg-l)/logg)**2) * 100  # distance [%]
            if d < dclose:
                iclose = i
                jclose = j
                dclose = d
    
    return z[iclose,jclose,:], v[iclose,jclose,:]

def get_axes():
    '''
    Get coordinates of logg and teff axes from files in data folder.

    Parameters
        None

    Returns
        arr_teff (np.ndarray): teff points
        arr_logg (np.ndarray): logg points
    '''

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

def get_spec_from_npy(teff, logg, xmax=1.0e9):
    '''
    Read stellar spectrum from npy file, locating it according to teff and logg.

    Parameters
        teff (float): effective temperature [K]
        logg (float): log surface gravity [log cm/s^2]

    Returns
        wl (np.ndarray): wavelengths [nm]
        fl (np.ndarray): spectral fluxes [erg s-1 cm-2 nm-1]
    '''

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

    # return spectrum
    return wl,fl

def extend_planck(teff, wl, fl, xmax=1.0e9):
    dx = wl[-1]*0.1
    print(dx)
    xp = np.arange(wl[-1]+dx, xmax, dx)
    yp = np.zeros(np.shape(xp))
    for i,x in enumerate(xp):
        lam = x * 1.0e-9  # nm -> m

        # Calculate planck function value [W m-2 sr-1 m-1]
        # http://spiff.rit.edu/classes/phys317/lectures/planck.html
        yp[i] = 2.0 * utils.h_pl * utils.c_vac**2.0 / lam**5.0   *   1.0 / ( np.exp(utils.h_pl * utils.c_vac / (lam * utils.k_B * teff)) - 1.0)

        # Integrate solid angle (hemisphere), convert units
        yp[i] = yp[i] * np.pi * 1.0e-9 # [W m-2 nm-1]
        yp[i] = yp[i] * 1000.0 # [erg s-1 cm-2 nm-1]

    wl = np.concatenate((wl, xp))
    fl = np.concatenate((fl, yp))

    return wl,fl
